"""Optional cache of the filtered pileup partitions (`--pileup_cache FILE`).

The pileup pass (coverage, minimum-modification and adjacency filters, split by bin and mod type) is the
same for every run with the same pileup, contig bins and thresholds. With a cache file, its result is kept
in one lz4 Arrow IPC file: one or more record batches per (bin, mod type), in any order. The last batch is
a one-row trailer holding the manifest (what the pass was made from) and the index (which batches belong
to which bin and mod type), so a later run reads only its own batches, memory-mapped, and no other file
is needed. (Not footer metadata: pyarrow < 20 can neither write nor read it.) nanomotif-rs writes the file during the pileup pass; without it, the in-memory partitions are
written with pyarrow.
"""
import hashlib
import json
import logging as log
import os

import polars as pl
import pyarrow as pa
import pyarrow.ipc

from nanomotif._version import __version__

SCHEMA = pa.schema([
    ("contig", pa.utf8()),
    ("position", pa.int64()),
    ("mod_type", pa.utf8()),
    ("strand", pa.utf8()),
    ("fraction_mod", pa.float64()),
    ("Nvalid_cov", pa.int64()),
    ("bin", pa.utf8()),
])
# trailer batch: contig TRAILER_TAG, mod_type the manifest (JSON), strand the index (JSON)
TRAILER_TAG = "nanomotif_cache"
# rows per record batch when writing in-memory partitions (the Rust partitioner writes up to 500k rows)
BATCH_ROWS = 500_000


def make_manifest(pileup_path: str, contig_bin: dict, min_coverage: int, threshold: float) -> dict:
    """What the partitions depend on. A cache is reused only when all of it is equal."""
    stat = os.stat(pileup_path)
    contig_bin_hash = hashlib.sha256(
        "\n".join(f"{contig}\t{bin_name}" for contig, bin_name in sorted(contig_bin.items())).encode()
    ).hexdigest()
    return {
        "nanomotif_version": __version__,
        "pileup": os.path.abspath(pileup_path),
        "pileup_size": stat.st_size,
        "pileup_mtime_ns": stat.st_mtime_ns,
        "contig_bin_sha256": contig_bin_hash,
        "threshold_valid_coverage": min_coverage,
        "methylation_threshold_high": threshold,
    }


def load(cache_path: str, manifest: dict):
    """The cache's partitions as {(bin, mod_type): (cache_path, [batch numbers])}, or None if it can't be used."""
    if not os.path.exists(cache_path):
        log.info(f"Pileup cache {cache_path} does not exist yet; it is written during the pileup pass")
        return None
    try:
        with pa.memory_map(cache_path) as source:
            reader = pa.ipc.open_file(source)
            trailer = reader.get_batch(reader.num_record_batches - 1).to_pydict()
        if trailer["contig"] != [TRAILER_TAG]:
            raise ValueError("no nanomotif trailer")
        cached_manifest = json.loads(trailer["mod_type"][0])
        index = json.loads(trailer["strand"][0])
    except (OSError, KeyError, ValueError, pa.ArrowException) as error:
        log.warning(f"Pileup cache {cache_path} is not readable ({error}); rebuilding it")
        return None
    differences = [key for key in manifest if cached_manifest.get(key) != manifest[key]]
    if differences:
        log.info(f"Pileup cache {cache_path} was made with a different {', '.join(differences)}; rebuilding it")
        return None
    log.info(f"Reusing the pileup partitions in {cache_path} (skipping the pileup pass)")
    return _partitions(cache_path, index)


def _partitions(cache_path, index):
    return {(bin_name, mod_type): (cache_path, batches) for bin_name, mod_type, batches in index}


def _finish(tmp_path, cache_path):
    # only a completely written cache replaces an older one
    os.replace(tmp_path, cache_path)
    log.info(f"Wrote the pileup partitions to {cache_path} ({os.path.getsize(cache_path) / 1e9:.2f} GB)")


def _tmp_path(cache_path):
    os.makedirs(os.path.dirname(os.path.abspath(cache_path)), exist_ok=True)
    return cache_path + ".tmp"


def write_with_rust(nanomotif_rs, cache_path: str, manifest: dict, pileup_path: str, contig_bin: dict,
                    min_coverage: int, threshold: float, threads: int) -> dict:
    """The pileup pass of nanomotif-rs, writing its partitions straight into the cache (no copy)."""
    tmp_path = _tmp_path(cache_path)
    batches = nanomotif_rs.partition_pileup_to_file(
        pileup_path, contig_bin, tmp_path, json.dumps(manifest),
        min_coverage=min_coverage, threshold=threshold, threads=threads,
    )
    _finish(tmp_path, cache_path)
    return {key: (cache_path, value) for key, value in batches.items()}


def write_from_frames(cache_path: str, manifest: dict, partitions: dict) -> dict:
    """Write in-memory partitions ({(bin, mod_type): polars DataFrame}) into the cache."""
    keys = sorted(partitions)
    index = []
    n_batches = 0
    for key in keys:
        count = -(-len(partitions[key]) // BATCH_ROWS)
        index.append([key[0], key[1], list(range(n_batches, n_batches + count))])
        n_batches += count
    trailer = pa.RecordBatch.from_pydict({
        "contig": [TRAILER_TAG], "position": [0], "mod_type": [json.dumps(manifest)], "strand": [json.dumps(index)],
        "fraction_mod": [0.0], "Nvalid_cov": [0], "bin": [""],
    }, schema=SCHEMA)
    tmp_path = _tmp_path(cache_path)
    options = pa.ipc.IpcWriteOptions(compression="lz4")
    with pa.OSFile(tmp_path, "wb") as sink, pa.ipc.new_file(sink, SCHEMA, options=options) as writer:
        for key in keys:
            frame = partitions[key]
            for offset in range(0, len(frame), BATCH_ROWS):
                # oldest compat level: plain/large strings, which every supported pyarrow can cast
                table = frame.slice(offset, BATCH_ROWS).to_arrow(compat_level=pl.CompatLevel.oldest())
                writer.write_batch(table.select(SCHEMA.names).cast(SCHEMA).combine_chunks().to_batches()[0])
        writer.write_batch(trailer)
    _finish(tmp_path, cache_path)
    return _partitions(cache_path, index)


def read_partition(cache_path: str, batches: list) -> pl.DataFrame:
    """One (bin, mod type) partition from the cache: only its batches are read."""
    with pa.memory_map(cache_path) as source:
        reader = pa.ipc.open_file(source)
        table = pa.Table.from_batches([reader.get_batch(i) for i in batches], schema=reader.schema)
    return pl.from_arrow(table)
