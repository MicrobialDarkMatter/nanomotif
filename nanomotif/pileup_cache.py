"""Optional cache of the filtered pileup partitions (`--pileup_cache FILE`).

The pileup pass (coverage, minimum-modification and adjacency filters, split by bin and mod type) is the
same for every run with the same pileup, contig bins and thresholds. With a cache file, its result is kept
in one lz4 Arrow IPC file: one or more record batches per (bin, mod type). The file's schema metadata holds
the manifest (what the pass was made from) and the index (which batches belong to which bin and mod type),
so a later run reads only its own batches, memory-mapped, and no other file is needed.
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
MANIFEST_KEY = b"nanomotif_manifest"
INDEX_KEY = b"nanomotif_index"
# rows per record batch when splitting in-memory partitions (the Rust partitioner writes 500k-row batches)
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
    """The cache's partitions as {(bin, mod_type): (cache_path, start, stop)}, or None if it can't be used."""
    if not os.path.exists(cache_path):
        log.info(f"Pileup cache {cache_path} does not exist yet; it is written after the pileup pass")
        return None
    try:
        with pa.memory_map(cache_path) as source:
            metadata = pa.ipc.open_file(source).schema.metadata or {}
        cached_manifest = json.loads(metadata[MANIFEST_KEY])
        index = json.loads(metadata[INDEX_KEY])
    except (OSError, KeyError, ValueError, pa.ArrowException) as error:
        log.warning(f"Pileup cache {cache_path} is not readable ({error}); rebuilding it")
        return None
    differences = [key for key in manifest if cached_manifest.get(key) != manifest[key]]
    if differences:
        log.info(f"Pileup cache {cache_path} was made with a different {', '.join(differences)}; rebuilding it")
        return None
    log.info(f"Reusing the pileup partitions in {cache_path} (skipping the pileup pass)")
    return {(bin_name, mod_type): (cache_path, start, stop) for bin_name, mod_type, start, stop in index}


def _batches_of_file(path):
    # plain reads, not a memory map: mapped pages of all partition files would add to the resident memory
    with pa.OSFile(path) as source:
        reader = pa.ipc.open_file(source)
        for i in range(reader.num_record_batches):
            yield reader.get_batch(i)


def _n_batches_of_file(path):
    with pa.OSFile(path) as source:
        return pa.ipc.open_file(source).num_record_batches


def _to_cache_schema(batch):
    table = pa.Table.from_batches([batch]).select(SCHEMA.names).cast(SCHEMA)
    return table.combine_chunks().to_batches()[0] if table.num_rows else None


def write_from_files(cache_path: str, manifest: dict, partitions: dict) -> dict:
    """Copy partition files ({(bin, mod_type): [IPC files]}) into the cache, batch by batch, deleting each
    file once copied (so the disk holds the data about once). Returns the cache's partitions."""
    keys = sorted(partitions)
    counts = {key: sum(_n_batches_of_file(path) for path in partitions[key]) for key in keys}

    def batches(key):
        for path in partitions[key]:
            yield from _batches_of_file(path)
            os.remove(path)

    return _write(cache_path, manifest, keys, counts, batches)


def write_from_frames(cache_path: str, manifest: dict, partitions: dict) -> dict:
    """Write in-memory partitions ({(bin, mod_type): polars DataFrame}) into the cache."""
    keys = sorted(partitions)
    counts = {key: -(-len(partitions[key]) // BATCH_ROWS) for key in keys}

    def batches(key):
        frame = partitions[key]
        for offset in range(0, len(frame), BATCH_ROWS):
            # oldest compat level: plain/large strings, which every supported pyarrow can cast
            table = frame.slice(offset, BATCH_ROWS).to_arrow(compat_level=pl.CompatLevel.oldest())
            yield table.combine_chunks().to_batches()[0]

    return _write(cache_path, manifest, keys, counts, batches)


def _copy_pool():
    # mimalloc (pyarrow's default on Linux) keeps the memory of freed batches: copying ZymoOral's
    # partitions raised the process from 1.1 to 3.3 GB resident, with at most 0.2 GB allocated. jemalloc
    # returns it (no rise above the pileup pass's peak); the system allocator is the fallback.
    try:
        return pa.jemalloc_memory_pool()
    except NotImplementedError:
        return pa.system_memory_pool()


def _write(cache_path, manifest, keys, counts, batches):
    previous_pool = pa.default_memory_pool()
    pa.set_memory_pool(_copy_pool())
    try:
        return _write_batches(cache_path, manifest, keys, counts, batches)
    finally:
        pa.set_memory_pool(previous_pool)


def _write_batches(cache_path, manifest, keys, counts, batches):
    index = []
    start = 0
    for key in keys:
        index.append([key[0], key[1], start, start + counts[key]])
        start += counts[key]
    schema = SCHEMA.with_metadata({MANIFEST_KEY: json.dumps(manifest), INDEX_KEY: json.dumps(index)})
    tmp_path = cache_path + ".tmp"
    os.makedirs(os.path.dirname(os.path.abspath(cache_path)), exist_ok=True)
    options = pa.ipc.IpcWriteOptions(compression="lz4")
    with pa.OSFile(tmp_path, "wb") as sink, pa.ipc.new_file(sink, schema, options=options) as writer:
        for key in keys:
            written = 0
            for batch in batches(key):
                batch = _to_cache_schema(batch)
                if batch is None:
                    # an empty batch still counts in the index; keep the count right
                    batch = pa.RecordBatch.from_pylist([], schema=SCHEMA)
                writer.write_batch(batch)
                written += 1
            assert written == counts[key], f"expected {counts[key]} batches for {key}, wrote {written}"
    # only a completely written cache replaces an older one
    os.replace(tmp_path, cache_path)
    log.info(f"Wrote the pileup partitions to {cache_path} ({os.path.getsize(cache_path) / 1e9:.2f} GB)")
    return {(bin_name, mod_type): (cache_path, s, e) for bin_name, mod_type, s, e in index}


def read_partition(cache_path: str, start: int, stop: int) -> pl.DataFrame:
    """One (bin, mod type) partition from the cache: only its batches are read."""
    with pa.memory_map(cache_path) as source:
        reader = pa.ipc.open_file(source)
        table = pa.Table.from_batches([reader.get_batch(i) for i in range(start, stop)], schema=SCHEMA)
    return pl.from_arrow(table)
