import numpy as np
# Data loading functionalities
import polars as pl

from nanomotif.seq import Assembly
import io
import sys
import pysam
import pyfastx
import nanomotif as nm
import logging as log
import time
import random


PILEUP_SCHEMA = {
    "column_1": pl.Utf8,
    "column_2": pl.Int64,
    "column_3": pl.Int64,
    "column_4": pl.Utf8,
    "column_5": pl.Int64,
    "column_6": pl.Utf8,
    "column_7": pl.Int64,
    "column_8": pl.Int64,
    "column_9": pl.Utf8,
    "column_10": pl.Int64,
    "column_11": pl.Float64,
    "column_12": pl.Int64,
    "column_13": pl.Int64,
    "column_14": pl.Int64,
    "column_15": pl.Int64,
    "column_16": pl.Int64,
    "column_17": pl.Int64,
    "column_18": pl.Int64
}

def load_fasta(path, trim_names=False, trim_character=" ") -> dict:
    """
    Reads a fasta file and returns a dictionary with the contig names as 
    keys and the sequences as values
    """
    with open(path, 'r') as f:
        lines = f.readlines()
    data = {}
    active_sequence_name = "no_header"
    for line in lines:
        line = line.rstrip()
        if line[0] == '>':
            active_sequence_name = line[1:]
            if trim_names:
                active_sequence_name = active_sequence_name.split(trim_character)[0]
            if active_sequence_name not in data:
                data[active_sequence_name] = ''
        else:
            data[active_sequence_name] += line
    return data

def load_fasta_fastx(path, trim_names=False, trim_character=" ") -> dict:
    """
    Reads a fasta file and returns a dictionary with the contig names as 
    keys and the sequences as values
    """
    fx = pyfastx.Fasta(path, build_index=False)
    data = {}
    for name, seq in fx.items():
        if trim_names:
            name = name.split(trim_character)[0]
        data[name] = seq
    assembly = nm.seq.Assembly(data)
    return assembly


def load_pileup(path: str, min_coverage: int = None, contigs: list = None):
    """
    Load pileup file from path to pileup.bed output of modkit pileup

    min_coverage and contigs (optional) are applied while reading, so rows that would be
    filtered out later are never held in memory.
    """
    pileup = (
        pl.scan_csv(
            path, 
            separator = "\t", 
            has_header = False,
            schema = PILEUP_SCHEMA,
            null_values=["NA", "null"]
        )
        .select(["column_1", "column_2","column_4", "column_6", "column_11", "column_10"])
    )
    if min_coverage is not None:
        pileup = pileup.filter(pl.col("column_10") >= min_coverage)
    if contigs is not None:
        pileup = pileup.filter(pl.col("column_1").is_in(pl.Series(list(contigs), dtype=pl.Utf8).implode()))
    pileup = pileup.with_columns(pl.col("column_11") / 100).collect()

    if pileup.is_empty() and min_coverage is None and contigs is None:
        print("Pileup is empty after initial load")
        sys.exit(1)
    pileup = pileup.rename({
        "column_1":"contig", 
        "column_2": "position", 
        "column_4": "mod_type", 
        "column_6": "strand", 
        "column_11": "fraction_mod", 
        "column_10": "Nvalid_cov"
    })
    return pileup

def load_contigs_pileup_bgzip(path: str, contigs: list[str]):
    """
    Rows of the given contigs from a bgzipped, tabix-indexed modkit pileup (pysam), with the columns
    of load_pileup.
    """
    log.debug(f"Querying pileup for {len(contigs)} contigs")
    with pysam.TabixFile(path) as tbx:
        indexed = set(tbx.contigs)
        lines = [line for contig in contigs if contig in indexed for line in tbx.fetch(contig)]
    if not lines:
        log.warning(f"No pileup data found for contigs: {contigs}")
        return pl.DataFrame(schema={
            "contig": pl.String,
            "position": pl.Int64,
            "mod_type": pl.Utf8,
            "strand": pl.Utf8,
            "fraction_mod": pl.Float64,
            "Nvalid_cov": pl.Int64,
        })
    pileup = pl.read_csv(
        io.StringIO("\n".join(lines)),
        separator="\t",
        has_header=False,
        schema=PILEUP_SCHEMA,
        null_values=["NA", "null"],
    )
    return pileup \
        .select(["column_1", "column_2", "column_4", "column_6", "column_11", "column_10"]) \
        .with_columns(pl.col("column_11") / 100) \
        .rename({
            "column_1": "contig",
            "column_2": "position",
            "column_4": "mod_type",
            "column_6": "strand",
            "column_11": "fraction_mod",
            "column_10": "Nvalid_cov",
        })



def load_low_coverage_positions(path_pileup: str, contig_mods_to_load: list[str], min_coverage: float = 5):
    """
    Load pileup file from path to pileup.bed output of modkit pileup
    """
    pileup = (
        pl.scan_csv(
            path_pileup, 
            separator = "\t", 
            has_header = False,
            # Schema overrides to match the expected columns
            schema = PILEUP_SCHEMA,
            null_values=["NA", "null"]
        )
        .filter(pl.col("column_10") <= min_coverage)
        .filter(pl.col("column_10") / (pl.col("column_10") + pl.col("column_17")) > 0.3)
        .select(["column_1", "column_2","column_4", "column_6", "column_10"])
        .with_columns((pl.col("column_1") + "_" + pl.col("column_4")).alias("contig_mod"))
        .filter(pl.col("contig_mod").is_in(contig_mods_to_load))
        .group_by("column_1", "column_4")
        .agg(
            pl.col("column_2"),
            pl.col("column_6")
        )
        .collect()
    )
    pileup = pileup.rename({"column_1":"contig", "column_2": "position", "column_4": "mod_type", "column_6": "strand"})
    return pileup


def load_assembly(path: str):
    """
    Load assembly from path to fasta file
    """
    return Assembly(load_fasta_fastx(path))

def filter_pileup(
        pileup: pl.DataFrame,
        min_coverage: int = 5,
): 
    """
    Filter pileup to only include positions with at least min_coverage valid reads (Nvalid_cov)
    """
    pileup = pileup.filter(pl.col("Nvalid_cov") >= min_coverage)
    return pileup

def filter_pileup_minimummod_frequency(
        pileup: pl.DataFrame, 
        methylation_threshold: float = 0.7,
        min_mod_frequency: int = 0.0001,
        min_mods_pr_contig: int = 50
):
    """
    Filter pileup to only include contig_mods with sufficient modification frequency and number of modified positions
    """
    contigs_with_mods = pileup \
        .group_by("contig", "mod_type") \
        .agg(
            n_positions=pl.len(),
            n_mod_positions=(pl.col("fraction_mod") >= methylation_threshold).sum()
        ) \
        .filter((pl.col("n_mod_positions") / pl.col("n_positions")) > min_mod_frequency) \
        .filter(pl.col("n_mod_positions") > min_mods_pr_contig) \
        .select("contig", "mod_type")

    return pileup.join(contigs_with_mods, on=["contig", "mod_type"], how="semi", maintain_order="left")

def filter_pileup_adjacency_filter(
        pileup: pl.DataFrame,
        methylation_threshold: float = 0.7,
        adjacency_distance: int = 8
):
    """
    Filter pileup to only include positions that are not within adjacency_distance of another modified position
    """
    # Rows below the threshold are always kept, and a row at or above it can only be removed by a
    # higher neighbour, which is then above it too: only the high rows need comparing. A high row is
    # kept when no other high row of its contig and strand within adjacency_distance is higher. With
    # them sorted by position, compare every row with its k-th neighbour for k = 1, 2, ... until no
    # k-th neighbour is within the distance.
    high = pileup.with_row_index("_row") \
        .filter(pl.col("fraction_mod") >= methylation_threshold) \
        .sort(["contig", "strand", "position"])
    group = high.select(pl.struct("contig", "strand").rle_id()).to_series().to_numpy()
    position = high.get_column("position").cast(pl.Int64).to_numpy()
    fraction = high.get_column("fraction_mod").to_numpy()
    is_max = np.ones(len(high), dtype=bool)
    k = 1
    while k < len(high):
        near = (group[k:] == group[:-k]) & (position[k:] - position[:-k] <= adjacency_distance)
        if not near.any():
            break
        is_max[:-k] &= ~(near & (fraction[k:] > fraction[:-k]))
        is_max[k:] &= ~(near & (fraction[:-k] > fraction[k:]))
        k += 1
    removed = high.get_column("_row").filter(pl.Series(~is_max))
    pileup = pileup.with_row_index("_row").filter(~pl.col("_row").is_in(removed.implode())).drop("_row")
    return pileup
