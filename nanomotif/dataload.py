# Data loading functionalities
import polars as pl
from .seq import Assembly
from epymetheus import query_pileup_records, PileupColumn
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


def load_pileup(path: str):
    """
    Load pileup file from path to pileup.bed output of modkit pileup
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
        .with_columns(pl.col("column_11") / 100)
        .collect()
    )

    if pileup.is_empty():
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
    Load pileup file from path to pileup.bed output of modkit pileup.
    Retries opening the tabix file and reading records if initial attempts fail.
    """
    # Step 2: query pileup from Rust
    log.debug(f"Querying pileup for {len(contigs)} contigs")
    pileup = query_pileup_records(
        path,
        contigs,
        columns = [
            PileupColumn.Contig,
            PileupColumn.Start,
            PileupColumn.ModType,
            PileupColumn.Strand,
            PileupColumn.FractionModified,
            PileupColumn.NValidCov
        ]
    )
    log.debug(f"Renaming and removing unnecessary columns")
    pileup = pileup.rename({
        "contig": "contig",
        "start": "position",
        "mod_type": "mod_type",
        "strand": "strand",
        "fraction_modified": "fraction_mod",
        "n_valid_cov": "Nvalid_cov",
    }) \
    .with_columns(pl.col("fraction_mod") / 100)

    return pileup



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
        min_modtype_fraction: float = 0.3,
        min_coverage: int = 5,
): 
    """
    Filter pileup to only include positions with sufficient modification fraction and coverage
    """
    pileup = pileup.filter(pl.col("Nvalid_cov") > min_coverage)
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
    pileup = pileup.with_columns([
        (pl.col("contig") + "_" + pl.col("mod_type")).alias("contig_mod")
    ])
    contigs_with_mods = pileup \
        .group_by("contig_mod") \
        .agg(
            n_positions=pl.count().alias("n_positions"),
            n_mod_positions=(pl.col("fraction_mod") > methylation_threshold).sum()
        ) \
        .filter((pl.col("n_mod_positions") / pl.col("n_positions")) > min_mod_frequency) \
        .filter(pl.col("n_mod_positions") > min_mods_pr_contig) \
        .get_column("contig_mod").unique().to_list()

    pileup = pileup.filter(pl.col("contig_mod").is_in(contigs_with_mods)) \
        .drop("contig_mod")
    return pileup
