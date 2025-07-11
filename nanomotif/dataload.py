# Data loading functionalities
import polars as pl
from .seq import Assembly
from .feature import Pileup
import sys

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

def load_pileup(path: str, min_coverage: int, min_fraction: float = 0):
    """
    Load pileup file from path to pileup.bed output of modkit pileup
    """
    pileup = (
        pl.scan_csv(
            path, 
            separator = "\t", 
            has_header = False,
            schema = PILEUP_SCHEMA
        )

        .filter(pl.col("column_10") > min_coverage)
        .filter(pl.col("column_11") > min_fraction*100)
        .filter(pl.col("column_10") / (pl.col("column_10") + pl.col("column_17")) > 0.3)
        .select(["column_1", "column_2","column_4", "column_6", "column_11", "column_10"])
        .with_columns(pl.col("column_11") / 100)
        .collect()
    )

    if pileup.is_empty():
        # TODO: import logger for this
        print("Pileup is empty after initial load")
        sys.exit(1)
    pileup = pileup.rename({"column_1":"contig", "column_2": "position", "column_4": "mod_type", "column_6": "strand", "column_11": "fraction_mod", "column_10":"Nvalid_cov"})
    return Pileup(pileup)

def extract_contig_mods_with_sufficient_information(pileup: Pileup, assembly: Assembly, min_mods_pr_contig: int, min_mod_frequency: int):
    contigs_in_assembly = list(assembly.keys())
    pileup = pileup.filter(pl.col("contig").is_in(contigs_in_assembly))

    if pileup.is_empty():
        #TODO: 
        print("Pileup empty after filtering contigs in assembly. Check pileup and assembly mismatch!")
        sys.exit(1)

    contigmods_with_more_than_min_mods = pileup.group_by("contig_mod").count().filter(
            pl.col("count") > min_mods_pr_contig
        ).get_column("contig_mod").to_list()

    assm_lengths = pl.DataFrame({
        "contig": list(assembly.keys()),
        "length": [len(contig) for contig in assembly.values()]
    })
    contigmods_with_sufficient_mod_frequency = pileup \
        .group_by(["contig", "mod_type"]) \
        .agg(pl.count()) \
        .join(assm_lengths, on = "contig") \
        .filter(pl.col("count") > pl.col("length")/min_mod_frequency) \
        .with_columns([
            (pl.col("contig") + "_" + pl.col("mod_type")).alias("contig_mod")
        ]) \
        .get_column("contig_mod").unique().to_list()
    
    # Get contig_mods to keep and to remove
    contig_mods_to_keep = list(set(contigmods_with_more_than_min_mods) & set(contigmods_with_sufficient_mod_frequency))
    contig_mods_to_remove = list(set(pileup.get_column("contig_mod").unique().to_list()) - set(contig_mods_to_keep))
    return contig_mods_to_keep, contig_mods_to_remove


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
            schema = PILEUP_SCHEMA
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
    return Assembly(load_fasta(path))
