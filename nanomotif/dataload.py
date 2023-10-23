# Data loading functionalities
import polars as pl
from .seq import Assembly
from .feature import Pileup

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

def load_pileup(path: str):
    """
    Load pileup file from path to pileup.bed output of modkit pileup
    """
    pileup = pl.read_csv(path, separator = "\t", has_header = False)
    pileup = pileup.filter(pl.col("column_10") / (pl.col("column_10") + pl.col("column_17")) > 0.1)
    pileup = pileup.select(["column_1", "column_2","column_4", "column_6", "column_11", "column_10"]) \
        .rename({"column_1":"contig", "column_2": "position", "column_4": "mod_type", "column_6": "strand", "column_11": "fraction_mod", "column_10":"Nvalid_cov"}) \
        .with_columns(pl.col("fraction_mod") / 100) \
        .sort("position")
    return Pileup(pileup)

def load_assembly(path: str):
    """
    Load assembly from path to fasta file
    """
    return Assembly(load_fasta(path))
