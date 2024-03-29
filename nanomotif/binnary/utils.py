import polars as pl
import re

def split_bin_contig(df):
    df = df \
        .with_columns([
            pl.col("bin_compare").str.extract(r"(.+)_contig_\d+", 1).alias("contig_bin"),
            pl.col("bin_compare").str.extract(r"contig_(\d+)", 1).alias("contig_number")
        ]) \
        .with_columns(
            ("contig_" + pl.col("contig_number")).alias("contig")
        ) \
        .drop("contig_number")

    return df


def get_bin(bin_contig):
    return re.search(r"(.+)_contig_\d+", bin_contig).group(1)