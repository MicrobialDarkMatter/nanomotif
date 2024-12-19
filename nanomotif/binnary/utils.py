import subprocess
import os
import nanomotif
import platform
import polars as pl
import re
import logging

def split_bin_contig(df):
    def custom_replace(bin_compare, contig_bin):
        if bin_compare is not None and contig_bin is not None:
            pattern = re.escape(contig_bin) + "_"
            return re.sub(pattern, "", bin_compare, count=1)
        return bin_compare

    df = df.with_columns(
        pl.when(pl.col("bin_compare").is_not_null()).then(
            pl.struct(["bin_compare", "contig_bin"]).map_elements(
                lambda x: custom_replace(x["bin_compare"], x["contig_bin"])
            , return_dtype = pl.String)
        ).alias("contig")
    )
    return df


def get_bin(bin_contig):
    return re.search(r"(.+)_contig_\d+", bin_contig).group(1)


def add_compare_df(df, bin_contig):
    df = df.filter(pl.col("bin_compare") == bin_contig)\
            .rename({"bin": "contig_bin"})
    
    return df

