import polars as pl
import re

def split_bin_contig(df):
    def custom_replace(bin_compare, contig_bin):
        if bin_compare is not None and contig_bin is not None:
            pattern = re.escape(contig_bin) + "_"
            return re.sub(pattern, "", bin_compare)
        return bin_compare

    df = df.with_columns(
        pl.when(pl.col("bin_compare").is_not_null()).then(
            pl.struct(["bin_compare", "contig_bin"]).apply(
                lambda x: custom_replace(x["bin_compare"], x["contig_bin"])
            )
        ).alias("contig")
    )
    
        # .with_columns([
        #     pl.col("bin_compare").str.replace(pl.col("bin") + "_", "").alias("contig_bin")
        # ])

    return df


def get_bin(bin_contig):
    return re.search(r"(.+)_contig_\d+", bin_contig).group(1)


def add_compare_df(df, bin_contig):
    df = df.filter(pl.col("bin_compare") == bin_contig)\
            .rename({"bin": "contig_bin"})
    
    return df