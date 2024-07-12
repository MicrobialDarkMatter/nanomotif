import polars as pl
import re

def split_bin_contig(df):
    df = df \
        .with_columns(
            pl.when(pl.col("bin_compare").is_not_null()).then(
                pl.col("bin_compare").str.replace(pl.col("contig_bin").fill_null("").first() + "_", "", literal = True)
            ).over("bin_compare").alias("contig")
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