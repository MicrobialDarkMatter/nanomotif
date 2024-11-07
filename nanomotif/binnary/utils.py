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
            return re.sub(pattern, "", bin_compare)
        return bin_compare

    df = df.with_columns(
        pl.when(pl.col("bin_compare").is_not_null()).then(
            pl.struct(["bin_compare", "contig_bin"]).map_elements(
                lambda x: custom_replace(x["bin_compare"], x["contig_bin"])
            , return_dtype = pl.String)
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


# def run_methylation_utils(
#     pileup,
#     assembly,
#     motifs,
#     threads,
#     min_valid_read_coverage,
#     output
# ):
#     logger = logging.getLogger(__name__)
#     system = platform.system()
#     tool = f"methylation_utils"
#     if system == "Windows":
#         tool += ".exe"
        
#     env = os.environ.copy()
#     env["POLARS_MAX_THREADS"] = str(threads)

#     bin_path = os.path.join(os.path.dirname(nanomotif.__file__), "bin", tool)
#     logger.info("Running methylation_utils")
#     try:


#         cmd_args = [
#             "--pileup", pileup,
#             "--assembly", assembly,
#             "--motifs", *motifs,
#             "--threads", str(threads),
#             "--min-valid-read-coverage", str(min_valid_read_coverage),
#             "--output", os.path.join(output, "motifs-scored-read-methylation.tsv")
#         ]
        
#         subprocess.run([bin_path] + cmd_args, check = True, env = env)
#     except subprocess.CalledProcessError as e:
#         print(f"Error: Command '{e.cmd}' failed with return code {e.returncode}")
#         print(f"Output: {e.output}")
#         return e.returncode
