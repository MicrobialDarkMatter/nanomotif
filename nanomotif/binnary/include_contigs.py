import pandas as pd
import polars as pl
import numpy as np
from nanomotif.binnary import data_processing as dp 
from nanomotif.binnary import utils as ut
import logging
import hdbscan



def include_contigs(contig_methylation, contig_lengths):
    logger = logging.getLogger(__name__)
    logger.info("Starting include_contigs analysis...")
    
    binned_contig_methylation = dp.impute_contig_methylation_within_bin(
        contig_methylation
    )

    bins_w_methylation = binned_contig_methylation\
        .group_by(["bin", "motif_mod"])\
        .agg(
            pl.col("median").mean().alias("median")
        )\
        .with_columns(
            (pl.col("median") > 0.5).cast(pl.Int8).alias("is_methylated")
        )\
        .filter(
            pl.col("is_methylated") == 1
        )\
        .get_column("bin").unique()

    binned_contig_methylation = binned_contig_methylation\
        .filter(pl.col("bin").is_in(bins_w_methylation))

    unbinned_contig_methylation = dp.impute_unbinned_contigs(contig_methylation)

    contig_names, matrix = dp.create_matrix(contig_methylation)

    logger.info("Performing HDBSCAN")
    clusterer = hdbscan.HDBSCAN(min_samples=10, min_cluster_size=2, metric = "euclidean", allow_single_cluster=False)
    labels = clusterer.fit_predict(matrix)

    contig_bin = contig_methylation\
        .select(["bin", "contig"])\
        .unique() 

    results = pl.DataFrame({
                   "contig": contig_names,
                   "cluster": labels
              })\
                .join(contig_bin, on="contig")\

    bin_size = contig_bin\
        .filter(pl.col("bin") != "unbinned")\
        .join(contig_lengths, on = "contig")\
        .group_by("bin")\
        .agg(
            pl.col("length").sum().alias("bin_length"),
            pl.col("contig").count().alias("n_contigs_bin")
        )

    bin_cluster_sizes = results\
        .filter(pl.col("bin") != "unbinned")\
        .join(contig_lengths, on="contig")\
        .group_by(["bin","cluster"])\
        .agg(
            pl.col("contig").count().alias("n_contigs"),
            pl.col("length").sum().alias("cluster_length")
        )\
        .join(bin_size, on = "bin")\
        .with_columns(
            (pl.col("n_contigs") / pl.col("n_contigs_bin")).alias("fraction_contigs"),
            (pl.col("cluster_length") / pl.col("bin_length")).alias("fraction_length")
        )

    assigned_cluster = bin_cluster_sizes\
        .group_by(["bin"])\
        .agg(
            pl.col("cluster_length").max().alias("cluster_length")
        )\
        .join(bin_cluster_sizes, on = ["bin", "cluster_length"], how = "left")\
        .rename({
                    "bin": "assigned_bin"
                })\
        .drop(["cluster_length", "n_contigs"])

    unbinned_contigs_assigned_to_bin_cluster = results\
        .filter(pl.col("bin") == "unbinned")\
        .join(assigned_cluster, on = ["cluster"])

    uniquely_assigned_contigs = unbinned_contigs_assigned_to_bin_cluster\
        .group_by("contig")\
        .agg(
            pl.col("bin").count().alias("n_assigned_bins")
        )\
        .filter(pl.col("n_assigned_bins") == 1)\
        .get_column("contig")

    unbinned_contigs_assigned_to_bin_cluster = unbinned_contigs_assigned_to_bin_cluster\
        .with_columns(
            pl.when(pl.col("contig").is_in(uniquely_assigned_contigs))
                .then(True)
                .otherwise(False)
                .alias("assignment_is_unique")
        )
    return unbinned_contigs_assigned_to_bin_cluster
