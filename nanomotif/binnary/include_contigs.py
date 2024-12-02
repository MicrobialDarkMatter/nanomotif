import polars as pl
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier
from nanomotif.binnary import data_processing as dp 
from nanomotif.binnary import utils as ut
import logging



def include_contigs(contig_methylation, contig_lengths):
    logger = logging.getLogger(__name__)
    logger.info("Starting include_contigs analysis...")
    
    contig_bin = contig_methylation\
        .select(["bin", "contig"])\
        .unique()
    
    binned_contig_methylation = dp.impute_contig_methylation_within_bin(
        contig_methylation
    )

    bins_w_methylation = binned_contig_methylation\
        .group_by(["bin", "motif_mod"])\
        .agg(
            pl.col("median").mean().alias("median")
        )\
        .with_columns(
            (pl.col("median") > 0.25).cast(pl.Int8).alias("is_methylated")
        )\
        .filter(
            pl.col("is_methylated") == 1
        )\
        .get_column("bin").unique()

    binned_contig_methylation = binned_contig_methylation\
        .filter(pl.col("bin").is_in(bins_w_methylation))

    binned_contig_names, binned_matrix = dp.create_matrix(binned_contig_methylation)
    bins = pl.DataFrame({
                            "contig": binned_contig_names,
                        })\
            .join(contig_bin, on= "contig", how = "left")


    unbinned_contig_methylation = dp.impute_unbinned_contigs(contig_methylation)
    unbinned_contig_names, unbinned_matrix = dp.create_matrix(unbinned_contig_methylation)

    lda = LinearDiscriminantAnalysis()
    knn = KNeighborsClassifier(n_neighbors = 3)
    rf = RandomForestClassifier(n_estimators=20, random_state=42)

    lda_pred = lda.fit(binned_matrix, bins.get_column("bin"))
    lda_labels = lda_pred.predict(unbinned_matrix)

    knn_pred = knn.fit(binned_matrix, bins.get_column("bin"))
    knn_labels = knn_pred.predict(unbinned_matrix)
    
    rf_pred = rf.fit(binned_matrix, bins.get_column("bin"))
    rf_labels = rf_pred.predict(unbinned_matrix)
    
    results = pl.DataFrame({
                    "contig": unbinned_contig_names,
                    "lda": lda_labels,
                    "knn": knn_labels,
                    "rf": rf_labels,
                })\
                .join(contig_bin, on="contig")\
                .melt(
                    id_vars = ["contig", "bin"],
                    value_vars=["lda", "knn", "rf"],
                    value_name="cluster",
                    variable_name="method"
                )

    n_assigned_bins = results\
        .group_by(["contig"])\
        .agg(
            pl.col("cluster").n_unique().alias("n_assigned_bins")
        )

    high_confidence_assignements = results\
        .filter(pl.col("contig").is_in(n_assigned_bins.filter(pl.col("n_assigned_bins") == 1).get_column("contig")))\
        .with_columns(
            pl.lit("high_confidence").alias("confidence"),
            pl.col("cluster").alias("assigned_bin")
        )\
        .sort(["contig"])



    medium_confidence = results\
        .filter(pl.col("contig").is_in(n_assigned_bins.filter(pl.col("n_assigned_bins") == 2).get_column("contig")))\
        .group_by(["contig", "cluster"])\
        .agg(
            pl.col("cluster").count().alias("n_assigned_bins"),
        )\
        .filter(pl.col("n_assigned_bins") == 2)\
        .with_columns(
            pl.lit("medium_confidence").alias("confidence"),
            pl.col("cluster").alias("assigned_bin")
        )\
        .select(["contig", "cluster", "confidence", "assigned_bin"])


    medium_confidence_assignements = results\
        .join(medium_confidence, on = ["contig", "cluster"], how = "left")\
        .filter(pl.col("confidence").is_not_null())


    assigned_contigs = pl.concat([high_confidence_assignements, medium_confidence_assignements])
    return assigned_contigs
