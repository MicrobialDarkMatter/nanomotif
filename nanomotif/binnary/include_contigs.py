import polars as pl
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.decomposition import PCA
from nanomotif.binnary import data_processing as dp 
from nanomotif.binnary import utils as ut
import logging



def include_contigs(contig_methylation, contig_lengths, mean_probability):
    logger = logging.getLogger(__name__)
    logger.info("Starting include_contigs analysis...")
    
    contig_bin = contig_methylation\
        .select(["bin", "contig"])\
        .unique()
    
    binned_contig_methylation = dp.impute_contig_methylation_within_bin(
        contig_methylation
    )


    binned_contig_names, binned_matrix = dp.create_matrix(binned_contig_methylation)

    pca_variance=0.90
    pca = PCA(n_components=pca_variance, svd_solver = "full")
    binned_matrix = pca.fit_transform(binned_matrix)
    
    
    bins = pl.DataFrame({
                            "contig": binned_contig_names,
                        })\
            .join(contig_bin, on= "contig", how = "left")


    unbinned_contig_methylation = dp.impute_unbinned_contigs(contig_methylation)
    unbinned_contig_methylation_feature_completions = pl.DataFrame({
                                                           "contig": unbinned_contig_methylation.get_column("contig").unique()
                                                       })\
            .join(
                pl.DataFrame({"motif_mod": binned_contig_methylation.get_column("motif_mod").unique()}), how = "cross"
            )\
            .join(unbinned_contig_methylation, on = ["contig", "motif_mod"], how = "left")
            

    
    unbinned_contig_names, unbinned_matrix = dp.create_matrix(unbinned_contig_methylation_feature_completions)
    unbinned_matrix = pca.transform(unbinned_matrix)

    knn = KNeighborsClassifier(n_neighbors = 3)
    rf = RandomForestClassifier(n_estimators=100, random_state=42)
    lda = LinearDiscriminantAnalysis()

    lda_pred = lda.fit(binned_matrix, bins.get_column("bin"))
    lda_labels = lda_pred.predict(unbinned_matrix)
    lda_prob = lda_pred.predict_proba(unbinned_matrix).max(axis = 1)
    
    knn_pred = knn.fit(binned_matrix, bins.get_column("bin"))
    knn_labels = knn_pred.predict(unbinned_matrix)
    knn_prob = knn_pred.predict_proba(unbinned_matrix).max(axis = 1)
    
    rf_pred = rf.fit(binned_matrix, bins.get_column("bin"))
    rf_labels = rf_pred.predict(unbinned_matrix)
    rf_prob = rf_pred.predict_proba(unbinned_matrix).max(axis = 1)
    
    rf_df = pl.DataFrame({
        "contig": unbinned_contig_names,
        "method": "rf",
        "pred": rf_labels,
        "prob": rf_prob
    })

    knn_df = pl.DataFrame({
        "contig": unbinned_contig_names,
        "method": "knn",
        "pred": knn_labels,
        "prob": knn_prob
    })

    lda_df = pl.DataFrame({
        "contig": unbinned_contig_names,
        "method": "lda",
        "pred": lda_labels,
        "prob": lda_prob
    })

    prob_df = pl.concat([rf_df, knn_df, lda_df])
    
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
                    value_name="assigned_bin",
                    variable_name="method"
                )
    
    columns = ["contig", "bin", "assigned_bin", "method", "prob", "mean_prob", "confidence"]
    n_assigned_bins = results\
        .group_by(["contig"])\
        .agg(
            pl.col("assigned_bin").n_unique().alias("n_assigned_bins")
        )

    high_confidence_assignements = results\
    .filter(pl.col("contig").is_in(n_assigned_bins.filter(pl.col("n_assigned_bins") == 1).get_column("contig")))


    high_mean_prob = high_confidence_assignements\
        .join(prob_df, on = ["contig", "method"])\
        .group_by(["contig"])\
        .agg(
            pl.col("prob").mean().alias("mean_prob")
        )


    high_confidence_assignements = high_confidence_assignements\
        .join(prob_df, on = ["contig", "method"])\
        .join(high_mean_prob, on = ["contig"])\
        .with_columns(
            pl.when(pl.col("mean_prob") >= mean_probability).then(pl.lit("high_confidence")).otherwise(pl.lit("medium_confidence")).alias("confidence")
        )\
        .sort(["contig"])\
        .select(columns)

    low_confidence_assignments = results\
    .filter(pl.col("contig").is_in(n_assigned_bins.filter(pl.col("n_assigned_bins") == 2).get_column("contig")))

    
    low_confidence_classification = low_confidence_assignments\
        .group_by(["contig", "assigned_bin"])\
        .agg(
            pl.col("assigned_bin").count().alias("n_assigned_bins"),
        )\
        .filter(pl.col("n_assigned_bins") == 2)\
        .drop(["n_assigned_bins"])


    low_confidence_assignments = low_confidence_classification\
        .join(low_confidence_assignments, on = ["contig", "assigned_bin"])\
        .join(prob_df, on = ["contig", "method"])
        

    low_mean_prob = low_confidence_assignments\
        .group_by(["contig"])\
        .agg(
            pl.col("prob").mean().alias("mean_prob")
        )

    low_confidence_assignments = low_confidence_assignments\
        .join(low_mean_prob, on = ["contig"])\
        .with_columns(
            pl.lit("low_confidence").alias("confidence")
        )\
        .sort(["contig"])\
        .select(columns)
    assigned_contigs = pl.concat([high_confidence_assignements, low_confidence_assignments])\
        .sort(["confidence", "contig"])
    return assigned_contigs
