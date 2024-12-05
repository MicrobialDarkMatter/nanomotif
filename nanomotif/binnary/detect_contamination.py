import hdbscan
from sklearn.cluster import AgglomerativeClustering, SpectralClustering
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler
import polars as pl
from nanomotif.binnary import data_processing as dp
import logging

def detect_contamination(contig_methylation, contig_lengths, num_consensus = 4):
    logger = logging.getLogger(__name__)
    logger.info("Starting contamination detection analysis...")

    # Create the bin_consensus dataframe for scoring
    logger.info("Creating bin_consensus dataframe for scoring...")
    contig_methylation = dp.impute_contig_methylation_within_bin(
        contig_methylation
    )
    # INFO: It does not makes sense to find contamination in bins with a single contig
    single_contig_bins = contig_methylation\
        .select(["bin", "contig"])\
        .unique()\
        .group_by("bin")\
        .agg(pl.col("contig").count().alias("n_contigs"))\
        .filter(pl.col("n_contigs") == 1)\
        .get_column("bin")

    contig_methylation = contig_methylation\
        .filter(~pl.col("bin").is_in(single_contig_bins))

    contig_names, matrix = dp.create_matrix(contig_methylation)

    scaler = StandardScaler()
    matrix = scaler.fit_transform(matrix)

    
    n_bins = len(contig_methylation.get_column("bin").unique())

    spectral = SpectralClustering(n_clusters=n_bins, affinity = 'nearest_neighbors', random_state=42)
    spectral_labels = spectral.fit_predict(matrix)

    agg = AgglomerativeClustering(n_clusters = n_bins)
    agg_labels = agg.fit(matrix).labels_
    

    dscan = hdbscan.HDBSCAN(min_samples=2, min_cluster_size=2, metric = "euclidean", allow_single_cluster=False)
    dscan_labels = dscan.fit_predict(matrix)

    gmm = GaussianMixture(n_components=n_bins, covariance_type='full', random_state=42)
    gmm.fit(matrix)
    gmm_labels = gmm.predict(matrix)


    contig_bin = contig_methylation\
        .select(["bin", "contig"])\
        .unique() 

    results = pl.DataFrame({
                    "contig": contig_names,
                    "spectral": spectral_labels,
                    "agg": agg_labels,
                    "hdbscan": dscan_labels,
                    "gmm": gmm_labels
                })\
                .join(contig_bin, on="contig")\
                .melt(
                    id_vars = ["contig", "bin"],
                    value_vars=["spectral", "agg", "hdbscan", "gmm"],
                    value_name="cluster",
                    variable_name="method"
                )




    bin_size = contig_bin\
        .join(contig_lengths, on = "contig")\
        .group_by("bin")\
        .agg(
            pl.col("length").sum().alias("bin_length"),
            pl.col("contig").count().alias("n_contigs_bin")
        )

    cluster_sizes = results\
        .join(contig_lengths, on="contig")\
        .group_by(["bin","method", "cluster"])\
        .agg(
            pl.col("contig").count().alias("n_contigs"),
            pl.col("length").sum().alias("cluster_length")
        )\
        .join(bin_size, on = "bin")\
        .with_columns(
            (pl.col("n_contigs") / pl.col("n_contigs_bin")).alias("fraction_contigs"),
            (pl.col("cluster_length") / pl.col("bin_length")).alias("fraction_length")
        )

    assigned_cluster = cluster_sizes\
        .group_by(["method","bin"])\
        .agg(
            pl.col("cluster_length").max().alias("cluster_length")
        )\
        .join(cluster_sizes, on = ["bin", "method","cluster_length"], how = "left")\
        .rename({
                    "cluster": "bin_cluster"
                })\
        .drop(["cluster_length", "n_contigs"])\
        .filter(pl.col("fraction_length") >= 0.85)



    results = results\
        .join(assigned_cluster, on = ["bin", "method"])\
        .sort(["bin", "contig"])


    confident_contaminants = results\
        .with_columns(
           (pl.col("bin_cluster") != pl.col("cluster")).cast(pl.Int8).alias("is_contaminant")
        )\
        .group_by(["contig"])\
        .agg(
            pl.col("is_contaminant").sum().alias("sum_predictions"),
        )\
        .filter(pl.col("sum_predictions") >= num_consensus)\
        .get_column("contig")

    contamination_contigs = results\
        .filter(pl.col("contig").is_in(confident_contaminants))
    
    return contamination_contigs
