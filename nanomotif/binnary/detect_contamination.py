# import hdbscan
from sklearn.ensemble import IsolationForest
import polars as pl
from nanomotif.binnary import data_processing as dp
import numpy as np
import logging

def detect_contamination(contig_methylation, contig_lengths, min_bin_size=3):
    logger = logging.getLogger(__name__)
    logger.info("Starting contamination detection analysis...")

    # Create the bin_consensus dataframe for scoring
    logger.info("Creating bin_consensus dataframe for scoring...")
    contig_methylation = dp.impute_contig_methylation_within_bin(
        contig_methylation
    )
    # INFO: It does not makes sense to find contamination in bins with a single contig
    small_bins = contig_methylation\
        .select(["bin", "contig"])\
        .unique()\
        .group_by("bin")\
        .agg(pl.col("contig").count().alias("n_contigs"))\
        .filter(pl.col("n_contigs") < 3)\
        .get_column("bin").to_list()

    # INFO: Remove bins with few contigs as contamination detection may become unstable
    contig_methylation = contig_methylation\
        .filter(~pl.col("bin").is_in(small_bins))
    
    multi_contig_bins = contig_methylation.get_column("bin").unique()

    results = []
    for bin in multi_contig_bins:
        bin_methylation = contig_methylation.filter(pl.col("bin") == bin)

        contig_names, matrix = dp.create_matrix(bin_methylation)
        matrix = dp.filter_motifs(matrix, min_var=0.1)
        logger.debug(f"Bin {bin} has shape after filtering: {matrix.shape}")
        from scipy.spatial.distance import pdist, squareform

        if matrix.shape[1] == 0:
            logger.debug(f"Bin {bin} does not have any variable motifs for contamination detection")
            continue
        
        # Centroid weighted by length of contig
        contig_order = pl.DataFrame({"contig": contig_names})
        lengths = contig_order.join(contig_lengths, how="left", on = "contig").get_column("length").to_numpy()
        weights = lengths / lengths.sum()
        centroid = (matrix.T @ weights).T
        distances_to_centroid = np.linalg.norm(matrix - centroid, axis=1)

        # median_dist = np.median(distances_to_centroid)
        # MAD = np.median(abs(distances_to_centroid - median_dist))
        # threshold = median_dist + 3 * MAD

        sorted_indices = np.argsort(distances_to_centroid)
        sorted_distances = distances_to_centroid[sorted_indices]
        sorted_weights = weights[sorted_indices]
        cumsum_weights = np.cumsum(sorted_weights)
        median_idx = np.searchsorted(cumsum_weights, 0.5)
        weighted_median = sorted_distances[median_idx]

        # Weighted MAD
        abs_deviations = np.abs(distances_to_centroid - weighted_median)
        sorted_dev_indices = np.argsort(abs_deviations)
        sorted_devs = abs_deviations[sorted_dev_indices]
        sorted_dev_weights = weights[sorted_dev_indices]
        cumsum_dev_weights = np.cumsum(sorted_dev_weights)
        mad_idx = np.searchsorted(cumsum_dev_weights, 0.5)
        weighted_MAD = sorted_devs[mad_idx]

        threshold = max(weighted_median + 3 * weighted_MAD, 0.5)


        # clf = IsolationForest(contamination = "auto", random_state=42)
        # predictions = clf.fit_predict(matrix)
        # scores = clf.score_samples(matrix)

        # clusterer = hdbscan.HDBSCAN(min_samples=2, min_cluster_size=2, metric = "euclidean", allow_single_cluster=True)
        # clusterer.fit(matrix)
        # labels = clusterer.labels_
        # probabilities = clusterer.probabilities_
        # persistence_dict = {i: p for i, p in enumerate(clusterer.cluster_persistence_)}

        bin_result = pl.DataFrame({
            "bin": bin,
            "contig": contig_names,
            "centroid_dist": distances_to_centroid,
            "weight": weights,
            }).join(contig_lengths, on = "contig", how="left")\
            .with_columns(
                pl.when(pl.col("centroid_dist") > threshold).then(-1).otherwise(0).alias("prediction")
                # pl.col("label").map_elements(lambda x: persistence_dict.get(x, None), return_dtype=pl.Float64).alias("cluster_persistence")
            )

        results.append(bin_result)

    results = pl.concat(results)
    
    contamination = results#.filter(pl.col("prediction") == -1)


    return contamination





    

    # 
    # original_features = matrix.shape[1]
    # print(f"Original number of features: {original_features}")

    # print("Applying PCA for feature decorrelation")
    # pca_variance=0.90
    # pca = PCA(n_components=pca_variance, svd_solver = "full")
    # matrix = pca.fit_transform(matrix)

    # reduced_components = matrix.shape[1]
    # print(f"PCA reduced the feature space from {original_features} to {reduced_components}")
    # print(f"Total explained variance by pca: {pca.explained_variance_ratio_.sum():.2f}")
    
    # n_bins = len(contig_methylation.get_column("bin").unique())

    # print("Running Spectral")
    # spectral = SpectralClustering(n_clusters=n_bins, affinity = 'nearest_neighbors', random_state=42, n_jobs = threads, n_neighbors = spectral_n_neighbors)
    # spectral_labels = spectral.fit_predict(matrix)

    # print("Running Agglomerative Clustering")
    # agg = AgglomerativeClustering(n_clusters = n_bins)
    # agg_labels = agg.fit(matrix).labels_
    
    # print("Running HDBSCAN")
    # dscan = hdbscan.HDBSCAN(min_samples=3, min_cluster_size=2, metric = "euclidean", allow_single_cluster=False)
    # dscan_labels = dscan.fit_predict(matrix)

    # print("Running Gaussian Mixture Model")
    # gmm = GaussianMixture(n_components=n_bins, covariance_type='full', random_state=42)
    # gmm.fit(matrix)
    # gmm_labels = gmm.predict(matrix)


    # contig_bin = contig_methylation\
    #     .select(["bin", "contig"])\
    #     .unique() 

    # results = pl.DataFrame({
    #                 "contig": contig_names,
    #                 "spectral": spectral_labels,
    #                 "agg": agg_labels,
    #                 "hdbscan": dscan_labels,
    #                 "gmm": gmm_labels
    #             })\
    #             .join(contig_bin, on="contig")\
    #             .unpivot(
    #                 index = ["contig", "bin"],
    #                 on=["spectral", "agg", "hdbscan", "gmm"],
    #                 value_name="cluster",
    #                 variable_name="method"
    #             )


    # bin_size = contig_bin\
    #     .join(contig_lengths, on = "contig")\
    #     .group_by("bin")\
    #     .agg(
    #         pl.col("length").sum().alias("bin_length"),
    #         pl.col("contig").count().alias("n_contigs_bin")
    #     )

    # cluster_sizes = results\
    #     .join(contig_lengths, on="contig")\
    #     .group_by(["bin","method", "cluster"])\
    #     .agg(
    #         pl.col("contig").count().alias("n_contigs"),
    #         pl.col("length").sum().alias("cluster_length")
    #     )\
    #     .join(bin_size, on = "bin")\
    #     .with_columns(
    #         (pl.col("n_contigs") / pl.col("n_contigs_bin")).alias("fraction_contigs"),
    #         (pl.col("cluster_length") / pl.col("bin_length")).alias("fraction_length")
    #     )

    # assigned_cluster = cluster_sizes\
    #     .group_by(["method","bin"])\
    #     .agg(
    #         pl.col("cluster_length").max().alias("cluster_length")
    #     )\
    #     .join(cluster_sizes, on = ["bin", "method","cluster_length"], how = "left")\
    #     .rename({
    #                 "cluster": "bin_cluster"
    #             })\
    #     .drop(["cluster_length", "n_contigs"])\
    #     .filter(pl.col("fraction_length") >= 0.85)



    # results = results\
    #     .join(assigned_cluster, on = ["bin", "method"])\
    #     .sort(["bin", "contig"])


    # confident_contaminants = results\
    #     .with_columns(
    #        (pl.col("bin_cluster") != pl.col("cluster")).cast(pl.Int8).alias("is_contaminant")
    #     )\
    #     .group_by(["contig"])\
    #     .agg(
    #         pl.col("is_contaminant").sum().alias("sum_predictions"),
    #     )\
    #     .filter(pl.col("sum_predictions") >= num_consensus)\
    #     .get_column("contig")

    # contamination_contigs = results\
    #     .filter(pl.col("contig").is_in(confident_contaminants))
    
