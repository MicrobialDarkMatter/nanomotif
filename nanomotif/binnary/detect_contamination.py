from sklearn.ensemble import IsolationForest
import polars as pl
from nanomotif.binnary import data_processing as dp
import numpy as np
import logging

def detect_contamination(contig_methylation, contig_lengths, min_motif_var=0.1, mad_multiplier=3, min_threshold=0.5):
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
        matrix = dp.filter_motifs(matrix, min_var=min_motif_var)
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

        threshold = max(weighted_median + mad_multiplier * weighted_MAD, min_threshold)

        contam_fraction = weights[distances_to_centroid > threshold].sum()
        # If > 30% is contamination we set a warning
        warning = ""
        if contam_fraction > 0.3:
            logger.warning(f"Bin {bin}: Large contamination fraction ({contam_fraction:.1%}). "
                           f"Consider bin splitting.")
            warning = f"contam_fraction > {contam_fraction:.1%}, consider splitting"

        bin_result = pl.DataFrame({
            "bin": bin,
            "contig": contig_names,
            "centroid_dist": distances_to_centroid,
            "weight": weights,
            "warning": warning,
            }).join(contig_lengths, on = "contig", how="left")\
            .with_columns(
                pl.when(pl.col("centroid_dist") > threshold).then(-1).otherwise(0).alias("prediction")
            )

        results.append(bin_result)

    results = pl.concat(results)
    
    return results

