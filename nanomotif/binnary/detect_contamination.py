from sklearn.ensemble import IsolationForest
import polars as pl
from nanomotif.binnary import data_processing as dp
import numpy as np
import logging

def detect_contamination(contig_methylation, contig_lengths, min_motif_var=0.1, mad_multiplier=3, min_threshold=0.1):
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

        contig_names, motif_mods, matrix = dp.create_matrix(bin_methylation)
        matrix, motif_mods = dp.filter_motifs(matrix, motif_mods, min_var=min_motif_var)
        logger.debug(f"Bin {bin} has shape after filtering: {matrix.shape}")

        if matrix.shape[1] == 0:
            logger.debug(f"Bin {bin} does not have any variable motifs for contamination detection")
            continue
        
        # TODO: We need to scale the matrix by motif occurences. More frequent motifs are more important than others.
        motif_counts = bin_methylation\
                .filter(pl.col("motif_mod").is_in(motif_mods))\
                .group_by(pl.col("motif_mod"))\
                .agg(
                        pl.col("motif_occurences_total").sum().alias("count")
                        )
        motif_counts = pl.DataFrame({"motif_mod": motif_mods}).join(motif_counts, on="motif_mod", how="left").get_column("count").to_numpy()
        print(motif_counts)
        motif_weights = motif_counts / motif_counts.sum()
        print(motif_weights)
        print(motif_mods)
        scaled_matrix = matrix * np.sqrt(motif_weights)


        # Centroid weighted by length of contig
        contig_order = pl.DataFrame({"contig": contig_names})
        lengths = contig_order.join(contig_lengths, how="left", on = "contig").get_column("length").to_numpy()
        weights = lengths / lengths.sum()
        centroid = scaled_matrix.T @ weights
        distances_to_centroid = np.linalg.norm(scaled_matrix - centroid, axis=1)
        print("Centroid: ", centroid)

        print(contig_names)
        print(lengths)
        print(distances_to_centroid)

        # BUG: Why are we not just basing it on the longest contig's distance?

        # sorted_indices = np.argsort(distances_to_centroid)
        # sorted_distances = distances_to_centroid[sorted_indices]
        # sorted_weights = weights[sorted_indices]
        # cumsum_weights = np.cumsum(sorted_weights)
        # median_idx = np.searchsorted(cumsum_weights, 0.5)
        longest_contig_idx = np.argmax(lengths)
        weighted_median = distances_to_centroid[longest_contig_idx]
        print("Weighted median: ", weighted_median)

        # Weighted MAD
        abs_deviations = np.abs(distances_to_centroid - weighted_median)
        sorted_dev_indices = np.argsort(abs_deviations)
        sorted_devs = abs_deviations[sorted_dev_indices]
        sorted_dev_weights = weights[sorted_dev_indices]
        cumsum_dev_weights = np.cumsum(sorted_dev_weights)
        mad_idx = np.searchsorted(cumsum_dev_weights, 0.5)
        weighted_MAD = sorted_devs[mad_idx]
        print("MAD: ", weighted_MAD)

        threshold = max(weighted_median + mad_multiplier * weighted_MAD, min_threshold)
        print(threshold)

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

    # results = [res for res in results if not res]
    print(results)
    results = pl.concat(results)
    
    return results

