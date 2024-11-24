from sklearn.manifold import TSNE
import polars as pl
from nanomotif.binnary import data_processing as dp
import logging

def detect_contamination(contig_methylation, contig_lengths, args):
    logger = logging.getLogger(__name__)
    logger.info("Starting contamination detection analysis...")

    contig_names, matrix = dp.create_matrix(contig_methylation)

    logger.info("Starting tSNE")
    tsne = TSNE(n_components = 3, perplexity = 30, random_state=42, n_jobs=args.threads)
    tsne_embedding = tsne.fit_transform(matrix)

    contig_weight = contig_lengths\
        .filter(pl.col("contig").is_in(contig_names))\
        .with_columns(
            pl.col("length").log(base=10).alias("weight")
        )\
        .get_column("weight")

    from sklearn.cluster import DBSCAN
    logger.info("Performing DBSCAN")
    dbscan = DBSCAN(eps=0.5, min_samples=5, n_jobs=args.threads)
    dbscan.fit(tsne_embedding, sample_weight = contig_weight)
    labels = dbscan.labels_

    results = pl.DataFrame({
                               "contig": contig_names,
                               "cluster": labels
                           })

    results = contig_methylation\
        .select(["bin", "contig"])\
        .unique()\
        .join(results, on = "contig", how = "left")

    logger.info("Finding contamination in bins")
    
    # Filter contig_bin == bin and contig_bin_comparison_score > 0
    contamination_contigs = contig_bin_comparison_score \
        .filter(
            (pl.col("bin") == pl.col("contig_bin")) &
            (pl.col("binary_methylation_mismatch_score") > 0)
        )
    
    contamination_contigs = contamination_contigs \
        .drop("contig_bin") \
            .rename(
                {
                    "bin_compare": "bin_contig_compare"
                }
            ) \
        .sort("bin", "bin_contig_compare")

    
    logger.info("Contamination detection complete")
    
    return contamination_contigs
