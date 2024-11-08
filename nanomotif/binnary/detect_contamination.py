import polars as pl
import pandas as pd
from nanomotif.binnary import data_processing as dp
from nanomotif.binnary import scoring as sc
from nanomotif.binnary.utils import split_bin_contig
import logging

def detect_contamination(contig_methylation, bin_methylation, args):
    logger = logging.getLogger(__name__)
    logger.info("Starting contamination detection analysis...")
    contig_methylation_wo_unbinned = contig_methylation \
        .filter(~pl.col("bin_contig").str.contains("unbinned"))
    
    
    contig_bin_comparison_score, contigs_w_no_methylation = sc.compare_methylation_pattern_multiprocessed(
        contig_methylation=contig_methylation_wo_unbinned,
        bin_methylation=bin_methylation,
        mode="contamination",
        args=args,
        num_processes=args.threads
    )
    
    logger.info("Finding contamination in bins")
    contig_bin_comparison_score = split_bin_contig(contig_bin_comparison_score)
    
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
