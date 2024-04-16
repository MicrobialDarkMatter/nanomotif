import pandas as pd
import polars as pl
import numpy as np
from nanomotif.binnary import data_processing as dp 
from nanomotif.binnary import scoring as sc
from nanomotif.binnary import utils as ut
import logging



def include_contigs(motifs_scored_in_bins, bin_consensus, contamination, args):
    """
    Takes the motifs_scored_in_bins and motifs_of_interest DataFrames and finds unbinned contigs with an exact methylation pattern as the bin.
    
    params:
        motifs_scored_in_bins: pd.DataFrame - DataFrame containing the motifs scored in each bin
        motifs_of_interest: list - List of motifs to be considered from bin_motif_binary
        args: argparse.Namespace - Namespace containing the arguments passed to the script
    
    """
    logger = logging.getLogger(__name__)
    logger.info("Starting include_contigs analysis...")
    
    
    # Remove bins with no methylation in consensus
    bins_w_no_methylation = bin_consensus \
        .group_by("bin") \
        .agg(
            pl.sum("methylation_binary").alias("binary_sum")    
        ) \
        .filter(pl.col("binary_sum") == 0) \
        .select("bin") \
        .unique()["bin"]
    

    bin_consensus = bin_consensus \
        .filter(~pl.col("bin").is_in(bins_w_no_methylation))
    
    # Retain only unbinned contigs or contigs in the contamination file
    contigs_for_comparison = motifs_scored_in_bins \
        .filter(
            (pl.col("bin_contig").str.contains("unbinned")) |
            (pl.col("bin_contig").is_in(contamination["bin_contig_compare"]))
        )
    

    contig_bin_comparison_score, contigs_w_no_methylation = sc.compare_methylation_pattern_multiprocessed(
        motifs_scored_in_bins=contigs_for_comparison,
        bin_consensus=bin_consensus,
        mode="include",
        args=args,
        num_processes=args.threads
    )
    
    contig_bin_comparison_score = ut.split_bin_contig(contig_bin_comparison_score)
    
    logger.info("Assigning contigs to bins...")
    # Filter contigs where motif comparisons are less than args.min_motif_comparisons
    # TODO: looking for 0 comparisons is now redundant. Also remove the column.
    contigs_of_interest = contig_bin_comparison_score \
        .filter(
            pl.col("non_na_comparisons") >= args.min_motif_comparisons,
            (~pl.col("bin_compare").is_in(contigs_w_no_methylation)),  # Remove contigs with no methylation
            (pl.col("binary_methylation_missmatch_score") == 0)        # Retain contigs with no methylation missmatch 
            
        ) \
        .sort("bin","bin_compare")
    
    
    single_contigs = contigs_of_interest \
        .group_by("contig") \
        .agg(
            pl.count("contig").alias("contig_count")
        ) \
        .filter(pl.col("contig_count") == 1)
    
    contigs_of_interest = contigs_of_interest \
        .filter(pl.col("contig").is_in(single_contigs["contig"]))
      

    logger.info("Finished include_contigs analysis.")
    return contigs_of_interest