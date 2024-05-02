import polars as pl
import pandas as pd
from nanomotif.binnary import data_processing as dp
from nanomotif.binnary import scoring as sc
from nanomotif.binnary.utils import split_bin_contig
import logging

def detect_contamination(motifs_scored_in_bins, bin_consensus, args):
    """
    Takes the bin_motif_binary and motifs_scored_in_bins DataFrames and performs the contamination detection analysis.
    Firstly bin_motif_binary is used to create a binary representation of the methylation status of each motif in each bin.
    This is the bases for the comparison. Only motifs that has a mean methylation value of at least 0.75 in bin_motif_binary are considered.
    
    Motifs_scored_in_bins is then used to create a binary representation of the methylation status of each motif in each contig.
    This is then compared to the binary representation of the bin to identify any mismatches.
    Only motifs that are observed more than 6 times are considered and must have a mean methylation value of at least 0.75-0.15.
    
    
    params:
        motifs_scored_in_bins: pd.DataFrame - DataFrame containing the motifs scored in each bin
        motifs_of_interest: list - List of motifs to be considered from bin_motif_binary
        args: argparse.Namespace - Namespace containing the arguments passed to the script
    """
    logger = logging.getLogger(__name__)
    logger.info("Starting contamination detection analysis...")
    motifs_scored_in_bins_wo_unbinned = motifs_scored_in_bins \
        .filter(~pl.col("bin_contig").str.contains("unbinned"))
    

    contig_bin_comparison_score, contigs_w_no_methylation = sc.compare_methylation_pattern_multiprocessed(
        motifs_scored_in_bins=motifs_scored_in_bins_wo_unbinned,
        bin_consensus=bin_consensus,
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
            (pl.col("binary_methylation_missmatch_score") > 0)
        )
    
    
    # contamination_contigs = contig_bin_comparison_score[
    #     # NOTE: This line also removes all contigs from bins with no methylation
    #     (contig_bin_comparison_score["bin"] == contig_bin_comparison_score["contig_bin"]) &
    #     (contig_bin_comparison_score["binary_methylation_missmatch_score"] > 0)
    # ]

    # logger.info("Finding alternative bin for contamination contigs")
    # # Find alternative bin for contamination contigs
    # ## Must have a perfect match
    # contamination_contigs_alternative_bin = contig_bin_comparison_score \
    #     .filter(
    #         (pl.col("bin") != pl.col("contig_bin")) &
    #         (pl.col("binary_methylation_missmatch_score") == 0) &
    #         (~pl.col("bin_compare").is_in(contigs_w_no_methylation))
    #     ) \
    #     .select(["contig", "bin", "binary_methylation_missmatch_score"]) \
    #     .rename(
    #         {
    #             "bin": "alternative_bin",
    #             "binary_methylation_missmatch_score": "alternative_bin_binary_methylation_missmatch_score"
    #         }
    #     )
    

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
