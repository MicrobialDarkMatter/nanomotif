from nanomotif.binnary import data_processing as dp
from nanomotif.binnary import utils as ut
from nanomotif.parallel import update_progress_bar
import os
import sys
import polars as pl
import multiprocessing
from multiprocessing import get_context
import logging
from logging.handlers import QueueHandler, QueueListener


def define_mean_methylation_thresholds(motif_binary_compare):
    """
    Define mean methylation thresholds for bin and contig motifs
    """
    # Calculate the mean methylation value for each motif in each bin
    motif_binary_compare = motif_binary_compare.with_columns([
        pl.when(pl.col("methylation_binary") == 1)
        .then(
            pl.when(pl.col("mean_bin_median") - 4 * pl.col("std_bin_median_filtered") > 0.25)
            .then(pl.col("mean_bin_median") - 4 * pl.col("std_bin_median_filtered"))
            .otherwise(0.25)
        )
        .otherwise(pl.lit(None))
        .alias("methylation_mean_threshold")
    ])

    # Calculate the binary methylation value for each motif in each bin where the bin consensus is 1
    motif_binary_compare = motif_binary_compare.with_columns([
        pl.when((pl.col("methylation_binary") == 1) & 
                ((pl.col("median") >= pl.col("methylation_mean_threshold")) | 
                (pl.col("median") > 0.4)))
        .then(1)
        .when((pl.col("methylation_binary") == 1) & pl.col("median").is_not_null())
        .then(0)
        .otherwise(pl.lit(None))
        .alias("methylation_binary_compare")
    ])

    # Calculate score for bin consensus is 0
    motif_binary_compare = motif_binary_compare.with_columns([
        pl.when(pl.col("methylation_binary") == 0)
        .then(0.25)
        .otherwise(pl.col("methylation_mean_threshold"))
        .alias("methylation_mean_threshold"),

        pl.when(pl.col("methylation_binary") == 0)
        .then((pl.col("median") >= 0.25).cast(pl.Int32))
        .otherwise(pl.col("methylation_binary_compare"))
        .alias("methylation_binary_compare")
    ])
    
    return motif_binary_compare



def compare_methylation_pattern(motif_binary_compare):
    """
    Compares the methylation pattern between bin and contig motifs using Polars and calculates the motif_comparison_score.
    """
    motif_comparison_score = (
        pl.when((motif_binary_compare['methylation_binary'] == 1) & (motif_binary_compare['methylation_binary_compare'] == 1))
        .then(0)
        .when((motif_binary_compare['methylation_binary'] == 1) & (motif_binary_compare['methylation_binary_compare'] == 0))
        .then(1)
        .when((motif_binary_compare['methylation_binary'] == 0) & (motif_binary_compare['methylation_binary_compare'] == 1))
        .then(1)
        .when((motif_binary_compare['methylation_binary'] == 0) & (motif_binary_compare['methylation_binary_compare'] == 0))
        .then(0)
        .when((motif_binary_compare['methylation_binary'] == 1) & (motif_binary_compare['methylation_binary_compare'].is_null()))
        .then(1)
        .when((motif_binary_compare['methylation_binary'] == 0) & (motif_binary_compare['methylation_binary_compare'].is_null()))
        .then(0)
        .otherwise(pl.lit(None))
    )

    # Add the 'motif_comparison_score' column to the DataFrame
    motif_binary_compare = motif_binary_compare.with_columns(motif_comparison_score.alias("motif_comparison_score"))
    
    # Group by bin and bin_compare and calculate the sum of the motif_comparison_score and the count of non-NA values
    contig_bin_comparison_score = motif_binary_compare \
        .group_by(["bin", "bin_compare", "contig_bin"]) \
        .agg([
            pl.sum("motif_comparison_score").alias("binary_methylation_missmatch_score"),
            pl.count("motif_comparison_score").alias("non_na_comparisons")
        ])
        
    return contig_bin_comparison_score


def process_bin_contig(
    task, 
    bin_motifs_from_motifs_scored_in_bins, 
    motifs_scored_in_contigs, 
    mode, 
    args,
    counter,
    lock
):
    """
    This function processes a single contig and compares the methylation pattern between the contig and the bin.
    Depending on the mode, the function will either look for contamination or include contigs.
    
    Parameters:
    - bin_contig: str
        The contig to process
    - bin_motifs_from_motifs_scored_in_bins: polars.DataFrame
        The motifs scored in the bin
    - motifs_scored_in_contigs: polars.DataFrame
        The motifs scored in the contigs
    - mode: str
        The mode to run the comparison in. Either "contamination" or "include"
    - args: argparse.Namespace
        The arguments passed to the script
    
    Returns:
    - contig_bin_comparison_score: polars.DataFrame
        The comparison score for the contig
    """
    
    bin, bin_contig = task
    
    if mode == "contamination":
        motif_binary_compare = bin_motifs_from_motifs_scored_in_bins \
            .filter(pl.col("bin") == bin) \
            .join(
                ut.add_compare_df(motifs_scored_in_contigs, bin_contig),
                on="motif_mod"
            )
    if mode == "include":
        motif_binary_compare = bin_motifs_from_motifs_scored_in_bins \
            .filter(pl.col("bin") != bin) \
            .join(
                ut.add_compare_df(motifs_scored_in_contigs, bin_contig),
                on="motif_mod",
                how="left"
            )\
            .with_columns(
                pl.col("contig_bin").fill_null(bin),
                pl.col("bin_compare").fill_null(bin_contig)
            )

    # Define methylation thresholds
    motif_binary_compare = define_mean_methylation_thresholds(motif_binary_compare)

    if args.save_scores:
        path = os.path.join(args.out, "scores", mode, "binary_compare" + bin_contig + ".csv")
        motif_binary_compare.write_csv(path)

    # Calculate the comparison score regardless of methylation presence
    contig_bin_comparison_score = compare_methylation_pattern(motif_binary_compare)
    
    # Check if the contig has no methylation and note it, but do not exclude it from further processing
    contigHasNMethylation = motif_binary_compare.filter(pl.col("methylation_binary_compare") == 1).height
    
    
    if mode == "include":
        contig_bin_comparison_score = contig_bin_comparison_score \
            .filter(
                pl.col("binary_methylation_missmatch_score") == 0
            )
    if mode == "contamination":
        contig_bin_comparison_score = contig_bin_comparison_score \
            .filter(
                pl.col("binary_methylation_missmatch_score") > 0
            )
            
    with lock:
        counter.value += 1
        
    if contigHasNMethylation == 0:
        return contig_bin_comparison_score, bin_contig   
    return contig_bin_comparison_score, None

def create_dir_if_not_exists(path):
    if not os.path.exists(path):
        os.makedirs(path)

def compare_methylation_pattern_multiprocessed(motifs_scored_in_bins, bin_consensus, mode, args, num_processes=1):
    logger = logging.getLogger(__name__)
    logger.info("Starting comparison of methylation patterns")
    
    motifs_scored_in_contigs = motifs_scored_in_bins \
        .filter(pl.col("N_motif_obs") >= args.n_motif_contig_cutoff)
    
    unique_tasks = motifs_scored_in_contigs.select(["bin", "bin_contig"]).unique()
    
    tasks = [(row['bin'], row['bin_contig']) for row in unique_tasks.iter_rows(named = True)]
    
    motifs_scored_in_contigs = motifs_scored_in_contigs \
        .select(["bin", "bin_contig", "motif_mod", "median"]) \
        .rename({"bin_contig": "bin_compare"})
    
    # Create a progress manager
    manager = multiprocessing.Manager()
    counter = manager.Value('i', 0)
    lock = manager.Lock()

    # Create a pool of workers
    pool = get_context("spawn").Pool(processes=num_processes)

    # Create a process for the progress bar
    progress_bar_process = multiprocessing.Process(target=update_progress_bar, args=(counter, len(tasks), True))
    progress_bar_process.start()
    
    if args.save_scores:
        create_dir_if_not_exists(os.path.join(args.out, "scores"))
        create_dir_if_not_exists(os.path.join(args.out, "scores", "pre-scoring"))
        create_dir_if_not_exists(os.path.join(args.out, "scores", mode))

    # Put them workers to work
    results = pool.starmap(process_bin_contig, [(
        task, 
        bin_consensus, 
        motifs_scored_in_contigs, 
        mode, 
        args,
        counter,
        lock
        ) for task in tasks])
    results = [result for result in results if result is not None] #TODO: Check if this is necessary

    # Close the pool
    pool.close()
    pool.join()

    # Close the progress bar
    progress_bar_process.join()
    
    comparison_score = pl.DataFrame()
    contigs_w_no_methylation = []
    
    for result, no_methylation in results:
        if result is not None:
            comparison_score = pl.concat([comparison_score, result])
        if no_methylation is not None:
            contigs_w_no_methylation.append(no_methylation)
    
    if comparison_score.shape[0] == 0:
        logger.warning(f"{args.command}: No contigs were found with the specified criteria.")
    
    return comparison_score, contigs_w_no_methylation

