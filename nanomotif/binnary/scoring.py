from nanomotif.binnary import data_processing as dp
from nanomotif.binnary import utils as ut
import os
import polars as pl
from multiprocessing import Pool, Queue, Process, get_context
import logging
from logging.handlers import QueueHandler, QueueListener

# Set up logging
## Create queue for logging
# log_queue = Queue()

# ## Set up a listener to handle logs from the queue
# def setup_logging_queue(queue):
#     while True:
#         record = queue.get()
#         if record is None:  # Use None as a sentinel to stop the listener
#             break
#         logger = logging.getLogger(record.name)
#         logger.handle(record)

# def worker_setup_logging(queue):
#     q_handler = QueueHandler(queue)
#     logger = logging.getLogger()
#     logger.setLevel(logging.INFO)
#     logger.addHandler(q_handler)

def define_mean_methylation_thresholds(motif_binary_compare):
    """
    Define mean methylation thresholds for bin and contig motifs
    """
    # Calculate the mean methylation value for each motif in each bin
    motif_binary_compare = motif_binary_compare.with_columns([
        pl.when(pl.col("methylation_binary") == 1)
        .then(
            pl.when(pl.col("mean_methylation") - 4 * pl.col("std_methylation_bin") > 0.1)
            .then(pl.col("mean_methylation") - 4 * pl.col("std_methylation_bin"))
            .otherwise(0.1)
        )
        .otherwise(pl.lit(None))
        .alias("methylation_mean_threshold")
    ])

    # Calculate the binary methylation value for each motif in each bin where the bin consensus is 1
    motif_binary_compare = motif_binary_compare.with_columns([
        pl.when((pl.col("methylation_binary") == 1) & 
                ((pl.col("mean") >= pl.col("methylation_mean_threshold")) | 
                (pl.col("mean") > 0.4)))
        .then(1)
        .when(pl.col("methylation_binary") == 1)
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
        .then((pl.col("mean") >= 0.25).cast(pl.Int32))
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
        .then(0)
        .when((motif_binary_compare['methylation_binary'] == 0) & (motif_binary_compare['methylation_binary_compare'].is_null()))
        .then(0)
        .otherwise(pl.lit(None))
    )

    # Add the 'motif_comparison_score' column to the DataFrame
    motif_binary_compare = motif_binary_compare.with_columns(motif_comparison_score.alias("motif_comparison_score"))
    
    # Group by bin and bin_compare and calculate the sum of the motif_comparison_score and the count of non-NA values
    contig_bin_comparison_score = motif_binary_compare \
        .group_by(["bin", "bin_compare"]) \
        .agg([
            pl.sum("motif_comparison_score").alias("binary_methylation_missmatch_score"),
            pl.count("motif_comparison_score").alias("non_na_comparisons")
        ])
        
    return contig_bin_comparison_score


def process_bin_contig(
    bin_contig, 
    bin_motifs_from_motifs_scored_in_bins, 
    motifs_scored_in_contigs, 
    mode, 
    args
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
    # worker_setup_logging(log_queue)
    # logger = logging.getLogger(__name__)
    # logger.info(f"Processing {bin_contig}")
    try:
        bin = ut.get_bin(bin_contig)
        
        if mode == "contamination":
            motif_binary_compare = bin_motifs_from_motifs_scored_in_bins \
                .filter(pl.col("bin") == bin) \
                .join(
                    motifs_scored_in_contigs.filter(pl.col("bin_compare") == bin_contig),
                    on="motif_mod"
                )
        if mode == "include":
            motif_binary_compare = bin_motifs_from_motifs_scored_in_bins \
                .filter(pl.col("bin") != bin) \
                .join(
                    motifs_scored_in_contigs.filter(pl.col("bin_compare") == bin_contig),
                    on="motif_mod"
                )
        
        # Define methylation thresholds
        motif_binary_compare = define_mean_methylation_thresholds(motif_binary_compare)

        if args.save_scores:
            path = os.path.join(args.out, "scores", args.command, "binary_compare" + bin_contig + ".csv")
            motif_binary_compare.write_csv(path)

        # Calculate the comparison score regardless of methylation presence
        contig_bin_comparison_score = compare_methylation_pattern(motif_binary_compare)
        
        # Check if the contig has no methylation and note it, but do not exclude it from further processing
        contigHasNMethylation = motif_binary_compare.filter(pl.col("methylation_binary_compare") == 1).height
        # logger.info(f"Finished processing {bin_contig}. Contig has {contigHasNMethylation} positive methylation comparisons.")
        
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
                
            
        if contigHasNMethylation == 0:
            return contig_bin_comparison_score, bin_contig
        
        return contig_bin_comparison_score, None
    except Exception as e:
        error_log_path = os.path.join(args.out, "logs", f"{args.command}_{bin_contig}_error.err")
        
        with open(error_log_path, 'w') as error_file:
            error_file.write(f"Error processing {bin_contig}: {str(e)}\n")
        return None, None

def compare_methylation_pattern_multiprocessed(motifs_scored_in_bins, bin_consensus, mode, args, num_processes=1):
    logger = logging.getLogger(__name__)
    logger.info("Starting comparison of methylation patterns")
    
    motifs_scored_in_contigs = motifs_scored_in_bins \
        .filter(pl.col("n_motifs") >= args.n_motif_contig_cutoff) \
        .select(["bin_contig", "motif_mod", "mean"]) \
        .rename({"bin_contig": "bin_compare"})
    
    if args.save_scores:
        dir = os.path.join(args.out, "scores", args.command)
        if not os.path.exists(dir):
            os.makedirs(dir)

    comparison_score = pl.DataFrame()
    contigs_w_no_methylation = []

    with get_context("spawn").Pool(processes=num_processes) as pool:
        results = pool.starmap(
            process_bin_contig,
            [
                (bin_contig, bin_consensus, motifs_scored_in_contigs, mode, args)
                for bin_contig in motifs_scored_in_contigs.select("bin_compare").unique()["bin_compare"]
            ]
        )
        
    for result, no_methylation in results:
        if result is not None:
            comparison_score = pl.concat([comparison_score, result])
        if no_methylation is not None:
            contigs_w_no_methylation.append(no_methylation)
        
    return comparison_score, contigs_w_no_methylation

