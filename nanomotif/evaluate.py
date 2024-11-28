import numpy as np
import polars as pl
import os
from polars import col
import logging as log
import itertools
from itertools import takewhile
from scipy.stats import entropy
import multiprocessing
from multiprocessing import get_context
import progressbar
from nanomotif.constants import *
from nanomotif.model import BetaBernoulliModel
from nanomotif.utils import subseq_indices, calculate_match_length
from nanomotif.seq import EqualLengthDNASet, DNAsequence, DNAarray
from nanomotif.candidate import Motif, MotifTree
from nanomotif.logger import configure_logger
from nanomotif.parallel import update_progress_bar
from nanomotif.seed import set_seed
import heapq as hq
import networkx as nx
import warnings
import time
from typing import Generator, Optional, List, Tuple

##########################################
# Motif candidate scoring
##########################################

def methylated_motif_occourances(
        motif, 
        sequence, 
        methylated_positions
    ) -> tuple:
    """
    Get occourances of a motif in a contig

    Parameters:
    - motif (str): The motif to search for.
    - seq (str): The contig to search in.
    - contig_meth_positions (list): The positions of the methylation sites in the contig.

    Returns:
    - tuple: A tuple of two numpy arrays, the first containing the positions of the methylated motifs, the second containing the positions of the non-methylated motifs.
    """
    assert len(motif) > 0, "Motif is empty"
    assert len(sequence) > 0, "Sequence is empty"
    assert type(motif) == Motif, "Motif is not a Motif type"
    # Get the index of the methylation in the motif in the contig
    motif_index = subseq_indices(motif.string, sequence) + motif.mod_position

    # Methylated motif positions
    meth_occurences = np.intersect1d(methylated_positions, motif_index)

    # Non-methylated motif positions
    nonmeth_occurences =  np.setdiff1d(motif_index, methylated_positions)

    return meth_occurences, nonmeth_occurences

def motif_model_contig(
        pileup, 
        contig: str, 
        motif, 
        save_motif_positions=False,
        na_positions = None
    ):
    """
    Get the posterior for a single motif. Uses number of methylated motifs as methylation count.

    Parameters:
    - pileup (Pileup): The pileup to be processed.
    - contig (str): The contig to be processed.
    - motif (str): The motif to be processed.
    - motif_meth_index (int): The index of the methylation site in the motif.

    Returns:
    - BetaBernoulliModel: The model for the candidate.
    """
    # Strip motif of periods    # Strip motif of periods
    motif_stripped = motif.new_stripped_motif()

    meth_positions_fwd = pileup.filter(pl.col("strand") == "+")["position"].to_numpy()
    meth_positions_rev = pileup.filter(pl.col("strand") == "-")["position"].to_numpy()

    index_meth_fwd, index_nonmeth_fwd = methylated_motif_occourances(motif_stripped, contig, meth_positions_fwd)
    index_meth_rev, index_nonmeth_rev = methylated_motif_occourances(motif_stripped.reverse_compliment(), contig, meth_positions_rev)

    if na_positions is not None:
        na_positions_fwd = na_positions["fwd"]
        index_meth_fwd = np.setdiff1d(index_meth_fwd, na_positions_fwd)
        index_nonmeth_fwd = np.setdiff1d(index_nonmeth_fwd, na_positions_fwd)
        na_positions_rev = na_positions["rev"]
        index_meth_rev = np.setdiff1d(index_meth_rev, na_positions_rev)
        index_nonmeth_rev = np.setdiff1d(index_nonmeth_rev, na_positions_rev)

    model = BetaBernoulliModel()
    model.update(len(index_meth_fwd) + len(index_meth_rev), len(index_nonmeth_fwd) + len(index_nonmeth_rev))
    
    if save_motif_positions:
        motif_data = {
            'index_meth_fwd': index_meth_fwd,
            'index_nonmeth_fwd': index_nonmeth_fwd,
            'index_meth_rev': index_meth_rev,
            'index_nonmeth_rev': index_nonmeth_rev
        }
        return model, motif_data
    else:
        return model


def methylated_reads_counts(
        pileup: pl.DataFrame, 
        sequence: str, 
        motif: Motif
    ) -> tuple:
    """
    Get the number of methylated and non-methylated reads for a motif.

    Parameters:
    - pileup (Pileup): The pileup to be processed. 
    - sequence (str): The sequence to be processed.
    - motif (Motif): The motif to be processed.

    Returns:
    - tuple: A tuple of two integers, the first containing the number of methylated reads, the second containing the number of non-methylated reads.
    """
    motif_index = subseq_indices(motif.string, sequence) + motif.mod_position

    pileup = pileup.filter(pl.col("position").is_in(motif_index))
    counts = pileup \
        .select(["n_mod", "n_nonmod"]) \
        .sum()
    n_mod = counts.get_column("n_mod")[0]
    n_nonmod = counts.get_column("n_nonmod")[0]
    return n_mod, n_nonmod


def motif_model_read(pileup: pl.DataFrame, contig: str, motif: Motif) -> BetaBernoulliModel:
    """
    Get the posterior for a single motif. Uses number of methylated reads as methylation count.

    Parameters:
    - pileup (Pileup): The pileup to be processed.
    - read (str): The read to be processed.
    - motif (str): The motif to be processed.

    Returns:
    - BetaBernoulliModel: The model for the candidate.
    """
    motif_stripped = motif.new_stripped_motif()
    pileup_counts = pileup.with_columns((pl.col("fraction_mod")*pl.col("Nvalid_cov")).alias("n_mod")) \
        .with_columns(pl.col("n_mod").round().cast(pl.Int32)) \
        .with_columns((pl.col("Nvalid_cov") - pl.col("n_mod")).alias("n_nonmod"))
    
    pileup_plus = pileup_counts.filter(pl.col("strand") == "+")
    pileup_minus = pileup_counts.filter(pl.col("strand") == "-")

    n_meth_plus, n_nonmeth_plus = methylated_reads_counts(pileup_plus, contig, motif_stripped)
    n_meth_minus, n_nonmeth_minus = methylated_reads_counts(pileup_minus, contig, motif_stripped.reverse_compliment())

    model = BetaBernoulliModel()
    model.update(n_meth_plus + n_meth_minus, n_nonmeth_plus + n_nonmeth_minus)
    return model


##########################################
# Motif candidate state space search
##########################################

def worker_function(
        args
    ):
    """
    Process a single subpileup for one contig and one modtype

    Parameters:
    - args (tuple): The arguments to the function: contig, modtype, subpileup
    - counter (multiprocessing.Value): The progress counter
    - lock (multiprocessing.Lock): The lock for the progress counter
    - assembly (Assembly): The assembly to be processed.
    - min_kl_divergence (float): Early stopping criteria, if max KL-divergence falls below, stops building motif.
    - padding (int): The padding to use for the motif.
    """
    (contig, modtype, subpileup, low_coverage_positions, counter, lock, assembly, 
     min_kl_divergence, padding, minimum_methylation_fraction_confident, 
     read_level_methylation, log_dir, verbose, seed) = args
    set_seed(seed=seed)
    warnings.filterwarnings("ignore")
    process_id = os.getpid()
    if log_dir is not None:
        log_file = f"{log_dir}/find-motifs.{process_id}.log"
        configure_logger(log_file, verbose=verbose)
    log.info(f"Worker {process_id} started")

    if low_coverage_positions is not None:
        low_coverage_positions = low_coverage_positions.explode("position", "strand")
        low_coverage_positions = {
            "fwd": low_coverage_positions.filter(pl.col("strand") == "+")["position"].to_numpy(),
            "rev": low_coverage_positions.filter(pl.col("strand") == "-")["position"].to_numpy()
        }
    

    try:
        result = process_subpileup(
            contig, 
            modtype, 
            subpileup, 
            assembly, 
            min_kl_divergence, 
            padding, 
            minimum_methylation_fraction_confident, 
            read_level_methylation,
            na_positions = low_coverage_positions
        )
        with lock:
            counter.value += 1
            log.debug(f"Task for ({contig}, {modtype}) complete; progress: {counter.value}")
        return result
    except Exception as e:
        log.error(f"Error processing {contig}, {modtype}: {e}")
        with lock:
            counter.value += 1
        return None


def process_subpileup(
        contig, 
        modtype, 
        contig_pileup, 
        assembly, 
        min_kl_divergence, 
        padding, 
        minimum_methylation_fraction_confident,
        read_level_methylation,
        na_positions: dict = None
    ):
    """
    Process a single subpileup for one contig and one modtype

    Parameters:
    - args (tuple): The arguments to the function: contig, modtype, subpileup
    - counter (multiprocessing.Value): The progress counter
    - lock (multiprocessing.Lock): The lock for the progress counter
    - assembly (Assembly): The assembly to be processed.
    - min_kl_divergence (float): Early stopping criteria, if max KL-divergence falls below, stops building motif.
    - padding (int): The padding to use for the motif.
    """
    log.info(f"Processing {contig} {modtype}")
    assert contig_pileup is not None, "Subpileup is None"
    assert len(contig_pileup) > 0, "Subpileup is empty"
    assert assembly is not None, "Assembly is None"
    assert min_kl_divergence >= 0, "min_kl_divergence must be greater than 0"
    assert contig_pileup.get_column("mod_type").unique().to_list() == [modtype], "subpileup modtype does not match modtype"

    contig_sequence = assembly.assembly[contig]

    motif_graph, best_candidates = find_best_candidates(
        contig_pileup, 
        contig_sequence, 
        modtype,
        minimum_methylation_fraction_confident,
        padding,
        na_positions = na_positions,
        read_level_methylation = read_level_methylation,
        min_kl = min_kl_divergence,
        max_dead_ends = 25,
        max_rounds_since_new_best = 15
    )
    identified_motifs = nxgraph_to_dataframe(motif_graph) \
        .filter(col("sequence").is_in(best_candidates))
    
    if len(identified_motifs) == 0:
        log.info("No motifs found")
        return None
    else:
        identified_motifs = identified_motifs.with_columns(
            pl.lit(contig).alias("contig"),
            pl.lit(modtype).alias("mod_type")
        ).with_columns([
            pl.col("model").map_elements(lambda x: x._alpha).alias("alpha"),
            pl.col("model").map_elements(lambda x: x._beta).alias("beta")
        ]).drop("model")
        return identified_motifs

def set_polars_env():
    os.environ["POLARS_MAX_THREADS"] = "1"



def process_sample_parallel(
        assembly, pileup, 
        low_coverage_positions = None,
        read_level_methylation = False,
        threads = 2,
        search_frame_size = 40,
        threshold_methylation_confident = 0.8,
        threshold_methylation_general = 0.7,
        threshold_valid_coverage = 5,
        minimum_kl_divergence = 0.05,
        verbose = False,
        log_dir = None,
        seed = None
    ):
    """
    Process a sample
    
    Parameters:
    - assembly (Assembly): The assembly of all contigs.
    - pileup (Pileup): The pileup to be processed.
    - max_candidate_size (int): The maximum size of the candidate motifs.
    - min_read_methylation_fraction (float): The minimum fraction of reads that must be methylated for a position to be considered methylated.
    - min_kl_divergence (float): Early stopping criteria, if max KL-divergence falls below, stops building motif.
    - min_cdf_score (float): Minimum score of 1 - cdf(cdf_position) for a motif to be considered valid.
    - cdf_position (float): The position to evaluate the cdf at.
    - min_motif_frequency (int): Used to get minimum number of sequences to evaluate motif at.
    """
    assert pileup is not None, "Pileup is None"
    assert len(pileup) > 0, "Pileup is empty" 
    assert assembly is not None, "Assembly is None"
    assert search_frame_size > 0, "search_frame_size must be greater than 0"
    assert threshold_methylation_confident >= 0 and threshold_methylation_confident <= 1, "threshold_methylation_confident must be between 0 and 1"
    assert threshold_methylation_general >= 0 and threshold_methylation_general <= 1, "threshold_methylation_general must be between 0 and 1"
    assert minimum_kl_divergence >= 0, "mininum_kl_divergence must be greater than 0"

    # Infer padding size from candidate_size
    padding = search_frame_size // 2
    
    log.debug("Partitining pileups to dict")
    pileup_dict = pileup.partition_by(["contig", "mod_type"], as_dict=True)
    low_coverage_positions_dict = low_coverage_positions.partition_by(["contig", "mod_type"], as_dict=True) if low_coverage_positions is not None else None

    def task_generator(pileup_dict, low_coverage_positions_dict):
        for (contig, modtype) in pileup_dict:
            subpileup = pileup_dict[(contig, modtype)]
            
            if low_coverage_positions is not None:
                if (contig, modtype) in low_coverage_positions_dict.keys():
                    low_coverage_positions_filtered = low_coverage_positions_dict[(contig, modtype)]
                else:
                    low_coverage_positions_filtered = None
            else:
                low_coverage_positions_filtered = None

            yield (contig, modtype, subpileup, low_coverage_positions_filtered, counter, lock, assembly, 
                       minimum_kl_divergence, search_frame_size // 2, threshold_methylation_confident, 
                       read_level_methylation, log_dir, verbose, seed)
                
    # Create a progress manager
    manager = multiprocessing.Manager()
    counter = manager.Value('i', 0)
    lock = manager.Lock()

    # Create a pool of workers
    pool = get_context("spawn").Pool(processes=threads, initializer=set_polars_env)

    # Create a process for the progress bar
    log.debug("Counting number of tasks")
    total_tasks = len(pileup.select("contig", "mod_type").unique())

    log.debug("Starting progress bar")
    progress_bar_process = multiprocessing.Process(target=update_progress_bar, args=(counter, total_tasks, True))
    progress_bar_process.start()

    # Put them workers to work
    log.debug("Starting workers")
    chunksize = max(min(100, total_tasks // threads), 1)
    results = pool.imap(worker_function, task_generator(pileup_dict, low_coverage_positions_dict), chunksize=chunksize)
    results = [result for result in results if result is not None]
    log.debug("Joining results")

    log.debug("Closing pool")
    pool.close()
    pool.join()

    log.debug("Joining progress bar")
    progress_bar_process.join()

    log.debug("Creating final dataframe")
    if len(results) == 0:
        log.debug("No motifs found")
        return None
    motifs = pl.concat(results, rechunk=True, parallel=False)

    model_col = []
    for a, b in zip(motifs.get_column("alpha").to_list(), motifs.get_column("beta").to_list()):
        model_col.append(BetaBernoulliModel(a, b))


    motifs = motifs.with_columns([
            pl.Series(model_col).alias("model"),
            pl.col("sequence").map_elements(lambda motif: motif[count_periods_at_start(motif):len(motif)-count_periods_at_end(motif)]).alias("motif"),
            pl.col("sequence").map_elements(lambda motif: padding - count_periods_at_start(motif)).alias("mod_position")
        ]).drop(["alpha", "beta"])

    return motifs




#########################################################################
# Motif candidate state space 

def find_best_candidates(        
        contig_pileup: pl.DataFrame, 
        contig_sequence: DNAsequence, 
        mod_type: str,
        minimum_methylation_fraction_confident: float,
        padding: int,
        na_positions = None,
        read_level_methylation: bool = False,
        min_kl: float = 0.2, 
        max_dead_ends: int = 25, 
        max_rounds_since_new_best: int = 15,
        score_threshold: float = 0.2,
        remaining_sequences_threshold: float = 0.01
    ) -> tuple[MotifTree, list[Motif]]:
    """
    Find the best motif candidates in a sequence.

    Parameters:
    - contig_pileup (Pileup): The pileup to be processed.
    - contig_sequence (DNAsequence): The sequence to be processed.
    - mod_type (str): The modtype to be processed.
    - minimum_methylation_fraction_confident (float): The minimum fraction of reads that must be methylated for a position to be considered confidently methylated and used in search.
    - padding (int): The padding to use for the motif.
    - min_kl (float): The minimum KL-divergence for a motif to be considered valid.
    - max_dead_ends (int): The maximum number of low scoring candidates before stopping the search.
    - max_rounds_since_new_best (int): The maximum number of rounds since a new best candidate was found before stopping the search.
    - score_threshold (float): The minimum score for a candidate to be considered valid.
    - remaining_sequences_threshold (float): The minimum fraction of sequences remaining before stopping the search.
    """
    subpileup_confident = contig_pileup.filter(pl.col("fraction_mod") >= minimum_methylation_fraction_confident)

    # Extract the sequences for confidently methylated positions
    index_plus = subpileup_confident.filter(pl.col("strand") == "+").get_column("position").to_list()
    index_minus = subpileup_confident.filter(pl.col("strand") == "-").get_column("position").to_list()
    
    # Sample sequences in contig to get background for KL-divergence
    contig_sequences_sample = contig_sequence.sample_n_subsequences(
        padding * 2 + 1, 10000
    )
    contig_pssm = contig_sequences_sample.pssm()

    methylation_sequences_string = []
    if len(index_plus) >= 1:
        methylation_sequences_string += contig_sequence.sample_at_indices(index_plus, padding).sequences
    if len(index_minus) >= 1:
        methylation_sequences_string += contig_sequence.sample_at_indices(index_minus, padding).reverse_compliment().sequences
    if len(methylation_sequences_string) == 0:
        log.info("No methylation sequences found")
        return None
    methylation_sequences = EqualLengthDNASet(methylation_sequences_string).convert_to_DNAarray()

    total_sequences = methylation_sequences.shape[0]
    
    root_motif = Motif("." * padding + MOD_TYPE_TO_CANONICAL[mod_type] + "." * padding, padding)
    methylation_sequences_clone = methylation_sequences.copy()
    best_candidates = []
    continue_search = True
    dead_ends = 0
    motif_graph = None

    while continue_search:
        if dead_ends >= max_dead_ends:
            log.debug("Stopping search, too many low scoring candidates")
            break
        # Find the initial guess within the tree
        searcher = MotifSearcher(
            root_motif, 
            contig_sequence, 
            contig_pssm,
            contig_pileup, 
            methylation_sequences_clone, 
            padding,
            read_level_methylation = read_level_methylation,
            na_positions = na_positions,
            motif_graph = motif_graph,
            min_kl = min_kl,
            max_rounds_since_new_best = max_rounds_since_new_best
        )
        motif_graph, naive_guess = searcher.run()

        # If there is no naive guess, stop the search
        if naive_guess == root_motif:
            log.debug("No naive guess found, stopping search")
            break

        # Remove new candidate from methylation sequences
        seq_before = methylation_sequences.shape[0]
        methylation_sequences_clone = methylation_sequences_clone.filter_sequence_matches(naive_guess.one_hot(), keep_matches = False)
        if methylation_sequences_clone is None:
            log.debug("No more sequences left")
            break

        # Check if we should continue the search
        seq_remaining = methylation_sequences_clone.shape[0]
        seq_remaining_percent = seq_remaining/total_sequences

        if motif_graph.nodes[naive_guess]["score"] < score_threshold:
            dead_ends += 1
            log.debug(f"Candidate has low score, {naive_guess}. {dead_ends} of {max_dead_ends} before temination")
            continue

        log.info(f"Keeping {naive_guess}, represented in {seq_before-seq_remaining} seqs. model: {motif_graph.nodes[naive_guess]['model']}. ({100*seq_remaining_percent:.1f} % of sequences remaining)")

        best_candidates.append(naive_guess)
        
        if (seq_remaining/total_sequences) < remaining_sequences_threshold:
            log.debug("Stopping search, too few sequences remaining")
            break
        log.debug("Continuing search")
    return motif_graph, best_candidates



class MotifSearcher:
    """
    Class to perform motif search for identifying enriched motifs in methylated sequences.
    """

    def __init__(
        self,
        root_motif: Motif,
        contig_sequence: DNAsequence,
        contig_pssm: DNAarray,
        contig_pileup: pl.DataFrame,
        methylation_sequences: DNAarray,
        padding: int,
        read_level_methylation: bool = False,
        na_positions = None,
        motif_graph: Optional[MotifTree] = None,
        min_kl: float = 0.1,
        freq_threshold: float = 0.25,
        max_rounds_since_new_best: int = 10,
        max_motif_length: int = 18
    ):
        """
        Initialize the MotifSearcher class.

        Parameters:
            root_motif (Motif): The initial motif to start the search from.
            contig_sequence: The contig sequence object with necessary methods.
            contig_pssm: The position-specific scoring matrix for the contig background.
            contig_pileup: The pileup data associated with the contig.
            methylation_sequences: The methylation sequences to be analyzed.
            padding (int): The padding size for sequence sampling.
            read_level_methylation (bool): Flag to use read-level methylation.
            na_positions (Optional[List[int]]): Positions with 'N/A' or missing data.
            motif_graph (Optional[MotifTree]): An optional pre-initialized motif graph.
            min_kl (float): Minimum Kullback-Leibler divergence threshold.
            freq_threshold (float): Frequency threshold for methylation.
            max_rounds_since_new_best (int): Max iterations without improvement.
            max_motif_length (int): Maximum allowed motif length.
        """
        # Input validation
        if contig_pileup.is_empty():
            raise ValueError("contig_pileup cannot be empty.")
        if padding < 0:
            raise ValueError("padding must be non-negative.")

        self.root_motif = root_motif
        self.contig_sequence = contig_sequence
        self.contig_pssm = contig_pssm
        self.contig_pileup = contig_pileup
        self.methylation_sequences = methylation_sequences
        self.padding = padding
        self.read_level_methylation = read_level_methylation
        self.na_positions = na_positions
        self.motif_graph = motif_graph or MotifTree()
        self.min_kl = min_kl
        self.freq_threshold = freq_threshold
        self.max_rounds_since_new_best = max_rounds_since_new_best
        self.max_motif_length = max_motif_length

    def _priority_function(self, next_model, root_model) -> float:
        """
        Calculate the priority for a motif based on the change in model parameters.

        Parameters:
            next_model: The model of the next motif.
            root_model: The model of the root motif.

        Returns:
            float: The calculated priority.
        """
        try:
            d_alpha = 1 - (next_model._alpha / root_model._alpha)
        except ZeroDivisionError:
            d_alpha = 1
        try:
            d_beta = next_model._beta / root_model._beta
        except ZeroDivisionError:
            d_beta = 1
        priority = d_alpha * d_beta
        return priority

    def _scoring_function(self, next_model, current_model) -> float:
        """
        Calculate the score for a motif based on its model.

        Parameters:
            next_model: The model of the next motif.
            current_model: The model of the current motif.

        Returns:
            float: The calculated score.
        """
        mean_diff = next_model.mean() - current_model.mean()
        try:
            score = (next_model.mean() * -np.log10(next_model.standard_deviation())) * mean_diff
        except (ValueError, ZeroDivisionError):
            score = 0
        return score

    def _motif_child_nodes_kl_dist_max(
        self,
        motif: Motif,
        meth_pssm: np.ndarray,
        contig_pssm: np.ndarray
    ) -> Generator[Motif, None, None]:
        """
        Generate child motifs based on maximum KL divergence.

        Parameters:
            motif (Motif): The current motif.
            meth_pssm (np.ndarray): Position-specific scoring matrix for methylation.
            contig_pssm (np.ndarray): Position-specific scoring matrix for the contig.

        Yields:
            Motif: The next motif in the search.
        """
        kl_divergence = entropy(meth_pssm, contig_pssm)
        split_motif = motif.split()

        # Positions to evaluate (positions with '.')
        evaluated_positions = np.array([i for i, base in enumerate(split_motif) if base == "."])
        if evaluated_positions.size == 0:
            return

        kl_divergence_masked = kl_divergence.copy()
        kl_divergence_masked[~np.isin(np.arange(len(split_motif)), evaluated_positions)] = 0

        max_kl = np.max(kl_divergence_masked)
        if max_kl < self.min_kl:
            return

        pos = int(np.argmax(kl_divergence_masked))

        # Methylation frequency must be above contig frequency
        index_meth_freq_higher = meth_pssm[:, pos] > contig_pssm[:, pos]

        # Methylation frequency must be above a threshold
        index_meth_freq_above_thresh = meth_pssm[:, pos] > self.freq_threshold

        # Combine the two filters
        index_position_filter = np.logical_and(index_meth_freq_higher, index_meth_freq_above_thresh)
        bases_indices = np.argwhere(index_position_filter).reshape(-1)
        bases_filtered = [BASES[int(i)] for i in bases_indices]

        if not bases_filtered:
            return

        # All combinations of the bases
        max_combination_length = min(len(bases_filtered), 4)
        for i in range(1, max_combination_length + 1):
            for base_tuple in itertools.combinations(bases_filtered, i):
                if len(base_tuple) > 1:
                    base_str = "[" + "".join(base_tuple) + "]"
                else:
                    base_str = base_tuple[0]
                new_motif_sequence = split_motif.copy()
                new_motif_sequence[pos] = base_str
                new_motif_str = "".join(new_motif_sequence)
                yield Motif(new_motif_str, motif.mod_position)

    def run(self) -> tuple[MotifTree, Motif]:
        """
        Execute the motif search algorithm.

        Returns:
            Tuple[MotifTree, Motif]: The graph of motifs explored and the best motif found.
        """
        # Initialize variables
        best_guess = self.root_motif
        if self.read_level_methylation:
            root_model = motif_model_read(
                self.contig_pileup, self.contig_sequence.sequence, self.root_motif
            )
        else:
            root_model = motif_model_contig(
                self.contig_pileup,
                self.contig_sequence.sequence,
                self.root_motif,
                na_positions=self.na_positions
            )
        best_score = self._scoring_function(root_model, root_model)
        rounds_since_new_best = 0
        visited_nodes: set[Motif] = set()
        
        # Initialize the search tree
        self.motif_graph.add_node(
            self.root_motif,
            model=root_model,
            motif=self.root_motif,
            visited=False,
            score=best_score
        )

        # Initialize priority queue
        priority_queue: list[tuple[float, Motif]] = []
        hq.heappush(priority_queue, (0, self.root_motif))

        # Search loop
        while priority_queue:
            # Get the current best candidate
            _, current_motif = hq.heappop(priority_queue)
            if current_motif in visited_nodes:
                continue
            if len(current_motif.strip()) > self.max_motif_length:
                log.debug(
                    f"{current_motif.string}, Skipping due to length > {self.max_motif_length}"
                )
                continue

            current_model = self.motif_graph.nodes[current_motif]["model"]
            visited_nodes.add(current_motif)
            log.debug(
                f"{current_motif.string} | Model: {current_model} | "
                f"Score: {self.motif_graph.nodes[current_motif]['score']:.2f} | "
                f"Queue size: {len(priority_queue)}"
            )
            self.motif_graph.nodes[current_motif]["visited"] = True
            rounds_since_new_best += 1

            active_methylation_sequences = self.methylation_sequences.copy().filter_sequence_matches(
                current_motif.one_hot(), keep_matches=True
            )
            if active_methylation_sequences is None:
                log.debug("No more sequences left after filtering")
                continue

            # Generate and evaluate neighbor motifs
            neighbors = list(self._motif_child_nodes_kl_dist_max(
                current_motif,
                active_methylation_sequences.pssm(),
                self.contig_pssm
            ))

            # Add neighbors to graph
            for next_motif in neighbors:
                if next_motif in self.motif_graph.nodes:
                    # Add only edge if motif already visited
                    next_model = self.motif_graph.nodes[next_motif]["model"]
                    score = self.motif_graph.nodes[next_motif]["score"]
                    self.motif_graph.add_edge(current_motif, next_motif)
                else:
                    if self.read_level_methylation:
                        next_model = motif_model_read(
                            self.contig_pileup,
                            self.contig_sequence.sequence,
                            next_motif
                        )
                    else:
                        next_model = motif_model_contig(
                            self.contig_pileup,
                            self.contig_sequence.sequence,
                            next_motif,
                            na_positions=self.na_positions
                        )

                    # Add neighbor to graph
                    self.motif_graph.add_node(
                        next_motif,
                        model=next_model,
                        motif=next_motif,
                        visited=False
                    )
                    self.motif_graph.add_edge(current_motif, next_motif)

                    score = self._scoring_function(next_model, current_model)
                    self.motif_graph.nodes[next_motif]["score"] = score

                # Add neighbor to priority queue if not visited
                if next_motif not in visited_nodes:
                    priority = self._priority_function(next_model, current_model)
                    hq.heappush(priority_queue, (priority, next_motif))

                if score > best_score:
                    best_score = score
                    best_guess = next_motif
                    rounds_since_new_best = 0

            # Stopping criteria
            if rounds_since_new_best >= self.max_rounds_since_new_best:
                break

        return self.motif_graph, best_guess

                    # Add neighbor to graph
                    self.motif_graph.add_node(
                        next_motif,
                        model=next_model,
                        motif=next_motif,
                        visited=False
                    )
                    self.motif_graph.add_edge(current_motif, next_motif)

                    score = self._scoring_function(next_model, current_model)
                    self.motif_graph.nodes[next_motif]["score"] = score

                # Add neighbor to priority queue if not visited
                if next_motif not in visited_nodes:
                    priority = self._priority_function(next_model, current_model)
                    hq.heappush(priority_queue, (priority, next_motif))

                if score > best_score:
                    best_score = score
                    best_guess = next_motif
                    rounds_since_new_best = 0

            # Stopping criteria
            if rounds_since_new_best >= self.max_rounds_since_new_best:
                break

        return self.motif_graph, best_guess

def count_periods_at_start(s):
    count = 0
    for char in s:
        if char == '.':
            count += 1
        else:
            break
    return count
def count_periods_at_end(s):
    s = s[::-1]
    count = 0
    for char in s:
        if char == '.':
            count += 1
        else:
            break
    return count
def nxgraph_to_dataframe(graph):
    return pl.DataFrame({
        "sequence":[i for i in graph.nodes],
        "model":[graph.nodes[i]["model"] for i in graph.nodes],
        "score":[graph.nodes[i]["score"] for i in graph.nodes]
    }).sort("score", descending=True)


