import numpy as np
import polars as pl
import os
from polars import col
import logging as log
import itertools
from itertools import count
from itertools import takewhile
from scipy.stats import entropy
import re
import multiprocessing
from multiprocessing import get_context
import progressbar
import nanomotif as nm
from nanomotif.constants import *
from nanomotif.model import BetaBernoulliModel
from nanomotif.utils import subseq_indices, calculate_match_length
from nanomotif.seq import EqualLengthDNASet, DNAsequence, DNAarray
from nanomotif.motif import Motif, MotifTree
from nanomotif.logger import configure_logger
from nanomotif.parallel import update_progress_bar
from nanomotif.seed import set_seed
import heapq as hq
import networkx as nx
import warnings
import time
from typing import Generator, Optional, List, Tuple


def set_polars_env():
    os.environ["POLARS_MAX_THREADS"] = "1"



def process_binned_sample_parallel(
        assembly, pileup, bin_contig,
        threads = 2,
        search_frame_size = 40,
        methylation_threshold_low = 0.3,
        methylation_threshold_high = 0.7,
        minimum_kl_divergence = 0.05,
        verbose = False,
        log_dir = None,
        seed = None,
        score_threshold = 2
    ):
    """
    Process a sample
    
    Parameters:
    - assembly (Assembly): The assembly of all bins.
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
    assert methylation_threshold_high >= 0 and methylation_threshold_high <= 1, "threshold_methylation_confident must be between 0 and 1"
    assert minimum_kl_divergence >= 0, "mininum_kl_divergence must be greater than 0"

    # Infer padding size from candidate_size
    padding = search_frame_size // 2
    
    pileup = pileup.join(bin_contig, on="contig", how = "inner")

    log.debug("Partitining pileups to dict")
    pileup_dict = pileup.partition_by(["bin", "mod_type"], as_dict=True)

    def task_generator(pileup_dict):
        for (bin_name, modtype) in pileup_dict:
            subpileup = pileup_dict[(bin_name, modtype)]
            bin_contig_dict = {
                bin_name: subpileup.get_column("contig").unique().to_list()
            }
            yield (bin_contig_dict, modtype, subpileup, counter, lock, assembly, 
                       minimum_kl_divergence, search_frame_size // 2,
                       log_dir, verbose, seed,
                       methylation_threshold_low, methylation_threshold_high,
                       score_threshold)

    # Create a progress manager
    manager = multiprocessing.Manager()
    counter = manager.Value('i', 0)
    lock = manager.Lock()

    # Create a pool of workers
    pool = get_context("spawn").Pool(processes=threads, initializer=set_polars_env)

    # Create a process for the progress bar
    log.debug("Counting number of tasks")
    total_tasks = len(pileup.select("bin", "mod_type").unique())

    log.debug("Starting progress bar")
    progress_bar_process = multiprocessing.Process(target=update_progress_bar, args=(counter, total_tasks, True))
    progress_bar_process.start()

    # Put them workers to work
    log.debug("Starting workers")
    chunksize = max(min(100, total_tasks // threads), 1)
    results = pool.imap(worker_function, task_generator(pileup_dict), chunksize=chunksize)
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


##########################################
# Motif candidate state space search
##########################################

def worker_function(
        args
    ):
    """
    Process a single subpileup for one bin and one modtype

    Parameters:
    - args (tuple): The arguments to the function: bin, modtype, subpileup
    - counter (multiprocessing.Value): The progress counter
    - lock (multiprocessing.Lock): The lock for the progress counter
    - assembly (Assembly): The assembly to be processed.
    - min_kl_divergence (float): Early stopping criteria, if max KL-divergence falls below, stops building motif.
    - padding (int): The padding to use for the motif.
    """
    (
        bin_contig, modtype, subpileup, counter, lock, assembly, 
        min_kl_divergence, padding, 
        log_dir, verbose, seed, 
        low_meth_threshold, high_meth_threshold,
        score_threshold
    ) = args
    bin_name = list(bin_contig.keys())[0]
    set_seed(seed=seed)
    warnings.filterwarnings("ignore")
    process_id = os.getpid()
    if log_dir is not None:
        log_file = f"{log_dir}/find-motifs.{process_id}.log"
        configure_logger(log_file, verbose=verbose)
    log.info(f"Worker {process_id} started")


    result = process_subpileup(
        bin_contig, modtype, subpileup, assembly, 
        min_kl_divergence, padding, 
        low_meth_threshold, high_meth_threshold,
        score_threshold
    )
    with lock:
        counter.value += 1
        log.debug(f"Task for ({bin_name}, {modtype}) complete; progress: {counter.value}")
    return result



def process_subpileup(
        bin_contig: dict, modtype, bin_pileup, assembly, 
        min_kl_divergence, padding, 
        low_meth_threshold, high_meth_threshold,
        score_threshold
    ):
    """
    Process a single subpileup for one bin and one modtype
    """
    bin_name = list(bin_contig.keys())[0]
    log.info(f"Processing {bin_name} {modtype}")
    log.debug(f"Contigs in bin: {list(bin_contig[bin_name])}")
    assert bin_pileup is not None, "Subpileup is None"
    assert len(bin_pileup) > 0, "Subpileup is empty"
    assert assembly is not None, "Assembly is None"
    assert min_kl_divergence >= 0, "min_kl_divergence must be greater than 0"
    assert bin_pileup.get_column("mod_type").unique().to_list() == [modtype], "subpileup modtype does not match modtype"

    log.debug(f"Subpileup for {bin_name} {modtype}: {bin_pileup}")
    bin_sequences = {contig: assembly[contig] for contig in bin_contig[bin_name]}

    motif_graph, best_candidates = find_best_candidates(
        bin_pileup, 
        bin_sequences, 
        modtype,
        low_meth_threshold=low_meth_threshold,
        high_meth_threshold=high_meth_threshold,
        padding=padding,
        min_kl=min_kl_divergence,
        max_dead_ends=25,
        max_rounds_since_new_best=15,
        score_threshold=score_threshold
    )
    identified_motifs = nxgraph_to_dataframe(motif_graph) \
        .filter(col("sequence").is_in(best_candidates))
    
    if len(identified_motifs) == 0:
        log.info("No motifs found")
        return None
    else:
        identified_motifs = identified_motifs.with_columns(
            pl.lit(bin_name).alias("bin"),
            pl.lit(modtype).alias("mod_type")
        ).with_columns([
            pl.col("model").map_elements(lambda x: x._alpha, return_dtype=pl.Float64).alias("alpha"),
            pl.col("model").map_elements(lambda x: x._beta, return_dtype=pl.Float64).alias("beta")
        ]).drop("model")
        return identified_motifs






#########################################################################
# Motif candidate state space 

def find_best_candidates(        
        bin_pileup: pl.DataFrame, 
        bin_sequences: dict[str, DNAsequence],
        mod_type: str,
        low_meth_threshold: float,
        high_meth_threshold: float,
        padding: int,
        min_kl: float = 0.2, 
        max_dead_ends: int = 25, 
        max_rounds_since_new_best: int = 15,
        score_threshold: float = 0.2,
        remaining_sequences_threshold: float = 0.01
    ) -> tuple[MotifTree, list[Motif]]:
    """
    Find the best motif candidates in a sequence.
    """
    subpileup_confident = bin_pileup.filter(pl.col("fraction_mod") >= high_meth_threshold)
    methylation_sequences = None
    background_sequences = None
    for contig in bin_pileup.get_column("contig").unique():
        subpileup_confident_contig = subpileup_confident.filter(pl.col("contig") == contig)
        contig_sequence = bin_sequences[contig]
        contig_length = len(contig_sequence)
        n_samples = int(max(contig_length/400, 100))
        # Extract the sequences for confidently methylated positions
        index_plus = subpileup_confident_contig.filter(pl.col("strand") == "+").get_column("position").to_list()
        index_minus = subpileup_confident_contig.filter(pl.col("strand") == "-").get_column("position").to_list()
        
        # Sample sequences in bin to get background for KL-divergence
        log.debug(f"Sampling {n_samples} sequences from {contig}")
        contig_sequences_sample = contig_sequence.sample_n_subsequences(
            padding * 2 + 1, n_samples
        )
        if background_sequences is None:
            background_sequences = contig_sequences_sample
        else:
            background_sequences = background_sequences + contig_sequences_sample.sequences

        log.debug(f"Sampling methylation sequences from {contig}")
        methylation_sequences_string = []
        if len(index_plus) >= 1:
            methylation_sequences_string += contig_sequence.sample_at_indices(index_plus, padding).sequences
        if len(index_minus) >= 1:
            methylation_sequences_string += contig_sequence.sample_at_indices(index_minus, padding).reverse_compliment().sequences
        if len(methylation_sequences_string) == 0:
            log.info("No methylation sequences found")
            return None
        methylation_sequences_contig = EqualLengthDNASet(methylation_sequences_string)
        if methylation_sequences is None:
            methylation_sequences = methylation_sequences_contig
        else:
            methylation_sequences = methylation_sequences + methylation_sequences_contig
    methylation_sequences = methylation_sequences.convert_to_DNAarray()
    bin_pssm = background_sequences.pssm()

    total_sequences = methylation_sequences.shape[0]
    root_motif = Motif("." * padding + MOD_TYPE_TO_CANONICAL[mod_type] + "." * padding, padding)
    methylation_sequences_clone = methylation_sequences.copy()
    best_candidates = []
    continue_search = True
    dead_ends = 0
    motif_graph = None

    while continue_search:
        log.debug(f"Searching for best candidates, {len(best_candidates)} found so far")
        if dead_ends >= max_dead_ends:
            log.debug("Stopping search, too many low scoring candidates")
            break
        # Find the initial guess within the tree
        searcher = MotifSearcher(
            root_motif, 
            bin_sequences, 
            bin_pssm,
            bin_pileup, 
            methylation_sequences_clone, 
            padding,
            high_meth_threshold,
            low_meth_threshold,
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
            log.debug(f"Candidate has low score ({motif_graph.nodes[naive_guess]['score']}), {naive_guess}. {dead_ends} of {max_dead_ends} before termination. Represented in {seq_before-seq_remaining} seqs. model: {motif_graph.nodes[naive_guess]['model']}. ({100*seq_remaining_percent:.1f} % of sequences remaining)")
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
        bin_sequence: DNAsequence,
        bin_pssm: DNAarray,
        bin_pileup: pl.DataFrame,
        methylation_sequences: DNAarray,
        padding: int,
        high_meth_threshold: float,
        low_meth_threshold: float,
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
            bin_sequence: The bin sequence object with necessary methods.
            bin_pssm: The position-specific scoring matrix for the bin background.
            bin_pileup: The pileup data associated with the bin.
            methylation_sequences: The methylation sequences to be analyzed.
            padding (int): The padding size for sequence sampling.
            read_level_methylation (bool): Flag to use read-level methylation.
            motif_graph (Optional[MotifTree]): An optional pre-initialized motif graph.
            min_kl (float): Minimum Kullback-Leibler divergence threshold.
            freq_threshold (float): Frequency threshold for methylation.
            max_rounds_since_new_best (int): Max iterations without improvement.
            max_motif_length (int): Maximum allowed motif length.
        """
        # Input validation
        if bin_pileup.is_empty():
            raise ValueError("bin_pileup cannot be empty.")
        if padding < 0:
            raise ValueError("padding must be non-negative.")

        self.root_motif = root_motif
        self.bin_sequence = bin_sequence
        self.bin_pssm = bin_pssm
        self.bin_pileup = bin_pileup
        self.methylation_sequences = methylation_sequences
        self.padding = padding
        self.motif_graph = motif_graph or MotifTree()
        self.min_kl = min_kl
        self.freq_threshold = freq_threshold
        self.max_rounds_since_new_best = max_rounds_since_new_best
        self.max_motif_length = max_motif_length
        self.low_meth_threshold = low_meth_threshold
        self.high_meth_threshold = high_meth_threshold

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

    def _scoring_function_predictive_evaluation(self, next_model, current_model) -> float:
        """
        Calculate the score for a motif based on its predictive model.

        Parameters:
            next_model: The model of the next motif.
            current_model: The model of the current motif.

        Returns:
            float: The calculated score.
        """
        extra_positive, extra_negative = current_model._alpha - next_model._alpha, current_model._beta - next_model._beta

        log_likelihood_next_per_obs = next_model.log_likelihood_per_obs()
        log_likelihood_extra_per_obs = next_model.posterior_predictive_per_obs(extra_positive, extra_negative)
        mean_ratio = next_model.mean() / current_model.mean() if current_model.mean() > 0 else next_model.mean()

        return mean_ratio * (log_likelihood_next_per_obs - log_likelihood_extra_per_obs)

    def _motif_child_nodes_kl_dist_max(
        self,
        motif: Motif,
        meth_pssm: np.ndarray,
        bin_pssm: np.ndarray
    ) -> Generator[Motif, None, None]:
        """
        Generate child motifs based on maximum KL divergence.

        Parameters:
            motif (Motif): The current motif.
            meth_pssm (np.ndarray): Position-specific scoring matrix for methylation.
            bin_pssm (np.ndarray): Position-specific scoring matrix for the bin.

        Yields:
            Motif: The next motif in the search.
        """
        kl_divergence = entropy(meth_pssm, bin_pssm)
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

        # Methylation frequency must be above bin frequency
        index_meth_freq_higher = meth_pssm[:, pos] > bin_pssm[:, pos] * 0.5

        # Methylation frequency must be above a threshold
        index_meth_freq_above_thresh = meth_pssm[:, pos] > self.freq_threshold

        # Combine the two filters
        index_position_filter = np.logical_and(index_meth_freq_higher, index_meth_freq_above_thresh)
        bases_indices = np.argwhere(index_position_filter).reshape(-1)
        bases_filtered = [BASES[int(i)] for i in bases_indices]

        if not bases_filtered:
            return

        # All combinations of the bases
        #max_combination_length = min(len(bases_filtered), 4)
        #for i in range(1, max_combination_length + 1):
        #    for base_tuple in itertools.combinations(bases_filtered, i):
        #        if len(base_tuple) > 1:
        #            base_str = "[" + "".join(base_tuple) + "]"
        #        else:
        #            base_str = base_tuple[0]
        #        new_motif_sequence = split_motif.copy()
        #        new_motif_sequence[pos] = base_str
        #        new_motif_str = "".join(new_motif_sequence)
        #        yield Motif(new_motif_str, motif.mod_position)
        for base in bases_filtered:
            new_motif_sequence = split_motif.copy()
            new_motif_sequence[pos] = base
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
        root_model = motif_model_bin(
            pileup=self.bin_pileup,
            contigs=self.bin_sequence,
            motif=self.root_motif,
            model=BetaBernoulliModel(),
            low_meth_threshold=self.low_meth_threshold,
            high_meth_threshold=self.high_meth_threshold,
        )

        best_score = self._scoring_function_predictive_evaluation(root_model, root_model)
        rounds_since_new_best = 0
        visited_nodes: set[Motif] = set()

        root_depth = 0
        self.motif_graph.add_node(
            self.root_motif,
            model=root_model,
            motif=self.root_motif,
            visited=False,
            score=best_score,
            priority=0,
            depth=root_depth
        )

        # Initialize priority queue 
        priority_queue: list[tuple] = []
        hq.heappush(
            priority_queue,
            (root_depth, 0,  self.root_motif)
        )

        # Search loop
        while priority_queue:
            # Get the current best candidate
            _, _, current_motif = hq.heappop(priority_queue)
            if current_motif in visited_nodes:
                continue

            current_attrs = self.motif_graph.nodes[current_motif]
            current_model = current_attrs["model"] 
            current_depth = current_attrs.get("depth", 0)

            # Enforce motif length safely (current_model is defined)
            if len(current_motif.strip()) > self.max_motif_length:
                log.debug(
                    f"{current_motif.string} expansion skipped due to length | Model: {current_model} | "
                    f"Score: {current_attrs['score']:.2f} | "
                    f"Queue size: {len(priority_queue)} | "
                    f"Came from {', '.join(str(node) for node in list(self.motif_graph.predecessors(current_motif)))}"
                )
                continue

            visited_nodes.add(current_motif)
            log.debug(
                f"{current_motif.string} expanded | Model: {current_model} | "
                f"Score: {current_attrs['score']:.2f} | "
                f"Queue size: {len(priority_queue)} | "
                f"Came from {', '.join(str(node) for node in list(self.motif_graph.predecessors(current_motif)))}"
            )
            self.motif_graph.nodes[current_motif]["visited"] = True
            rounds_since_new_best += 1

            # Filter sequences for this motif
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
                self.bin_pssm
            ))

            # Add neighbors to graph
            for next_motif in neighbors:
                next_depth = current_depth + 1

                # Build model for the next motif
                next_model = motif_model_bin(
                    pileup=self.bin_pileup,
                    contigs=self.bin_sequence,
                    motif=next_motif,
                    model=BetaBernoulliModel(),
                    low_meth_threshold=self.low_meth_threshold,
                    high_meth_threshold=self.high_meth_threshold,
                )

                # Score relative to current
                score = self._scoring_function_predictive_evaluation(next_model, current_model)
                priority = self._priority_function(next_model, root_model)

                if next_motif in self.motif_graph.nodes:
                    # Update existing node if new score is better
                    if self.motif_graph.nodes[next_motif]["score"] < score:
                        self.motif_graph.nodes[next_motif]["score"] = score
                else:
                    self.motif_graph.add_node(
                        next_motif,
                        model=next_model,
                        motif=next_motif,
                        visited=False,
                        score=score,
                        priority=priority,
                        depth=next_depth
                    )
                self.motif_graph.add_edge(current_motif, next_motif)

                # Add neighbor to priority queue if not visited
                if next_motif not in visited_nodes:
                    attrs = self.motif_graph.nodes[next_motif]
                    hq.heappush(
                        priority_queue,
                        (attrs["depth"], -attrs["priority"], next_motif)
                    )

                # Track best
                if score > best_score:
                    best_score = score
                    best_guess = next_motif
                    rounds_since_new_best = 0

            # Stopping criteria
            if rounds_since_new_best >= self.max_rounds_since_new_best:
                break

        return self.motif_graph, best_guess





def count_periods_at_start(string: str) -> int:
    """
    Count the number of periods at the start of a string

    Parameters:
    - string (str): The string to be processed.

    Returns:
    - int: The number of periods at the start of the string.

    Example:
    >>> count_periods_at_start("...ACGT")
    3
    """
    count = 0
    for char in string:
        if char == '.':
            count += 1
        else:
            break
    return count
def count_periods_at_end(string):
    string = string[::-1]
    count = 0
    for char in string:
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






##########################################
# Motif candidate scoring
##########################################

def methylated_motif_occourances(
        motif, 
        sequence, 
        methylated_positions,
        non_methylated_positions
    ) -> tuple:
    """
    Get occourances of a motif in a bin

    Parameters:
    - motif (str): The motif to search for.
    - seq (str): The bin to search in.
    - bin_meth_positions (list): The positions of the methylation sites in the bin.

    Returns:
    - tuple: A tuple of two numpy arrays, the first containing the positions of the methylated motifs, the second containing the positions of the non-methylated motifs.
    """
    assert len(motif) > 0, "Motif is empty"
    assert len(sequence) > 0, "Sequence is empty"
    assert type(motif) == Motif, "Motif is not a Motif type"
    # Get the index of the methylation in the motif in the bin
    motif_index = subseq_indices(motif.string, sequence) + motif.mod_position

    # Methylated motif positions
    meth_occurences = np.intersect1d(methylated_positions, motif_index)

    # Non-methylated motif positions
    nonmeth_occurences =  np.intersect1d(non_methylated_positions, motif_index)

    return meth_occurences, nonmeth_occurences

def motif_model_bin(
        pileup,
        contigs: dict[str, DNAsequence],
        motif: str,
        model,
        low_meth_threshold,
        high_meth_threshold,
):
    for contig, contig_sequence in contigs.items():
        pileup_contig = pileup.filter(pl.col("contig") == contig)
        model = motif_model_contig(
            pileup_contig,
            contig_sequence.sequence,
            model,
            motif,
            low_meth_threshold=low_meth_threshold,
            high_meth_threshold=high_meth_threshold
        )
    return model

def motif_model_contig(
        pileup, 
        contig: str, 
        prior: BetaBernoulliModel,
        motif, 
        low_meth_threshold=0.3,
        high_meth_threshold=0.7,
        save_motif_positions=False,
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
    pileup_high_meth = pileup.filter(pl.col("fraction_mod") >= high_meth_threshold)
    pileup_low_meth = pileup.filter(pl.col("fraction_mod") <= low_meth_threshold)

    meth_positions_fwd = pileup_high_meth.filter(pl.col("strand") == "+")["position"].to_numpy()
    non_meth_positions_fwd = pileup_low_meth.filter(pl.col("strand") == "+")["position"].to_numpy()
    meth_positions_rev = pileup_high_meth.filter(pl.col("strand") == "-")["position"].to_numpy()
    non_meth_positions_rev = pileup_low_meth.filter(pl.col("strand") == "-")["position"].to_numpy()

    index_meth_fwd, index_nonmeth_fwd = methylated_motif_occourances(motif_stripped, contig, meth_positions_fwd, non_meth_positions_fwd)
    index_meth_rev, index_nonmeth_rev = methylated_motif_occourances(motif_stripped.reverse_compliment(), contig, meth_positions_rev, non_meth_positions_rev)

    model = prior
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



def merge_motifs_in_df(motif_df, pileup, assembly, contig_bin, low_meth_threshold=0.3, high_meth_threshold=0.7, merge_threshold = -0.2):
    new_df = []
    for (bin_name, mod_type), df in motif_df.group_by("bin", "mod_type"):
        contigs = contig_bin.filter(col("bin") == bin_name).get_column("contig").unique().to_list()
        log.debug(f"Merging motifs of bin: {bin_name}, modtype: {mod_type}")
        
        # Create dictionary of bin contigs to sequence
        bin_sequences = {contig: assembly[contig] for contig in contigs}


        # Get list of motifs
        motif_seq = df["motif"].to_list()
        motif_pos = df["mod_position"].to_list()
        motifs = [nm.motif.Motif(seq, pos) for seq, pos in zip(motif_seq, motif_pos)]

        # Merge motifs
        merged_motifs = nm.motif.merge_motifs(motifs)
        all_merged_motifs = []
        all_premerge_motifs = []
        # Check mean shift of premerge motifs to merged motif is high enough
        for _, motifs in merged_motifs.items():
            merged_motif = motifs[0]
            premerge_motifs = motifs[1]

            merge_model = motif_model_bin(
                pileup.filter((col("contig").is_in(contigs)) & (col("mod_type") == mod_type)), 
                bin_sequences,
                merged_motif,
                BetaBernoulliModel(),
                low_meth_threshold=low_meth_threshold,
                high_meth_threshold=high_meth_threshold
            )

            pre_merge_model = BetaBernoulliModel()
            for pre_merge_motif in premerge_motifs:
                pre_merge_model = motif_model_bin(
                    pileup.filter((col("contig").is_in(contigs)) & (col("mod_type") == mod_type)),
                    bin_sequences,
                    pre_merge_motif,
                    pre_merge_model,
                    low_meth_threshold=low_meth_threshold,
                    high_meth_threshold=high_meth_threshold
                )
            mean_diff = merge_model.mean() - pre_merge_model.mean()
            if mean_diff > merge_threshold:
                log.debug(f"Merging motifs {premerge_motifs} to {merged_motif}, mean diff: {mean_diff}")
            else:
                log.debug(f"Not merging motifs {premerge_motifs} to {merged_motif}, mean diff: {mean_diff}")
                continue

            all_merged_motifs.append(merged_motif)
            all_premerge_motifs.extend(premerge_motifs)
        log.info(f"Motif merge for bin {bin_name} modtype {mod_type} complete")
        # Create a new dataframe with the non merged motifs
        if  len(all_premerge_motifs) == 0:
            # No motifs were merged
            new_df.append(df)
            log.info(f"No motifs were merged for bin {bin_name} modtype {mod_type}")
            continue
        log.info(f"New merged motifs: {all_merged_motifs}")
        single_df = df.filter(col("motif").is_in(all_premerge_motifs).not_())
        new_df.append(single_df)
        merged_df = []
        for motif in all_merged_motifs:
            merged_model = motif_model_bin(
                pileup.filter((col("contig").is_in(contigs)) & (col("mod_type") == mod_type)), 
                bin_sequences,
                motif,
                BetaBernoulliModel(),
                low_meth_threshold=low_meth_threshold,
                high_meth_threshold=high_meth_threshold
            )
            merged_df.append(pl.DataFrame({
                "sequence": motif.string,
                "score": -1.,
                "bin": bin_name,
                "mod_type": mod_type,
                "model": merged_model,
                "motif": motif.string,
                "mod_position": motif.mod_position
            }).cast({'bin': pl.String}))

        merged_df_con = pl.concat(merged_df)
        new_df.append(merged_df_con)
    new_df = pl.concat(new_df)
    return new_df
