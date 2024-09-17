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
from nanomotif.seq import EqualLengthDNASet, DNAsequence
from nanomotif.candidate import Motif, MotifTree
from nanomotif.logger import configure_logger
from nanomotif.parallel import update_progress_bar
from nanomotif.seed import set_seed
import heapq as hq
import networkx as nx
import warnings
import time


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
        na_positions: dict=None
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
        args, 
        counter, 
        lock, 
        assembly, 
        min_kl_divergence, 
        padding, 
        minimum_methylation_fraction_confident, 
        read_level_methylation,
        log_dir, 
        verbose, 
        seed
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
    set_seed(seed=seed)
    warnings.filterwarnings("ignore")
    contig, modtype, subpileup, low_coverage_positions = args

    if low_coverage_positions is not None:
        low_coverage_positions = {
            "fwd": low_coverage_positions.filter(pl.col("strand") == "+")["position"].to_numpy(),
            "rev": low_coverage_positions.filter(pl.col("strand") == "-")["position"].to_numpy()
        }
    
    process_id = os.getpid()
    if log_dir is not None:
        log_file = f"{log_dir}/find-motifs.{process_id}.log"
        configure_logger(log_file, verbose=verbose)

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
        return result
    except:
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
    - threshold_valid_coverage (int): The minimum number of reads that must cover a position for it to be considered valid.
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
    assert threshold_valid_coverage >= 0, "min_valid_coverage must be greater than 0"
    assert minimum_kl_divergence >= 0, "mininum_kl_divergence must be greater than 0"

    # Infer padding size from candidate_size
    padding = search_frame_size // 2

    # Filter pileup
    pileup = pileup \
            .filter(pl.col("Nvalid_cov") > threshold_valid_coverage) \
            .sort("contig") 
    
    # Create a list of tasks (TODO: not have a list of all data)
    if low_coverage_positions is not None:
        tasks = [(contig, modtype, subpileup, low_coverage_positions.filter((pl.col("contig")==contig) & (pl.col("mod_type") == modtype))) for (contig, modtype), subpileup in pileup.group_by(["contig", "mod_type"])]
    else:
        tasks = [(contig, modtype, subpileup, None) for (contig, modtype), subpileup in pileup.group_by(["contig", "mod_type"])]
    # Create a progress manager
    manager = multiprocessing.Manager()
    counter = manager.Value('i', 0)
    lock = manager.Lock()

    # Create a pool of workers
    pool = get_context("spawn").Pool(processes=threads)

    # Create a process for the progress bar
    progress_bar_process = multiprocessing.Process(target=update_progress_bar, args=(counter, len(tasks), True))
    progress_bar_process.start()

    # Put them workers to work
    results = pool.starmap(worker_function, [(
        task, 
        counter, lock, 
        assembly, minimum_kl_divergence, padding, threshold_methylation_confident, read_level_methylation,
        log_dir, verbose, seed
        ) for task in tasks])
    results = [result for result in results if result is not None]

    # Close the pool
    pool.close()
    pool.join()

    # Close the progress bar
    progress_bar_process.join()

    if len(results) == 0:
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
        na_positions: dict = None,
        read_level_methylation: bool = False,
        min_kl: float = 0.2, 
        max_dead_ends: int = 25, 
        max_rounds_since_new_best: int = 15,
        score_threshold: float = 0.2,
        remaining_sequences_threshold: float = 0.01
    ) -> tuple:
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
        motif_graph, naive_guess = motif_search(
            root_motif, 
            contig_sequence, 
            contig_pileup, 
            methylation_sequences_clone, 
            padding,
            read_level_methylation = read_level_methylation,
            na_positions = na_positions,
            motif_graph = motif_graph,
            min_kl = min_kl,
            max_rounds_since_new_best = max_rounds_since_new_best
        )

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


 
def motif_child_nodes_kl_dist_max(
        motif, 
        meth_pssm, 
        contig_pssm, 
        freq_threshold=0.25, 
        min_kl=0.1
    ):
    kl_divergence = entropy(meth_pssm, contig_pssm)
    split_motif = motif.split()

    evaluated = np.array([i for i, base in enumerate(split_motif) if base != "."])
    kl_divergence[evaluated] = 0
    if np.max(kl_divergence) < min_kl:
        return
    
    pos = np.where(kl_divergence == np.max(kl_divergence))[0][0]

    # Methylation frequency most be above contig frequency
    index_meth_frequncies_highest = meth_pssm[:, pos] > contig_pssm[:, pos]

    # Methylation frequency most be above a threshold
    index_meth_frequncies_above_threshold = meth_pssm[:, pos] > freq_threshold

    # Combine the two filters
    index_position_filt = np.logical_and(index_meth_frequncies_highest, index_meth_frequncies_above_threshold)
    bases_index = np.argwhere(index_position_filt).reshape(-1)
    bases_filt = [BASES[int(i)] for i in bases_index]
    
    # All combination of the bases
    combinations = []
    for i in range(1, min(len(bases_filt)+1, 4)):
        combinations += list(itertools.combinations(bases_filt, i))
    
    for base in combinations:
        if len(base) > 1:
            base = "[" + "".join(list(base)) + "]"
        else:
            base = "".join(list(base))
        new_motif = split_motif[:pos] + [base] + split_motif[pos+1:]
        yield Motif("".join(new_motif), motif.mod_position)

 
def motif_child_nodes_kl_dist_prune(motif, meth_pssm, contig_pssm, min_kl=0.1, freq_threshold=0.35):
    kl_divergence = entropy(meth_pssm, contig_pssm)
    split_motif = motif.split()
    
    for pos in np.where(kl_divergence > min_kl)[0]:
        if split_motif[pos] != ".":
            continue

        # Methylation frequency most be above contig frequency
        index_meth_frequncies_highest = meth_pssm[:, pos] > contig_pssm[:, pos]

        # Methylation frequency most be above a threshold
        index_meth_frequncies_above_threshold = meth_pssm[:, pos] > freq_threshold

        # Combine the two filters
        index_position_filt = np.logical_and(index_meth_frequncies_highest, index_meth_frequncies_above_threshold)
        bases_index = np.argwhere(index_position_filt).reshape(-1)
        bases_filt = [BASES[int(i)] for i in bases_index]
        
        # All combination of the bases
        combinations = []
        for i in range(1, min(len(bases_filt)+1, 4)):
            combinations += list(itertools.combinations(bases_filt, i))
        
        for base in combinations:
            if len(base) > 1:
                base = "[" + "".join(list(base)) + "]"
            else:
                base = "".join(list(base))
            new_motif = split_motif[:pos] + [base] + split_motif[pos+1:]
            yield Motif("".join(new_motif), motif.mod_position)

def motif_search(
        root_motif, 
        contig_sequence, 
        contig_pileup, 
        methylation_sequences, 
        padding,
        read_level_methylation = False,
        na_positions = None,
        motif_graph = None, 
        min_kl = 0.1, 
        freq_threshold = 0.25,
        max_rounds_since_new_best = 10, 
        max_motif_length = 18
    ) -> tuple:
    """
    Search algorithm for identifying enriched motif in  methylated sequences

    Parameters:
    - root_motif (list): The root motif to start the search from.
    - contig_sequence: The contig to be processed.
    - contig_pileup (Pileup): The pileup to be processed.
    - methylation_sequences (DNAarray): The methylation sequences to be processed.
    """
    
    # Function for calculating the priority of a motif
    def priority_function(next_model, root_model):
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


    # Function for scoring a motif
    def scoring_function(next_model, current_model):
        mean_diff =  next_model.mean() - current_model.mean()
        return (next_model.mean() * -np.log10(next_model.standard_deviation())) * mean_diff

    # Search setup
    best_guess = root_motif
    if read_level_methylation:
        root_model = motif_model_read(contig_pileup, contig_sequence.sequence, root_motif)
    else:
        root_model = motif_model_contig(contig_pileup, contig_sequence.sequence, root_motif, na_positions=na_positions)
    best_score = scoring_function(root_model, root_model)
    rounds_since_new_best = 0
    visisted_nodes = []

    # Sample sequence in contig to get background for KL-divergence
    contig_sequences = contig_sequence.sample_n_subsequences(padding*2 + 1, 10000)
    contig_pssm = contig_sequences.pssm()

    # Initialize the search tree
    if motif_graph is None:
        motif_graph = MotifTree()
    motif_graph.add_node(root_motif, model=root_model, motif=root_motif, visited=False, score=scoring_function(root_model, root_model))

    # Initialize priority que
    priority_que = []
    hq.heappush(priority_que, (0, root_motif))

    # Search loop
    while len(priority_que) > 0:
        # Get the current best candidate
        current = hq.heappop(priority_que)[1]
        if current in visisted_nodes:
            continue
        if len(current.strip()) > max_motif_length:
            log.debug(f"{current}, Skipping scoring and expansion due to length (>{max_motif_length})")
            continue


        current_model = motif_graph.nodes[current]["model"]
        visisted_nodes.append(current)
        log.debug(f"{''.join(current)} | {current_model} | {motif_graph.nodes[current]['score']:.2f} | {len(priority_que)}")
        motif_graph.nodes[current]["visited"] = True
        rounds_since_new_best += 1
        active_methylation_sequences = methylation_sequences.copy().filter_sequence_matches(current.one_hot(), keep_matches = True)
        if active_methylation_sequences is None:
            log.debug("No more sequences left")
            continue
        # Prune the neibhors of the current candidate
        neighbors = list(motif_child_nodes_kl_dist_max(
            current, 
            active_methylation_sequences.pssm(), 
            contig_pssm,
            freq_threshold=freq_threshold,
            min_kl=min_kl

        ))

        # Add neighbors to graph
        for next in neighbors:
            if next in motif_graph.nodes:
                # Add only edge if motif alredy visited
                next_model = motif_graph.nodes[next]["model"]
                score = motif_graph.nodes[next]["score"]
                motif_graph.add_edge(current, next)
            else:
                if read_level_methylation:
                    next_model = motif_model_read(contig_pileup, contig_sequence.sequence, next)
                else:
                    next_model = motif_model_contig(contig_pileup, contig_sequence.sequence, next, na_positions=na_positions)

                # Add neighbor to graph
                motif_graph.add_node(next, model=next_model, motif=next, visited=False)
                motif_graph.add_edge(current, next)

                score = scoring_function(next_model, current_model)
                motif_graph.nodes[next]["score"] = score

            # Add neighbor to priority que if has not been visited in this search
            if next not in visisted_nodes:
                priority =  priority_function(next_model, current_model)
                hq.heappush(priority_que, (priority, next))
            
            if score > best_score:
                best_score = score
                best_guess = next
                rounds_since_new_best = 0

        # Stopping criteria
        if rounds_since_new_best >= max_rounds_since_new_best:
            break
    return motif_graph, best_guess


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

def count_periods_at_end(string: str) -> int:
    """
    Count the number of periods at the end of a string

    Parameters:
    - string (str): The string to be processed.

    Returns:
    - int: The number of periods at the end of the string.

    Example:
    >>> count_periods_at_end("ACGT...")
    3
    """
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


