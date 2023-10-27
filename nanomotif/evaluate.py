import numpy as np
import polars as pl
from polars import col
import logging as log
import itertools
from itertools import takewhile
from scipy.stats import entropy
from nanomotif.constants import *
from nanomotif.model import BetaBernoulliModel
from nanomotif.utils import subseq_indices, calculate_match_length
from nanomotif.seq import EqualLengthDNASet, DNAsequence
from nanomotif.candidate import Motif, MotifTree
import heapq as hq
import networkx as nx

##########################################
# Motif candidate scoring
##########################################

def methylated_motif_occourances(motif: str, sequence: str, methylated_positions) -> tuple:
    """
    Get occourances of a motif in a contig

    Parameters:
    - motif (str): The motif to search for.
    - seq (str): The contig to search in.
    - contig_meth_positions (list): The positions of the methylation sites in the contig.

    Returns:
    - tuple: A tuple of two numpy arrays, the first containing the positions of the methylated motifs, the second containing the positions of the non-methylated motifs.
    """


    # Get the index of the methylation in the motif in the contig
    motif_index = subseq_indices(motif.string, sequence) + motif.mod_position

    # Methylated motif positions
    meth_occurences = np.intersect1d(methylated_positions, motif_index)

    # Non-methylated motif positions
    nonmeth_occurences =  np.setdiff1d(motif_index, methylated_positions)

    return meth_occurences, nonmeth_occurences

def motif_model_contig(pileup, contig: str, motif):
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

    model = BetaBernoulliModel()
    model.update(len(index_meth_fwd) + len(index_meth_rev), len(index_nonmeth_fwd) + len(index_nonmeth_rev))

    return model


def methylated_reads_counts(pileup: list, sequence: str, motif: str) -> tuple:
    motif_index = subseq_indices(motif.string, sequence) + motif.mod_position

    pileup = pileup.filter(pl.col("position").is_in(motif_index))
    counts = pileup \
        .select(["n_mod", "n_nonmod"]) \
        .sum()
    n_mod = counts.get_column("n_mod")[0]
    n_nonmod = counts.get_column("n_nonmod")[0]
    return n_mod, n_nonmod


def motif_model_read(pileup, contig: str, motif):
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

def process_sample(assembly, pileup, 
                   max_candidate_size = 40,
                   min_read_methylation_fraction = 0.8,
                   min_valid_coverage = 10,
                   min_kl_divergence = 0.1
                   ):
    """
    Process a single sample
    
    Parameters:
    - assembly (Assembly): The assembly to be processed.
    - pileup (Pileup): The pileup to be processed.
    - max_candidate_size (int): The maximum size of the candidate motifs.
    - min_read_methylation_fraction (float): The minimum fraction of reads that must be methylated for a position to be considered methylated.
    - min_valid_coverage (int): The minimum number of reads that must cover a position for it to be considered valid.
    - min_kl_divergence (float): Early stopping criteria, if max KL-divergence falls below, stops building motif.
    - min_cdf_score (float): Minimum score of 1 - cdf(cdf_position) for a motif to be considered valid.
    - cdf_position (float): The position to evaluate the cdf at.
    - min_motif_frequency (int): Used to get minimum number of sequences to evaluate motif at.
    """
    padding = max_candidate_size // 2
    result = []
    pileup = pileup \
            .filter(pl.col("Nvalid_cov") > min_valid_coverage) \
            .filter(pl.col("fraction_mod") >= min_read_methylation_fraction) \
            .sort("contig") \
            .groupby(["contig", "mod_type"])
    
    for (contig, modtype), subpileup in pileup:
        log.info(f"Processing {contig}")
        log.info(f"Processing {modtype}")

        contig_sequence = assembly.assembly[contig]
        index_plus = subpileup.filter(pl.col("fraction_mod") > min_read_methylation_fraction) \
            .filter(pl.col("Nvalid_cov") > min_valid_coverage) \
            .filter(pl.col("strand") == "+").get_column("position").to_list()
        index_minus = subpileup.filter(pl.col("fraction_mod") > min_read_methylation_fraction) \
            .filter(pl.col("Nvalid_cov") > min_valid_coverage) \
            .filter(pl.col("strand") == "-").get_column("position").to_list()
        if len(index_minus) <= 1 or len(index_plus) <= 1:
            continue

        sequences_plus = contig_sequence.sample_at_indices(index_plus, padding)
        sequences_minus = contig_sequence.sample_at_indices(index_minus, padding).reverse_compliment()
        sequences = sequences_plus + sequences_minus
        sequences_array = sequences.convert_to_DNAarray()

        motif_graph, best_candidates = find_best_candidates(
            sequences_array, 
            contig_sequence, 
            subpileup, 
            min_kl = min_kl_divergence
        )
        identified_motifs = nxgraph_to_dataframe(motif_graph) \
            .filter(col("sequence").is_in(best_candidates))
        
        if len(identified_motifs) == 0:
            continue
        else:
            result.append(identified_motifs.with_columns(
                pl.lit(contig).alias("contig"),
                pl.lit(modtype).alias("mod_type")
            ))
    if len(result) == 0:
        return None
    motifs = pl.concat(result) \
        .with_columns([
            pl.col("sequence").apply(lambda motif: motif[count_periods_at_start(motif):len(motif)-count_periods_at_end(motif)]).alias("motif"),
            pl.col("sequence").apply(lambda motif: padding - count_periods_at_start(motif)).alias("mod_position")
        ])
    return motifs


#########################################################################
# Motif candidate state space 

def find_best_candidates(methylation_sequences, sequence, pileup, min_kl = 0.2, max_dead_ends = 3):
    padding = methylation_sequences.shape[1] // 2
    total_sequences = methylation_sequences.shape[0]
    mod_type = pileup.get_column("mod_type").unique().to_list()[0]
    
    root_motif = Motif("." * padding + MOD_TYPE_TO_CANONICAL[mod_type] + "." * padding, padding) # Represent all possible motifs
    methylation_sequences_subset = methylation_sequences.copy()
    best_candidates = []
    continue_search = True
    dead_ends = 0
    motif_graph = None

    while continue_search:
        if dead_ends >= max_dead_ends:
            log.debug("Stopping search, too many low scoring candidates")
            break
        # Find the initial guess within the tree
        motif_graph, naive_guess = a_star_search(
            root_motif, 
            sequence, 
            pileup, 
            methylation_sequences_subset, 
            motif_graph = motif_graph,
            min_kl = min_kl,
            max_rounds_since_new_best = 30
        )

        # If there is no naive guess, stop the search
        if naive_guess == root_motif:
            log.info("No naive guess found, stopping search")
            break
        guess = naive_guess
        # next_guess = naive_guess
        # while guess != next_guess:
        #     # Find the best guess within the subtree of the naive guess
        #     guess = next_guess
        #     motif_graph, next_guess = a_star_search(
        #         guess, 
        #         sequence, 
        #         pileup, 
        #         methylation_sequences_subset, 
        #         motif_graph = motif_graph,
        #         min_kl = min_kl,
        #         max_rounds_since_new_best = 5
        #     )

        # Remove new candidate from methylation sequences
        seq_before = methylation_sequences_subset.shape[0]
        methylation_sequences_subset = methylation_sequences_subset.filter_sequence_matches(naive_guess.one_hot(), keep_matches = False)
        if methylation_sequences_subset is None:
            log.info("No more sequences left")
            break

        # Check if we should continue the search
        seq_remaining = methylation_sequences_subset.shape[0]
        seq_remaining_percent = seq_remaining/total_sequences
        log.info(f"{naive_guess}, {seq_before-seq_remaining} seqs. model: {motif_graph.nodes[naive_guess]['model']}. ({100*seq_remaining_percent:.1f} % left)")

        if motif_graph.nodes[naive_guess]["model"].mean() < 0.3:
            dead_ends += 1
            log.info(f"Candidate has low score. {dead_ends} of {max_dead_ends} before temination")
            continue

        best_candidates.append(naive_guess)
        
        if (seq_remaining/total_sequences) < 0.01:
            log.info("Stopping search, too few sequences left")
            break
        log.info("Continuing search")
    return motif_graph, best_candidates


 
def motif_child_nodes_kl_dist_max(motif, meth_pssm, contig_pssm, freq_threshold=0.25, min_kl=0.1):
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

def a_star_search(root_motif, contig, pileup, methylation_sequences, 
                  motif_graph = None, min_kl = 0.1, max_rounds_since_new_best = 10):
    """
    A* search for methylation motifs

    Parameters:
    - root_motif (list): The root motif to start the search from.
    - contig (str): The contig to be processed.
    - pileup (Pileup): The pileup to be processed.
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
    padding = methylation_sequences.shape[1] // 2
    total_sequences = methylation_sequences.shape[0]
    best_guess = root_motif
    root_model = motif_model_contig(pileup, contig.sequence, root_motif)
    best_score = scoring_function(root_model, root_model)
    rounds_since_new_best = 0
    visisted_nodes = []

    # Sample sequence in contig to get background for KL-divergence
    contig_sequences = contig.sample_n_subsequences(padding*2 + 1, 10000)
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
        if len(current.strip()) > 18:
            log.debug(f"{current}, Skipping due to length")
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
            contig_pssm
        ))

        # Add neighbors to graph
        for next in neighbors:
            if next in motif_graph.nodes:
                # Add only edge if motif -> next
                next_model = motif_graph.nodes[next]["model"]
                score = motif_graph.nodes[next]["score"]
                motif_graph.add_edge(current, next)
            else:
                next_model = motif_model_contig(pileup, contig.sequence, next)

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



if __name__ == "__main__":
    from nanomotif.dataload import load_assembly, load_pileup
    assembly = load_assembly("data/ecoli/assembly.polished.fasta")
    ecoli = load_pileup("data/ecoli/modkit.pileup.bed")
    result = process_sample(assembly, ecoli.pileup, min_read_methylation_fraction = 0.80)
