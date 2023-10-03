import numpy as np
import polars as pl
import logging as log
from itertools import takewhile
from scipy.stats import entropy
from nanomotif.model import BetaBernoulliModel
from nanomotif.utils import subseq_indices, calculate_match_length
from nanomotif.seq import reverse_complement, EqualLengthDNASet, DNAsequence
import heapq as hq
import networkx as nx
log.basicConfig(encoding='utf-8', level=log.DEBUG, format='%(levelname)s: %(message)s')

##########################################
# Motif candidate scoring
##########################################

def methylated_motif_occourances(motif: str, motif_meth_index: int, seq: str, contig_meth_positions: int) -> tuple:
    """
    Get occourances of a motif in a contig

    Parameters:
    - motif (str): The motif to search for.
    - motif_meth_index (int): The index of the methylation site in the motif.
    - seq (str): The contig to search in.
    - contig_meth_positions (list): The positions of the methylation sites in the contig.

    Returns:
    - tuple: A tuple of two numpy arrays, the first containing the positions of the methylated motifs, the second containing the positions of the non-methylated motifs.
    """
    # Remove leading and trailing dots from motif
    motif_stripped = motif.lstrip(".").rstrip(".")

    # Get the index of the methylation site in the motif, after removing variable positions e.g. '.'
    motif_meth_index = motif_meth_index - len(list(takewhile(lambda x: x == '.', motif)))

    # Get the index of the methylation in the motif in the contig
    motif_index = subseq_indices(motif_stripped, seq) + motif_meth_index

    # Methylated motif positions
    meth_occurences = np.intersect1d(contig_meth_positions, motif_index)

    # Non-methylated motif positions
    nonmeth_occurences =  np.setdiff1d(motif_index, contig_meth_positions)

    return meth_occurences, nonmeth_occurences

def score_candidate(pileup, contig: str, motif: str, motif_meth_index: int):
    """
    Get the posterior for a single candidate

    Parameters:
    - pileup (Pileup): The pileup to be processed.
    - contig (str): The contig to be processed.
    - motif (str): The motif to be processed.
    - motif_meth_index (int): The index of the methylation site in the motif.

    Returns:
    - BetaBernoulliModel: The model for the candidate.
    """
    modtypes = pileup["mod_type"].unique()
    meth_positions_fwd = pileup.filter(pl.col("strand") == "+")["position"].to_numpy()
    meth_positions_rev = pileup.filter(pl.col("strand") == "-")["position"].to_numpy()

    motif_length = calculate_match_length(motif)
    index_meth_fwd, index_nonmeth_fwd = methylated_motif_occourances(motif, motif_meth_index, contig, meth_positions_fwd)
    index_meth_rev, index_nonmeth_rev = methylated_motif_occourances(reverse_complement(motif), motif_length - 1 - motif_meth_index, contig, meth_positions_rev)

    model = BetaBernoulliModel()
    model.update(len(index_meth_fwd) + len(index_meth_rev), len(index_nonmeth_fwd) + len(index_nonmeth_rev))

    return model



##########################################
# Motif candidate state space searcg
##########################################


def identify_motifs(methylation_sequences, contig, pileup, 
                    min_sequences = 50, min_cdf_score = 0.5, cdf_limit = 0.55, 
                    n_contig_samples = 10000, canonical_base = "A", min_kl_divergence = 0.1):
    """
    Identify candidate methylation motifs in a contig

    Parameters:
    - methylation_sequences (EqualLengthDNASet): The methylation sequences to be processed.
    - contig (str): The contig to be processed.
    - pileup (Pileup): The pileup to be processed.
    - min_sequences (int): The minimum number of methylation sequences.
    - min_cdf_score (float): The minimum score of 1 - cdf(cdf_limit) for a motif to be considered valid.
    - cdf_limit (float): The position to evaluate the cdf at.
    - n_contig_samples (int): The number of samples to take from the contig.
    - canonical_base (str): The canonical base of methylation (6mA = A).
    - min_kl_divergence (float): Early stopping criteria, if max KL-divergence falls below, stops building motif.

    Returns:
    - tuple: A tuple of three lists, the first containing the candidate motifs, the second containing the models for the candidates, the third containing the scores for the candidates.
    """
    # Infer padding from sequence length
    padding = methylation_sequences.sequence_length // 2

    # Sample sequence in contig to get background for KL-divergence
    contig_sequences = contig.sample_n_at_subsequence(padding, canonical_base, n_contig_samples)
    contig_pssm = contig_sequences.pssm()

    # Initialize candidate lists
    candidate_motifs = []
    failed_canidates = []
    candidate_models = []
    candidate_scores = []

    # Bases in the DNA sequences (To make sure the order is the same as used internally in the class)
    bases = methylation_sequences.sequences[0].bases

    # Convert DNA sequences to numpy array
    methylation_sequences = convert_seqeunces_to_numpy(methylation_sequences.sequences)

    while True:
        active_methylation_sequences = methylation_sequences.copy()

        if len(candidate_motifs) > 0:
            # Remove all previously identified candidates from sequence set
            for motif in candidate_motifs:
                active_methylation_sequences = subset_DNA_array(active_methylation_sequences, motif, remove_matches = True)
        if len(failed_canidates) > 0:
            # Remove all previously failed candidates from sequence set
            for motif in failed_canidates:
                active_methylation_sequences = subset_DNA_array(active_methylation_sequences, motif, remove_matches = True)

        if active_methylation_sequences.shape[0] < min_sequences:
            break
        
        # Initialize the active candidate as blank sequence
        active_candidate = ["."] * padding + [canonical_base] + ["."] * padding
        evaluated_positions = []
        candidates = []
        best_candidate = None

        while True:
            # Calculate KL-divergence
            methylation_pssm = calculate_pssm(active_methylation_sequences)
            kl_divergence = entropy(methylation_pssm, contig_pssm)

            # Set KL-divergence to 0 for already evaluated positions
            kl_divergence[evaluated_positions] = 0

            # Get the index of the maximum KL-divergence
            max_kl_index = kl_divergence.argmax()

            # Update the active candidate
            new_base = [bases[i] for i in np.argwhere(methylation_pssm[:, max_kl_index] > 0.45)[:, 0]]
            new_base = "".join(new_base)

            # If multiple bases have a high probability, use a regex expression to represent: [bases]
            if len(new_base) == 0:
                new_base = bases[methylation_pssm[:, max_kl_index].argmax()]
            if len(new_base) > 1:
                new_base = "[" + new_base + "]"
            active_candidate[max_kl_index] = new_base

            model = score_candidate(pileup, contig.sequence, "".join(active_candidate), padding)
            score = 1 - model.cdf(cdf_limit)

            evaluated_positions.append(max_kl_index)
            log.debug(f"{''.join(active_candidate)} | cdf score: {score:.3f} | mean: {model.mean():.3f} | {model} | max kl: {kl_divergence[max_kl_index]:.3f}")
            
            active_methylation_sequences = subset_DNA_array(active_methylation_sequences, active_candidate, remove_matches = False)

            candidates.append((active_candidate, model, score))

            # Success criteria
            if score >= min_cdf_score:
                failed = False
                best_candidate = len(candidates) - 1
                break

            # Failure criteria
            if active_methylation_sequences.shape[0] < min_sequences:
                failed = True
                log.debug("Too few sequences left")
                break
            if len(evaluated_positions) >= padding * 2 + 1:
                failed = True
                log.debug("Too many positions evaluated")
                break
            if kl_divergence[max_kl_index] < min_kl_divergence:
                failed = True
                log.debug("Low KL divergence")
                break

        if failed:
            failed_canidates.append(active_candidate)
        else:
            log.debug("Saving candidate")
            candidate_motifs.append(candidates[best_candidate][0])
            candidate_models.append(candidates[best_candidate][1])
            candidate_scores.append(candidates[best_candidate][2])
    candidate_motifs_str = ["".join(motif) for motif in candidate_motifs]
    return candidate_motifs_str, candidate_models, candidate_scores


def process_sample(assembly, pileup, 
                   max_candidate_size = 40,
                   min_read_methylation_fraction = 0.8,
                   min_valid_coverage = 5,
                   min_kl_divergence = 0.1,
                   min_cdf_score = 0.8,
                   cdf_position = 0.55,
                   min_motif_frequency = 20000
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
    mod_dict = {"a":"A", "m":"C", "h":"C", "c":"C"}
    pileup = pileup \
            .filter(pl.col("Nvalid_cov") > min_valid_coverage) \
            .filter(pl.col("fraction_mod") > min_read_methylation_fraction) \
            .sort("contig") \
            .groupby(["contig", "mod_type"])
    for (contig, modtype), subpileup in pileup:
        log.info(f"Processing {contig}")
        log.info(f"Processing {modtype}")

        sequence = assembly.assembly[contig]

        methylation_index_fwd = subpileup \
            .filter(pl.col("strand") == "+")  \
            .get_column("position").to_list()
        if len(methylation_index_fwd) <= 1:
            continue
        methylation_index_rev = subpileup \
            .filter(pl.col("strand") == "-") \
            .get_column("position").to_list()
        if len(methylation_index_rev) <= 1:
            continue

        methylation_sequences_fwd = EqualLengthDNASet(
            [DNAsequence(sequence[(i - padding):(i + padding + 1)]) for i in methylation_index_fwd if (i > padding) and (i < (len(sequence) - padding))]
        )
        methylation_sequences_rev = EqualLengthDNASet(
            [DNAsequence(sequence[(i - padding):(i + padding + 1)]).reverse_complement() for i in methylation_index_rev if (i > padding) and (i < (len(sequence) - padding))]
        )
        methylation_sequences = methylation_sequences_fwd + methylation_sequences_rev
        
        min_sequences = max(min(len(sequence) // min_motif_frequency, 200), 20)

        identified_motifs = pl.DataFrame(identify_motifs(
            methylation_sequences, 
            sequence, 
            subpileup, 
            min_kl_divergence = min_kl_divergence, 
            min_sequences = min_sequences, 
            cdf_limit = cdf_position,
            min_cdf_score = min_cdf_score,
            canonical_base = mod_dict[modtype]
        ))
        if len(identified_motifs) == 0:
            continue
        else:
            result.append(identified_motifs.with_columns(
                pl.lit(contig).alias("contig"),
                pl.lit(modtype).alias("mod_type")
            ))
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
    if len(result) == 0:
        return pl.DataFrame({
            "motif": [],
            "mod_position": [],
            "mod_type": [],
            "cdf_score": []
        })
    motifs = pl.concat(result) \
        .rename({"column_0": "padded_motif", "column_1": "model", "column_2":"cdf_score"}) \
        .with_columns([
        pl.col("padded_motif").apply(lambda motif: motif[count_periods_at_start(motif):-max(count_periods_at_end(motif), 1)]).alias("motif"),
            pl.col("padded_motif").apply(lambda motif: padding - count_periods_at_start(motif)).alias("mod_position")
        ])
    return motifs



def convert_seqeunces_to_numpy(sequences):
    # Mapping the characters to indices for easy lookup
    # Index [sequence, position, base]
    base_to_int = {'A': [1,0,0,0], 'T': [0,1,0,0], 'G': [0,0,1,0], 'C': [0,0,0,1]}
    DNAarray = np.array([[base_to_int[base] for base in sequence] for sequence in sequences])
    return DNAarray

def calculate_pssm(DNAarray):
    return DNAarray.sum(axis = 0).transpose() / DNAarray.shape[0]

def convert_nucletotide_string_to_numpy(sequence):
    base_to_int = {'A': [1,0,0,0], 'T': [0,1,0,0], 'G': [0,0,1,0], 'C': [0,0,0,1], ".": [1,1,1,1]}
    DNAarray = np.array([[base_to_int[base] for base in sequence]])
    return DNAarray

def convert_regex_to_numpy(regex):
    base_to_int = {'A': [1,0,0,0], 'T': [0,1,0,0], 'G': [0,0,1,0], 'C': [0,0,0,1], ".": [1,1,1,1]}
    DNAarray = np.zeros((1, len(regex), 4), dtype = np.int32)
    for base, arr in base_to_int.items():
        for i, seq_base in enumerate(regex):
            if base in seq_base:
                DNAarray[0, i, :] += arr
    return DNAarray

def subset_DNA_array(DNAarray, sequence: list, remove_matches = True):
    # Convert the query sequence into a numpy array, with NaN wherever there's a dot
    pattern = convert_regex_to_numpy(sequence)

    if remove_matches:
        return DNAarray[np.any(DNAarray > pattern, axis=(1,2)), :]
    else:
        return DNAarray[np.all(DNAarray <= pattern, axis=(1,2)), :]



def find_motif_neighbors_kl_distance(motif, meth_pssm, contig_pssm, min_kl=0.05):
    # Calculate KL-divergence
    kl_divergence = entropy(meth_pssm, contig_pssm)
    bases = ["A", "T", "G", "C"]
    # Update the active candidate
    for pos in np.where(kl_divergence > min_kl)[0]:
        for base in [bases[i] for i in np.argwhere(np.logical_and(meth_pssm[:, pos] > contig_pssm[:, pos], meth_pssm[:, pos] > 0.35))[:, 0]]:
            if motif[pos] == ".":
                yield motif[:pos] + base + motif[pos+1:]

    
def a_star_search(root_motif, contig, pileup, methylation_sequences, min_score = 0.6, extra_rounds = 2, min_kl = 0.1):
    # Search setup
    count_down = extra_rounds
    initiate_stopping = False
    padding = methylation_sequences.sequence_length // 2
    best_score = 0

    # Sample sequence in contig to get background for KL-divergence
    contig_sequences = contig.sample_n_at_subsequence(padding, "G", 10000)
    contig_pssm = contig_sequences.pssm()

    # Convert DNA sequences to numpy array
    methylation_sequences = methylation_sequences.convert_to_DNAarray()

    # Initialize the root motif
    root_model = score_candidate(pileup, contig.sequence, root_motif, padding)
    motif_graph = nx.DiGraph()
    motif_graph.add_node(root_motif, model=root_model)

    # Initialize priority que
    priority_que = []
    hq.heappush(priority_que, (0, root_motif))
    
    # Search loop
    while len(priority_que) > 0:
        # Get the next candidate
        current = hq.heappop(priority_que)
        current = current[1]

        current_model = motif_graph.nodes[current]["model"]
        log.debug(f"{''.join(current)} | {current_model} ")

        # Prune the neibhors of the current candidate
        active_methylation_sequences = methylation_sequences.copy().filter_sequence_matches(current, keep_matches = True)
        meth_pssm = active_methylation_sequences.pssm()
        neighbors = list(find_motif_neighbors_kl_distance(
            current, 
            active_methylation_sequences.pssm(), 
            contig_pssm, 
            min_kl=min_kl
        ))

        # Add neighbors to graph
        for next in neighbors:
            # Skip if already in graph, only add new edge
            if next in motif_graph.nodes:
                motif_graph.add_edge(current, next, d_alpha=motif_graph.nodes[next]["model"]._alpha - current_model._alpha, d_beta=motif_graph.nodes[next]["model"]._beta - current_model._beta)
                continue

            next_model = score_candidate(pileup, contig.sequence, next, padding)

            # Add neighbor to graph
            motif_graph.add_node(next, model=next_model)
            motif_graph.add_edge(current, next)

            # Add neighbor to priority que
            distance = 1 - (next_model._alpha / root_model._alpha)
            heuristic = 1 - (next_model._beta / root_model._beta)
            priority =  distance + heuristic
            hq.heappush(priority_que, (priority, next))

            score = 1 - next_model.cdf(0.55)
            if score > best_score:
                best_score = score

        # Early stopping criteria
        if best_score > min_score:
            initiate_stopping = True
        if initiate_stopping:
            count_down -= 1
            if count_down == 0:
                break
            
    
    return motif_graph


            

if __name__ == "__main__":
    from nanomotif.dataload import load_assembly, load_pileup
    assembly = load_assembly("../data/ecoli/assembly.polished.fasta")
    ecoli = load_pileup("../data/ecoli/modkit.pileup.bed")
    result = process_sample(assembly, ecoli.pileup, min_cdf_score = 0.6, min_read_methylation_fraction = 0.85)
