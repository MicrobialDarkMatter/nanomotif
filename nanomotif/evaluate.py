import numpy as np
import polars as pl
import logging as log
from scipy.stats import entropy
from nanomotif.model import BetaBernoulliModel
from nanomotif.utils import subseq_indices, calculate_match_length
from nanomotif.seq import reverse_complement, EqualLengthDNASet, DNAsequence
log.basicConfig(encoding='utf-8', level=log.DEBUG, format='%(levelname)s: %(message)s')

def methylated_motif_occourances(motif: str, motif_meth_index: int, seq: str, contig_meth_positions) -> tuple:
    """
    Get occourances of a motif in a contig
    """
    # Motif match index
    motif_index = subseq_indices(motif, seq) + motif_meth_index
    # Methylated motif positions
    meth_occurences = np.intersect1d(contig_meth_positions, motif_index)
    # Nonmethylated motid positions
    nonmeth_occurences =  np.setdiff1d(motif_index, contig_meth_positions)
    return meth_occurences, nonmeth_occurences

def score_candidate(pileup, contig: str, motif: str, motif_meth_index: int):
    """
    Score a single candidate
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

def score_candidates(pileup, contig: str, candidates):
    """
    Score candidates for a single contig
    """
    modtypes = pileup["mod_type"].unique()
    meth_positions_fwd = pileup.filter(pl.col("strand") == "+")["position"].to_numpy()
    meth_positions_rev = pileup.filter(pl.col("strand") == "-")["position"].to_numpy()

    score_canidates = {
        "motif": [],
        "motif_meth_index":[],
        "mod_type": [],
        "posterior": []
    }
    for motif, motif_meth_index, mod_type in candidates:
        motif_length = calculate_match_length(motif)
        index_meth_fwd, index_nonmeth_fwd = methylated_motif_occourances(motif, motif_meth_index, contig, meth_positions_fwd)
        index_meth_rev, index_nonmeth_rev = methylated_motif_occourances(reverse_complement(motif), motif_length - 1 - motif_meth_index, contig, meth_positions_rev)

        model = BetaBernoulliModel()
        model.update(len(index_meth_fwd) + len(index_meth_rev), len(index_nonmeth_fwd) + len(index_nonmeth_rev))

        score_canidates["motif"].append(motif)
        score_canidates["motif_meth_index"].append(motif_meth_index)
        score_canidates["mod_type"].append(mod_type)
        score_canidates["posterior"].append(model)
    return pl.DataFrame(score_canidates)

def identify_motifs(methylation_sequences, contig, pileup, 
                    min_sequences = 50, min_cdf_score = 0.5, cdf_limit = 0.55, 
                    n_contig_samples = 10000, canonical_base = "A", min_kl_divergence = 0.1):
    padding = methylation_sequences.sequence_length // 2
    contig_sequences = contig.sample_n_at_subsequence(padding, canonical_base, n_contig_samples)
    contig_pssm = contig_sequences.pssm()
    candidate_motifs = []
    failed_canidates = []
    candidate_models = []
    candidate_scores = []
    bases = methylation_sequences.sequences[0].bases

    # Convert DNA sequences to numpy array
    methylation_sequences = convert_seqeunce_to_numpy(methylation_sequences.sequences)

    while True:
        active_methylation_sequences = methylation_sequences.copy()
        if len(candidate_motifs) > 0:
            for motif in candidate_motifs:
                active_methylation_sequences = subset_DNA_array(active_methylation_sequences, motif, remove_matches = True)
        if len(failed_canidates) > 0:
            for motif in failed_canidates:
                active_methylation_sequences = subset_DNA_array(active_methylation_sequences, motif, remove_matches = True)
        if active_methylation_sequences.shape[0] < min_sequences:
            break

        active_candidate = "." * (padding * 2 + 1)
        evaluated_positions = []
        previous_score = 0
        candidates = []
        while True:
            sequence_pssm = calculate_pssm(active_methylation_sequences)
            kl_divergence = entropy(sequence_pssm, contig_pssm)
            kl_divergence[evaluated_positions] = -1


            new_base_index = kl_divergence.argmax()
            new_base = bases[sequence_pssm[:, new_base_index].argmax()]
            active_candidate = active_candidate[:new_base_index] + new_base + active_candidate[new_base_index + 1:]

            model = score_candidate(pileup, contig.sequence, active_candidate, padding)
            score = 1 - model.cdf(cdf_limit)

            # n = 10
            # score_weighted_cdf = 0
            # previous_cdf = 0
            # for i in range(n):
            #     current_cdf = model.cdf((i + 1)/ n) 
            #     score_weighted_cdf += ((i + 1)/ n) * (current_cdf - previous_cdf)
            #     previous_cdf = current_cdf

            evaluated_positions.append(new_base_index)
            log.debug(f"{active_candidate} | cdf score: {score:.3f} | mean: {model.mean():.3f} | n seqs: {active_methylation_sequences.shape[0]:8} | max kl: {kl_divergence[new_base_index]:.3f}")
            
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
            if kl_divergence[new_base_index] < min_kl_divergence:
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
    return candidate_motifs, candidate_models, candidate_scores


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
    methylation_types = pileup["mod_type"].unique()
    result = []
    for modtype in methylation_types:
        log.info(f"Processing {modtype}")
        for contig, sequence in assembly.assembly.items():
            log.info(f"Processing {contig}")
            pileup_meth = pileup \
                .filter(pl.col("contig") == contig) \
                .filter(pl.col("mod_type") == modtype) \
                .filter(pl.col("fraction_mod") > min_read_methylation_fraction) 

            methylation_index_fwd = pileup_meth \
                .filter(pl.col("strand") == "+") \
                .filter(pl.col("Nvalid_cov") > min_valid_coverage) \
                .get_column("position").to_list()
            if len(methylation_index_fwd) <= 1:
                continue
            methylation_index_rev = pileup_meth \
                .filter(pl.col("strand") == "-") \
                .filter(pl.col("Nvalid_cov") > min_valid_coverage) \
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
            
            min_sequences = max(min(len(sequence) // min_motif_frequency, 200), 5)

            identified_motifs = pl.DataFrame(identify_motifs(
                methylation_sequences, 
                sequence, 
                pileup_meth, 
                min_kl_divergence = min_kl_divergence, 
                min_sequences = min_sequences, 
                cdf_limit = cdf_position,
                min_cdf_score = min_cdf_score
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



def convert_seqeunce_to_numpy(sequences):
    # Mapping the characters to indices for easy lookup
    base_to_int = {'A': 0, 'T': 1, 'G': 2, 'C': 3}
    array = np.array([[base_to_int[base] for base in sequence] for sequence in sequences])
    return array

def calculate_pssm(DNAarray):
    # Get the number of sequences
    return np.array([
        (DNAarray == 0).sum(axis=0),
        (DNAarray == 1).sum(axis=0),
        (DNAarray == 2).sum(axis=0),
        (DNAarray == 3).sum(axis=0)
    ]) / DNAarray.shape[0]

def subset_DNA_array(DNAarray, sequence, remove_matches = True):
    # Convert the query sequence into a numpy array, with NaN wherever there's a dot
    pattern = np.array([{'A': 0, 'T': 1, 'G': 2, 'C': 3, '.': np.nan}[base] for base in sequence])
    mask = np.logical_not(np.isnan(pattern))
    if remove_matches:
        return DNAarray[np.any(DNAarray[:, mask] != pattern[mask], axis=1), :]
    else:
        return DNAarray[np.all(DNAarray[:, mask] == pattern[mask], axis=1), :]


if __name__ == "__main__":
    from nanomotif.dataload import load_assembly, load_pileup
    assembly = load_assembly("../data/ecoli/assembly.polished.fasta")
    ecoli = load_pileup("../data/ecoli/modkit.pileup.bed")
    result = process_sample(assembly, ecoli.pileup, min_cdf_score = 0.6, min_read_methylation_fraction = 0.85)
