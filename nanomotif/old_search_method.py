

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

