import numpy as np
import polars as pl
from .model import BetaBernoulliModel
from .utils import subseq_indices 
from .seq import reverse_complement

def methylated_motif_occourances(motif: str, motif_meth_index: int, seq: str, contig_meth_positions) -> tuple:
    """
    Get occourances of a motif in a contig
    """
    motif_index = subseq_indices(motif, seq) + motif_meth_index
    meth_occurences = np.intersect1d(contig_meth_positions, motif_index)
    nonmeth_occurences =  np.setdiff1d(motif_index, contig_meth_positions)
    return meth_occurences, nonmeth_occurences

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
        index_meth_fwd, index_nonmeth_fwd = methylated_motif_occourances(motif, motif_meth_index, contig, meth_positions_fwd)
        index_meth_rev, index_nonmeth_rev = methylated_motif_occourances(reverse_complement(motif), len(motif) - 1 - motif_meth_index, contig, meth_positions_rev)

        model = BetaBernoulliModel()
        model.update(len(index_meth_fwd) + len(index_meth_rev), len(index_nonmeth_fwd) + len(index_nonmeth_rev))

        score_canidates["motif"].append(motif)
        score_canidates["motif_meth_index"].append(motif_meth_index)
        score_canidates["mod_type"].append(mod_type)
        score_canidates["posterior"].append(model)
    return pl.DataFrame(score_canidates)


if __name__ == "__main__":
    from .candidate import generate_random_candidates
    assembly = nm.load_assembly("../data/ecoli/assembly.polished.fasta")
    ecoli = nm.load_pileup("../data/ecoli/modkit.pileup.bed")
    ecoli_6mA_80p = ecoli.pileup.filter(pl.col("mod_type") == "a").filter(pl.col("fraction_mod") > 0.8)
    motif_candidates = generate_random_candidates(4, "a")
    scored_candidates = score_candidates(ecoli_6mA_80p, assembly["contig_3"], motif_candidates)