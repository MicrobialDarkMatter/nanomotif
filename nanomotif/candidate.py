"""
Motif candidate funcitonalities
"""
from .utils import generate_kmers

class MotifCandidates():
    """

    """
    def __init__(self, motif: list, mod_position: list, mod_type: list):
        self.motif = motif
        self.mod_position = mod_position
        self.mod_type = mod_type

    def __iter__(self):
        return iter(zip(self.motif, self.mod_position, self.mod_type))

    def __len__(self):
        return len(self.motif)
    
    def __repr__(self):
        return f"MotifCandidates with {len(self.motif)} candidates"
    
    def add_candidates(self, motif: list, mod_position: list, mod_type: list):
        self.motif += motif
        self.mod_position += mod_position
        self.mod_type += mod_type


def generate_random_candidates(k, modtype):
    """
    Generate a list of candidate kmers and their respective indices for a given k and modification type.

    Parameters:
    - k (int): Max length of the kmer sequences to be generated.
    - modtype (str): The type of modification to be applied to the kmers (a or m).

    Returns:
    - MotifCandidates: An object that encapsulates the generated candidate kmers, their indices, 
      and a list with the same modification type repeated for each kmer.

    Notes:
    - The indices for each kmer represent positions in the kmer string. Therefore, for a kmer of length 'k', 
      there will be 'k' indices, resulting in each kmer being repeated 'k' times in the candidate_kmers list.
    """
    canonical = modtype_canonical = {"a": "A", "m": "C"}
    kmers = generate_kmers(k, alphabet=["A", "C", "G", "T"])

    candidate_index = [i for motif in kmers for i in range(len(motif))]
    candidate_kmers = [motif for motif in kmers for i in range(len(motif))]

    return MotifCandidates(candidate_kmers, candidate_index, [modtype]*len(candidate_kmers))



