"""
Utility functions
"""
from itertools import product
import numpy as np
import re
import numpy as  np
import random

np.random.seed(1)

def generate_kmers(k: int, alphabet = ["A", "T", "C", "G"]):
    """
    Generate all possible k-mers of a given alphabet of sizes 1 to k
    """
    return [''.join(p) for i in range(1, k + 1) for p in product(alphabet, repeat=i)]
    


def subseq_indices(subseq, seq):
    """
    Find occourance indeces of a subseq in a seq
    """
    index = [match.start() for match in re.finditer(subseq, seq)]
    return np.array(index)

def calculate_match_length(regex):
    # Remove non-capturing groups for simplicity
    regex = re.sub(r'\[.*\]', '.', regex)

    return len(regex)