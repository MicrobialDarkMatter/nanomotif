"""
Utility functions
"""
from itertools import product
import numpy as np
import re
import random
import nanomotif as nm
np.random.seed(1)

def has_n_character_stretches(sequence, n, character):
    """
    Check if the given sequence has a segment of three or more consecutive dots.

    Parameters:
    sequence (str): The sequence to be checked.

    Returns:
    bool: True if there is a segment of three or more consecutive dots, False otherwise.
    """
    count = 0
    previous_char = ""
    for char in sequence:
        if char == character:
            if previous_char ==character:
                previous_char = char
                continue
            else:
                previous_char = char
                count += 1
        else:
            previous_char = char
    if count >= n:
        return True
    else:
        return False
def motif_type(motif_str):
    if has_n_character_stretches(motif_str, 2, "N"):
        return "ambiguous"
    elif re.search(r"(N){3,}", motif_str):
        return "bipartite"
    elif nm.seq.reverse_compliment(motif_str) == motif_str:
        return "palindrome"
    else:
        return "non-palindrome"

def generate_kmers(k: int, alphabet = ["A", "T", "C", "G"]):
    """
    Generate all possible k-mers of a given alphabet of sizes 1 to k
    """
    return [''.join(p) for i in range(1, k + 1) for p in product(alphabet, repeat=i)]
    


def subseq_indices(subseq, seq):
    """
    Find occourance indeces of a subseq in a seq
    """
    compiled_subseq = re.compile(subseq)
    index = [match.start() for match in re.finditer(compiled_subseq, seq)]
    return np.array(index)

def calculate_match_length(regex):
    # Remove non-capturing groups for simplicity
    regex = re.sub(r'\[.*\]', '.', regex)

    return len(regex)


def all_lengths_equal(iterator):
    """
    Checks whether the lengths of all elements in an iterable are equal.

    The function will return True even if the iterable is empty. It requires that the elements of the iterable
    also be iterable, such as strings, lists, and tuples.

    Args:
        iterator (Iterable): An iterable object containing other iterable elements.

    Returns:
        bool: True if all lengths are equal or if iterable is empty, False otherwise.

    Examples:
    >>> all_lengths_equal(['abc', 'def', 'ghi'])
    True
    >>> all_lengths_equal([[1, 2, 3], [4, 5, 6]])
    True
    >>> all_lengths_equal([])
    True
    >>> all_lengths_equal(['abc', 'de'])
    False
    """
    iterator = iter(iterator)
    try:
        first = len(next(iterator))
    except StopIteration:
        return True
    return all(first == len(x) for x in iterator)
