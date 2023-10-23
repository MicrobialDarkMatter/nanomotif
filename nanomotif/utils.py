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