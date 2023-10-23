import numpy as np
import re
import itertools
import networkx as nx
from scipy.stats import entropy
from nanomotif.constants import *

class Motif(str):
    def __new__(cls, motif_string, *args, **kwargs):
        return str.__new__(cls, motif_string)
    
    def __init__(self, _, mod_position):
        self.mod_position = mod_position
        self.string = self.__str__()

    def child_of(self, other_motif):
        """
        Check if a motif is in another motif.
        """
        assert isinstance(other_motif, Motif)
        # Strip for redundant posititons (flanking .)
        self_stripped = self.new_stripped_motif()
        other_stripped = other_motif.new_stripped_motif()

        if self_stripped.length() < other_stripped.length():
            return False

        # Split into list of bases
        self_split = self_stripped.split()
        other_split = other_stripped.split()

        index_offset = other_stripped.mod_position - self_stripped.mod_position
        if index_offset > 0:
            return False
        for i, base in enumerate(self_split):
            if i + index_offset < 0:
                continue
            try:
                if set(base).intersection(set(other_split[i + index_offset])) != set(base) and other_split[i + index_offset] != ".": 
                    # Using set to make sure nucleotide order doesn't matter
                    return False
            except IndexError:
                return True
        return True

    def length(self) -> int:
        """
        Return the length of the motif.
        """
        return len(self.split())
    
    def strip(self, character = ".") -> str:
        """
        Return the motif string with single characters surrounded by dots removed.
        """
        new_motif = self
        return self.lstrip(character).rstrip(character)
    
    def new_stripped_motif(self, character = "."):
        """
        Return the motif string with single characters surrounded by dots removed.
        """
        if re.search("[^\.]", self.string) is None:
            return self
        
        new_position = self.mod_position - re.search("[^\.]", self.string).start()
        new_motif = self.strip(character)

        return Motif(new_motif, new_position)

    def split(self) -> list:
        """
        Split a string into a list of characters, keeping bracketed characters together.

        Returns:
        - list: A list of bases, where bracketed characters are kept together.
        """
        result = []
        i = 0
        while i < len(self):
            if self[i] == "[":
                j = self.find("]", i)
                result.append(self[i:j+1])
                i = j + 1
            else:
                result.append(self[i])
                i += 1
        return result
    
    def one_hot(self) -> np.ndarray:
        """
        Convert the motif to a numpy array.
        """
        split_motif = self.split()
        squence_array = np.zeros((len(split_motif), 4), dtype=int)
        for position, bases in enumerate(split_motif):

            for base in bases:
                if base in BASE_TO_ONE_HOT.keys():
                    squence_array[position, :] += BASE_TO_ONE_HOT[base]
        return squence_array
    
    def reverse_compliment(self):
        """
        Reverse compliment a motif.
        """
        reversed_motif = "".join([COMPLEMENT[base] for base in reversed(self.string)])
        new_position = self.length() - self.mod_position - 1
        return Motif(reversed_motif, new_position)

def remove_child_motifs(motifs):
    """
    Remove motifs that are children of other motifs in the list.
    """
    new_motifs = []
    for i, motif in enumerate(motifs):
        is_child = False
        for j, other_motif in enumerate(motifs):
            if motif.child_of(other_motif) and i != j:
                is_child = True
        if not is_child:
            new_motifs.append(motif)
    return new_motifs

class MotifTree(nx.DiGraph):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    