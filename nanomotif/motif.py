import numpy as np
import re
import itertools
import networkx as nx
import logging as log
from scipy.stats import entropy
from nanomotif.constants import *
from nanomotif.seq import regex_to_iupac, iupac_to_regex, reverse_compliment
class Motif(str):
    def __new__(cls, motif_string, *args, **kwargs):
        return str.__new__(cls, motif_string)
    
    def __init__(self, _, mod_position):
        self.mod_position = mod_position
        self.string = self.__str__()
    def from_iupac(self):
        """
        Create a motif from an IUPAC string.
        """
        regex = iupac_to_regex(self.string)
        return Motif(regex, self.mod_position)
    def iupac(self):
        """
        Get IUPAC sequence of motif
        """
        return regex_to_iupac(self.string)
    def identical(self, other_motif):
        """
        Check if two motifs are identical.
        """
        assert isinstance(other_motif, Motif)
        str_equal = self.string == other_motif.string
        pos_equal = self.mod_position == other_motif.mod_position
        return str_equal and pos_equal
    def sub_motif_of(self, other_motif):
        """
        Check if a motif is in another motif.
        """
        assert isinstance(other_motif, Motif)
        # Strip for redundant posititons (flanking .)
        if self.string == other_motif.string:
            return False
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

    def sub_string_of(self, other_motif):
        """
        Check if a motif is in another motif.
        """
        assert isinstance(other_motif, Motif)
        # Strip for redundant posititons (flanking .)
        self_stripped = self.new_stripped_motif()
        other_stripped = other_motif.new_stripped_motif()

        if self_stripped.length() < other_stripped.length():
            return False
        if self_stripped.string == other_stripped.string:
            return False
        size_difference = self_stripped.length() - other_stripped.length()

        # Split into list of bases
        self_split = self_stripped.split()
        other_split = other_stripped.split()
        
        for i in range (size_difference+1):
            match = True
            for j, base in enumerate(other_split):
                if not set(self_split[j+i]).issubset(set(base)) and base != ".": 
                    # Using set to make sure nucleotide order doesn't matter
                    match = False
            if match:
                return True
        return False 
    
    def distance(self, other_motif):
        assert isinstance(other_motif, Motif)
        self_split = self.split()
        other_motif_split = other_motif.split()

        index_offset = self.mod_position - other_motif.mod_position

        start = (-self.mod_position, -other_motif.mod_position)
        end = (len(self_split) - self.mod_position, len(other_motif_split) - other_motif.mod_position)
        search_size = range(min(start), max(end))
    
        
        distance = 0
        for i in search_size:
            if i < start[0]:
                if other_motif_split[i - start[1]] != ".":
                    distance += 1
            elif i < start[1]:
                if self_split[i - start[0]] != ".":
                    distance += 1
            elif i >= end[0]:
                if other_motif_split[i - start[1]] != ".":
                    distance += 1
            elif i >= end[1]:
                if self_split[i - start[0]] != ".":
                    distance += 1
            else:
                if set(self_split[i - start[0]]) == set(other_motif_split[i - start[1]]):
                    continue
                else:
                    distance += 1
        return distance

    def have_isolated_bases(self, isolation_size=2):
        """
        Check if a motif has isolated bases.
        """
        motif_split = self.split()
        isolated = False
        for pos in range(len(motif_split)):
            if motif_split[pos] == ".":
                continue
            index_start = max(pos - isolation_size, 0)
            index_end = min(pos + isolation_size + 1, len(motif_split) - 1)
            # If all surrounding positions are ".", it is isolated
            if set(motif_split[index_start:pos] + motif_split[pos+1:index_end]) == set(["."]):
                isolated = True
            if set(motif_split[index_start:pos] + motif_split[pos+1:index_end]) == set(["N"]):
                isolated = True
        return isolated
        

    def length(self) -> int:
        """
        Return the length of the motif.
        """
        return len(self.split())
    def trimmed_length(self) -> int:
        """
        Return the length of the motif without dots.
        """
        return len(self.split()) - self.count(".")
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
                if j == -1:
                    raise ValueError("Unmatched bracket")
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
    
    def merge(self, other):
        assert isinstance(other, Motif)
        self_trim = self.new_stripped_motif()
        other_trim = other.new_stripped_motif()
        self_seq = self_trim.split()
        other_seq = other_trim.split()

        offset = self_trim.mod_position - other_trim.mod_position
        if offset > 0:
            self_seq = self_seq[offset:]
        elif offset < 0:
            other_seq = other_seq[-offset:]

        final_length = min(len(self_seq), len(other_seq))
        self_seq = self_seq[:final_length]
        other_seq = other_seq[:final_length]
        
        # Merge the motifs base by base, keeping the modification position aligned
        merged_motif_string = ''.join([self_trim.merge_bases(base, base_other) for base, base_other in zip(self_seq, other_seq)])
        merged_mod_position = min(self_trim.mod_position, other_trim.mod_position)

        # Create and return a new Motif object with the merged motif string and modification position
        merged_motif = Motif(merged_motif_string, merged_mod_position)
        merged_motif = merged_motif.new_stripped_motif(".")
        return merged_motif
        
    def merge_bases(self, base1, base2):
        canonical_bases = {'A', 'C', 'G', 'T'}
        # Extract bases from brackets if present
        base1_set = set(base1[1:-1]) if '[' in base1 else set(base1)
        base2_set = set(base2[1:-1]) if '[' in base2 else set(base2)

        # Merge the sets of bases
        merged_set = base1_set.union(base2_set)

        # If the merged set has more than one base, enclose them in brackets
        if merged_set == canonical_bases:
            return '.'
        elif "." in merged_set:
            return '.'
        elif len(merged_set.intersection(canonical_bases)) > 1:
            return f"[{''.join(sorted(merged_set.intersection(canonical_bases)))}]"
        else:
            return ''.join(merged_set)  # return as a string without brackets if only one base
    
    def explode_motif(self):
        """
        Explode a motif into all possible combinations of bases.
        """
        exploded_strings = explode_sequence(self.string)
        return [Motif(motif, self.mod_position) for motif in exploded_strings]


def explode_sequence(seq):
    # Splitting the sequence into parts
    parts = []
    temp = ''
    in_brackets = False

    for char in seq:
        if char == '[':
            if temp:
                parts.append([temp])
                temp = ''
            in_brackets = True
            temp += char
        elif char == ']':
            temp += char
            parts.append([temp])
            temp = ''
            in_brackets = False
        else:
            if in_brackets:
                temp += char
            else:
                parts.append([char])

    # Handling the case where the sequence ends with a non-bracketed part
    if temp:
        parts.append([temp])

    # Processing bracketed parts to extract options
    for i, part in enumerate(parts):
        if part[0].startswith('['):
            # Extracting characters between brackets and splitting them into separate options
            parts[i] = list(part[0][1:-1])

    # Generating all combinations
    exploded_sequences = [''.join(combination) for combination in itertools.product(*parts)]

    return exploded_sequences

def check_all_nodes_connected(graph):
        n = len(graph.nodes()) 
        expected_edges = n * (n - 1) / 2
        is_fully_connected = len(graph.edges()) == expected_edges
        return is_fully_connected

class MotifDistanceGraph(nx.Graph):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
    def add_motifs(self, motifs, min_distance = 2):
        for motif in motifs:
            self.add_motif(motif, min_distance)
    def add_motif(self, motif, min_distance = 2):
        self.add_node(motif)
        for other_motif in self.nodes():
            if not motif.identical(other_motif):
                distance = motif.distance(other_motif)
                if distance <= min_distance:
                    self.add_edge(motif, other_motif, dist=distance)
    def get_connected_clusters(self):
        return list(nx.connected_components(self))
    
    def get_fully_connected_clusters(self):
        clusters = self.get_connected_clusters()
        return [cluster for cluster in clusters if check_all_nodes_connected(self.subgraph(cluster).copy())]


def merge_motifs(motifs, connectivity_dist=2, min_length=4):
    """
    Merges a Motifs that are closely related. 
    
    Args:
    - motifs (list): List of Motifs to merge.
    - connectivity_dist (int): Maximum edit distance between motifs to be merged.
    - min_length (int): Minimum length of motifs to be merged. 4 Based on estimate from REBASE motifs

    Returns:
    - list: List of merged motifs.
    """

    # Subset to motif greater than minimum length
    motifs = [motif for motif in motifs if motif.trimmed_length() > min_length]
    distance_graph = MotifDistanceGraph()
    distance_graph.add_motifs(motifs, connectivity_dist)
    clusters = distance_graph.get_fully_connected_clusters()


    # Now 'clusters' is a list of sets, where each set is a cluster of nodes
    motif_merged = {}

    for i, cluster in enumerate(clusters, start=1):
        if len(cluster) == 1:
            continue
        log.debug(f'Cluster {i}: {cluster}')
        cluster_motifs = []
        merged_motif = None
        for node in cluster:
            cluster_motifs.append(node)
            # Merge motifs
            if merged_motif is None:
                merged_motif = node
            else:
                merged_motif = merged_motif.merge(node)


        log.debug(f'Merged motif: {merged_motif}')
        motif_merged[i] = [merged_motif, cluster_motifs]
    return motif_merged


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


if __name__ == "__main__":
    motif1 = Motif("AT.G", 0)
    motif2 = Motif("ATCG", 0)
    motif1.sub_string_of(motif2)        