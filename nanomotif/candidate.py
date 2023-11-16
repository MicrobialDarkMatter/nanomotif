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

    def from_iupac(iupac_motif):
        """
        Generate Motif from IUPAC seqeunce
        """
        

    def sub_motif_of(self, other_motif):
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
            index_end = min(pos + isolation_size, len(motif_split) - 1)
            # If all surrounding positions are ".", it is isolated
            if set(motif_split[index_start:pos] + motif_split[pos+1:index_end]) == set(["."]):
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
        mod_pos = self.mod_position
        mod_pos_other = other.mod_position

        # Determine the offset required to align the modification positions
        offset =  mod_pos - mod_pos_other

        # Split the motifs into lists of bases
        motif_split = self.split()
        other_split = other.split()

        motif_right = len(motif_split) - mod_pos 
        other_right = len(other_split) - mod_pos_other
        # Pad the shorter motif with '.' to ensure equal lengths for merging
        if mod_pos < mod_pos_other:
            motif_split = ['.'] * abs(offset) + motif_split
        else:
            other_split = ['.'] * abs(offset) + other_split

        if motif_right < other_right:
            motif_split = motif_split + ['.'] * (other_right - motif_right)
        else:
            other_split = other_split + ['.'] * (motif_right - other_right)
        # Merge the motifs base by base, keeping the modification position aligned
        merged_motif_string = ''.join([self.merge_bases(base, base_other) for base, base_other in zip(motif_split, other_split)])

        # The modification position in the merged motif is the same as in the original motifs
        merged_mod_position = mod_pos + abs(offset)

        # Create and return a new Motif object with the merged motif string and modification position
        return Motif(merged_motif_string, merged_mod_position)
        
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

def merge_motifs(motifs, connectivity_dist=2):
    # Create network where each node is a motif and the edge weight is the distance between them
    G = nx.Graph()
    for i, motif in enumerate(motifs):
        G.add_node(i, motif=motif)
        # Used Rebase motifs to estimate this threashold 
        # Motif shorter than5 generally dont have multi base positions
        if motif.trimmed_length() > 4:
            for j, other in enumerate(motifs):
                if i != j:
                    distance = motif.distance(other)
                    if distance <= connectivity_dist:
                        G.add_edge(i, j, weight=distance)

    # Step 1: Create a copy of the graph
    G_filtered = G.copy()

    # Step 2: Remove edges below a certain threshold
    threshold = connectivity_dist  # Set your threshold value here
    edges_to_remove = [(u, v) for u, v, d in G_filtered.edges(data=True) if d['weight'] > threshold]
    G_filtered.remove_edges_from(edges_to_remove)


    # Step 3: Find connected components (clusters)
    clusters = list(nx.connected_components(G_filtered))

    # Now 'clusters' is a list of sets, where each set is a cluster of nodes
    new_motifs = []
    single_motifs = []
    for i, cluster in enumerate(clusters, start=1):
        if len(cluster) == 1:
            single_motifs.append(motifs[list(cluster)[0]])
            continue
        print(f'Cluster {i}: {cluster}')
        merged_motif = None
        for node in cluster:
            print(f' - Node {node}: {G.nodes[node]["motif"]}')
            # Merge motifs
            if merged_motif is None:
                merged_motif = motifs[node]
            else:
                merged_motif = merged_motif.merge(G.nodes[node]["motif"])

        print(f' - Merged motif: {merged_motif}')
        new_motifs.append(merged_motif)
    return new_motifs, single_motifs


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