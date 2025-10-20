import numpy as np
import re
import itertools
from functools import reduce
from itertools import product
import networkx as nx
import logging as log
from typing import overload, Union, Sequence, Any
from scipy.stats import entropy
from nanomotif.constants import *
from nanomotif.seq import regex_to_iupac, iupac_to_regex, reverse_compliment
from nanomotif.model import BetaBernoulliModel
import nanomotif.utils as nm_utils
import polars as pl
from collections.abc import Mapping
class Motif(str):
    def __new__(cls, motif_string, *args, **kwargs):
        return str.__new__(cls, motif_string)
    
    def __init__(self, _, mod_position):
        self.mod_position = mod_position
        self.string = self.__str__()
    
    def __eq__(self, other):
        if isinstance(other, Motif):
            return self.string == other.string and self.mod_position == other.mod_position
        return False

    def __repr__(self):
        # What shows in lists or when inspecting in the console
        return f"Motif({super().__str__()!r}, pos={self.mod_position})"

    def __hash__(self):
        return hash((self.string, self.mod_position))

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

        if self_stripped.string == other_stripped.string:
            return False
        size_difference = self_stripped.length() - other_stripped.length()

        # Split into list of bases
        self_split = self_stripped.split()
        other_split = other_stripped.split()
        
        for i in range (size_difference+1):
            match = True
            for j, base in enumerate(other_split):
                if j + i >= len(self_split):
                    continue
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
    
    def count_isolated_bases(self, isolation_size=2):
        """
        Count the number of isolated bases in a motif.
        """
        motif_split = self.split()
        isolated_count = 0
        for pos in range(len(motif_split)):
            if motif_split[pos] == ".":
                continue
            index_start = max(pos - isolation_size, 0)
            index_end = min(pos + isolation_size + 1, len(motif_split) - 1)
            # If all surrounding positions are ".", it is isolated
            if set(motif_split[index_start:pos] + motif_split[pos+1:index_end]) == set(["."]):
                isolated_count += 1
            if set(motif_split[index_start:pos] + motif_split[pos+1:index_end]) == set(["N"]):
                isolated_count += 1
        return isolated_count

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
        matches =  re.search("[^.]", self.string)
        if matches is None:
            return self

        new_position = self.mod_position - matches.start()
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
    
    def merge_no_strip(self, other):
        assert isinstance(other, Motif)
        self_seq = self.split()
        other_seq = other.split()
        length_self = len(self_seq)
        length_other = len(other_seq)

        offset = self.mod_position - other.mod_position
        if offset > 0:
            self_seq = self_seq[offset:]
        elif offset < 0:
            other_seq = other_seq[-offset:]

        final_length = min(len(self_seq), len(other_seq))
        self_seq = self_seq[:final_length]
        other_seq = other_seq[:final_length]
        
        # Merge the motifs base by base, keeping the modification position aligned
        merged_motif_string = ''.join([self.merge_bases(base, base_other) for base, base_other in zip(self_seq, other_seq)])
        merged_mod_position = min(self.mod_position, other.mod_position)

        # Pad to ensure keeps all positions in both of the original motifs
        if self.mod_position < other.mod_position:
            padding = other.mod_position - self.mod_position
            merged_motif_string = "." * padding + merged_motif_string
            merged_mod_position = self.mod_position
        elif other.mod_position < self.mod_position:
            padding = self.mod_position - other.mod_position
            merged_motif_string = "." * padding + merged_motif_string
            merged_mod_position = other.mod_position
        if length_self - self.mod_position > length_other - other.mod_position:
            padding = (length_self - self.mod_position) - (length_other - other.mod_position)
            merged_motif_string = merged_motif_string + "." * padding
        elif length_other - other.mod_position > length_self - self.mod_position:
            padding = (length_other - other.mod_position) - (length_self - self.mod_position)
            merged_motif_string = merged_motif_string + "." * padding

        # Create and return a new Motif object with the merged motif string and modification position
        merged_motif = Motif(merged_motif_string, merged_mod_position)
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


def align_motifs(motifs):
    """
    Align motifs by their modification positions.
    Returns list of aligned motifs.
    """

    if not motifs:
        return []

    # Align left side
    min_mod_pos = min(motif.mod_position for motif in motifs)
    max_mod_pos = max(motif.mod_position for motif in motifs)
    aligned_motifs = []

    for motif in motifs:
        prefix_dots = max_mod_pos - motif.mod_position
        new_motif_string = '.' * prefix_dots + motif.string
        new_mod_position = max_mod_pos
        aligned_motifs.append(Motif(new_motif_string, new_mod_position))
    # Align right side
    max_length = max(motif.length() for motif in aligned_motifs)
    for i, motif in enumerate(aligned_motifs):
        suffix_dots = max_length - motif.length()
        new_motif_string = motif.string + '.' * suffix_dots
        aligned_motifs[i] = Motif(new_motif_string, motif.mod_position) 
    return aligned_motifs

def explode_with_mask(motif, mask):
    """
    Explode motif but only at positions given by mask (list of indices).
    Returns set of exploded strings, preserving '.' at unmasked positions.
    """
    segments = motif.split()

    # Build a list of possible bases for each masked position
    bases_options = []
    for i in mask:
        seg = segments[i]
        if seg == ".":
            bases_options.append(["A", "C", "G", "T"])
        elif seg.startswith("["):
            bases_options.append(list(seg[1:-1]))
        else:
            bases_options.append([seg])

    variants = set()
    for combo in product(*bases_options):
        new_seq = segments.copy()
        for idx, pos in enumerate(mask):
            new_seq[pos] = combo[idx]
        variants.add("".join(new_seq))

    # Convert to Motif objects
    variants = {Motif(v, motif.mod_position) for v in variants}  
    return variants



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
        return list(nx.find_cliques(self))

def merge_and_find_new_variants(motifs):
    """
    Merge motifs and return:
      - merged motif
      - pre-merge exploded variants (masked)
      - new variants that appear after merge
    """
    if not motifs:
        return None, set(), set()
    motifs = align_motifs(motifs)

    # Determine mask (positions that are not '.' in ANY motif)
    all_segments = [m.split() for m in motifs]
    n = len(all_segments[0])
    mask = [
        i for i in range(n)
        if any(seg[i] != "." for seg in all_segments)
    ]

    # Explode pre-merge variants
    pre_variants = set()
    for m in motifs:
        pre_variants |= explode_with_mask(m, mask)

    # Merge motifs
    merged = reduce(lambda a, b: a.merge_no_strip(b), motifs)

    # Explode post-merge variants using same mask
    post_variants = explode_with_mask(merged, mask)

    new_variants = post_variants - pre_variants

    new_variants = {v.new_stripped_motif() for v in new_variants}
    pre_variants = {v.new_stripped_motif() for v in pre_variants}
    merged = merged.new_stripped_motif()
    return merged, pre_variants, new_variants


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
        cluster_motifs = list(cluster)
        merged_motif, pre_variants, new_variants = merge_and_find_new_variants(cluster_motifs)
        if merged_motif is None:
            continue
        if merged_motif.trimmed_length() < min_length:
            continue



        log.debug(f'Merged motif: {merged_motif}')
        motif_merged[i] = [merged_motif, cluster_motifs, pre_variants, new_variants]
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


class MotifSearchResult(pl.DataFrame):

    REQUIRED_COLUMNS = {
        "reference": pl.Utf8, 
        "motif": pl.Utf8,
        "mod_type": pl.Utf8, 
        "mod_position": pl.Int64, 
        "model": pl.Object,
        "score": pl.Float64,
    }

    DERIVED_COLUMNS = {
        "n_mod": pl.Int64,
        "n_nomod": pl.Int64,
        "motif_iupac": pl.Utf8,
        "mod_position_iupac": pl.Int64
    }

    COMPLEMENTARY_COLUMNS = {
        "motif_complement": pl.Utf8,
        "mod_position_complement": pl.Int64,
        "score_complement": pl.Float64,
        "model_complement": pl.Object,
        "n_mod_complement": pl.Int64,
        "n_nomod_complement": pl.Int64,
        "motif_iupac_complement": pl.Utf8,
        "mod_position_iupac_complement": pl.Int64
    }


    def __init__(self, data=None, *args, **kwargs):
        if isinstance(data, pl.DataFrame):
            # Use Polars' internal constructor to wrap without copying
            super().__init__(pl.DataFrame._from_pydf(data._df))
        else:
            super().__init__(data, *args, **kwargs)

        self._normalize_model_column()
        self._ensure_required()
        self._ensure_derived()
        self._ensure_column_order()
        self._coerce_dtypes()

    def __getattribute__(self, name):
        # Ignore dunder methods
        if name.startswith("__"):
            return object.__getattribute__(self, name)
        attr = object.__getattribute__(self, name)

        if isinstance(getattr(self.__class__, name, None), property):
            return attr
        if callable(attr):
            def wrapper(*args, **kwargs):
                out = attr(*args, **kwargs)
                if isinstance(out, pl.DataFrame) and not isinstance(out, self.__class__):
                    obj = self.__class__.__new__(self.__class__)
                    obj._df = out._df
                    return obj
                return out
            return wrapper
        return attr


    def __getstate__(self):
        """Serialize to a pickle-friendly dict, converting model objects to state dicts."""
        serialized = {}
        for col in self.columns:
            values = self[col].to_list()
            if any(isinstance(v, BetaBernoulliModel) for v in values if v is not None):
                serialized[col] = [v.__getstate__() if v is not None else None for v in values]
            else:
                serialized[col] = values
        return {"data": serialized, "schema": self.schema}

    def __setstate__(self, state: dict):
        """Restore from pickle-friendly dict, rehydrating model objects."""
        restored = {}
        for col, values in state["data"].items():
            if isinstance(values, list) and values and isinstance(values[0], dict) and "_alpha" in values[0]:
                # Looks like a BetaBernoulliModel state dict
                restored[col] = [
                    self._restore_model_from_state(v) if v is not None else None
                    for v in values
                ]
            else:
                restored[col] = values

        df = pl.DataFrame(restored)
        self._df = df._df  # reinitialize Polars parent
        self._normalize_model_column()  # <-- this rebuilds model column as pl.Object
        self._ensure_column_order()

    @staticmethod
    def _restore_model_from_state(state_dict):
        model = BetaBernoulliModel.__new__(BetaBernoulliModel)
        model.__setstate__(state_dict)
        return model

    def _load_motifs(self):
        assert "motif" in self.columns
        assert "mod_position" in self.columns
        assert len(self) > 0
        motif_strings = self.get_column("motif").to_list()
        positions = self.get_column("mod_position").to_list()
        motifs = [Motif(motif_string, pos) for motif_string, pos in zip(motif_strings, positions)]
        return motifs

    def _ensure_required(self):
        # Ensure required columns are present
        # Convert contig or bin to reference if present
        if "contig" in self.columns:
            self._df = self.rename({"contig": "reference"})._df
        elif "bin" in self.columns:
            self._df = self.rename({"bin": "reference"})._df

        missing = [c for c in self.REQUIRED_COLUMNS if c not in self.columns]

        if missing:
            raise ValueError(f"Missing required columns: {missing}")
        
    def _ensure_derived(self):
        """Compute missing derived columns if possible."""
        to_add = {}
        for col in self.columns:
            if self[col].dtype == pl.Object:
                values = self[col].to_list()
                if not values:
                    continue
                if not any(isinstance(v, BetaBernoulliModel) for v in values if v is not None):
                    continue  # skip if no BetaBernoulliModel inside

                # Figure out suffix, e.g. model -> "", model_complement -> "_complement"
                suffix = col.replace("model", "", 1)

                n_mod_col = f"n_mod{suffix}"
                n_nomod_col = f"n_nomod{suffix}"

                if n_mod_col not in self.columns or n_nomod_col not in self.columns:
                    n_mod_vals, n_nomod_vals = [], []
                    for m in values:
                        if m is None:
                            n_mod_vals.append(None)
                            n_nomod_vals.append(None)
                        else:
                            n_mod_vals.append(m._alpha - m._alpha_prior)
                            n_nomod_vals.append(m._beta - m._beta_prior)
                    to_add[n_mod_col] = n_mod_vals
                    to_add[n_nomod_col] = n_nomod_vals
        if "motif_iupac" not in self.columns:
            if "motif" not in self.columns:
                raise ValueError(
                    "Cannot derive motif_iupac because 'motif' column is missing"
                )
            if "mod_position" not in self.columns:
                raise ValueError(
                    "Cannot derive motif_iupac because 'mod_position' column is missing"
                )
            motifs = self._load_motifs()
            motifs_trimmed = [m.new_stripped_motif() for m in motifs]
            to_add["mod_position_iupac"] = [m.mod_position for m in motifs_trimmed]
            to_add["motif_iupac"] = [m.iupac() for m in motifs_trimmed]
        if to_add:
            # Concatenate derived columns
            derived_df = pl.DataFrame(to_add)
            self._df = self.hstack(derived_df)._df 
    
    def _coerce_dtypes(self):
        """Coerce any existing columns to expected dtype if mismatched."""
        exprs = []
        columns_to_check = {**self.REQUIRED_COLUMNS, **self.DERIVED_COLUMNS, **self.COMPLEMENTARY_COLUMNS}
        for c in self.columns:
            expected = columns_to_check.get(c)

            if expected is None:
                # no explicit expectation, but handle Nulls
                if self[c].dtype == pl.Null:
                    # Default fallback: make Nulls Utf8 (safe for metadata) 
                    exprs.append(pl.col(c).cast(pl.Utf8))
            else:
                if self[c].dtype != expected:
                    # Cast both mismatches and Nulls
                    exprs.append(pl.col(c).cast(expected))

        if exprs:
            coerced = self.with_columns(exprs)
            self._df = coerced._df

    def _ensure_column_order(self):
        if self.COMPLEMENTARY_COLUMNS.keys() & set(self.columns):
            fixed_order_columns = (
                list(self.REQUIRED_COLUMNS.keys()) +
                list(self.DERIVED_COLUMNS.keys()) +
                list(self.COMPLEMENTARY_COLUMNS.keys())
            )
        else:
            fixed_order_columns = list(self.REQUIRED_COLUMNS.keys()) + list(self.DERIVED_COLUMNS.keys())
        ordered_columns = fixed_order_columns + [col for col in self.columns if col not in fixed_order_columns]
        self._df = self.select(ordered_columns)._df

    def _normalize_model_column(self):
        """Ensure `model` column is pl.Object with BetaBernoulliModel instances."""
        model_cols = [c for c in self.columns if c.startswith("model")]

        for col in model_cols:
            s = self[col]
            # If already correct (Object dtype with only BetaBernoulliModel or None), skip
            if s.dtype == pl.Object and all(
                (v is None or isinstance(v, BetaBernoulliModel)) for v in s.to_list()
            ):
                continue

            restored = []
            for v in s.to_list():
                if v is None:
                    restored.append(None)
                elif isinstance(v, BetaBernoulliModel):
                    restored.append(v)
                elif isinstance(v, dict) or hasattr(v, "keys"):
                    restored.append(self._restore_model_from_state(dict(v)))
                else:
                    # Last resort: try attribute access (if v looks like a model-ish object)
                    try:
                        restored.append(self._restore_model_from_state({
                            "_alpha": v._alpha,
                            "_beta": v._beta,
                            "_alpha_prior": v._alpha_prior,
                            "_beta_prior": v._beta_prior,
                        }))
                    except Exception:
                        restored.append(None)

            # Drop and reinsert column to guarantee dtype=pl.Object
            other_cols = [c for c in self.columns if c != col]
            df = self.select(other_cols)
            df = df.with_columns(pl.Series(col, restored, dtype=pl.Object))
            self._df = df._df

    def write_motifs(self, file_path):
        df = pl.DataFrame(self)
        df = df.sort(["reference", "mod_type", "motif"])
        for column in df.columns:
            if df[column].dtype == pl.Object:
                df = df.drop(column)
        df.write_csv(file_path, separator="\t")

    def write_motif_formatted(self, file_path):
        df = pl.DataFrame(self)
        df_formatted = df.select([
            "reference", "motif_iupac", "mod_position_iupac", "mod_type", "n_mod", "n_nomod"
        ])
        df_formatted = df_formatted.rename({
            "motif_iupac": "motif",
            "mod_position_iupac": "mod_position"
        })
        # Add motif type column
        df_formatted = df_formatted.with_columns([
            pl.col("motif").map_elements(lambda x: nm_utils.motif_type(x)).alias("motif_type")
        ])
        # Include complement columns if present
        complement_cols = [
            "motif_iupac_complement", "mod_position_iupac_complement", "n_mod_complement", "n_nomod_complement"
        ]
        if all(col in df.columns for col in complement_cols):
            df_complement = df.select([
                pl.col("motif_iupac_complement").alias("motif_complement"),
                pl.col("mod_position_iupac_complement").alias("mod_position_complement"),
                pl.col("n_mod_complement"),
                pl.col("n_nomod_complement")
            ])
            df_formatted = df_formatted.hstack(df_complement)

        df_formatted = df_formatted.sort(["reference", "mod_type", "motif"])
        df_formatted.write_csv(file_path, separator="\t")



if __name__ == "__main__":
    motif1 = Motif("AT.G", 0)
    motif2 = Motif("ATCG", 0)
    motif1.sub_string_of(motif2)        