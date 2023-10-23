from nanomotif.candidate import *
from nanomotif.constants import *
import pytest
from hypothesis import given, strategies as st
import itertools
import pandas as pd
import polars.testing

# Strategy for a single base
base = st.sampled_from("ATCG.")
# Strategy for a bracketed sequence
bracketed = st.text(min_size=2, max_size=3, alphabet="ATCG", ).map(lambda s: f"[{''.join(list(set(s)))}]")
# Strategy for a list with bracketed and unbracketed sequences
combined = st.lists(st.one_of(base, bracketed), min_size=1)

class TestMotif:
    @given(st.text(alphabet=["A", "C", "G", "T", "."], min_size=5, max_size=10))
    def test_strip(self, motif_string):
        motif = Motif(motif_string, 0)
        stripped = motif.strip()
        assert stripped == motif_string.lstrip(".").rstrip(".")
    
    @given(combined)
    def test_split(self, string_list):
        motif_string = "".join(string_list)
        motif = Motif(motif_string, 0)
        split_motif = motif.split()
        # Splitting of motif removes braces
        assert split_motif == string_list

    @given(combined)
    def test_array(self, motif_list):
        motif_string = "".join(motif_list)
        motif = Motif(motif_string, 0)
        array = motif.one_hot()
        assert array.shape == (len(motif_list), 4)
        assert array.dtype == int
        assert array.min() >= 0
        assert array.max() <= 1
    
    def test_reverse_compliment(self):
        motif = Motif("ATCG", 0)
        assert motif.reverse_compliment() == Motif("CGAT", 3)
        motif = Motif("AT[CG]G", 0)
        assert motif.reverse_compliment() == Motif("C[CG]AT", 3)
        motif = Motif("ATC.G.", 0)
        assert motif.reverse_compliment() == Motif(".C.GAT", 3)
        motif = Motif("ATAC.G.", 2)
        assert motif.reverse_compliment() == Motif(".C.GTAT", 4)
    
    def test_new_stripped_motif(self):
        motif = Motif("....ATCG", 4)
        assert motif.new_stripped_motif() == Motif("ATCG", 0)
        motif = Motif("ATCG....", 0)
        assert motif.new_stripped_motif() == Motif("ATCG", 0)
        motif = Motif("....AT..CG", 4)
        assert motif.new_stripped_motif() == Motif("AT..CG", 0)
    
    def test_strip_reverse(self):
        motif = Motif("....AT..CG..", 4)
        motif_strip = motif.new_stripped_motif()
        motif_strip_reverse = motif_strip.reverse_compliment()
        assert motif_strip_reverse == Motif("CG..AT", 5)
    
    def test_child_of(self):
        motif1 = Motif("ATCG", 0)
        motif2 = Motif("AT", 0)
        assert motif1.child_of(motif2) == True

        motif1 = Motif("A[TCG]CG", 0)
        motif2 = Motif("A.C", 0)
        assert motif1.child_of(motif2) == True

        motif1 = Motif("ATCG", 0)
        motif2 = Motif("AT", 0)
        assert motif1.child_of(motif2) == False

        motif1 = Motif("ATCG", 2)
        motif2 = Motif("CG", 0)
        assert motif1.child_of(motif2) == False

        motif1 = Motif("ATCG", 2)
        motif2 = Motif("CG", 2)
        assert motif1.child_of(motif2) == False

        motif1 = Motif("ATCG", 2)
        motif2 = Motif("CG", 0)
        assert motif1.child_of(motif2) == True

        motif1 = Motif("CG", 2)
        motif2 = Motif("ATCG", 2)
        assert motif1.child_of(motif2) == False



