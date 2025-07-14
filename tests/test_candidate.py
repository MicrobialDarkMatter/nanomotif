from nanomotif.motif import *
from nanomotif.constants import *
from nanomotif.utils import *
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
    
    def test_sub_motif_of(self):
        motif1 = Motif("ATCG", 2)
        motif2 = Motif("ATCG", 2)
        assert motif1.sub_motif_of(motif2) == False

        motif1 = Motif("ATCG", 0)
        motif2 = Motif("AT", 0)
        assert motif1.sub_motif_of(motif2) == True

        motif1 = Motif("A[TCG]CG", 0)
        motif2 = Motif("A.C", 0)
        assert motif1.sub_motif_of(motif2) == True

        motif1 = Motif("ATCG", 0)
        motif2 = Motif("AT", 0)
        assert motif1.sub_motif_of(motif2) == True

        motif1 = Motif("ATCG", 2)
        motif2 = Motif("CG", 0)
        assert motif1.sub_motif_of(motif2) == True

        motif1 = Motif("ATCG", 2)
        motif2 = Motif("CG", 2)
        assert motif1.sub_motif_of(motif2) == False

        motif1 = Motif("CG", 1)
        motif2 = Motif("ATCG", 2)
        assert motif1.sub_motif_of(motif2) == False

    def test_sub_string_of(self):
        # Equal motifs strings
        motif1 = Motif("ATCG", 2)
        motif2 = Motif("ATCG", 0)
        assert motif1.sub_string_of(motif2) == False

        # Single base difference
        motif1 = Motif("AGCG", 2)
        motif2 = Motif("ATCG", 0)
        assert motif1.sub_string_of(motif2) == False

        # Longer sub string of shorter
        motif1 = Motif("ATCG", 0)
        motif2 = Motif("AT", 0)
        assert motif1.sub_string_of(motif2) == True

        # Workrking with braces
        motif1 = Motif("ACC", 0)
        motif2 = Motif("A[CG]C", 2)
        assert motif1.sub_string_of(motif2) == True

        # Working with dots
        motif1 = Motif("ATCG", 0)
        motif2 = Motif("AT.G", 2)
        assert motif1.sub_string_of(motif2) == True

        # Works when braces in motif1
        motif1 = Motif("AT.G", 0)
        motif2 = Motif("ATCG", 2)
        assert motif1.sub_string_of(motif2) == False

        # Works when braces in motif1
        motif1 = Motif("AT[ACG]G", 0)
        motif2 = Motif("ATCG", 0)
        assert motif1.sub_string_of(motif2) == False

        # Works when motif1 is longer
        motif1 = Motif("T.....ATCG...C", 0)
        motif2 = Motif("ATCG", 2)
        assert motif1.sub_string_of(motif2) == True

        # Works when motif2 is longer
        motif1 = Motif("ATCG", 2)
        motif2 = Motif("T.....ATCG...C", 0)
        assert motif1.sub_string_of(motif2) == False

        # Dots and braces in both
        motif1 = Motif("[AT]..AT..CG[TA]C..C", 2)
        motif2 = Motif("[TA]T..CG[AT][GC]", 0)
        assert motif1.sub_string_of(motif2) == True

    def test_distance(self):
        motif1 = Motif("A[TCG]CG", 0)
        motif2 = Motif("A.C", 0)
        assert motif1.distance(motif2) == 2

        motif1 = Motif("ATCG", 0)
        motif2 = Motif("AT", 0)
        assert motif1.distance(motif2) == 2

        motif1 = Motif("ATCG", 2)
        motif2 = Motif("CG", 0)
        assert motif1.distance(motif2) == 2

        motif1 = Motif("CG", 0)
        motif2 = Motif("ATCG", 2)
        assert motif1.distance(motif2) == 2

        motif1 = Motif("C...ATCG", 6)
        motif2 = Motif("CG", 0)
        assert motif1.distance(motif2) == 3

        motif1 = Motif("CG", 0)
        motif2 = Motif("C...ATCG", 6)
        assert motif1.distance(motif2) == 3

        motif1 = Motif("C..CG...G", 3)
        motif2 = Motif("CG", 0)
        assert motif1.distance(motif2) == 2

        motif1 = Motif("CG", 0)
        motif2 = Motif("C..CG...G", 3)
        assert motif1.distance(motif2) == 2

        motif1 = Motif("CG", 0)
        motif2 = Motif("C..CG...[GC]", 3)
        assert motif1.distance(motif2) == 2

    def test_merge(self):
        motif1 = Motif("..A..", 2)
        motif2 = Motif("..C..", 2)
        assert motif1.merge(motif2) == Motif("[AC]", 0)

        motif1 = Motif(".ATGC.", 1)
        motif2 = Motif(".TAAG", 2)
        merged = motif1.merge(motif2)
        assert merged == Motif("A[AT]G", 0)
        assert merged.mod_position == 0

        motif1 = Motif("AGG", 0)
        motif2 = Motif("VAAG", 1)
        merged = motif1.merge(motif2)
        assert merged == Motif("A[AG]G", 1)
        assert merged.mod_position == 0

        motif1 = Motif("GATC...R", 1)
        motif2 = Motif("A...TATC", 5)
        merged = motif1.merge(motif2)
        assert merged == Motif("[GT]ATC", 1)


    
    def test_iupac(self):
        motif = Motif("A[TCG]CG", 0)
        assert motif.iupac() == "ABCG"
        motif = Motif("A.C", 0)
        assert motif.iupac() == "ANC"
        motif = Motif("ATCG", 0)
        assert motif.iupac() == "ATCG"
        motif = Motif("ATCG", 2)
        assert motif.iupac() == "ATCG"
        motif = Motif("CG", 0)
        assert motif.iupac() == "CG"
        motif = Motif("C...ATCG", 6)
        assert motif.iupac() == "CNNNATCG"
        motif = Motif("C..CG...G", 3)
        assert motif.iupac() == "CNNCGNNNG"
        motif = Motif("C..CG...[GC]", 3)
        assert motif.iupac() == "CNNCGNNNS"

        

def test_has_n_character_stretches_of_length_m():
    assert has_n_character_stretches_of_length_m("NNA", 1, 2, "N") == True
    assert has_n_character_stretches_of_length_m("NNA", 2, 1, "N") == False
    assert has_n_character_stretches_of_length_m("NNA", 2, 2, "N") == False
    assert has_n_character_stretches_of_length_m("NNANN", 2, 2, "N") == True
    assert has_n_character_stretches_of_length_m("NNANN", 1, 2, "N") == True
    assert has_n_character_stretches_of_length_m("NNANN", 2, 1, "N") == True
    assert has_n_character_stretches_of_length_m("NNANN", 1, 1, "N") == True
    assert has_n_character_stretches_of_length_m("NNANN", 3, 1, "N") == False
    assert has_n_character_stretches_of_length_m("NNANN", 4, 1, "N") == False
    assert has_n_character_stretches_of_length_m("NNANN", 3, 2, "N") == False
    assert has_n_character_stretches_of_length_m("NNANN", 4, 2, "N") == False
    assert has_n_character_stretches_of_length_m("NNANN", 3, 3, "N") == False
    assert has_n_character_stretches_of_length_m("NNANN", 4, 3, "N") == False
    assert has_n_character_stretches_of_length_m("NaNNaNNNaNNNN", 4, 1, "N") == True
    assert has_n_character_stretches_of_length_m("NaNNaNNNaNNNN", 4, 2, "N") == False
    assert has_n_character_stretches_of_length_m("NaNNaNNNaNNNN", 4, 3, "N") == False
    assert has_n_character_stretches_of_length_m("NaNNaNNNaNNNN", 4, 4, "N") == False
    assert has_n_character_stretches_of_length_m("NaNNaNNNaNNNN", 3, 1, "N") == True
    assert has_n_character_stretches_of_length_m("NaNNaNNNaNNNN", 3, 2, "N") == True
    assert has_n_character_stretches_of_length_m("NaNNaNNNaNNNN", 3, 3, "N") == False
    assert has_n_character_stretches_of_length_m("NaNNaNNNaNNNN", 3, 4, "N") == False
    assert has_n_character_stretches_of_length_m("NaNNaNNNaNNNN", 2, 1, "N") == True
    assert has_n_character_stretches_of_length_m("NaNNaNNNaNNNN", 2, 2, "N") == True
    assert has_n_character_stretches_of_length_m("NaNNaNNNaNNNN", 2, 3, "N") == True
    assert has_n_character_stretches_of_length_m("NaNNaNNNaNNNN", 2, 4, "N") == False


