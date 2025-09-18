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
        assert motif.reverse_compliment() == Motif(".C.GAT", 5)
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

        motif1 = Motif("AGG", 0)
        motif2 = Motif("VAAG", 1)
        merged = motif1.merge(motif2)
        assert merged == Motif("A[AG]G", 0)
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


import pytest

def test_merge_and_find_new_variants_no_generalization():
    # Both motifs are identical -> merged should equal original, no new variants
    m1 = Motif("G.GATC", 3)
    m2 = Motif("G.GATC", 3)
    merged, pre_variants, new_variants = merge_and_find_new_variants([m1, m2])
    assert str(merged) == "G.GATC", "Merged motif should equal original"
    assert len(new_variants) == 0
    # pre_variants should contain 4 variants (pos 1 has '.' so expands to 4)
    assert pre_variants == {Motif("G.GATC", 3)}


def test_merge_and_find_new_variants_simple_generalization():
    # Merge G.GATC and C.GATG -> merged motif should generalize pos 0 and pos 5
    m1 = Motif("G.GATC", 3)
    m2 = Motif("C.GATG", 3)
    merged, pre_variants, new_variants = merge_and_find_new_variants([m1, m2])

    # Check merged motif
    assert merged == Motif("[CG].GAT[CG]", 3)

    # There should be 2 pre-variants
    assert len(pre_variants) == 2

    # There should be 2 new variants (combinations not seen before)
    assert len(new_variants) == 2

    # Check a couple of expected new variants
    expected_new = {Motif("C.GATC", 3), Motif("G.GATG", 3)}
    assert expected_new.issubset(new_variants)


def test_merge_with_fully_general_motif():
    # ..GAT. already covers G.GATC, so merge should not create new variants
    m1 = Motif("G.GATC", 3)
    m2 = Motif("..GAT.", 3)
    merged, pre_variants, new_variants = merge_and_find_new_variants([m1, m2])

    assert merged == Motif("GAT", 1)  # merged is fully general
    assert len(new_variants) == 0   # no new variants created
    assert pre_variants == {
            Motif('A.GATA', 3),
            Motif('A.GATC', 3),
            Motif('A.GATG', 3),
            Motif('A.GATT', 3),
            Motif('C.GATA', 3),
            Motif('C.GATC', 3),
            Motif('C.GATG', 3),
            Motif('C.GATT', 3),
            Motif('G.GATA', 3),
            Motif('G.GATC', 3),
            Motif('G.GATG', 3),
            Motif('G.GATT', 3),
            Motif('T.GATA', 3),
            Motif('T.GATC', 3),
            Motif('T.GATG', 3),
            Motif('T.GATT', 3)
      }  # pre_variants is just the general motif

def test_merge_with_merged_motif_general():
    m1 = Motif("G.GATC", 3)
    m2 = Motif("C.GAT.", 3)
    m3 = Motif("A.GATC", 3)
    m4 = Motif("T.GAT.", 3)
    merged, pre_variants, new_variants = merge_and_find_new_variants([m1, m2, m3, m4])

    assert merged == Motif("GAT", 1) 
    assert new_variants == {
        Motif('A.GATA', 3),
        Motif('A.GATG', 3),
        Motif('A.GATT', 3),
        Motif('G.GATA', 3),
        Motif('G.GATG', 3),
        Motif('G.GATT', 3)
    }
    assert pre_variants == {
        Motif('A.GATC', 3),
        Motif('C.GATA', 3),
        Motif('C.GATC', 3),
        Motif('C.GATG', 3),
        Motif('C.GATT', 3),
        Motif('G.GATC', 3),
        Motif('T.GATA', 3),
        Motif('T.GATC', 3),
        Motif('T.GATG', 3),
        Motif('T.GATT', 3)
    } 

def test_merge_with_merged_motif_general_offset_length():
    m1 = Motif(".G.GATC", 4)
    m2 = Motif("C.GAT.", 3)
    m3 = Motif("A.GATC", 3)
    m4 = Motif("T.GAT.", 3)
    merged, pre_variants, new_variants = merge_and_find_new_variants([m1, m2, m3, m4])

    assert merged == Motif("GAT", 1) 
    assert new_variants == {
        Motif('A.GATA', 3),
        Motif('A.GATG', 3),
        Motif('A.GATT', 3),
        Motif('G.GATA', 3),
        Motif('G.GATG', 3),
        Motif('G.GATT', 3)
    }
    assert pre_variants == {
        Motif('A.GATC', 3),
        Motif('C.GATA', 3),
        Motif('C.GATC', 3),
        Motif('C.GATG', 3),
        Motif('C.GATT', 3),
        Motif('G.GATC', 3),
        Motif('T.GATA', 3),
        Motif('T.GATC', 3),
        Motif('T.GATG', 3),
        Motif('T.GATT', 3)
    } 

def test_merge_with_fully_general_motif_different_length():
    # ..GAT. already covers G.GATC, so merge should not create new variants
    m1 = Motif("G.GATC", 3)
    m2 = Motif("GAT.", 1)
    merged, pre_variants, new_variants = merge_and_find_new_variants([m1, m2])

    assert merged == Motif("GAT", 1)  # merged is fully general
    assert len(new_variants) == 0   # no new variants created

def test_merge_three_motifs_creates_more_variants():
    # Three motifs that cover different corners of variant space
    m1 = Motif("G.GATC", 3)
    m2 = Motif("C.GATC", 3)
    m3 = Motif("G.GATG", 3)

    merged, pre_variants, new_variants = merge_and_find_new_variants([m1, m2, m3])

    assert merged == Motif("[CG].GAT[CG]", 3)

    # There should be new variants where we combine C at pos0 and G at pos5
    assert Motif("C.GATG", 3) in new_variants
    # And G at pos0 and C at pos5 (which existed before) should not appear in new_variants
    assert Motif("G.GATC", 3) not in new_variants


def test_merge_with_only_dot_positions_returns_empty_mask():
    # All motifs have only '.', so mask should be empty and explode returns a single variant
    m1 = Motif(".....", 2)
    m2 = Motif(".....", 2)
    merged, pre_variants, new_variants = merge_and_find_new_variants([m1, m2])

    assert merged == Motif(".....", 2)
    assert pre_variants == {Motif(".....", 2)}  # no explosion since no explicit positions
    assert len(new_variants) == 0


def test_merge_no_strip():
    # Merging motifs without stripping should still work correctly
    m1 = Motif("..A..", 2)
    m2 = Motif("..C..", 2)
    m_merge = m1.merge_no_strip(m2)

    assert str(m_merge) == "..[AC].."

    m1 = Motif("G.A..", 2)
    m2 = Motif("..C..", 2)
    m_merge = m1.merge_no_strip(m2)

    assert str(m_merge) == "..[AC].."

    m1 = Motif(".A..", 1)
    m2 = Motif("..C..", 2)
    m_merge = m1.merge_no_strip(m2)

    assert str(m_merge) == "..[AC].."

    m1 = Motif(".A......", 1)
    m2 = Motif("..C..", 2)
    m_merge = m1.merge_no_strip(m2)

    assert str(m_merge) == "..[AC]......"

    m1 = Motif("G.GATC", 3)
    m2 = Motif("C.GATC", 3)
    m_merge = m1.merge_no_strip(m2)
    assert str(m_merge) == "[CG].GATC"

    m1 = Motif("[CG].GATC", 3)
    m2 = Motif("G.GATG", 3)
    m_merge = m1.merge_no_strip(m2)
    assert str(m_merge) == "[CG].GAT[CG]"


def test_align_motifs():
    m1 = Motif("G.GATC", 3)
    m2 = Motif("C.GATC", 3)
    m3 = Motif("G.GATG", 3)
    aligned = align_motifs([m1, m2, m3])
    assert aligned == [m1, m2, m3]

    m1 = Motif("..A..", 2)
    m2 = Motif("..C..", 2)
    aligned = align_motifs([m1, m2])
    assert aligned == [m1, m2]

    m1 = Motif("G.GATC", 3)
    m2 = Motif("GAT.", 1)
    aligned = align_motifs([m1, m2])
    assert aligned == [Motif("G.GATC", 3), Motif("..GAT.", 3)]

    m1 = Motif("G.GATC", 3)
    m2 = Motif("GAT..", 1)
    m3 = Motif(".A...", 1)
    aligned = align_motifs([m1, m2, m3])
    assert aligned == [Motif("G.GATC.", 3), Motif("..GAT..", 3), Motif("...A...", 3)]