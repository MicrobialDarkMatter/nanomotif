from nanomotif.constants import *
from nanomotif.motif import Motif
from nanomotif.seq import DNAsequence, DNAarray
import nanomotif as nm
import pytest
from hypothesis import given, strategies as st
import itertools
import pandas as pd
import polars as pl
import polars.testing
import numpy as np


class TestMethylatedMotifOccurrences:
    def test_known_motif_sequence(self):
        motif = nm.motif.Motif('ACG', 0)
        sequence = 'TACGGACGCCACG'
        methylated_positions = np.array([1, 5])
        non_methylated_positions = np.array([10])

        expected_methylated = np.array([1, 5])
        expected_nonmethylated = np.array([10])
        
        result = nm.find_motifs_bin.methylated_motif_occourances(motif, sequence, methylated_positions, non_methylated_positions)
        np.testing.assert_array_equal(result[0], expected_methylated)
        np.testing.assert_array_equal(result[1], expected_nonmethylated)

    def test_known_motif_sequence_no_methylated(self):
        motif = nm.motif.Motif('ACG', 0)
        sequence = 'TACGGACGCCACG'
        methylated_positions = np.array([])
        non_methylated_positions = np.array([1, 10])

        expected_methylated = np.array([])
        expected_nonmethylated = np.array([1, 10])
        
        result = nm.find_motifs_bin.methylated_motif_occourances(motif, sequence, methylated_positions, non_methylated_positions)
        np.testing.assert_array_equal(result[0], expected_methylated)
        np.testing.assert_array_equal(result[1], expected_nonmethylated)


def test_get_motif_parental_relationships():
    motif1 = Motif("ACGT", 0)
    motif2 = Motif("ACG", 0)
    motif3 = Motif("CGT", 1)
    motif4 = Motif("ACGTG", 0)
    motif5 = Motif("TGCA", 1)
    
    motifs = [motif1, motif2, motif3, motif4, motif5]
    expected_relationships = set([
        (motif2, motif1),
        (motif2, motif4),
        (motif3, motif1),
        (motif3, motif4),
        (motif1, motif4)
    ])
    relationships = set(nm.postprocess.get_motif_parental_relationship(motifs))
    assert relationships == expected_relationships

    
    motif6 = Motif("GGCA[AT]", 2)
    motif7 = Motif("GGCAAT", 2)
    motif8 = Motif("GGCAAT", 4)
    motif9 = Motif("AATTT", 0)
    motif10 = Motif("AATTT", 1)
    motif11 = Motif("AATTTT", 0)
    motifs = [motif6, motif7, motif8, motif9, motif10, motif11]
    expected_relationships = set([
        (motif6, motif7),
        (motif6, motif8),
        (motif9, motif11),
        (motif10, motif11)
    ])
    relationships = set(nm.postprocess.get_motif_parental_relationship(motifs))
    assert relationships == expected_relationships

    motif12 = Motif("TTAAGGAG", 6)
    motif13 = Motif("TTAA", 3)
    motifs = [motif12, motif13]
    expected_relationships = set([
        (motif13, motif12)
    ])
    relationships = set(nm.postprocess.get_motif_parental_relationship(motifs))
    assert relationships == expected_relationships