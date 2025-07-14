from nanomotif.constants import *
from nanomotif.motif import Motif
from nanomotif.seq import DNAsequence, DNAarray
from nanomotif.find_motifs import MotifSearcher
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
        
        result = nm.find_motifs.methylated_motif_occourances(motif, sequence, methylated_positions, non_methylated_positions)
        np.testing.assert_array_equal(result[0], expected_methylated)
        np.testing.assert_array_equal(result[1], expected_nonmethylated)

    def test_known_motif_sequence_no_methylated(self):
        motif = nm.motif.Motif('ACG', 0)
        sequence = 'TACGGACGCCACG'
        methylated_positions = np.array([])
        non_methylated_positions = np.array([1, 10])

        expected_methylated = np.array([])
        expected_nonmethylated = np.array([1, 10])
        
        result = nm.find_motifs.methylated_motif_occourances(motif, sequence, methylated_positions, non_methylated_positions)
        np.testing.assert_array_equal(result[0], expected_methylated)
        np.testing.assert_array_equal(result[1], expected_nonmethylated)


class TestMotifSearcher:
    def test_init_with_valid_inputs(self):
        root_motif = Motif("ACGT", 0)
        contig_sequence = DNAsequence("ACGTACGTACGT")
        contig_sequences_sample = contig_sequence.sample_n_subsequences(
            2 * 2 + 1, 10
        )
        contig_pssm = contig_sequences_sample.pssm()
        contig_pileup = pl.DataFrame({"position": [1,2,3], "value": [0.1, 0.2, 0.3]})
        methylation_sequences = DNAarray(["ACGT", "CGTA", "GTAC"])
        padding = 2

        # Create MotifSearcher instance
        searcher = MotifSearcher(
            root_motif,
            contig_sequence,
            contig_pssm,
            contig_pileup,
            methylation_sequences,
            padding,
            0.7,
            0.3
        )

        # Check that attributes are set correctly
        assert searcher.root_motif == root_motif
        assert searcher.contig_sequence == contig_sequence
        assert searcher.contig_pileup.equals(contig_pileup)
        assert np.array_equal(searcher.methylation_sequences, methylation_sequences)
        assert searcher.padding == padding

    def test_priority_function_zero_division_beta(self):
        root_motif = Motif("ACGT", 0)
        contig_sequence = DNAsequence("ACGTACGTACGT")
        contig_sequences_sample = contig_sequence.sample_n_subsequences(
            2 * 2 + 1, 10
        )
        contig_pssm = contig_sequences_sample.pssm()
        contig_pileup = pl.DataFrame({"position": [1,2,3], "value": [0.1, 0.2, 0.3]})
        methylation_sequences = DNAarray(["ACGT", "CGTA", "GTAC"])
        padding = 2

        # Create MotifSearcher instance
        searcher = MotifSearcher(
            root_motif,
            contig_sequence,
            contig_pssm,
            contig_pileup,
            methylation_sequences,
            padding,
            0.7,
            0.3
        )

        class MockModel:
            def __init__(self, alpha, beta):
                self._alpha = alpha
                self._beta = beta

        next_model = MockModel(alpha=2, beta=3)
        root_model = MockModel(alpha=4, beta=0)

        priority = searcher._priority_function(next_model, root_model)

        expected_d_alpha = 1 - (next_model._alpha / root_model._alpha)  # 1 - (2/4) = 0.5
        expected_d_beta = 1  # ZeroDivisionError caught
        expected_priority = expected_d_alpha * expected_d_beta  # 0.5 * 1 = 0.5

        assert priority == pytest.approx(expected_priority)

    def test_scoring_function(self):
        root_motif = Motif("ACGT", 0)
        contig_sequence = DNAsequence("ACGTACGTACGT")
        contig_sequences_sample = contig_sequence.sample_n_subsequences(
            2 * 2 + 1, 10
        )
        contig_pssm = contig_sequences_sample.pssm()
        contig_pileup = pl.DataFrame({"position": [1,2,3], "value": [0.1, 0.2, 0.3]})
        methylation_sequences = DNAarray(["ACGT", "CGTA", "GTAC"])
        padding = 2

        # Create MotifSearcher instance
        searcher = MotifSearcher(
            root_motif,
            contig_sequence,
            contig_pssm,
            contig_pileup,
            methylation_sequences,
            padding,
            0.7,
            0.3
        )

        class MockModel:
            def __init__(self, mean_value, std_dev):
                self._mean = mean_value
                self._std_dev = std_dev

            def mean(self):
                return self._mean

            def standard_deviation(self):
                return self._std_dev

        next_model = MockModel(mean_value=0.6, std_dev=0.1)
        current_model = MockModel(mean_value=0.4, std_dev=0.2)

        score = searcher._scoring_function(next_model, current_model)

        mean_diff = next_model.mean() - current_model.mean()  # 0.6 - 0.4 = 0.2
        expected_score = (next_model.mean() * -np.log10(next_model.standard_deviation())) * mean_diff
        expected_score = 0.6 * 1 * 0.2

        assert score == pytest.approx(expected_score)
