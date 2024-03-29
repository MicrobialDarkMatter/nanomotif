from nanomotif.constants import *
import nanomotif as nm
import pytest
from hypothesis import given, strategies as st
import itertools
import pandas as pd
import polars as pl
import polars.testing
import numpy as np


class TestMethylatedMotifOccourances:
    def test_known_motif_sequence(self):
        motif = nm.candidate.Motif('ACG', 0)
        sequence = 'TACGGACGCCACG'
        methylated_positions = np.array([1, 5])

        expected_methylated = np.array([1, 5])
        expected_nonmethylated = np.array([10])
        
        result = nm.evaluate.methylated_motif_occourances(motif, sequence, methylated_positions)
        np.testing.assert_array_equal(result[0], expected_methylated)
        np.testing.assert_array_equal(result[1], expected_nonmethylated)

class TestProcessSubpileup:
    def find_motifs(self):
        pileup_path = nm.datasets.geobacillus_plasmids_pileup_path()
        pileup = nm.dataload.load_pileup(pileup_path, min_fraction=0.5)
        modtype = "a"
        subpileup = pileup.pileup.filter(pl.col("mod_type") == modtype)

        assembly = nm.datasets.geobacillus_plasmids_assembly()
        contig = "contig_3"
        min_kl_divergence = 0.1
        padding = 20
        minimum_methylation_fraction_confident = 0.75

        result = nm.evaluate.process_subpileup(contig, modtype ,subpileup, assembly, min_kl_divergence, padding, minimum_methylation_fraction_confident)
        return result

    def test_no_errors(self):
        self.find_motifs()
    
    def test_correct_output(self):
        result = self.find_motifs()
        expected_output = pl.DataFrame({
            "sequence":[
                "...................GATC..................", 
                "................ACCCA....................", 
                "................CCAAAT...................", 
                "...............G[AG].GAAG[TC].................."],
            "score":[1.687542, 1.475938, 1.388421, 0.572978],
            "contig":["contig_3", "contig_3", "contig_3", "contig_3"],
            "mod_type":["a", "a", "a", "a"],
            "alpha":[706, 65, 80, 84],
            "beta":[3, 1, 1, 1]
        })
        polars.testing.assert_frame_equal(result, expected_output)


class TestProcessSampleParallel():
    def run_function(self):
        pileup_path = nm.datasets.geobacillus_plasmids_pileup_path()
        pileup = nm.dataload.load_pileup(pileup_path, min_fraction=0.5)
        pileup = pileup.pileup
        assembly = nm.datasets.geobacillus_plasmids_assembly()

        return nm.evaluate.process_sample_parallel(assembly, pileup, threads=2, seed=1)
    def test_no_errors(self):
        self.run_function()

    def test_reprocuible_result(self):
        res1 = self.run_function().drop("model").sort(["contig", "mod_type", "motif"])
        res2 = self.run_function().drop("model").sort(["contig", "mod_type", "motif"])
        polars.testing.assert_frame_equal(res1, res2)

