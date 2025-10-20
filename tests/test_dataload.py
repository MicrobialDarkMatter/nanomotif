from nanomotif.motif import *
from nanomotif.constants import *
import nanomotif as nm
import pytest
from hypothesis import given, strategies as st
import itertools
import pandas as pd
import polars as pl
import polars.testing


class TestPileupLoad:
    def test_no_error(self):
        pileup = nm.datasets.geobacillus_plasmids_pileup()
        assert isinstance(pileup, pl.DataFrame)

    def test_column_types(self):
        pileup = nm.datasets.geobacillus_plasmids_pileup()
        assert pileup["contig"].dtype == pl.Utf8
        assert pileup["position"].dtype == pl.Int64
        assert pileup["mod_type"].dtype == pl.Utf8
        assert pileup["strand"].dtype == pl.Utf8
        assert pileup["fraction_mod"].dtype == pl.Float64
        assert pileup["Nvalid_cov"].dtype == pl.Int64

    def test_exotic_seq_name_format(self):
        pileup = nm.load_pileup("nanomotif/datasets/test_pileup.bed")
        assert pileup["contig"].dtype == pl.Utf8
        assert pileup["position"].dtype == pl.Int64
        assert pileup["mod_type"].dtype == pl.Utf8
        assert pileup["strand"].dtype == pl.Utf8
        assert pileup["fraction_mod"].dtype == pl.Float64
        assert pileup["Nvalid_cov"].dtype == pl.Int64



def test_filter_pileup_adjacency_filter():
    pileup = pl.DataFrame({
        "contig": ["contig1"] * 10,
        "position": list(range(10)),
        "mod_type": ["m6A"] * 10,
        "strand": ["+"] * 10,
        "fraction_mod": [0.8, 
                         0.9, 0.1, 0.95, 
                         0.85, 
                         0.2, 
                         0.75, 
                         0.9, 0.05, 
                         0.8],
        "Nvalid_cov": [10] * 10
    })
    filtered = nm.dataload.filter_pileup_adjacency_filter(pileup, methylation_threshold=0.7, adjacency_distance=1)
    expected_positions = [1, 2, 3, 5, 7, 8, 9]
    assert filtered["position"].to_list() == expected_positions

def test_filter_pileup_adjacency_filter_more_groups():
    pileup = pl.DataFrame({ 
        "contig": ["contig1"] * 5 + ["contig2"] * 5,
        "position": list(range(5)) + list(range(5)),
        "mod_type": ["m6A", "5mC"] * 5,
        "strand": ["+"] * 5 + ["-"] * 5,
        "fraction_mod": [0.8, 0.9, 0.1, 0.95, 0.85, 
                         0.2, 0.75, 0.9, 0.05, 0.8],
        "Nvalid_cov": [10] * 10
    })
    filtered = nm.dataload.filter_pileup_adjacency_filter(pileup, methylation_threshold=0.7, adjacency_distance=1)
    expected_positions_contig1 = [1, 2, 3]
    expected_positions_contig2 = [0,2,3,4]
    assert filtered.filter(pl.col("contig") == "contig1")["position"].to_list() == expected_positions_contig1
    assert filtered.filter(pl.col("contig") == "contig2")["position"].to_list() == expected_positions_contig2