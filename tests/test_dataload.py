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

def test_filter_pileup_min_coverage_is_inclusive():
    pileup = pl.DataFrame({
        "contig": ["contig1"] * 4,
        "position": [0, 1, 2, 3],
        "mod_type": ["a"] * 4,
        "strand": ["+"] * 4,
        "fraction_mod": [0.9] * 4,
        "Nvalid_cov": [4, 5, 6, 20]
    })
    assert nm.dataload.filter_pileup(pileup, min_coverage=5)["position"].to_list() == [1, 2, 3]
    assert nm.dataload.filter_pileup(pileup, min_coverage=6)["position"].to_list() == [2, 3]



def test_filter_pileup_minimummod_frequency_uses_threshold():
    # 60 positions at 0.75: methylated at threshold 0.7 (inclusive), not at 0.8
    pileup = pl.DataFrame({
        "contig": ["contig1"] * 60,
        "position": list(range(60)),
        "mod_type": ["a"] * 60,
        "strand": ["+"] * 60,
        "fraction_mod": [0.75] * 60,
        "Nvalid_cov": [10] * 60
    })
    assert len(nm.dataload.filter_pileup_minimummod_frequency(pileup, methylation_threshold=0.75)) == 60
    assert len(nm.dataload.filter_pileup_minimummod_frequency(pileup, methylation_threshold=0.8)) == 0


def test_load_pileup_filters_while_reading():
    path = nm.datasets.geobacillus_plasmids_pileup_path()
    full = nm.dataload.load_pileup(path)
    expected = nm.dataload.filter_pileup(full, min_coverage=8).filter(pl.col("contig") == "contig_2")
    got = nm.dataload.load_pileup(path, min_coverage=8, contigs=["contig_2"])
    polars.testing.assert_frame_equal(got, expected)


def test_filter_pileup_adjacency_filter_matches_naive():
    import numpy as np
    rng = np.random.default_rng(5)
    n = 400
    pileup = pl.DataFrame({
        "contig": rng.choice(["c1", "c2"], n),
        "position": rng.integers(0, 150, n),
        "mod_type": rng.choice(["a", "m"], n),
        "strand": rng.choice(["+", "-"], n),
        "fraction_mod": rng.choice([0.1, 0.5, 0.7, 0.8, 0.9, 1.0], n),
        "Nvalid_cov": [10] * n,
    })
    rows = pileup.rows(named=True)
    expected = [
        r for r in rows
        if r["fraction_mod"] < 0.7 or r["fraction_mod"] == max(
            o["fraction_mod"] for o in rows
            if o["contig"] == r["contig"] and o["strand"] == r["strand"] and abs(o["position"] - r["position"]) <= 3
        )
    ]
    got = nm.dataload.filter_pileup_adjacency_filter(pileup, methylation_threshold=0.7, adjacency_distance=3)
    key = lambda r: (r["contig"], r["strand"], r["position"], r["mod_type"], r["fraction_mod"])
    assert sorted(got.rows(named=True), key=key) == sorted(expected, key=key)
