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



