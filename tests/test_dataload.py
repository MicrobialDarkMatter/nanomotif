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
        assert isinstance(pileup.pileup, pl.DataFrame)

    def test_column_types(self):
        pileup = nm.datasets.geobacillus_plasmids_pileup()
        assert pileup.pileup["contig"].dtype == pl.Utf8
        assert pileup.pileup["position"].dtype == pl.Int64
        assert pileup.pileup["mod_type"].dtype == pl.Utf8
        assert pileup.pileup["strand"].dtype == pl.Utf8
        assert pileup.pileup["fraction_mod"].dtype == pl.Float64
        assert pileup.pileup["Nvalid_cov"].dtype == pl.Int64

    def test_exotic_seq_name_format(self):
        pileup = nm.load_pileup("nanomotif/datasets/test_pileup.bed", min_coverage=5, min_fraction=0.5)
        assert pileup.pileup["contig"].dtype == pl.Utf8
        assert pileup.pileup["position"].dtype == pl.Int64
        assert pileup.pileup["mod_type"].dtype == pl.Utf8
        assert pileup.pileup["strand"].dtype == pl.Utf8
        assert pileup.pileup["fraction_mod"].dtype == pl.Float64
        assert pileup.pileup["Nvalid_cov"].dtype == pl.Int64

    def test_fraction_filter(self):
        pileup_path = nm.datasets.geobacillus_plasmids_pileup_path()
        pileup = nm.dataload.load_pileup(pileup_path, min_fraction=0.5, min_coverage=5)
        assert (pileup.pileup["fraction_mod"] > 0.5).all()
        # Assert that the filtered pileup is a subset of the complete pileup
        pileup_complete = nm.dataload.load_pileup(pileup_path, min_fraction=0.5, min_coverage=5)
        assert len(pileup.pileup.join(pileup_complete.pileup, on=pileup.pileup.columns, how="inner")) == len(pileup.pileup)


