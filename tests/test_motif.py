import io
import polars as pl
import pytest
from pathlib import Path
from nanomotif.motif import MotifSearchResult
import nanomotif as nm
from nanomotif.motif import Motif
from nanomotif.model import BetaBernoulliModel

import inspect
print("MSR is class? ", inspect.isclass(MotifSearchResult), "repr:", repr(MotifSearchResult))


def make_base_df():
    return pl.DataFrame({
        "reference": ["r1", "r2"],
        "motif": ["CCAGG", "CCTGG"],
        "mod_type": ["m", "m"],
        "mod_position": [1, 2],
        "model": [BetaBernoulliModel(50, 100), BetaBernoulliModel(500, 50)],
        "score": [1.0, 2.0],
    })


def test_initialization_with_valid_data():
    df = make_base_df()
    result = MotifSearchResult(df)

    # Should contain all required columns
    for col in MotifSearchResult.REQUIRED_COLUMNS:
        assert col in result.columns

    # Derived columns should be created
    for col in ["n_mod", "n_nomod", "motif_iupac", "mod_position_iupac"]:
        assert col in result.columns


def test_rename_contig_to_reference():
    df = make_base_df().rename({"reference": "contig"})
    result = MotifSearchResult(df)
    assert "reference" in result.columns
    assert "contig" not in result.columns


def test_rename_bin_to_reference():
    df = make_base_df().rename({"reference": "bin"})
    result = MotifSearchResult(df)
    assert "reference" in result.columns
    assert "bin" not in result.columns


def test_missing_required_column_raises():
    df = make_base_df().drop("motif")
    with pytest.raises(ValueError, match="Missing required columns"):
        MotifSearchResult(df)


def test_derived_column_values():
    df = make_base_df()
    result = MotifSearchResult(df)
    # Derived n_mod and n_nomod should match model attributes
    expected_n_mod = [m._alpha - m._alpha_prior for m in df["model"]]
    expected_n_nomod = [m._beta - m._beta_prior for m in df["model"]]
    assert result["n_mod"].to_list() == expected_n_mod
    assert result["n_nomod"].to_list() == expected_n_nomod


def test_load_motifs_returns_motif_objects():
    df = make_base_df()
    result = MotifSearchResult(df)
    motifs = result._load_motifs()
    assert len(motifs) == len(df)
    assert all(isinstance(m, Motif) for m in motifs)


def test_sort_returns_same_class():
    df = make_base_df()
    result = MotifSearchResult(df)
    sorted_result = result.sort("score")
    assert isinstance(sorted_result, MotifSearchResult)


def test_write_motifs(tmp_path: Path):
    df = make_base_df()
    result = MotifSearchResult(df)
    output_file = tmp_path / "motifs.tsv"
    result.write_motifs(output_file)

    # File should exist and contain tab-separated data
    content = output_file.read_text()
    assert "reference" in content
    assert "motif" in content


def test_hstack_does_not_break_class():
    df = make_base_df()
    result = MotifSearchResult(df)
    result = result.with_columns(pl.Series("extra", [1, 2]))
    assert isinstance(result, MotifSearchResult)
    assert "extra" in result.columns