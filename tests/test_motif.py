import polars as pl
import pytest
import inspect
import pickle
from pathlib import Path
import nanomotif as nm
from nanomotif.motif import MotifSearchResult, Motif
from nanomotif.model import BetaBernoulliModel

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




# --- Import the class under test ---
from nanomotif.motif import MotifSearchResult  # adjust path as needed


# ---------------- Tests ----------------

def make_minimal_df():
    """Helper to create a DataFrame with required columns and a model object."""
    model = BetaBernoulliModel(5, 5)
    model.update(50, 100)  # make sure model has non-prior values
    return pl.DataFrame({
        "reference": ["contig1"],
        "motif": ["GATC"],
        "mod_type": ["m6A"],
        "mod_position": [2],
        "model": [model],
        "score": [1.23],
    })

def make_minimal_df_complement():
    """Helper to create a DataFrame with required columns and a model object."""
    model = BetaBernoulliModel(5, 5)
    model.update(50, 100)  # make sure model has non-prior values
    model_complement = BetaBernoulliModel(5, 5)
    model_complement.update(30, 70)
    return pl.DataFrame({
        "reference": ["contig1"],
        "motif": ["GATC"],
        "mod_type": ["m6A"],
        "mod_position": [2],
        "model": [model],
        "score": [1.23],
        "n_mod_complement": [30],
        "n_nomod_complement": [70],
        "motif_complement": ["CTAG"],
        "motif_iupac_complement": ["CTAG"],
        "score_complement": [0.98],
        "mod_position_complement": [3],
        "mod_position_iupac_complement": [3],
        "model_complement": [model_complement],
    })

def test_init_from_dataframe_normalizes_model():
    df = make_minimal_df()
    # Polars may keep model as Object already
    msr = MotifSearchResult(df)
    assert isinstance(msr, MotifSearchResult)
    assert msr["model"].dtype == pl.Object
    assert isinstance(msr["model"].to_list()[0], BetaBernoulliModel)


def test_init_from_struct_like_model_dicts():
    # Simulate worker returning struct-like values (dicts) instead of objects
    df = pl.DataFrame({
        "reference": ["contig1"],
        "motif": ["GATC"],
        "mod_type": ["m6A"],
        "mod_position": [2],
        "model": [{"_alpha": 5, "_beta": 10, "_alpha_prior": 1, "_beta_prior": 1}],
        "score": [1.23],
    })
    msr = MotifSearchResult(df)
    val = msr["model"].to_list()[0]
    assert isinstance(val, BetaBernoulliModel)
    assert val._alpha == 5 and val._beta == 10


def test_pickle_and_unpickle_roundtrip(tmp_path):
    df = make_minimal_df()
    msr = MotifSearchResult(df)

    pkl = tmp_path / "msr.pkl"
    with open(pkl, "wb") as fh:
        pickle.dump(msr, fh)

    with open(pkl, "rb") as fh:
        msr2 = pickle.load(fh)

    assert isinstance(msr2, MotifSearchResult)
    assert msr2["model"].dtype == pl.Object
    restored_model = msr2["model"].to_list()[0]
    assert isinstance(restored_model, BetaBernoulliModel)
    assert restored_model._alpha == 55


def test_concat_with_mixed_model_types():
    # one with real object
    df1 = make_minimal_df()
    msr1 = MotifSearchResult(df1)

    # one with dict -> should normalize
    df2 = pl.DataFrame({
        "reference": ["contig2"],
        "motif": ["CCGG"],
        "mod_type": ["m5C"],
        "mod_position": [1],
        "model": [{"_alpha": 2, "_beta": 8, "_alpha_prior": 1, "_beta_prior": 1}],
        "score": [0.42],
    })
    msr2 = MotifSearchResult(df2)

    concatenated = pl.concat([msr1, msr2])
    assert concatenated["model"].dtype == pl.Object
    vals = concatenated["model"].to_list()
    assert all(isinstance(v, BetaBernoulliModel) for v in vals if v is not None)


def test_write_motifs_drops_object(tmp_path):
    df = make_minimal_df()
    msr = MotifSearchResult(df)

    out = tmp_path / "motifs.tsv"
    msr.write_motifs(out)

    content = out.read_text().strip().splitlines()
    header = content[0].split("\t")
    assert "model" not in header  # object column dropped


def test_derived_columns_computed():
    df = make_minimal_df()
    msr = MotifSearchResult(df)

    # n_mod and n_nomod derived from model priors
    assert "n_mod" in msr.columns
    assert "n_nomod" in msr.columns
    assert msr["n_mod"].to_list()[0] == 50
    assert msr["n_nomod"].to_list()[0] == 100


def test_column_order_is_respected():
    df = make_minimal_df()
    msr = MotifSearchResult(df)
    required_order = list(MotifSearchResult.REQUIRED_COLUMNS.keys()) + list(MotifSearchResult.DERIVED_COLUMNS.keys())
    assert all(col in msr.columns for col in required_order)
    # Check that required columns come first
    assert msr.columns[:len(required_order)] == required_order



def test_derive_complement_counts_from_model_complement():
    # Create DF with model_complement but no n_mod_complement/n_nomod_complement
    model = BetaBernoulliModel(alpha=10, beta=10)
    model.update(15, 25)
    model_complement = BetaBernoulliModel(alpha=10, beta=10)
    model_complement.update(20, 30)
    df = pl.DataFrame({
        "reference": ["Enterococcus_faecalis_complete_genome"],
        "motif": ["C[AG]AA......[AG]TTG"],
        "mod_type": ["a"],
        "mod_position": [3],
        "score": [16.22159133144396],
        "model": [model],
        # model_complement provided instead of n_mod_complement/n_nomod_complement
        "motif_iupac": ["CRAANNNNNNRTTG"],
        "mod_position_iupac": [3],
        "motif_complement": ["CAA[CT]......TT[CT]G"],
        "motif_iupac_complement": ["CTAG"],
        "mod_position_iupac_complement": [2],
        "mod_position_complement": [2],
        "model_complement": [model_complement],
        "score_complement": [16.137373428650985],
    })

    msr = MotifSearchResult(df)

    # Ensure complement-derived columns are created
    assert "n_mod_complement" in msr.columns
    assert "n_nomod_complement" in msr.columns

    model = msr["model_complement"].to_list()[0]
    expected_n_mod = model._alpha - model._alpha_prior
    expected_n_nomod = model._beta - model._beta_prior

    # Verify correct derivation
    assert msr["n_mod_complement"].to_list()[0] == expected_n_mod
    assert msr["n_nomod_complement"].to_list()[0] == expected_n_nomod

    # Also check non-complement data is intact
    assert msr["motif"].to_list()[0] == "C[AG]AA......[AG]TTG"
    assert msr["motif_complement"].to_list()[0] == "CAA[CT]......TT[CT]G"

def test_pickle_and_unpickle_complement_model():
    df = make_minimal_df_complement()
    msr = MotifSearchResult(df)

    pkl = "msr_complement.pkl"
    with open(pkl, "wb") as fh:
        pickle.dump(msr, fh)

    with open(pkl, "rb") as fh:
        msr2 = pickle.load(fh)

    assert isinstance(msr2, MotifSearchResult)
    assert msr2["model"].dtype == pl.Object
    assert msr2["model_complement"].dtype == pl.Object
    restored_model = msr2["model"].to_list()[0]
    restored_model_complement = msr2["model_complement"].to_list()[0]
    assert isinstance(restored_model, BetaBernoulliModel)
    assert isinstance(restored_model_complement, BetaBernoulliModel)
    assert restored_model._alpha == 55
    assert restored_model_complement._alpha == 35

def test_column_order_is_enforced():
    model = BetaBernoulliModel(5, 5)
    model.update(50, 100)  # make sure model has non-prior values
    df = pl.DataFrame({
        "score": [1.23],
        "motif": ["GATC"],
        "reference": ["contig1"],
        "model": [model],
        "mod_position": [2],
        "mod_type": ["m6A"],
    })
    msr = MotifSearchResult(df)
    required_order = list(MotifSearchResult.REQUIRED_COLUMNS.keys()) + list(MotifSearchResult.DERIVED_COLUMNS.keys())
    assert all(col in msr.columns for col in required_order)
    # Check that required columns come first
    assert msr.columns[:len(required_order)] == required_order

def test_column_order_is_enforced_complementary():
    model = BetaBernoulliModel(5, 5)
    model.update(50, 100)  # make sure model has non-prior values
    df = pl.DataFrame({
        "score": [1.23],
        "motif": ["GATC"],
        "reference": ["contig1"],
        "model": [model],
        "mod_position": [2],
        "mod_type": ["m6A"],
        "n_mod_complement": [30],
        "motif_complement": ["CTAG"],
        "motif_iupac_complement": ["CTAG"],
        "score_complement": [0.98],
        "mod_position_complement": [3],
        "model_complement": [model],
        "n_nomod_complement": [70],
        "mod_position_iupac_complement": [3],
    })
    msr = MotifSearchResult(df)
    required_order = list(MotifSearchResult.REQUIRED_COLUMNS.keys()) + list(MotifSearchResult.DERIVED_COLUMNS.keys()) + list(MotifSearchResult.COMPLEMENTARY_COLUMNS.keys())
    assert all(col in msr.columns for col in required_order)
    # Check that required columns come first
    assert msr.columns[:len(required_order)] == required_order
