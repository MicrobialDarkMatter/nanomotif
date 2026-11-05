import polars as pl
import pytest

# Monkeypatch nm.motif.reverse_compliment inside the test
import nanomotif as nm  # <-- change to actual module path

from nanomotif.motif import MotifSearchResult
from nanomotif.postprocess import join_motif_complements
from nanomotif.model import BetaBernoulliModel


def make_motif_df():
    # Build a tiny DataFrame with two motifs that are complements of each other
    return 


def test_join_motif_complements_basic():
    df = MotifSearchResult(pl.DataFrame({
        "reference": ["ref1", "ref1"],
        "motif": ["AAGGTT", "AACCTT"],
        "mod_type": ["m", "m"],
        "mod_position": [0, 0],
        "score": [10.0, 20.0],
        "model": [BetaBernoulliModel(), BetaBernoulliModel()]  # dummy objects just to satisfy schema
    }))

    result = join_motif_complements(df)
    print(result)

    assert isinstance(result, MotifSearchResult)
    assert "motif_complement" in result.columns
    assert any(c.endswith("_complement") for c in result.columns)
    assert len(result) >= 1
    assert result["motif"].to_list() == ["AAGGTT"]
    assert result["motif_complement"].to_list() == ["AACCTT"]

def test_join_motif_complements_no_complements():
    df = MotifSearchResult(pl.DataFrame({
        "reference": ["ref1", "ref1"],
        "motif": ["GATCC", "GATCG"],
        "mod_type": ["m", "m"],
        "mod_position": [0, 0],
        "score": [10.0, 20.0],
        "model": [BetaBernoulliModel(), BetaBernoulliModel()]  # dummy objects just to satisfy schema
    }))

    result = join_motif_complements(df)
    print(result)

    assert "motif_complement" in result.columns
    assert any(c.endswith("_complement") for c in result.columns)
    assert len(result) == 2
    assert set(result["motif"].to_list()) == {"GATCC", "GATCG"}
    assert set(result["motif_complement"].to_list()) == {None, None}

def test_join_motif_complements_palindrome():
    df = MotifSearchResult(pl.DataFrame({
        "reference": ["ref1", "ref1"],
        "motif": ["GCGC", "GATC"],
        "mod_type": ["m", "m"],
        "mod_position": [1, 2],
        "score": [10.0, 20.0],
        "model": [BetaBernoulliModel(), BetaBernoulliModel()]  # dummy objects just to satisfy schema
    }))

    result = join_motif_complements(df)
    print(result)

    assert "motif_complement" in result.columns
    assert any(c.endswith("_complement") for c in result.columns)
    assert len(result) == 2
    assert set(result["motif"].to_list()) == {"GCGC", "GATC"}
    assert set(result["motif_complement"].to_list()) == {"GCGC", "GATC"}

def test_join_motif_complements_two_palindromes():
    df = MotifSearchResult(pl.DataFrame({
        "reference": ["ref1", "ref1"],
        "motif": ["GCGC", "GCGC"],
        "mod_type": ["m", "m"],
        "mod_position": [1, 3],
        "score": [10.0, 20.0],
        "model": [BetaBernoulliModel(), BetaBernoulliModel()]  # dummy objects just to satisfy schema
    }))

    result = join_motif_complements(df)
    print(result)

    assert "motif_complement" in result.columns
    assert any(c.endswith("_complement") for c in result.columns)
    # GCGC methylated at position 1 and at position 3 are two motifs; a palindrome is its own
    # complement, so each is paired with itself (not also with the other: that gave 4 crossed rows)
    assert len(result) == 2
    assert result["motif"].to_list() == ["GCGC", "GCGC"]
    assert result["motif_complement"].to_list() == ["GCGC", "GCGC"]
    assert result["mod_position_iupac"].to_list() == result["mod_position_iupac_complement"].to_list()


def test_join_motif_complements_mruber():
    df = MotifSearchResult(pl.DataFrame({
        "reference": ["Meiothermus-ruber_DSMZ1279"] * 9,
        "motif": [
            "AATT",
            "GATC",
            "CCA......TGCC",
            "CAGACG..G",
            "GGCA......TGG",
            "GGGAGC",
            "TTAA",
            "CTCGAG",
            "GCAGATG"
        ],
        "mod_type": ["a"] * 9,
        "mod_position": [1, 1, 2, 3, 3, 3, 3, 4, 4],
        "score": [36.11927388312779, 3.969322870370931, 17.394075780998488, 16.787817553222727,
                  15.420804062089653, 23.557242011998742, 49.514246510667256, 8.33360785017712,
                  5.768720177536614],
        "model": [BetaBernoulliModel(5944 + 5, 1 + 5), BetaBernoulliModel(14643 + 5, 3620 + 5),
                  BetaBernoulliModel(1203 + 5, 0 + 5), BetaBernoulliModel(38 + 5, 4 + 5),
                  BetaBernoulliModel(1240 + 5, 0 + 5), BetaBernoulliModel(1816 + 5, 0 + 5),
                  BetaBernoulliModel(4439 + 5, 2 + 5), BetaBernoulliModel(14089 + 5, 30 + 5),
                  BetaBernoulliModel(150 + 5, 39 + 5)]
    }))

    result = join_motif_complements(df)
    print(result)

    assert isinstance(result, MotifSearchResult)
    assert "motif_complement" in result.columns
    assert any(c.endswith("_complement") for c in result.columns)
    assert len(result) >= 8
    assert result["motif"].to_list() == [
        "AATT",
        "GATC",
        "CAGACG..G",
        "GGCA......TGG",
        "GGGAGC",
        "TTAA",
        "CTCGAG",
        "GCAGATG"
    ]
    assert result["motif_complement"].to_list() == [
        "AATT",
        "GATC",
        None,
        "CCA......TGCC",
        None,
        "TTAA",
        "CTCGAG",
        None
    ]
def test_remove_noisy_motifs():
    df = MotifSearchResult(pl.DataFrame({
        "reference": ["ref1", "ref1", "ref1"],
        "motif": ["AAGGTT", "AACCTT", "GATCC"],
        "mod_type": ["m", "m", "m"],
        "mod_position": [0, 0, 0],
        "score": [10.0, 20.0, 5.0],
        "model": [BetaBernoulliModel(50, 100), BetaBernoulliModel(60, 90), BetaBernoulliModel(5, 200)]
    }))

    result = nm.postprocess.remove_noisy_motifs(df)
    print(result)

    assert len(result) == 3
    assert isinstance(result, MotifSearchResult)
    assert set(result["motif"].to_list()) == {"AAGGTT", "AACCTT", "GATCC"}

def test_bin_motifs_has_no_duplicates_or_crossed_palindrome_complements(tmp_path):
    """A palindrome methylated at two positions (TTAA at either A) is two motifs, each its own
    complement; motifs listed twice in the search result are written once."""
    from nanomotif.model import BetaBernoulliModel
    pad = "." * 5
    motifs = nm.motif.MotifSearchResult(pl.DataFrame({
        "motif": [pad + "TTAA" + pad, pad + "TTAA" + pad, pad + "GATC" + pad, pad + "GATC" + pad],
        "mod_position": [7, 8, 6, 6],
        "reference": ["bin1"] * 4,
        "mod_type": ["a"] * 4,
        "model": [BetaBernoulliModel(), BetaBernoulliModel(), BetaBernoulliModel(), BetaBernoulliModel()],
        "score": [2.0, 2.0, 2.0, 2.0],
    }))
    out = tmp_path / "bin-motifs.tsv"
    nm.motif.MotifSearchResult(nm.postprocess.join_motif_complements(motifs)).write_motif_formatted(str(out))
    written = pl.read_csv(out, separator="\t")
    rows = sorted(written.select("motif", "mod_position", "motif_complement", "mod_position_complement").rows())
    assert rows == [("GATC", 1, "GATC", 1), ("TTAA", 2, "TTAA", 2), ("TTAA", 3, "TTAA", 3)]
