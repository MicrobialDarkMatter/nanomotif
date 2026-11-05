import os
import shutil

import polars as pl
import pytest

import nanomotif as nm
import nanomotif.find_motifs_bin as fmb


def _config(tmp_path, pileup, cache, **kwargs):
    return fmb.ProcessorConfig(
        assembly=nm.datasets.geobacillus_plasmids_assembly(),
        pileup_path=pileup,
        bin_contig=nm.datasets.geobacillus_plasmids_bin(), threads=2, search_frame_size=40,
        methylation_threshold_low=0.3, methylation_threshold_high=0.7,
        minimum_kl_divergence=0.05, score_threshold=1.5, log_dir=None, seed=1, output_dir=str(tmp_path / "out"),
        pileup_cache=cache, **kwargs,
    )


def _partitions(strategy):
    """{(bin, mod_type): rows} of a strategy, whichever way they are stored."""
    out = {}
    for key, value in strategy.partitioned_pileup.items():
        if isinstance(value, tuple):
            frame = nm.pileup_cache.read_partition(*value)
        else:
            frame = pl.concat([pl.read_ipc(path) for path in value])
        out[key] = frame.sort("contig", "strand", "position")
    return out


@pytest.fixture
def pileup(tmp_path):
    # a copy, so the test can change it
    path = tmp_path / "pileup.bed"
    shutil.copy(nm.datasets.geobacillus_plasmids_pileup_path(), path)
    return str(path)


@pytest.mark.parametrize("rust", [True, False])
def test_cache_holds_the_same_partitions(tmp_path, pileup, monkeypatch, rust):
    if rust:
        pytest.importorskip("nanomotif_rs")
    else:
        monkeypatch.setattr(fmb, "nanomotif_rs", None)
    reference = fmb.BinMotifProcessorBuilder(_config(tmp_path, pileup, None)).build().task_strategy
    expected = _partitions(reference)
    reference.cleanup()
    cache = str(tmp_path / "cache.arrow")
    got = _partitions(fmb.BinMotifProcessorBuilder(_config(tmp_path, pileup, cache)).build().task_strategy)
    assert set(got) == set(expected)
    for key in expected:
        assert got[key].equals(expected[key].select(got[key].columns).cast(got[key].schema))
    # one file, no partition files left behind
    assert os.listdir(tmp_path / "out" / "temp" / "pileup_partitions") == []


def test_cache_is_reused_and_kept(tmp_path, pileup, monkeypatch):
    cache = str(tmp_path / "cache.arrow")
    first = fmb.BinMotifProcessorBuilder(_config(tmp_path, pileup, cache)).build().task_strategy
    first.cleanup()
    assert os.path.exists(cache)

    def no_pass(*args, **kwargs):
        raise AssertionError("pileup pass should be skipped")
    monkeypatch.setattr(fmb, "nanomotif_rs", None)
    monkeypatch.setattr(fmb.BinMotifProcessorBuilder, "_partition_plain_python", no_pass)
    second = fmb.BinMotifProcessorBuilder(_config(tmp_path, pileup, cache)).build().task_strategy
    assert second.partitioned_pileup == first.partitioned_pileup


@pytest.mark.parametrize("change", ["threshold", "coverage", "pileup", "bins"])
def test_cache_is_rebuilt_when_its_inputs_change(tmp_path, pileup, change, caplog):
    cache = str(tmp_path / "cache.arrow")
    fmb.BinMotifProcessorBuilder(_config(tmp_path, pileup, cache)).build()
    config = _config(tmp_path, pileup, cache)
    if change == "threshold":
        config.methylation_threshold_high = 0.8
    elif change == "coverage":
        config.min_valid_coverage = 10
    elif change == "pileup":
        with open(pileup, "a") as f:
            f.write("")
        os.utime(pileup, ns=(0, 0))
    else:
        config.bin_contig = config.bin_contig.with_columns(pl.lit("other_bin").alias("bin"))
    caplog.set_level("INFO")
    strategy = fmb.BinMotifProcessorBuilder(config).build().task_strategy
    assert "rebuilding it" in caplog.text
    expected = _partitions(fmb.BinMotifProcessorBuilder(_config(tmp_path, pileup, None)).build().task_strategy) \
        if change == "pileup" else None
    if expected is not None:
        assert _partitions(strategy).keys() == expected.keys()


def test_unreadable_cache_is_rebuilt(tmp_path, pileup, caplog):
    cache = tmp_path / "cache.arrow"
    cache.write_text("not an arrow file")
    caplog.set_level("INFO")
    strategy = fmb.BinMotifProcessorBuilder(_config(tmp_path, pileup, str(cache))).build().task_strategy
    assert "not readable" in caplog.text
    assert len(_partitions(strategy)) > 0
