import pytest
import polars as pl
import nanomotif as nm

nanomotif_rs = pytest.importorskip("nanomotif_rs")


def test_rust_partitioning_matches_python(tmp_path):
    """nanomotif_rs.partition_pileup gives the same rows as the pure-Python plain-pileup path."""
    contig_bin = nm.datasets.geobacillus_plasmids_bin()
    config = nm.find_motifs_bin.ProcessorConfig(
        assembly=nm.datasets.geobacillus_plasmids_assembly(),
        pileup_path=nm.datasets.geobacillus_plasmids_pileup_path(),
        bin_contig=contig_bin, threads=2, search_frame_size=40,
        methylation_threshold_low=0.3, methylation_threshold_high=0.7,
        minimum_kl_divergence=0.05, score_threshold=1.5, log_dir=None, seed=1, output_dir=str(tmp_path),
    )
    builder = nm.find_motifs_bin.BinMotifProcessorBuilder(config)
    (tmp_path / "py").mkdir()
    (tmp_path / "rs").mkdir()
    expected = builder._partition_plain_python(config.pileup_path, str(tmp_path / "py"))
    got = nanomotif_rs.partition_pileup(
        config.pileup_path, dict(contig_bin.iter_rows()), str(tmp_path / "rs"), min_coverage=5, threshold=0.7, threads=2
    )
    assert set(got) == set(expected)
    key = ["contig", "strand", "position", "mod_type"]
    for part in expected:
        a = pl.concat([pl.read_ipc(f, memory_map=False) for f in expected[part]]).sort(key)
        b = pl.concat([pl.read_ipc(f, memory_map=False) for f in got[part]]).sort(key)
        assert a.drop("fraction_mod").equals(b.drop("fraction_mod"))
        assert (a["fraction_mod"] - b["fraction_mod"]).abs().max() < 1e-12
