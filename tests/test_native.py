import pytest
import polars as pl
import polars.testing
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
        a = pl.concat([pl.read_ipc(f) for f in expected[part]]).sort(key)
        b = pl.concat([pl.read_ipc(f) for f in got[part]]).sort(key)
        assert a.drop("fraction_mod").equals(b.drop("fraction_mod"))
        assert (a["fraction_mod"] - b["fraction_mod"]).abs().max() < 1e-12


def test_rust_partitioning_of_gzip_matches_plain(tmp_path):
    """A bgzipped pileup is partitioned to the same rows as the plain one."""
    contig_bin = dict(nm.datasets.geobacillus_plasmids_bin().iter_rows())
    plain = nm.datasets.geobacillus_plasmids_pileup_path()
    (tmp_path / "plain").mkdir()
    (tmp_path / "gz").mkdir()
    a = nanomotif_rs.partition_pileup(plain, contig_bin, str(tmp_path / "plain"), min_coverage=5, threshold=0.7, threads=2)
    b = nanomotif_rs.partition_pileup(plain + ".gz", contig_bin, str(tmp_path / "gz"), min_coverage=5, threshold=0.7, threads=2)
    assert set(a) == set(b)
    key = ["contig", "strand", "position", "mod_type"]
    for part in a:
        x = pl.concat([pl.read_ipc(f) for f in a[part]]).sort(key)
        y = pl.concat([pl.read_ipc(f) for f in b[part]]).sort(key)
        assert x.equals(y)


def test_bgzip_reader_without_rust_matches_plain():
    """The pysam reader (used without nanomotif_rs) returns the rows of load_pileup."""
    path = nm.datasets.geobacillus_plasmids_pileup_path()
    expected = nm.dataload.load_pileup(path).filter(pl.col("contig") == "contig_2")
    got = nm.dataload.load_contigs_pileup_bgzip(path + ".gz", ["contig_2"])
    polars.testing.assert_frame_equal(got.sort("position", "mod_type", "strand"), expected.sort("position", "mod_type", "strand"))


def test_column_bitsets_match_packbits_layout():
    import numpy as np
    rng = np.random.default_rng(1)
    windows = rng.choice(np.frombuffer(b"ACGTN", np.uint8), size=(1003, 41))
    windows[5, 3] = 0
    codes, bits = nanomotif_rs.column_bitsets(windows)
    for column in range(41):
        for i, code in enumerate(codes):
            assert np.array_equal(bits[column, i], nm.find_motifs_bin._pack_bits(windows[:, column] == code))
