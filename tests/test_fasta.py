import os
import tempfile
import pytest
from pathlib import Path

import polars as pl
from nanomotif.fasta import generate_contig_bin, FastaReader


class DummyArgs:
    def __init__(self, contig_bin=None, files=None, directory=None, extension=None, out=None):
        self.contig_bin = contig_bin
        self.files = files
        self.directory = directory
        self.extension = extension
        self.out = out


# -------------------------------
# Test: Input from contig_bin file
# -------------------------------
def test_generate_contig_bin_from_tsv(tmp_path):
    input_path = tmp_path / "input.tsv"
    input_path.write_text("contigA\tbin1\ncontigB\tbin2\n")

    args = DummyArgs(contig_bin=str(input_path))
    df = generate_contig_bin(args)

    assert df.shape == (2, 2)
    assert set(df.columns) == {"contig", "bin"}
    assert "contigA" in df["contig"].to_list()
    assert "bin1" in df["bin"].to_list()


# -------------------------------
# Test: Input from list of files
# -------------------------------
def test_generate_contig_bin_from_files(tmp_path):
    tmp_path.mkdir(exist_ok=True)
    fasta1 = tmp_path / "a.fasta"
    fasta1.write_text(">c1\nATGC\n>c2\nGGCC\n")

    fastq = tmp_path / "b.fastq"
    fastq.write_text(">c3\nAATT")

    out_dir = tmp_path / "out"
    out_dir.mkdir()

    args = DummyArgs(files=[str(fasta1), str(fastq)], out=str(out_dir))
    df = generate_contig_bin(args)

    # Verify file output
    output_path = out_dir / "temp" / "contig_bin.tsv"
    assert output_path.exists()

    with open(output_path) as f:
        lines = f.read().splitlines()
    assert "c1\ta" in lines
    assert "c3\tb" in lines

    # Verify DataFrame
    assert df.shape == (3, 2)
    assert set(df.columns) == {"contig", "bin"}
    assert "c2" in df["contig"].to_list()


# -------------------------------
# Test: Input from directory with extension
# -------------------------------
def test_generate_contig_bin_from_directory(tmp_path):
    tmp_path.mkdir(exist_ok=True)
    (tmp_path / "one.fasta").write_text(">alpha\nATGCGC\n>beta\nTTAAGG\n")
    (tmp_path / "two.fasta").write_text("@gamma\nGGTT\n+\nIIII\n")

    out_dir = tmp_path / "results"
    out_dir.mkdir()

    args = DummyArgs(directory=str(tmp_path), extension=".fasta", out=str(out_dir))
    df = generate_contig_bin(args)

    output_file = out_dir / "temp" / "contig_bin.tsv"
    assert output_file.exists()

    content = output_file.read_text()
    assert "alpha" in content
    assert "beta" in content
    assert "gamma" in content

    assert df.shape == (3, 2)
    assert set(df.columns) == {"contig", "bin"}
    assert "one" in df["bin"].to_list()
    assert "two" in df["bin"].to_list()
