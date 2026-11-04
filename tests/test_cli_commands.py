import pytest
import subprocess
import pandas as pd
import os
import glob
import shutil



def test_motif_discovery():
    """
    """
    outdir = "tests/cli_test_motif_discovery"
    
    cmd = [
        "nanomotif", "motif_discovery",
        "-t", "1",
        "nanomotif/datasets/geobacillus-plasmids.assembly.fasta",
        "nanomotif/datasets/geobacillus-plasmids.pileup.bed",
        "-c", "nanomotif/datasets/geobacillus-contig-bin.tsv",
        "--out", outdir
    ]
    
    result = subprocess.run(cmd)
    shutil.rmtree(outdir)
    
    # Check that the CLI tool executed successfully
    assert result.returncode == 0, "CLI tool did not exit successfully"

def test_motif_discovery_gzip():
    """
    """
    outdir = "tests/cli_test_motif_discovery"
    
    cmd = [
        "nanomotif", "motif_discovery",
        "-t", "1",
        "nanomotif/datasets/geobacillus-plasmids.assembly.fasta",
        "nanomotif/datasets/geobacillus-plasmids.pileup.bed.gz",
        "-c", "nanomotif/datasets/geobacillus-contig-bin.tsv",
        "--out", outdir
    ]
    
    result = subprocess.run(cmd)
    shutil.rmtree(outdir)
    
    # Check that the CLI tool executed successfully
    assert result.returncode == 0, "CLI tool did not exit successfully"

def test_check_installation():
    """
    """
    
    cmd = [
        "nanomotif", "check_installation"
    ]
    
    result = subprocess.run(cmd)
    assert result.returncode == 0, "CLI tool did not exit successfully"

def test_version():
    """
    """
    
    cmd = [
        "nanomotif", "--version"
    ]
    
    result = subprocess.run(cmd)
    assert result.returncode == 0, "CLI tool did not exit successfully"


def test_help():
    """
    """
    
    cmd = [
        "nanomotif", "--help"
    ]
    
    result = subprocess.run(cmd)
    assert result.returncode == 0, "CLI tool did not exit successfully"

@pytest.mark.parametrize("command", ["detect_contamination", "include_contigs"])
def test_removed_commands(command):
    """
    Bin contamination and inclusion were removed; the commands must be rejected.
    """
    result = subprocess.run(["nanomotif", command], capture_output=True)
    assert result.returncode != 0, f"{command} should no longer be accepted"


def _run_motif_discovery(outdir, pileup="nanomotif/datasets/geobacillus-plasmids.pileup.bed", extra=()):
    cmd = [
        "nanomotif", "motif_discovery",
        "-t", "1",
        "nanomotif/datasets/geobacillus-plasmids.assembly.fasta",
        pileup,
        "-c", "nanomotif/datasets/geobacillus-contig-bin.tsv",
        "--out", str(outdir),
        *extra
    ]
    result = subprocess.run(cmd)
    assert result.returncode == 0, "CLI tool did not exit successfully"
    return pd.read_csv(os.path.join(outdir, "bin-motifs.tsv"), sep="\t")


def test_motif_discovery_finds_motifs(tmp_path):
    motifs = _run_motif_discovery(tmp_path / "out")
    assert len(motifs) > 0


def test_threshold_valid_coverage_is_used(tmp_path):
    # no position in the test pileup has this coverage, so nothing can be found
    motifs = _run_motif_discovery(tmp_path / "out", extra=["--threshold_valid_coverage", "100000"])
    assert len(motifs) == 0


def test_motif_discovery_is_reproducible(tmp_path):
    """Same seed, same intermediate results: the background sample depends on contig order."""
    outputs = []
    for i in range(3):
        outdir = tmp_path / f"out{i}"
        _run_motif_discovery(outdir, extra=["-t", "2"])
        outputs.append(sorted(
            (p.relative_to(outdir).as_posix(), p.read_text())
            for p in outdir.rglob("*") if p.suffix in (".txt", ".tsv") and "logs" not in p.parts
        ))
    assert outputs[0] == outputs[1] == outputs[2]


def test_works_without_snakemake():
    """Only MTase-linker needs snakemake: without it, the CLI still loads, and MTase-linker exits with a
    message instead of a traceback."""
    import sys
    block = "import sys; sys.modules['snakemake'] = None; "
    loads = subprocess.run([sys.executable, "-c", block + "import nanomotif.main; nanomotif.main.nm.argparser.create_parser()"],
                           capture_output=True, text=True)
    assert loads.returncode == 0, loads.stderr
    linker = subprocess.run(
        [sys.executable, "-c", block + "import sys; sys.argv = ['nanomotif', 'MTase-linker', 'run', '--assembly', 'a', "
         "'--contig_bin', 'b', '--bin_motifs', 'c']; import nanomotif.main; nanomotif.main.main()"],
        capture_output=True, text=True)
    assert linker.returncode == 1
    assert "MTase-linker needs snakemake" in linker.stderr + linker.stdout
    assert "Traceback" not in linker.stderr
