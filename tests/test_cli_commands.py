import pytest
import subprocess
import pandas as pd
import os
import glob
import shutil


def test_bin_consensus():
    """
    """
    outdir = "tests/cli/test_bin_consensus"
    
    cmd = [
        "nanomotif", "bin_consensus",
        "-t", "1",
        "nanomotif/datasets/geobacillus-plasmids.assembly.fasta",
        "nanomotif/datasets/geobacillus-plasmids.pileup.bed",
        "nanomotif/datasets/geobacillus-plasmids.motifs.tsv",
        "nanomotif/datasets/geobacillus-contig-bin.tsv",
        "nanomotif/datasets/geobacillus-plasmids.motifs-scored.tsv",
        "--out", outdir
    ]
    
    result = subprocess.run(cmd)
    
    # Check that the CLI tool executed successfully
    assert result.returncode == 0, "CLI tool did not exit successfully"

    #TODO: remove outputdir
    #


def test_score_motifs():
    """
    """
    outdir = "tests/cli/test_score_motifs"
    
    cmd = [
        "nanomotif", "score_motifs",
        "-t", "1",
        "nanomotif/datasets/geobacillus-plasmids.assembly.fasta",
        "nanomotif/datasets/geobacillus-plasmids.pileup.bed",
        "nanomotif/datasets/geobacillus-plasmids.motifs.tsv",
        "--out", outdir
    ]
    
    result = subprocess.run(cmd)
    
    # Check that the CLI tool executed successfully
    assert result.returncode == 0, "CLI tool did not exit successfully"

    #TODO: remove outputdir

