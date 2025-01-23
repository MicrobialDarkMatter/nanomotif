import pytest
import subprocess
import pandas as pd
import os
import glob
import shutil

def test_find_motifs():
    """
    """
    outdir = "tests/cli_test_find_motifs"
    
    cmd = [
        "nanomotif", "find_motifs",
        "-t", "1",
        "nanomotif/datasets/geobacillus-plasmids.assembly.fasta",
        "nanomotif/datasets/geobacillus-plasmids.pileup.bed",
        "--out", outdir
    ]
    
    result = subprocess.run(cmd)
    
    # Check that the CLI tool executed successfully
    shutil.rmtree(outdir)
    assert result.returncode == 0, "CLI tool did not exit successfully"


def test_score_motifs():
    """
    """
    outdir = "tests/cli_test_score_motifs"
    
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
    shutil.rmtree(outdir)
    assert result.returncode == 0, "CLI tool did not exit successfully"


    #TODO: remove outputdir

def test_bin_consensus():
    """
    """
    outdir = "tests/cli_test_bin_consensus"
    
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
    shutil.rmtree(outdir)
    assert result.returncode == 0, "CLI tool did not exit successfully"


    #TODO: remove outputdir
    #

def test_motif_discovery():
    """
    """
    outdir = "tests/cli_test_motif_discovery"
    
    cmd = [
        "nanomotif", "motif_discovery",
        "-t", "1",
        "nanomotif/datasets/geobacillus-plasmids.assembly.fasta",
        "nanomotif/datasets/geobacillus-plasmids.pileup.bed",
        "nanomotif/datasets/geobacillus-contig-bin.tsv",
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

def test_detect_contamination():
    """
    """
    # TODO: Add test for detect_contamination
    # Require sample with unbinned contigs in datasets
    pass


def test_include_contigs():
    """
    """
    # TODO: Add test for include_contigs
    # Require sample with unbinned contigs in datasets
    pass