import pytest
import subprocess
import pandas as pd
import os
import glob
import shutil

def remove_dir_with_contents(path):
    if path in ["/", "C:\\", os.environ["HOME"]]:
        print("Dangerous deletion path, aborting.")
        return
    try:
        logs = glob.glob(f"{path}/logs/*.log")
        for log in logs:
            os.remove(log)
            
        os.rmdir(f"{path}/logs")

        tsvs = glob.glob(f"{path}/*.tsv")
        for tsv in tsvs:
            os.remove(tsv)
        
        bins = glob.glob(f"{path}/*_bins/*.fa")
        for bin in bins:
            os.remove(bin)
        
        bin_dirs = glob.glob(f"{path}/*_bins")
        for bin_dir in bin_dirs:
            os.rmdir(bin_dir)
            
        jsons = glob.glob(f"{path}/*.json")
        for json in jsons:
            os.remove(json)
        
        os.rmdir(path)
            
    except Exception as e:
        print(f"Error removing directory {path}: {e}")


def test_contamination_detection():
    """
    Test contamination detection.
    
    nanomotif detect_contamination \
        --motifs_scored datasets/binnary_testdata/motifs-scored.tsv \
        --bin_motifs datasets/binnary_testdata/bin-motifs.tsv \
        --contig_bins datasets/binnary_testdata/contig_bin.tsv \
        --write_bins --assembly_file datasets/binnary_testdata/assembly_file.fasta \
        -t 1 --out testrun_contamination
    
    """
    outdir = "tests/binnary/testrun_contamination_1"
    bin_contamination_file = f"{outdir}/bin_contamination.tsv"
    
    cmd = [
        "nanomotif", "detect_contamination",
        "--motifs_scored", "datasets/binnary_testdata/motifs-scored.tsv",
        "--bin_motifs", "datasets/binnary_testdata/bin-motifs.tsv",
        "--contig_bins", "datasets/binnary_testdata/contig_bin.tsv",
        "--write_bins",
        "--assembly_file", "datasets/binnary_testdata/assembly_file.fasta",
        "-t", "1",
        "--out", outdir
    ]
    
    result = subprocess.run(cmd)
    
    # Check that the CLI tool executed successfully
    assert result.returncode == 0, "CLI tool did not exit successfully"

    # Check that the output directory was created
    assert os.path.isdir(outdir), "Output directory was not created"

    # Check that the TSV file exists
    assert os.path.isfile(bin_contamination_file), "TSV file was not created"

    # Check content
    contamination = pd.read_csv(bin_contamination_file, sep="\t")
    
    # bin	bin_contig_compare	binary_methylation_missmatch_score	non_na_comparisons	contig
    # b3	b3_contig_12	2	4	contig_12
    # b3	b3_contig_13	2	4	contig_13
    # b3	b3_contig_6	3	4	contig_6

    
    assert contamination.shape[0] == 3, "Incorrect number of rows in contamination file"
    assert contamination["bin"].unique() == ["b3"], "Incorrect bin in contamination file"

    # Cleanup after test
    remove_dir_with_contents(outdir)
    
    
def test_contamination_detection_with_file():
    """
    Test contamination detection.
    
    nanomotif detect_contamination \
        --motifs_scored datasets/binnary_testdata/motifs-scored.tsv \
        --bin_motifs datasets/binnary_testdata/bin-motifs.tsv \
        --contig_bins datasets/binnary_testdata/contig_bin.tsv \
        --contamination_file datasets/binnary_testdata/bin_contamination.tsv \
        --write_bins --assembly_file datasets/binnary_testdata/assembly_file.fasta \
        -t 1 --out tests/binnary/testrun_contamination
    
    """
    outdir = "tests/binnary/testrun_contamination_2"
    
    cmd = [
        "nanomotif", "detect_contamination",
        "--motifs_scored", "datasets/binnary_testdata/motifs-scored.tsv",
        "--bin_motifs", "datasets/binnary_testdata/bin-motifs.tsv",
        "--contig_bins", "datasets/binnary_testdata/contig_bin.tsv",
        "--write_bins",
        "--assembly_file", "datasets/binnary_testdata/assembly_file.fasta",
        "--contamination_file", "datasets/binnary_testdata/bin_contamination.tsv",
        "-t", "1",
        "--out", outdir
    ]
    
    result = subprocess.run(cmd)
    
    # Check that the CLI tool executed successfully
    assert result.returncode == 0, "CLI tool did not exit successfully"

    # Check that the output directory was created
    assert os.path.isdir(outdir), "Output directory was not created"

    # Check that the detect_contamination_bins/*.fa exists
    bins_in_dir = glob.glob(f"{outdir}/detect_contamination_bins/*.fa")
    assert len(bins_in_dir) == 4, "No bins were written"
    
    
    remove_dir_with_contents(outdir)
    
    
       
def test_include_with_run_detect_contamination():
    """
    Test include contigs with run contamination flag.
    
    nanomotif include_contigs \
        --motifs_scored nanomotif/datasets/binnary_testdata/motifs-scored.tsv \
        --bin_motifs nanomotif/datasets/binnary_testdata/bin-motifs.tsv \
        --contig_bins nanomotif/datasets/binnary_testdata/contig_bin.tsv \
        --run_detect_contamination \
        --write_bins --assembly_file nanomotif/datasets/binnary_testdata/assembly_file.fasta \
        --min_motif_comparisons 2 \
        -t 1 --out tests/binnary/include_contigs
    
    """
    outdir = "tests/binnary/testrun_include_contigs_1"
    
    cmd = [
        "nanomotif", "include_contigs",
        "--motifs_scored", "datasets/binnary_testdata/motifs-scored.tsv",
        "--bin_motifs", "datasets/binnary_testdata/bin-motifs.tsv",
        "--contig_bins", "datasets/binnary_testdata/contig_bin.tsv",
        "--run_detect_contamination",
        "--write_bins",
        "--assembly_file", "datasets/binnary_testdata/assembly_file.fasta",
        "-t", "1",
        "--out", outdir
    ]
    
    result = subprocess.run(cmd)
    
    # Check that the CLI tool executed successfully
    assert result.returncode == 0, "CLI tool did not exit successfully"

    # Check that the output directory was created
    assert os.path.isdir(outdir), "Output directory was not created"

    # Check that the include_contigs_bins/*.fa exists
    bins_in_dir = glob.glob(f"{outdir}/include_contigs_bins/*.fa")
    assert len(bins_in_dir) == 4, "No bins were written"
    
    remove_dir_with_contents(outdir)
    
    
           
def test_include_with_contamination_file():
    """
    Test include contigs with contamination file.
    
    nanomotif include_contigs \
        --motifs_scored nanomotif/datasets/binnary_testdata/motifs-scored.tsv \
        --bin_motifs nanomotif/datasets/binnary_testdata/bin-motifs.tsv \
        --contig_bins nanomotif/datasets/binnary_testdata/contig_bin.tsv \
        --contamination_file datasets/binnary_testdata/bin_contamination.tsv \
        --write_bins --assembly_file nanomotif/datasets/binnary_testdata/assembly_file.fasta \
        --min_motif_comparisons 2 \
        -t 1 --out tests/binnary/include_contigs
    
    """
    outdir = "tests/binnary/testrun_include_contigs_2"
    
    cmd = [
        "nanomotif", "include_contigs",
        "--motifs_scored", "datasets/binnary_testdata/motifs-scored.tsv",
        "--bin_motifs", "datasets/binnary_testdata/bin-motifs.tsv",
        "--contig_bins", "datasets/binnary_testdata/contig_bin.tsv",
        "--contamination_file", "datasets/binnary_testdata/bin_contamination.tsv",
        "--write_bins",
        "--assembly_file", "datasets/binnary_testdata/assembly_file.fasta",
        "-t", "1",
        "--out", outdir
    ]
    
    result = subprocess.run(cmd)
    
    # Check that the CLI tool executed successfully
    assert result.returncode == 0, "CLI tool did not exit successfully"

    # Check that the output directory was created
    assert os.path.isdir(outdir), "Output directory was not created"
    
    bins_in_dir = glob.glob(f"{outdir}/include_contigs_bins/*.fa")
    assert len(bins_in_dir) == 4, "No bins were written"
    
    remove_dir_with_contents(outdir)
    
    
         
def test_include_no_contamination():
    """
    Test include contigs without running contamination.
    
    nanomotif include_contigs \
        --motifs_scored nanomotif/datasets/binnary_testdata/motifs-scored.tsv \
        --bin_motifs nanomotif/datasets/binnary_testdata/bin-motifs.tsv \
        --contig_bins nanomotif/datasets/binnary_testdata/contig_bin.tsv \
        --write_bins --assembly_file nanomotif/datasets/binnary_testdata/assembly_file.fasta \
        --min_motif_comparisons 2 \
        -t 1 --out tests/binnary/include_contigs
    
    """
    outdir = "tests/binnary/testrun_include_contigs_3"
    
    cmd = [
        "nanomotif", "include_contigs",
        "--motifs_scored", "datasets/binnary_testdata/motifs-scored.tsv",
        "--bin_motifs", "datasets/binnary_testdata/bin-motifs.tsv",
        "--contig_bins", "datasets/binnary_testdata/contig_bin.tsv",
        "--write_bins",
        "--assembly_file", "datasets/binnary_testdata/assembly_file.fasta",
        "-t", "1",
        "--out", outdir
    ]
    
    result = subprocess.run(cmd)
    
    # Check that the CLI tool executed successfully
    assert result.returncode == 0, "CLI tool did not exit successfully"

    # Check that the output directory was created
    assert os.path.isdir(outdir), "Output directory was not created"

    # Check that the detect_contamination_bins/*.fa does not exists
    assert os.path.isdir(f"{outdir}/detect_contamination_bins") == False, "Bins were written"
    
    bins_in_dir = glob.glob(f"{outdir}/include_contigs_bins/*.fa")
    assert len(bins_in_dir) == 4, "No bins were written"
    
    remove_dir_with_contents(outdir)
    
    
         
def test_contamination_with_output_scores():
    """
    Test detect_contamination with outputting scores.
    
    nanomotif detect_contamination \
        --motifs_scored nanomotif/datasets/binnary_testdata/motifs-scored.tsv \
        --bin_motifs nanomotif/datasets/binnary_testdata/bin-motifs.tsv \
        --contig_bins nanomotif/datasets/binnary_testdata/contig_bin.tsv \
        --min_motif_comparisons 2 \
        --save_scores \
        -t 1 --out <out>
    
    """
    outdir = "tests/binnary/testrun_savescores"
    
    cmd = [
        "nanomotif", "detect_contamination",
        "--motifs_scored", "datasets/binnary_testdata/motifs-scored.tsv",
        "--bin_motifs", "datasets/binnary_testdata/bin-motifs.tsv",
        "--contig_bins", "datasets/binnary_testdata/contig_bin.tsv",
        "--save_scores",
        "-t", "1",
        "--out", outdir
    ]
    
    result = subprocess.run(cmd)
    
    # Check that the CLI tool executed successfully
    assert result.returncode == 0, "CLI tool did not exit successfully"

    # Check that the output directory was created
    assert os.path.isdir(outdir), "Output directory was not created"
    
    bins_in_dir = glob.glob(f"{outdir}/scores/detect_contamination/*.csv")
    assert len(bins_in_dir) == 13, "No scores were written"
    
    csvs = glob.glob(f"{outdir}/scores/*/*.csv")
    for csv in csvs:
        os.remove(csv)
    
    csvdirs = glob.glob(f"{outdir}/scores/*")
    for csvdir in csvdirs:
        os.rmdir(csvdir)
    
    os.rmdir(f"{outdir}/scores")
    remove_dir_with_contents(outdir)