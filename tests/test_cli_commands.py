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
        "nanomotif/datasets/geobacillus-plasmids.assembly.fasta",
        "nanomotif/datasets/geobacillus-plasmids.pileup.bed",
        "-t", "1",
        "--out", outdir
    ]
    
    result = subprocess.run(cmd)
    
    # Check that the CLI tool executed successfully
    shutil.rmtree(outdir)
    assert result.returncode == 0, "CLI tool did not exit successfully"



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
    import polars as pl
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    infile = "nanomotif/datasets/geobacillus-plasmids.assembly.fasta"
    outfile_a = "nanomotif/datasets/geobacillus-plasmids.assembly.duplicated.fasta"

    # Read all records from the original file
    records = list(SeqIO.parse(infile, "fasta"))

    # We'll store our new records here
    new_records = []

    # For each record, if it's contig_2 or contig_3, create 5 duplicates
    # with IDs appended by _1.._5. Otherwise, keep it as is.
    for record in records:
        if record.id in ["contig_2", "contig_3"]:
            for i in range(1, 6):
                new_id = f"{record.id}_{i}"
                # Create a new SeqRecord with the same sequence
                new_record = SeqRecord(
                    record.seq,
                    id=new_id,
                    description=""
                )
                new_records.append(new_record)
        else:
            # For non-contig_2/3, simply keep the original record
            new_records.append(record)

    # Write out the new FASTA file
    SeqIO.write(new_records, outfile_a, "fasta")

    print(f"Duplicated contigs written to: {outfile_a}")


    infile = "nanomotif/datasets/geobacillus-plasmids.pileup.bed"
    outfile_p = "nanomotif/datasets/geobacillus-plasmids.pileup.duplicated.bed"
    
    p = pl.read_csv(infile, has_header = False, separator = "\t")
    
    p_dup = pl.DataFrame()
    for contig in ["contig_3", "contig_2"]:
        p_tmp = p.filter(pl.col("column_1") == contig)

        for i in range(1, 6):
            p_i = p_tmp.with_columns(
                (pl.col("column_1") + f"_{i}").alias("column_1")
            )

            p_dup = pl.concat([p_dup, p_i])

    p_dup.write_csv(outfile_p, separator = "\t", include_header = False)

    outfile_b = "nanomotif/datasets/geobacillus-plasmids.contig_bin.tmp.tsv"
    contig_bin = pl.DataFrame({
                                  "contig": [f"contig_{i}_{j}" for i in [2, 3] for j in range(1, 6)],
                                  "bin": ["bin_1"] * 10,
                              })

    contig_bin.write_csv(outfile_b, separator ="\t", include_header = False)
    
    outdir = "tests/cli_test_detect_contamination"

    cmd = [
        "nanomotif", "detect_contamination",
        "-t", "1",
        "--force",
        "--assembly", outfile_a,
        "--pileup", outfile_p,
        "--contig_bins", outfile_b,
        "--bin_motifs", "nanomotif/datasets/geobacillus-plasmids.bin-motifs.tsv",
        "--out", outdir
    ]
    result = subprocess.run(cmd)
    
    # Check that the CLI tool executed successfully
    shutil.rmtree(outdir)
    assert result.returncode == 0, "CLI tool did not exit successfully"
    for outfile in [outfile_a, outfile_p, outfile_b]:
        if os.path.exists(outfile):
            os.remove(outfile)
            print(f"Deleted: {outfile}")
        else:
            print(f"File not found: {outfile}")


def test_include_contigs():
    """
    """
    import polars as pl
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    infile = "nanomotif/datasets/geobacillus-plasmids.assembly.fasta"
    outfile_a = "nanomotif/datasets/geobacillus-plasmids.assembly.duplicated.fasta"

    # Read all records from the original file
    records = list(SeqIO.parse(infile, "fasta"))

    # We'll store our new records here
    new_records = []

    # For each record, if it's contig_2 or contig_3, create 5 duplicates
    # with IDs appended by _1.._5. Otherwise, keep it as is.
    for record in records:
        if record.id in ["contig_2", "contig_3"]:
            for i in range(1, 7):
                new_id = f"{record.id}_{i}"
                # Create a new SeqRecord with the same sequence
                new_record = SeqRecord(
                    record.seq,
                    id=new_id,
                    description=""
                )
                new_records.append(new_record)
        else:
            # For non-contig_2/3, simply keep the original record
            new_records.append(record)

    # Write out the new FASTA file
    SeqIO.write(new_records, outfile_a, "fasta")

    print(f"Duplicated contigs written to: {outfile_a}")


    infile = "nanomotif/datasets/geobacillus-plasmids.pileup.bed"
    outfile_p = "nanomotif/datasets/geobacillus-plasmids.pileup.duplicated.bed"
    
    p = pl.read_csv(infile, has_header = False, separator = "\t")
    
    p_dup = pl.DataFrame()
    for contig in ["contig_3", "contig_2"]:
        p_tmp = p.filter(pl.col("column_1") == contig)

        for i in range(1, 7):
            p_i = p_tmp.with_columns(
                (pl.col("column_1") + f"_{i}").alias("column_1")
            )

            p_dup = pl.concat([p_dup, p_i])

    p_dup.write_csv(outfile_p, separator = "\t", include_header = False)

    outfile_b = "nanomotif/datasets/geobacillus-plasmids.contig_bin.tmp.tsv"
    contig_bin = pl.DataFrame({
                                  "contig": [f"contig_{i}_{j}" for i in [2, 3] for j in range(1, 6)],
                                  "bin": ["bin_1"] * 5 + ["bin_2"] * 5,
                              })

    contig_bin.write_csv(outfile_b, separator ="\t", include_header = False)
    
    outdir = "tests/cli_test_include_contigs"

    cmd = [
        "nanomotif", "include_contigs",
        "-t", "1",
        "--force",
        "--run_detect_contamination",
        "--assembly", outfile_a,
        "--pileup", outfile_p,
        "--contig_bins", outfile_b,
        "--bin_motifs", "nanomotif/datasets/geobacillus-plasmids.bin-motifs.tsv",
        "--out", outdir
    ]
    result = subprocess.run(cmd)
    
    # Check that the CLI tool executed successfully
    shutil.rmtree(outdir)
    assert result.returncode == 0, "CLI tool did not exit successfully"
    for outfile in [outfile_a, outfile_p, outfile_b]:
        if os.path.exists(outfile):
            os.remove(outfile)
            print(f"Deleted: {outfile}")
        else:
            print(f"File not found: {outfile}")

