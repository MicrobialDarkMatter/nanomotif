import nanomotif.dataload as dl
import nanomotif.fasta as fasta
import polars as pl
import os


def geobacillus_plasmids_pileup_path():
    data_path = os.path.join(os.path.dirname(__file__), "datasets", "geobacillus-plasmids.pileup.bed")
    return data_path

def geobacillus_plasmids_pileup():
    data_path = geobacillus_plasmids_pileup_path()
    pileup = dl.load_pileup(data_path)
    return pileup

def geobacillus_plasmids_assembly_path():
    data_path = os.path.join(os.path.dirname(__file__), "datasets", "geobacillus-plasmids.assembly.fasta")
    return data_path
def geobacillus_plasmids_assembly():
    data_path = geobacillus_plasmids_assembly_path()
    assembly = fasta.load_fasta_fastx(data_path)
    return assembly

def geobacillus_plasmids_bin_path():
    data_path = os.path.join(os.path.dirname(__file__), "datasets", "geobacillus-contig-bin.tsv")
    return data_path
def geobacillus_plasmids_bin():
    data_path = geobacillus_plasmids_bin_path()
    contig_bin = pl.read_csv(data_path, separator="\t", has_header=False)
    contig_bin = contig_bin.rename(dict(zip(contig_bin.columns[:2], ["contig", "bin"])))
    return contig_bin
