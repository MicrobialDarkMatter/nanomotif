import nanomotif.dataload as dl
import polars as pl
import os


def geobacillus_plasmids_pileup_path():
    data_path = os.path.join(os.path.dirname(__file__), "datasets", "geobacillus-plasmids.pileup.bed")
    return data_path

def geobacillus_plasmids_pileup():
    data_path = geobacillus_plasmids_pileup_path()
    pileup = dl.load_pileup(data_path, min_fraction=0.3, min_coverage=5)
    return pileup

def geobacillus_plasmids_assembly_path():
    data_path = os.path.join(os.path.dirname(__file__), "datasets", "geobacillus-plasmids.assembly.fasta")
    return data_path
def geobacillus_plasmids_assembly():
    data_path = geobacillus_plasmids_assembly_path()
    assembly = dl.load_assembly(data_path)
    return assembly

def geobacillus_plasmids_bin_path():
    data_path = os.path.join(os.path.dirname(__file__), "datasets", "geobacillus-contig-bin.tsv")
    return data_path
def geobacillus_plasmids_bin():
    data_path = geobacillus_plasmids_bin_path()
    contig_bin = pl.read_csv(data_path, separator="\t", has_header=False) \
        .rename({"column_1":"contig", "column_2":"bin"})
    return contig_bin
