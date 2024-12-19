import pytest
from nanomotif.binnary import data_processing as dp
import pandas as pd
import glob
import os

def test_create_contig_bin_contamination(loaded_data):
    """
    Test create_contig_bin function.
    """
    contig_bin = loaded_data["contig_bins"]
    
    contamination = pd.DataFrame({
        'bin': ['b3', 'b3', 'b3'],
        'bin_contig_compare': ['b3_contig_12', 'b3_contig_13', 'b3_contig_6'],
        'binary_methylation_mismatch_score': [2, 2, 3],
        'non_na_comparisons': [4, 4, 4],
        'contig': ['contig_12', 'contig_13', 'contig_6']
    })
        
    new_contig_bin = dp.create_contig_bin_file(contig_bin.to_pandas(), include=None, contamination=contamination)

    assert new_contig_bin is not None, "Output is None"
    # Contigs in contamination may not be in the new contig_bin
    assert not new_contig_bin["contig"].isin(contamination["contig"]).all(), "Contigs in contamination are not in the new contig_bin"



def test_create_contig_bin_include_and_no_contamination(loaded_data):
    """
    Test create_contig_bin function.
    """
    contig_bin = loaded_data["contig_bins"]
    
    # Just for testing purposes it is not used in the function.
    contamination = pd.DataFrame({
        'bin': ['b3', 'b3', 'b3'],
        'bin_contig_compare': ['b3_contig_12', 'b3_contig_13', 'b3_contig_6'],
        'binary_methylation_mismatch_score': [2, 2, 3],
        'non_na_comparisons': [4, 4, 4],
        'contig': ['contig_12', 'contig_13', 'contig_6']
    })
    
    include = pd.DataFrame({
        'bin': ['b2', 'b2', 'b3'],
        'bin_compare': ['b3_contig_6', 'unbinned_contig_14', 'unbinned_contig_15'],
        'binary_methylation_mismatch_score': [0, 0, 0],
        'non_na_comparisons': [4, 4, 2],
        'contig_bin': ['b3', 'unbinned', 'unbinned'],
        'contig': ['contig_6', 'contig_14', 'contig_15']
    })
    contamination_dummy = pd.DataFrame({
        'bin': [],
        'bin_contig_compare': [],
        'binary_methylation_mismatch_score': [],
        'non_na_comparisons': [],
        'contig': []
    })

    new_contig_bin = dp.create_contig_bin_file(contig_bin.to_pandas(), include=include, contamination=contamination_dummy)

    assert new_contig_bin is not None, "Output is None"
    # Contigs in contamination may not be in the new contig_bin
    assert contamination["contig"].isin(new_contig_bin["contig"]).all(), "Contigs in contamination are not in the new contig_bin"
    assert include["contig"].isin(new_contig_bin["contig"]).all(), "Contigs in include are not in the new contig_bin"
    
  

def test_create_contig_bin_include_and_contamination(loaded_data):
    """
    Test create_contig_bin function.
    """
    contig_bin = loaded_data["contig_bins"]
    
    contamination = pd.DataFrame({
        'bin': ['b3', 'b3', 'b3'],
        'bin_contig_compare': ['b3_contig_12', 'b3_contig_13', 'b3_contig_6'],
        'binary_methylation_mismatch_score': [2, 2, 3],
        'non_na_comparisons': [4, 4, 4],
        'contig': ['contig_12', 'contig_13', 'contig_6']
    })
    
    include = pd.DataFrame({
        'bin': ['b2', 'b2', 'b3'],
        'bin_compare': ['b3_contig_6', 'unbinned_contig_14', 'unbinned_contig_15'],
        'binary_methylation_mismatch_score': [0, 0, 0],
        'non_na_comparisons': [4, 4, 2],
        'contig_bin': ['b3', 'unbinned', 'unbinned'],
        'contig': ['contig_6', 'contig_14', 'contig_15']
    })


    new_contig_bin = dp.create_contig_bin_file(contig_bin.to_pandas(), include=include, contamination=contamination)

    assert new_contig_bin is not None, "Output is None"
    # Contigs in contamination may not be in the new contig_bin
    assert not new_contig_bin["contig"].isin(contamination["contig"]).all(), "Contigs in contamination are in the new contig_bin"
    assert include["contig"].isin(new_contig_bin["contig"]).all(), "Contigs in include are not in the new contig_bin"
    
    
def test_write_bins(loaded_data):
    """
    Check if the bins are written correctly
    """
    contig_bin = loaded_data["contig_bins"]
    
    contamination = pd.DataFrame({
        'bin': ['b3', 'b3', 'b3'],
        'bin_contig_compare': ['b3_contig_12', 'b3_contig_13', 'b3_contig_6'],
        'binary_methylation_mismatch_score': [2, 2, 3],
        'non_na_comparisons': [4, 4, 4],
        'contig': ['contig_12', 'contig_13', 'contig_6']
    })
        
    new_contig_bin = dp.create_contig_bin_file(contig_bin.to_pandas(), include=None, contamination=contamination)
    assembly = dp.read_fasta("datasets/binnary_testdata/assembly_file.fasta")
    dp.write_bins_from_contigs(new_contig_bin, assembly, "tests/binnary/bins_test")
    
    bins_in_dir = glob.glob("tests/binnary/bins_test/*.fa")
    assert set(new_contig_bin["bin"]) == set([os.path.basename(x).split(".")[0] for x in bins_in_dir])
    
    contig_ids = []
    with open("tests/binnary/bins_test/b3.fa", 'r') as file:
        for line in file:
            if line.startswith('>'):
                # Extract the contig ID from the description line
                contig_id = line.split()[0][1:]  # Remove '>' and any extra info
                contig_ids.append(contig_id)
    print(contig_ids)
    assert set(contig_ids) == set(new_contig_bin[new_contig_bin["bin"] == "b3"]["contig"])
    
    
    # remove bins 
    for f in bins_in_dir:
        os.remove(f)
    os.rmdir("tests/binnary/bins_test")
