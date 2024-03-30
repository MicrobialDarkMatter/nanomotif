import pytest
from nanomotif.binnary import data_processing as dp
import pandas as pd


def test_create_contig_bin(loaded_data):
    """
    Test create_contig_bin function.
    """
    contig_bin = loaded_data["contig_bins"]
    
    contamination = pd.DataFrame({
        'bin': ['b3', 'b3', 'b3'],
        'bin_contig_compare': ['b3_contig_12', 'b3_contig_13', 'b3_contig_6'],
        'binary_methylation_missmatch_score': [2, 2, 3],
        'non_na_comparisons': [4, 4, 4],
        'contig': ['contig_12', 'contig_13', 'contig_6']
    })
        
    new_contig_bin = dp.create_contig_bin_file(contig_bin.to_pandas(), include=None, contamination=contamination)

    assert new_contig_bin is not None, "Output is None"
    # Contigs in contamination may not be in the new contig_bin
    assert not new_contig_bin["contig"].isin(contamination["contig"]).all(), "Contigs in contamination are not in the new contig_bin"
    