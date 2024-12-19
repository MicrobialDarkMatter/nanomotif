import pytest
from nanomotif.binnary.detect_contamination import detect_contamination
from nanomotif.binnary import data_processing as dp
from .conftest import MockArgs
import os
import polars as pl


def test_detect_contamination():
    contig_methylation_1 = pl.DataFrame({
                                          "contig": ["contig_1","contig_1","contig_1","contig_2","contig_2","contig_2","contig_3","contig_3","contig_3","contig_4","contig_4","contig_4"],
                                          "bin": ["bin1","bin1","bin1","bin1","bin1","bin1","bin1","bin1","bin1","bin1","bin1","bin1"],
                                          "median": [0.9, 0.85, 0.0, 0.9, 0.85, 0.001, 0.0, 0.001, 0.9, 0.99, 0.9, 0.05],
                                          "N_motif_obs": [1000, 500, 600, 270, 100, 50, 50, 100, 100, 1000, 500, 600,],
                                          "motif_mod": ["m1", "m2", "m3","m1", "m2", "m3","m1", "m2", "m3","m1", "m2", "m3"]
                                      }) 
    contig_methylation_2 = pl.DataFrame({
                                          "contig": ["contig_5","contig_5","contig_5","contig_6","contig_6","contig_6","contig_7","contig_7","contig_7","contig_8","contig_8","contig_8"],
                                          "bin": ["bin2","bin2","bin2","bin2","bin2","bin2","bin2","bin2","bin2","bin2","bin2","bin2"],
                                          "median": [0.04, 0.01, 0.9, 0.0, 0.0, 0.97, 0.1, 0.0, 0.9, 0.0, 0.01, 0.95],
                                          "N_motif_obs": [1000, 500, 600, 270, 100, 50, 50, 100, 100, 1000, 500, 600,],
                                          "motif_mod": ["m1", "m2", "m3","m1", "m2", "m3","m1", "m2", "m3","m1", "m2", "m3"]
                                      }) 
    contig_methylation_3 = pl.DataFrame({
                                          "contig": ["contig_9","contig_9","contig_9","contig_10","contig_10","contig_10","contig_11","contig_11","contig_11","contig_12","contig_12","contig_12"],
                                          "bin": ["bin3","bin3","bin3","bin3","bin3","bin3","bin3","bin3","bin3","bin3","bin3","bin3"],
                                          "median": [0.04, 0.9, 0.9, 0.0, 0.92, 0.97, 0.0, 0.99, 0.9, 0.01, 0.99, 0.95],
                                          "N_motif_obs": [1000, 500, 600, 270, 100, 50, 50, 100, 100, 1000, 500, 600,],
                                          "motif_mod": ["m1", "m2", "m3","m1", "m2", "m3","m1", "m2", "m3","m1", "m2", "m3"]
                                      }) 

    contig_methylation = pl.concat([contig_methylation_1, contig_methylation_2,contig_methylation_3])
    contig_lengths = pl.DataFrame({
                                      "contig": ["contig_1", "contig_2", "contig_3", "contig_4","contig_5", "contig_6", "contig_7", "contig_8","contig_9", "contig_10", "contig_11", "contig_12"],
                                      "length": [100000, 60000, 20000, 80000, 100000, 60000, 20000, 80000, 100000, 60000, 20000, 80000],
                                  })

    contamination = detect_contamination(contig_methylation, contig_lengths, 4, 1, 10)
    print(contamination)

    assert contamination.shape[0] == 4
    assert contamination.get_column("contig").unique().to_list() == ["contig_3"]
    
    
  
