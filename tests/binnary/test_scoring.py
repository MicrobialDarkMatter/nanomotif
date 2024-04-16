import pytest
from nanomotif.binnary import scoring as sc
from nanomotif.binnary import data_processing as dp
import polars as pl
from nanomotif.binnary.utils import split_bin_contig

from .conftest import MockArgs

def test_compare_methylation_pattern_contamination(loaded_data, motifs_scored_in_bins_and_bin_motifs):
    """
    
    """
    motifs_scored_in_bins = motifs_scored_in_bins_and_bin_motifs["motifs_scored_in_bins"]
    
    args = MockArgs()
    
    motifs_scored_in_bins = motifs_scored_in_bins \
        .filter(~pl.col("bin_contig").str.contains("unbinned"))
    
    bin_consensus = dp.construct_bin_consensus_from_motifs_scored_in_bins(
        motifs_scored_in_bins,
        args
    )
    
    # Define the corresponding choices for each condition
    # choices = [
    #     0,  # bin motif is methylated, contig motif is methylated
    #     1,  # bin motif is methylated, contig motif is not methylated
    #     1,  # bin motif is not methylated, contig motif is methylated
    #     0,  # bin motif is not methylated, contig motif is not methylated
    #     0,  # bin motif is methylated, contig motif is not observed
    #     0,  # bin motif is not methylated, contig motif is not observed
    # ]
    
    contig_bin_comparison_score, contigs_w_no_methylation = sc.compare_methylation_pattern_multiprocessed(
        motifs_scored_in_bins=motifs_scored_in_bins,
        bin_consensus=bin_consensus,
        mode="contamination",
        # choices=choices,
        args=args
    )
    contig_bin_comparison_score = split_bin_contig(contig_bin_comparison_score)
    print(contig_bin_comparison_score)
    print(contigs_w_no_methylation)
    
    contig_bin_comparison_score = contig_bin_comparison_score.to_pandas()
    
    assert contig_bin_comparison_score is not None
    assert set(contig_bin_comparison_score["bin"].unique()) == {'b3'}
    assert len(contig_bin_comparison_score["contig"].unique()) == 3
    # assert contig_bin_comparison_score[(contig_bin_comparison_score["bin"] == "b2") & (contig_bin_comparison_score["contig"] == "contig_4")]["binary_methylation_missmatch_score"].values[0] == 0
    

def test_compare_methylation_pattern_include(loaded_data, motifs_scored_in_bins_and_bin_motifs):
    """
    
    """
    motifs_scored_in_bins = motifs_scored_in_bins_and_bin_motifs["motifs_scored_in_bins"]
    
    args = MockArgs()
    
    motifs_scored_in_bins = motifs_scored_in_bins \
        .filter(~pl.col("bin_contig").str.contains("unbinned"))
    
    bin_consensus = dp.construct_bin_consensus_from_motifs_scored_in_bins(
        motifs_scored_in_bins,
        args
    )
    
    contig_bin_comparison_score, contigs_w_no_methylation = sc.compare_methylation_pattern_multiprocessed(
        motifs_scored_in_bins=motifs_scored_in_bins,
        bin_consensus=bin_consensus,
        mode="include",
        args=args
    )
    contig_bin_comparison_score = split_bin_contig(contig_bin_comparison_score)
    print(contig_bin_comparison_score)
    print(contigs_w_no_methylation)
    
    contig_bin_comparison_score = contig_bin_comparison_score.to_pandas()
    
    assert contig_bin_comparison_score is not None
    assert set(contig_bin_comparison_score["bin"].unique()) =={'b2', 'b3'} 
    assert len(contig_bin_comparison_score["contig"].unique()) == 2
    
