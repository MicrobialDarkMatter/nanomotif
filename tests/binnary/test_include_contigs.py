import pytest
from nanomotif.binnary import include_contigs as ic
from nanomotif.binnary import detect_contamination as dc
from nanomotif.binnary import data_processing as dp
from .conftest import MockArgs

def test_include_contigs(loaded_data, motifs_scored_in_bins_and_bin_motifs):
    """
    GIVEN loaded_data
    WHEN include_contigs is called
    THEN assert that the output contains only contig 3 and the expected columns
    """
    args = MockArgs()
    
    motifs_scored_in_bins = motifs_scored_in_bins_and_bin_motifs["motifs_scored_in_bins"]
    
    
    bin_motifs_from_motifs_scored_in_bins = dp.construct_bin_consensus_from_motifs_scored_in_bins(
        motifs_scored_in_bins,
        args
    )
    
    contamination = dc.detect_contamination(motifs_scored_in_bins, bin_motifs_from_motifs_scored_in_bins, args)
    
    include_contigs_df = ic.include_contigs(motifs_scored_in_bins, bin_motifs_from_motifs_scored_in_bins, contamination, args)
    
    include_contigs_df = include_contigs_df.to_pandas()
    
    assert include_contigs_df is not None
    assert include_contigs_df.shape[0] == 3
    assert include_contigs_df[include_contigs_df["bin"] == "b3"]["contig"].iloc[0] == "contig_15"
    
    
    
def test_include_with_too_high_min_motif_comparisons(loaded_data, motifs_scored_in_bins_and_bin_motifs):
    """
    GIVEN loaded_data
    WHEN include_contigs is called with a min_motif_comparisons of 10
    THEN assert that the output contains only contig 3 and the expected columns
    """
    args = MockArgs()
    args.min_motif_comparisons = 10
    
    
    
    motifs_scored_in_bins = motifs_scored_in_bins_and_bin_motifs["motifs_scored_in_bins"]
    bin_motifs_from_motifs_scored_in_bins = dp.construct_bin_consensus_from_motifs_scored_in_bins(
        motifs_scored_in_bins,
        args
    )
    
    contamination = dc.detect_contamination(motifs_scored_in_bins, bin_motifs_from_motifs_scored_in_bins, args)
    
    include = ic.include_contigs(motifs_scored_in_bins, bin_motifs_from_motifs_scored_in_bins, contamination, args)
    
    # Assert df is empty
    assert include.shape[0] == 0