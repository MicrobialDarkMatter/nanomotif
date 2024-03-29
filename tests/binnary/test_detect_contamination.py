import pytest
from nanomotif.binnary import detect_contamination
from nanomotif.binnary import data_processing as dp
from .conftest import MockArgs
import os


def test_detect_contamination(loaded_data, motifs_scored_in_bins_and_bin_motifs):
    """
    GIVEN loaded_data
    WHEN detect_contamination is called
    THEN assert that the output contains only contig 3 and the expected columns
    """
    args = MockArgs()
    
    motifs_scored_in_bins = motifs_scored_in_bins_and_bin_motifs["motifs_scored_in_bins"]
    
    
    bin_motifs_from_motifs_scored_in_bins = dp.construct_bin_consensus_from_motifs_scored_in_bins(
        motifs_scored_in_bins,
        args
    )
    
    contaminated_contigs = detect_contamination.detect_contamination(
        motifs_scored_in_bins,
        bin_motifs_from_motifs_scored_in_bins,
        args
    )
    
    contaminated_contigs = contaminated_contigs.to_pandas()
    
    assert contaminated_contigs is not None
    assert contaminated_contigs["bin"].unique().tolist() == ["b3"]
    assert sorted(contaminated_contigs["contig"].unique().tolist()) == ["contig_12", "contig_13", "contig_6"]
    assert contaminated_contigs[contaminated_contigs["contig"] == "contig_6"]["binary_methylation_missmatch_score"].values[0] == 3.0


def test_load_contamination_file(loaded_data, motifs_scored_in_bins_and_bin_motifs):
    """
    Test that the contamination file is loaded correctly.
    """
    # setup
    args = MockArgs()
    
    motifs_scored_in_bins = motifs_scored_in_bins_and_bin_motifs["motifs_scored_in_bins"]
    
    
    bin_motifs_from_motifs_scored_in_bins = dp.construct_bin_consensus_from_motifs_scored_in_bins(
        motifs_scored_in_bins,
        args
    )
    
    contaminated_contigs = detect_contamination.detect_contamination(
        motifs_scored_in_bins,
        bin_motifs_from_motifs_scored_in_bins,
        args
    )
    
    contaminated_contigs = contaminated_contigs.to_pandas().to_csv("./tests/contaminated_contigs.tsv", index=False,sep="\t")
    
    # test
    
    contamination_file = dp.load_contamination_file("./tests/contaminated_contigs.tsv")
    required_columns = ["bin", "bin_contig_compare", "binary_methylation_missmatch_score", "non_na_comparisons", "contig"]
    
    contamination_file = contamination_file.to_pandas()
    
    assert not contamination_file.empty, "The contamination file should not be empty."
    assert list(contamination_file.columns) == required_columns, "The contamination file columns are incorrect."
    
    os.remove("./tests/contaminated_contigs.tsv")
    
    