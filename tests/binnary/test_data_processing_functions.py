import pytest
from nanomotif.binnary import data_processing
from .conftest import MockArgs
import polars as pl

def test_feature_with_loaded_data(loaded_data):
    # Access the loaded data directly if returned as a dictionary
    motifs_scored = loaded_data["motifs_scored"]
    bin_motifs = loaded_data["bin_motifs"]
    contig_bins = loaded_data["contig_bins"]

    # Now you can use the data in your test assertions
    assert motifs_scored is not None

    # bin motifs
    assert bin_motifs is not None
    assert set(bin_motifs["bin"].unique()) == {"b1", "b2", "b3"}

    # contig bins
    assert contig_bins is not None
    assert set(contig_bins["bin"].unique()) == {"b1", "b2", "b3", "b4"}

    contig_set = {f"contig_{i}" for i in range(1, 14)} | {"contig_16"}
    assert set(contig_bins["contig"].unique()) == contig_set



def test_prepare_bin_consensus(loaded_data):
    """
    GIVEN loaded_data
    WHEN prepare_bin_consensus is called
    THEN assert that the output contains only the expected columns
    """
    args = MockArgs()
    
    bin_motif_binary = data_processing.prepare_bin_consensus(loaded_data["bin_motifs"], args)
    
    assert bin_motif_binary is not None
    assert bin_motif_binary.columns == ["bin", "motif_mod", "mean_methylation", "methylation_binary"]
    
    assert bin_motif_binary.filter((pl.col("bin") == "b3") & (pl.col("motif_mod") == "m6_a-1")).select("methylation_binary").item() == 1
    
    # Assert that there are 4 motifs in bin 1
    assert bin_motif_binary.filter(pl.col("bin") == "b3").height == 4
    assert set(bin_motif_binary.filter(pl.col("bin") == "b3").to_pandas()["motif_mod"].to_list()) == {"m1_a-1", "m2_a-1", "m3_a-1", "m6_a-1"}
    
    # Assert m7_a is filtered because of too few observations
    assert bin_motif_binary.filter(pl.col("motif_mod") == "m7_a-1").is_empty()
    


def test_motifs_scored_in_bins(loaded_data):
    # Access the loaded data directly if returned as a dictionary
    motifs_scored = loaded_data["motifs_scored"]
    bin_motifs = loaded_data["bin_motifs"]
    contig_bins = loaded_data["contig_bins"]

    # Step 1 create bin_motif_binary
    args = MockArgs()
    bin_motif_binary = data_processing.prepare_bin_consensus(bin_motifs, args)

    # Step 2: create motifs_scored_in_bins
    motifs_scored_in_bins = data_processing.prepare_motifs_scored_in_bins(
        motifs_scored,
        bin_motif_binary.select("motif_mod").unique()["motif_mod"],
        contig_bins
    )

    # motifs_scored_in_bins
    assert sorted(bin_motif_binary.select("motif_mod").unique()["motif_mod"]) == sorted(
        motifs_scored_in_bins.select("motif_mod").unique()["motif_mod"]
    )
    
    # Assert that the number of columns is 12
    assert len(motifs_scored_in_bins.columns) == 12  # number of columns
    # Assert contig 1 belongs to bin 1
    assert motifs_scored_in_bins.filter(pl.col("contig") == "contig_1").select("bin").unique()["bin"].item() == "b1"
    



def test_remove_ambiguous_motifs(loaded_data, motifs_scored_in_bins_and_bin_motifs):
    """
    GIVEN loaded_data and motifs_scored_in_bins_and_bin_motifs
    WHEN construct_bin_motifs_from_motifs_scored_in_bins is called
    THEN assert that the output contains only the expected columns
    """
    motifs_scored_in_bins = motifs_scored_in_bins_and_bin_motifs["motifs_scored_in_bins"]
    
    args = MockArgs()
    
    
    bins_w_no_ambiguous_motifs = data_processing.remove_ambiguous_motifs_from_bin_consensus(
        motifs_scored_in_bins,
        args
    )
    print(bins_w_no_ambiguous_motifs)
    
    print(bins_w_no_ambiguous_motifs.filter(pl.col("bin") == "b1").to_pandas())
    assert bins_w_no_ambiguous_motifs is not None
    assert bins_w_no_ambiguous_motifs.columns == ['bin', 'motif_mod']
    assert bins_w_no_ambiguous_motifs.filter((pl.col("bin") == "b1") & (pl.col("motif_mod") == "m1_a")).is_empty()
    
    


def test_bin_motifs_from_motifs_scored_in_bins(loaded_data, motifs_scored_in_bins_and_bin_motifs):
    """
    GIVEN loaded_data and motifs_scored_in_bins_and_bin_motifs
    WHEN construct_bin_motifs_from_motifs_scored_in_bins is called
    THEN assert that the output contains only the expected columns
    """
    motifs_scored_in_bins = motifs_scored_in_bins_and_bin_motifs["motifs_scored_in_bins"]
    
    args = MockArgs()
    print(motifs_scored_in_bins)
    bin_motifs_from_motifs_scored_in_bins = data_processing.construct_bin_consensus_from_motifs_scored_in_bins(
        motifs_scored_in_bins,
        args
    )
    
    assert bin_motifs_from_motifs_scored_in_bins is not None
    assert bin_motifs_from_motifs_scored_in_bins.columns == ['bin', 'motif_mod', 'n_mod', 'n_nomod', 'n_motifs_bin', 'mean_methylation', 'mean_methylation_bin', 'std_methylation_bin', 'n_contigs', 'methylation_binary']
    
    
# def test_calculate_binary_motif_comparison_matrix(loaded_data, motifs_scored_in_bins_and_bin_motifs):
#     """
#     GIVEN loaded_data and motifs_scored_in_bins_and_bin_motifs
#     WHEN calculate_binary_motif_comparison_matrix is called
#     THEN assert that the output contains only the expected columns
#     """
#     motifs_scored_in_bins = motifs_scored_in_bins_and_bin_motifs["motifs_scored_in_bins"]
    
#     args = MockArgs()
    
#     motifs_scored_in_bins_filtered = motifs_scored_in_bins \
#         .filter(~pl.col("bin_contig").str.contains("unbinned"))
    
#     motif_binary_compare = data_processing.calculate_binary_motif_comparison_matrix(
#         motifs_scored_in_bins_filtered,
#         args
#     )
    
#     assert motif_binary_compare is not None
#     # Assert that no bin_contig contains "unbinned"
#     assert motif_binary_compare.filter(pl.col("bin_compare").str.contains("unbinned")).is_empty()
    
#     # b3 = motif_binary_compare[(motif_binary_compare["bin"] == "b3") & (motif_binary_compare["bin_compare"].str.contains("b3"))]
#     b3 = motif_binary_compare.filter((pl.col("bin") == "b3") & (pl.col("bin_compare").str.contains("b3")))
    
#     assert set(b3.to_pandas()["motif_mod"].unique()) == set(["m1_a", "m2_a", "m3_a", "m6_a"])
#     assert b3.filter(pl.col("motif_mod") == "m6_a").select("methylation_binary").get_column("methylation_binary")[0] == 1
#     assert b3.filter(pl.col("motif_mod") == "m2_a").select("methylation_binary").get_column("methylation_binary")[0] == 0
    