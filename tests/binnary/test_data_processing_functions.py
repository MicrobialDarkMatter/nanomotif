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
    
    assert bin_motif_binary.filter((pl.col("bin") == "b3") & (pl.col("motif_mod") == "m6_a_1")).select("methylation_binary").item() == 1
    
    # Assert that there are 4 motifs in bin 1
    assert bin_motif_binary.filter(pl.col("bin") == "b3").height == 4
    assert set(bin_motif_binary.filter(pl.col("bin") == "b3").to_pandas()["motif_mod"].to_list()) == {"m1_a_1", "m2_a_1", "m3_a_1", "m6_a_1"}
    
    # Assert m7_a is filtered because of too few observations
    assert bin_motif_binary.filter(pl.col("motif_mod") == "m7_a_1").is_empty()
    


def test_add(loaded_data):
    # Access the loaded data directly if returned as a dictionary
    motifs_scored = loaded_data["motifs_scored"]
    bin_motifs = loaded_data["bin_motifs"]
    contig_bins = loaded_data["contig_bins"]

    # Step 1 create bin_motif_binary
    args = MockArgs()
    bin_motif_binary = data_processing.prepare_bin_consensus(bin_motifs, args)

    # Step 2: create motifs_scored_in_bins
    motifs_scored_in_bins = data_processing.add_bin(
        motifs_scored,
        contig_bins
    )
    print(motifs_scored_in_bins)
    # motifs_scored_in_bins
    assert sorted(bin_motif_binary.select("motif_mod").unique()["motif_mod"]) == sorted(
        motifs_scored_in_bins.select("motif_mod").unique()["motif_mod"]
    )
    
    # Assert that the number of columns is 12
    assert set(motifs_scored_in_bins.columns) == set(['contig', 'median', 'N_motif_obs', 'motif_mod', 'bin']) # number of columns
    # Assert contig 1 belongs to bin 1
    assert motifs_scored_in_bins.filter(pl.col("contig") == "contig_1").select("bin").unique()["bin"].item() == "b1"
    

   


def test_impute_contig_methylation_within_bin(loaded_data, motifs_scored_in_bins_and_bin_motifs):
    """
    GIVEN loaded_data and motifs_scored_in_bins_and_bin_motifs
    WHEN construct_bin_motifs_from_motifs_scored_in_bins is called
    THEN assert that the output contains only the expected columns
    """
    motifs_scored_in_bins = motifs_scored_in_bins_and_bin_motifs["motifs_scored_in_bins"]
    
    args = MockArgs()
    print(motifs_scored_in_bins)
    imputed_binned_contig_methylation = data_processing.impute_contig_methylation_within_bin(
        motifs_scored_in_bins
    )
    print(imputed_binned_contig_methylation)
    assert imputed_binned_contig_methylation is not None
    assert set(imputed_binned_contig_methylation.columns) == set(['contig', 'bin', 'motif_mod', 'mean_bin_median', 'median' ])
    assert None not in imputed_binned_contig_methylation.to_pandas()["median"] 
    
def test_impute_contig_methylation_within_bin2():
    # Sample input DataFrame
    contig_methylation = pl.DataFrame({
        "contig": ["contig_1", "contig_1", "contig_1", "contig_2", "contig_3"],
        "bin": ["bin1", "bin1", "bin1", "bin1", "bin2"],
        "motif_mod": ["mod1", "mod2", "mod3", "mod3", "mod2"],
        "median": [0.5, 0.0, 0.9, 0.5, 0.9],
        "N_motif_obs": [10, 5, 15, 20, 25]
    })

    # Expected output DataFrame after imputation
    expected_output = pl.DataFrame({
        "contig": ["contig_1", "contig_1", "contig_1", "contig_2", "contig_2", "contig_2", "contig_3"],
        "bin": ["bin1", "bin1", "bin1", "bin1", "bin1", "bin1", "bin2"],
        "motif_mod": ["mod1", "mod2","mod3", "mod1","mod2", "mod3","mod2"],
        "median": [0.5, 0.0, 0.9, 0.5, 0.0, 0.5, 0.9]
    })

    # Since args is not used in the function (commented out), we can pass None
    args = None

    # Call the function to test
    output = data_processing.impute_contig_methylation_within_bin(contig_methylation)


    # Select and sort columns for comparison
    output = output.select(expected_output.columns).sort(["bin", "contig", "motif_mod"])
    expected_output = expected_output.sort(["bin", "contig", "motif_mod"])
    print(output)
    print(expected_output)

    # Assert that the output matches the expected output
    assert output.equals(expected_output, null_equal=True), "The imputed DataFrame does not match the expected output."
    
def test_impute_unbinned_contigs():
    import numpy as np
    np.random.seed(42)
    # Mock input data
    contig_methylation = pl.DataFrame({
        "bin": ["unbinned", "binned", "unbinned", "binned"],
        "contig": ["contig_1", "contig_2", "contig_3", "contig_4"],
        "motif_mod": ["mod1", "mod2", "mod1", "mod3"],
        "median": [0.0, 0.3, 0.95, 0.5]
    })

    # Define the expected behavior of the pseudo-random number generator
    result = data_processing.impute_unbinned_contigs(contig_methylation)
    print(result)

    # Validate the resulting DataFrame structure
    assert "median" in result.columns
    assert len(result.get_column("contig")) == 6
    assert set(result.get_column("motif_mod").unique()) == set(["mod1", "mod2", "mod3"])
    assert None not in result.to_pandas()["median"]
    # Expected values for the "median" column (based on your output)
    expected_median = [0.0, 0.142607, 0.109799, 0.95, 0.023403, 0.023399]

    # Extract actual "median" values from the DataFrame
    actual_median = result.get_column("median").to_list()

    # Assert the two lists are equal with some tolerance for floating-point precision
    assert np.allclose(actual_median, expected_median, atol=1e-6), \
        f"Median values do not match: {actual_median} != {expected_median}"

    

