import pytest
from nanomotif.binnary import data_processing

class MockArgs:
    def __init__(self):
        self.motifs_scored = "datasets/binnary_testdata/motifs-scored.tsv"
        self.bin_motifs = "datasets/binnary_testdata/bin-motifs.tsv"
        self.contig_bins = "datasets/binnary_testdata/contig_bin.tsv"
        self.mean_methylation_cutoff = 0.25
        self.n_motif_bin_cutoff = 500
        self.n_motif_contig_cutoff = 10
        self.ambiguous_motif_percentage_cutoff = 0.40
        self.min_motif_comparisons = 2
        self.save_scores = False
        self.out = "tests/binnary/test_output"
        self.threads = 1

# Remove tests/test_output directory
import shutil
shutil.rmtree("tests/binnary/test_output", ignore_errors=True)


@pytest.fixture(scope="session")
def loaded_data():
    # Instantiate MockArgs
    mock_args = MockArgs()

    # Load the data using the mock_args object
    data = data_processing.load_data(mock_args)

    # Unpack the data tuple to individual variables if needed
    motifs_scored, bin_motifs, contig_bins = data

    # Return the data as a dictionary or as individual variables, depending on your preference
    return {
        "motifs_scored": motifs_scored,
        "bin_motifs": bin_motifs,
        "contig_bins": contig_bins
    }
    

@pytest.fixture(scope="session")
def motifs_scored_in_bins_and_bin_motifs(loaded_data):
    args = MockArgs()
    
    # Access the loaded data directly if returned as a dictionary
    motifs_scored = loaded_data["motifs_scored"]
    bin_motifs = loaded_data["bin_motifs"]
    contig_bins = loaded_data["contig_bins"]
    
    bin_motif_binary = data_processing.prepare_bin_consensus(bin_motifs, args)
    
    motifs_in_bins = bin_motif_binary.select("motif_mod").unique()["motif_mod"]
    
    # Step 2: create motifs_scored_in_bins
    motifs_scored_in_bins = data_processing.prepare_motifs_scored_in_bins(
        motifs_scored, motifs_in_bins, contig_bins 
    )
    
    return {
        "motifs_scored_in_bins": motifs_scored_in_bins,
        "bin_motif_binary": bin_motif_binary
    }