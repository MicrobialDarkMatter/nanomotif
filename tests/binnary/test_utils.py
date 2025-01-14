import polars as pl
import pandas as pd
import nanomotif.binnary.utils as ut
from pymethylation_utils.utils import run_epimetheus
from .conftest import MockArgs
import os
from pathlib import Path

def test_split_bin_contig():
    # Input DataFrame
    input_data = pl.DataFrame({
        "contig_bin": ["bin1", "bin1", "bin_c1"],
        "bin_compare": ["bin1_contig_23", "bin1_contig_40", "bin_c1_contig_23"]
    })

    # Expected output DataFrame
    expected_output = pd.DataFrame({
        "contig_bin": ["bin1", "bin1", "bin_c1"],
        "bin_compare": ["bin1_contig_23", "bin1_contig_40", "bin_c1_contig_23"],
        "contig": ["contig_23", "contig_40", "contig_23"]
    })

    # Run the function
    output = ut.split_bin_contig(input_data)
    output = output.to_pandas()
    print(output)
    
    # Verify the output
    pd.testing.assert_frame_equal(output, expected_output, check_like=True)

def test_methylation_utils():
    args = MockArgs()

    run_epimetheus(
        pileup = "nanomotif/datasets/geobacillus-plasmids.pileup.bed",
        assembly = "nanomotif/datasets/geobacillus-plasmids.assembly.fasta",
        motifs = ["GATC_m_3", "GATC_a_1"],
        threads = 1,
        min_valid_read_coverage = args.min_valid_read_coverage,
        output = os.path.join(args.out, "motifs-scored-read-methylation.tsv")
    )

    file = Path(os.path.join(args.out, "motifs-scored-read-methylation.tsv"))
    assert file.exists(), "motifs_scored-read-methylation.tsv does not exist"

    res = pl.read_csv(file, separator = "\t")
    print(res.columns)
    assert res.columns == ['contig', 'motif', 'mod_type', 'mod_position', 'median', 'mean_read_cov', 'N_motif_obs', 'motif_occurences_total']
    assert res.shape == (4, 8), "Shape does not match"

    
