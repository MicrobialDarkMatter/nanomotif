import polars as pl
import pandas as pd
import nanomotif.binnary.utils as ut
from epymetheus.epymetheus import methylation_pattern, MethylationOutput
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
    cwd = os.getcwd()
    print(cwd)
    print(args.out)
    
    # Ensure output directory exists
    os.makedirs(args.out, exist_ok=True)

    methylation_pattern(
        pileup = "./nanomotif/datasets/geobacillus-plasmids.pileup.bed",
        assembly = "./nanomotif/datasets/geobacillus-plasmids.assembly.fasta",
        motifs = ["GATC_m_3", "GATC_a_1"],
        threads = 1,
        min_valid_read_coverage = args.min_valid_read_coverage,
        output = os.path.join(args.out, "motifs-scored-read-methylation.tsv"),
        batch_size=1000,
        min_valid_cov_to_diff_fraction=0.8,
        allow_assembly_pileup_mismatch=False,
        output_type = MethylationOutput.Median
    )

    file = Path(os.path.join(args.out, "motifs-scored-read-methylation.tsv"))
    assert file.exists(), "motifs_scored-read-methylation.tsv does not exist"

    res = pl.read_csv(file, separator = "\t")
    print(res.columns)
    assert res.columns == ['contig', 'motif', 'mod_type', 'mod_position', 'methylation_value', 'mean_read_cov', 'n_motif_obs']
    assert res.shape == (4, 7), "Shape does not match"

    
