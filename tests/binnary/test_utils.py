import polars as pl
import pandas as pd
from nanomotif.binnary import utils as ut



def test_split_bin_contig():
    # Input DataFrame
    input_data = pl.DataFrame({
        "bin_compare": ["bin1_contig_23", "bin1_contig_40", "bin_c1_contig_23"]
    })

    # Expected output DataFrame
    expected_output = pd.DataFrame({
        "bin_compare": ["bin1_contig_23", "bin1_contig_40", "bin_c1_contig_23"],
        "contig_bin": ["bin1", "bin1", "bin_c1"],
        "contig": ["contig_23", "contig_40", "contig_23"]
    })

    # Run the function
    output = ut.split_bin_contig(input_data)
    output = output.to_pandas()
    print(output)
    
    # Verify the output
    pd.testing.assert_frame_equal(output, expected_output, check_like=True)