import os
import pytest
import pandas as pd
from nanomotif.binnary.data_processing import generate_output
from pathlib import Path

def test_generate_output_new_directory(tmp_path):
    """
    Test generate_output when the output directory does not exist.
    It should create the directory and the output file.
    """
    df = pd.DataFrame({'col1': [1, 2], 'col2': [3, 4]})
    
    output_dir = tmp_path / "new_directory"  # Use Path object for directory
    output_path = output_dir / "output.tsv"  # Full path including filename
    
    generate_output(df, str(output_dir), "output.tsv")
    
    assert output_dir.is_dir(), "Output directory was not created."
    assert output_path.is_file(), "Output file was not created."
    
    result_df = pd.read_csv(output_path, sep="\t")
    pd.testing.assert_frame_equal(result_df, df, check_dtype=False)

def test_generate_output_existing_directory(tmp_path):
    """
    Test generate_output when the output directory already exists.
    """
    df = pd.DataFrame({'col1': [5, 6], 'col2': [7, 8]})
    
    output_path = tmp_path / "output.tsv"  # Path object for the file
    
    generate_output(df, str(tmp_path), "output.tsv")
    
    assert output_path.is_file(), "Output file was not created in the existing directory."
    
    result_df = pd.read_csv(output_path, sep="\t")
    pd.testing.assert_frame_equal(result_df, df, check_dtype=False)