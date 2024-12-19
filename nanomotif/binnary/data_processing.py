import polars as pl
import pandas as pd
import numpy as np
np.random.seed(1)

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import gzip
import os
import logging

def load_data(args):
    """
    Load data for the analysis.
    """
    logger = logging.getLogger(__name__)
    logger.info("Loading data.")
    # Load the data from the provided files
    # Load the data from the provided files using Polars
    bin_motifs = pl.read_csv(args.bin_motifs, separator="\t")
    
    # Polars automatically infers headers; if args.contig_bins file doesn't have headers, you need to specify it
    contig_bins = pl.read_csv(args.contig_bins, separator="\t", has_header=False, infer_schema_length=10000) \
        .rename({
            "column_1": "contig",
            "column_2": "bin"
        })

    logger.info("Data loaded.")
    return bin_motifs, contig_bins

def read_fasta(file_path):
    # Check if the file exists
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")    
    
    # Check if the file has a valid FASTA extension
    valid_extensions = ['.fasta', '.fa', '.fna', '.gz']
    if not any(file_path.endswith(ext) for ext in valid_extensions):
        raise ValueError(f"Unsupported file extension. Please provide a FASTA file with one of the following extensions: {', '.join(valid_extensions)}")

    # Check if the file is a gzipped FASTA file
    if file_path.endswith('.gz'):
        with gzip.open(file_path, "rt") as handle:  # "rt" mode for reading as text
            return {record.id: str(record.seq) for record in SeqIO.parse(handle, "fasta")}
    else:
        # Read a regular (uncompressed) FASTA file
        with open(file_path, "r") as handle:
            return {record.id: str(record.seq) for record in SeqIO.parse(handle, "fasta")}

def find_contig_lengths(assembly):
    # Parse the assembly file using Biopython
    contigs = []
    lengths = []
    
    for id, seq in assembly.items():
        contigs.append(id)       # Contig name (FASTA header)
        lengths.append(len(seq)) # Contig length
    
    # Create a Polars DataFrame
    df = pl.DataFrame({
        "contig": contigs,
        "length": lengths
    })
    
    return df

def write_bins_from_contigs(new_contig_bins, assembly_dict, output_dir):
    """
    Creates new bin files based on contig assignments in new_contig_bins DataFrame.

    Args:
    - new_contig_bins (DataFrame): DataFrame with contig and bin assignments.
    - assembly_dict (dict): Dictionary with contig IDs as keys and sequences as values, from the assembly file.
    - output_dir (str): Path to the directory where the new bin files will be saved.
    """
    logger = logging.getLogger(__name__)
    
    # check if the output directory exists
    logger.info(f"Writing bins to {output_dir}")
    output_bins_dir = os.path.join(output_dir)
    if not os.path.exists(output_bins_dir):
        os.makedirs(output_bins_dir)
    
    # Group the DataFrame by 'bin'
    grouped = new_contig_bins.groupby('bin')

    for bin_name, group in grouped:
        # Initialize a list to hold SeqRecord objects for the current bin
        bin_records = []

        for contig in group['contig']:
            # Create a SeqRecord object for each contig in the bin, if it exists in the assembly dictionary
            if contig in assembly_dict:
                seq_record = SeqRecord(Seq(assembly_dict[contig]), id=contig, description="")
                bin_records.append(seq_record)

        # Define the output file name for the current bin
        output_file = f"{bin_name}.fa"
        output_path = os.path.join(output_bins_dir, output_file)

        # Write the SeqRecord objects to a FASTA file
        with open(output_path, "w") as output_handle:
            SeqIO.write(bin_records, output_handle, "fasta")

        logger.info(f"Written {len(bin_records)} contigs to {output_file}")



def generate_output(output_df, outdir, filename, header=True):
    """
    Generate the output files for the analysis.
    """
    # If the output directory does not exist, create it
    if not os.path.exists(outdir) and outdir != "":
        os.makedirs(outdir)
    
    file_path = os.path.join(outdir, filename)
    
    # Generate the output files
    output_df.to_csv(file_path, sep="\t", index=False, header=header)


def prepare_bin_consensus(bin_motifs, args):
    """
    Prepares the bin_consensus_from_bin_motifs DataFrame by calculating the mean methylation per bin and motif_mod and converting it to binary.    
    """
    # Combine 'motif' and 'mod_type' into 'motif_mod'
    bin_motifs = bin_motifs.with_columns((pl.col("motif") + "_" + pl.col("mod_type") + "_" + pl.col("mod_position").cast(pl.Utf8)).alias("motif_mod"))

    # Calculate total motifs and mean methylation
    bin_motifs = bin_motifs.with_columns([
        (pl.col("n_mod_bin") + pl.col("n_nomod_bin")).alias("n_motifs"),
        (pl.col("n_mod_bin") / (pl.col("n_mod_bin") + pl.col("n_nomod_bin"))).alias("mean_methylation")
    ])

    # Create 'methylation_binary' column based on mean methylation cutoff
    bin_motifs = bin_motifs.with_columns((pl.col("mean_methylation") >= args.mean_methylation_cutoff).cast(int).alias("methylation_binary"))

    # Filter motifs that are not observed more than n_motif_bin_cutoff times
    bin_motifs = bin_motifs.filter(pl.col("n_motifs") >= args.n_motif_bin_cutoff)

    # Filter for rows where 'methylation_binary' is 1 and select relevant columns
    bin_motif_binary = bin_motifs.filter(pl.col("methylation_binary") == 1)[["bin", "motif_mod", "mean_methylation", "methylation_binary"]]

    # Create a binary index for bin and motif_mod
    # Assuming unique_bins and unique_motif_mods are lists of unique values
    unique_bins = bin_motif_binary["bin"].unique().to_list()
    unique_motif_mods = bin_motif_binary["motif_mod"].unique().to_list()

    # Create a DataFrame with all combinations of "bin" and "motif_mod"
    bin_motif_combinations = pl.DataFrame({
        "bin": [bin for bin in unique_bins for _ in unique_motif_mods],
        "motif_mod": [motif_mod for _ in unique_bins for motif_mod in unique_motif_mods]
    })

    # Join the combinations DataFrame with bin_motif_binary
    bin_motif_binary_enriched = bin_motif_combinations.join(
        bin_motif_binary, on=["bin", "motif_mod"], how="left"
    )

    # Fill missing values if necessary and adjust columns as needed
    bin_motif_binary_enriched = bin_motif_binary_enriched.with_columns([
        pl.col("mean_methylation").fill_null(0).alias("mean_methylation"),
        pl.col("methylation_binary").fill_null(0).alias("methylation_binary")
    ])
    

    return bin_motif_binary_enriched


def add_bin(contig_methylation, contig_bins):
    # Filter and enhance contig_methylation 
    contig_methylation_w_bin = contig_methylation \
        .with_columns(
            (pl.col("motif") + "_" + pl.col("mod_type") + "_" + pl.col("mod_position").cast(pl.Utf8)).alias("motif_mod")
        )
    
    # Merge with contig_bins
    contig_methylation_w_bin = contig_methylation_w_bin.join(contig_bins, on="contig", how="left") \
        .with_columns(pl.col("bin").fill_null("unbinned")) \
        .drop(["mod_position", "mod_type", "motif"])
    
    return contig_methylation_w_bin


def impute_contig_methylation_within_bin(contig_methylation):
    contig_methylation = contig_methylation\
        .filter(pl.col("bin") != "unbinned")
    
    bin_methylation = contig_methylation \
        .group_by(["bin", "motif_mod"]) \
        .agg(
            (
                (pl.col("median") * pl.col("N_motif_obs")).sum() / pl.col("N_motif_obs").sum()
            ).alias("mean_bin_median")
        )

    contig_bin_cross = contig_methylation\
        .unique(subset=["contig", "bin"], maintain_order = True)\
        .select(["contig", "bin"])

    contig_methylation_imputed = contig_bin_cross\
        .join(bin_methylation, on = ["bin"], how = "left")\
        .join(contig_methylation, on = ["bin", "contig", "motif_mod"], how = "left")\
        .sort(["bin", "contig", "motif_mod"])\
        .with_columns(
            pl.when(pl.col("median").is_null()).then(pl.col("mean_bin_median")).otherwise(pl.col("median")).alias("median")
        )\
        .drop("N_motif_obs")
        
    return contig_methylation_imputed

def impute_unbinned_contigs(contig_methylation):
    unbinned_contig_methylation = contig_methylation\
        .filter(pl.col("bin") == "unbinned")

    unbinned_contigs = unbinned_contig_methylation\
        .select(["contig"])\
        .unique()
        
    cross_motif_mod_contig = contig_methylation\
        .select(["motif_mod"])\
        .unique()\
        .join(unbinned_contigs, how = "cross")\
        .with_columns(
            pl.lit("unbinned").alias("bin")
        )\
        .sort(["contig", "motif_mod"])

    num_rows = cross_motif_mod_contig.height
    random_values = np.random.uniform(0.0, 0.15, size = num_rows)
    
    cross_motif_mod_contig = cross_motif_mod_contig\
        .with_columns(
            pl.Series("pseudo_median", random_values)
        )
    
    imputed_unbinned_contig_methylation = cross_motif_mod_contig\
        .join(unbinned_contig_methylation, on=["contig", "motif_mod", "bin"], how="left")\
        .with_columns(
            pl.when(pl.col("median").is_null())
                .then(pl.col("pseudo_median"))
                .otherwise(pl.col("median"))
                .alias("median")
        )\
        .drop("pseudo_median")\
        .sort(["contig", "motif_mod"])

    return imputed_unbinned_contig_methylation
    

def create_matrix(contig_methylation):
    matrix_df = contig_methylation\
        .select(["contig", "motif_mod", "median"])\
        .pivot(
            values = "median",
            index = "contig",
            columns="motif_mod",
            aggregate_function = None,
            sort_columns = True
        )\
        .fill_null(0)

    contig_names = matrix_df.get_column("contig")
    matrix = matrix_df.drop("contig").to_numpy()
    return contig_names, matrix

    
def load_contamination_file(contamination_file):
    """
    Load the contamination file from the provided path.
    """
    contamination = pl.read_csv(contamination_file, separator= "\t")
    
    # Check if the file contains the required columns
    required_columns = ["contig"]
    if not all(column in contamination.columns for column in required_columns):
        raise ValueError("The contamination file does not contain the required columns.")
    
    return contamination


def create_contig_bin_file(contig_bins, contamination, include=None):
    """
    Create a new contig_bin file based on the analysis results and contamination file.
    """
    # Remove contigs in the contamination file from the contig_bins
    contig_bins = contig_bins[~contig_bins["contig"].isin(contamination["contig"])]
    
    # Add the contigs in the include DataFrame to the contig_bins
    if include is not None:
        contig_bins = pd.concat([contig_bins, include[["contig", "bin"]]], ignore_index=True)
    
    # Sort the contig_bins by bin and contig
    contig_bins = contig_bins.sort_values(by=["bin", "contig"])
    
    return contig_bins
