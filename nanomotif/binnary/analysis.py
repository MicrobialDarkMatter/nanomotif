import pandas as pd
import numpy as np
import itertools
from scipy.spatial.distance import euclidean

from .data_processing import prepare_motifs_scored_in_bins


def perform_analysis(motifs_scored, bin_motifs, contig_bins, assembly_stats, assembly_file, args):
    """
    Performs the core analysis of the tool. This function is called from the main entry point.
    """
    # Find relevant motifs in bins and contigs
    motifs_scored_in_bins = prepare_motifs_scored_in_bins(motifs_scored, bin_motifs, contig_bins, assembly_stats)
    
    # Calculate belonging score
    belonging_score = calculate_belonging_score(motifs_scored_in_bins, args)
    
    # Determine contamination
    belonging_score = determine_contamination(belonging_score)
    
    # Calculate kmer distance for contigs that can belong to multuple bins
    contamination_kmer_distance_dict = calculate_contamination_kmer_distance(belonging_score, assembly_file, contig_bins, args)
    # Initialize a list to hold your formatted data
    formatted_data = []

    # Iterate over the dictionary to format the string
    for contig, bins in contamination_kmer_distance_dict.items():
        # Format the bin distances into a single string
        matching_bins = '|'.join([f"{bin}:{distance}" for bin, distance in bins.items()])
        # Append the contig and the formatted string to your data list
        formatted_data.append([contig, matching_bins])

    # Create a DataFrame from the formatted data
    contamination_kmer_distance_df = pd.DataFrame(formatted_data, columns=['contig', 'matching_bins'])

    
    # Dataframe for output
    output_df = belonging_score[["bin", "contig", "bin_id", "contamination_group"]]\
        .rename(columns={"bin": "matching_bins", "bin_id": "current_bin"})
    
    ## Merge contamination_df with contamination_kmer_distance
    contamination_matching_multiple_bins = output_df[
        (output_df["contamination_group"] == "multiple_contamination") & (output_df["current_bin"] != "unbinned")
    ].drop(columns=["matching_bins"])\
    .drop_duplicates(subset='contig')\
    .merge(contamination_kmer_distance_df, on="contig", how="left")\
    .dropna()
    
    ## 
    unbinned_matching_multiple_bins = output_df[
        (output_df["contamination_group"] == "multiple_contamination") & (output_df["current_bin"] == "unbinned")
    ]
    unbinned_matching_multiple_bins = unbinned_matching_multiple_bins.groupby('contig')['matching_bins'].apply(lambda x: '|'.join(x)).reset_index()
    unbinned_matching_multiple_bins.rename(columns={'matching_bins': 'matching_bins'}, inplace=True)
    unbinned_matching_multiple_bins["current_bin"] = "unbinned"
    
    correct_assignment_df = output_df[output_df["contamination_group"] == "correct"]
    
    ## Concat all dataframes
    output_final = pd.concat([correct_assignment_df, contamination_matching_multiple_bins, unbinned_matching_multiple_bins])\
        .drop(columns=["contamination_group"])
    
    return output_final





def calculate_belonging_score(motifs_scored_in_bins, args):
    """
    Finds the closest match between contig and bin based on methylation pattern.
    """
    # Step 1: Group by bin and motif_mod and calculate mean methylation
    bin_motif_binary = motifs_scored_in_bins[motifs_scored_in_bins['bin'] != "unbinned"].groupby(['bin', 'motif_mod'])['mean'].mean().reset_index(name='mean_methylation')
    
    # Step 2: Convert mean methylation values to binary
    bin_motif_binary['methylation_binary'] = (bin_motif_binary['mean_methylation'] >= args.mean_methylation_cutoff).astype(int)

    
    # Create contig_motif_binary
    ## Filter motifs that are not observed in bins
    contig_motif_binary = motifs_scored_in_bins[motifs_scored_in_bins["motif_mod"].isin(bin_motif_binary["motif_mod"])]
    
    ## Filter motifs that are not observed more than n_motif_cutoff times
    contig_motif_binary = contig_motif_binary[contig_motif_binary["n_motifs"] >= args.n_motif_cutoff]
    
    ## Convert mean methylation values to binary
    contig_motif_binary["methylation_binary"] = (contig_motif_binary["mean"] >= args.mean_methylation_cutoff).astype(int)
    
    ## Rename bin_contig to bin
    contig_motif_binary = contig_motif_binary[["bin_contig", "motif_mod", "methylation_binary"]]
    contig_motif_binary.rename(columns={"bin_contig": "bin"}, inplace=True)
    
    # Pivot the DataFrame
    contig_motif_binary_pivoted = contig_motif_binary.pivot_table(
        index='bin', columns='motif_mod', values='methylation_binary', fill_value=None
    )
    # Unpivot the DataFrame back to long format
    contig_motif_binary = contig_motif_binary_pivoted.reset_index().melt(
        id_vars=['bin'], var_name='motif_mod', value_name='methylation_binary'
    )
    
    # Check if there is a bin with no methylated motifs
    bin_with_no_methylations_exists = bin_motif_binary.groupby("bin")["methylation_binary"].sum().min() == 0
    
    if bin_with_no_methylations_exists == False:
        print("No bin with no methylated motifs exists. Removing contigs with no methylated motifs...")
        # Remove contigs with no methylated motifs
        contig_motif_binary = contig_motif_binary[contig_motif_binary.groupby("bin")["methylation_binary"].transform("sum") > 0]

    
    # Combine bin_motif_binary and contig_motif_binary
    contig_motif_binary = contig_motif_binary.rename(columns={
        'bin': 'bin_compare',
        'methylation_binary': 'methylation_binary_compare'
    })
    
    motif_binary_compare = pd.merge(bin_motif_binary, contig_motif_binary, on='motif_mod')
    
    # Define the conditions
    conditions = [
        (motif_binary_compare['methylation_binary'] == 1) & (motif_binary_compare['methylation_binary_compare'] == 1),
        (motif_binary_compare['methylation_binary'] == 1) & (motif_binary_compare['methylation_binary_compare'] == 0),
        (motif_binary_compare['methylation_binary'] == 0) & (motif_binary_compare['methylation_binary_compare'] == 1),
        (motif_binary_compare['methylation_binary'] == 0) & (motif_binary_compare['methylation_binary_compare'] == 0),
        (motif_binary_compare['methylation_binary'] == 1) & (motif_binary_compare['methylation_binary_compare'].isna()),
        (motif_binary_compare['methylation_binary'] == 0) & (motif_binary_compare['methylation_binary_compare'].isna())
    ]

    # Define the corresponding choices for each condition
    choices = [1, -1, -1, 0, 0, 0]

    # Use numpy.select to apply these conditions and choices
    motif_binary_compare['motif_comparison_score'] = np.select(conditions, choices, default=np.nan)
    
    
    # Calculate belonging score
    ## Calculate the sum of motif_comparison_score for each combination of bin and bin_compare
    belonging_score = motif_binary_compare.groupby(["bin", "bin_compare"])["motif_comparison_score"].sum().reset_index(name="belonging_score")
    
    # Find the max belonging_score for each bin_compare
    max_scores = belonging_score.groupby("bin_compare")["belonging_score"].max().reset_index(name="max_belonging_score")

    # Merge the max_scores back with the original belonging_score DataFrame to get the corresponding bin values
    # This merge operation ensures that we only retain the rows with the max belonging_score for each bin_compare
    belonging_score_with_max = pd.merge(belonging_score, max_scores, how='inner', left_on=["bin_compare", "belonging_score"], right_on=["bin_compare", "max_belonging_score"])

    # Count the number of best matches for each bin_compare
    belonging_score_with_max["belonging_bins"] = belonging_score_with_max.groupby("bin_compare")["bin_compare"].transform("count")
    
    
    # Drop the max_belonging_score column as it's redundant now
    belonging_score_final = belonging_score_with_max.drop(columns=["max_belonging_score"])

    separated_columns = belonging_score_final['bin_compare'].str.split('_', expand=True)
    belonging_score_final[['bin_id', 'contig', 'contig_number', 'length']] = separated_columns

    # Step 2: Modify 'contig' column
    belonging_score_final['contig'] = 'contig_' + belonging_score_final['contig_number']

    # Step 3: Convert 'length' column to numeric
    belonging_score_final['length'] = pd.to_numeric(belonging_score_final['length'])

    # Step 4: Remove 'contig_number' column
    belonging_score_final.drop(columns=['contig_number'], inplace=True)
    
    
    return belonging_score_final


def determine_contamination(belonging_score):
    """
    Determines whether a contig is contaminated based on belonging_score.
    """
    # Define the conditions
    conditions = [
        (belonging_score["bin"] == belonging_score["bin_id"]) & (belonging_score["belonging_bins"] == 1),
        (belonging_score["bin"] != belonging_score["bin_id"]) & (belonging_score["belonging_bins"] == 1),
        (belonging_score["belonging_bins"] > 1)
    ]

    # Define the corresponding choices for each condition
    choices = [
        "correct", # pattern between bin and contig matches and there is only one best match 
        "single_contamination", # contig pattern does not match bin pattern and there is only one best match
        "multiple_contamination" # contig pattern does not match bin pattern and there are multiple best matches
    ]

    belonging_score["contamination_group"] = np.select(conditions, choices, default=np.nan)
    
    return belonging_score
    


def calculate_contamination_kmer_distance(belonging_score, assembly_file, contig_bins, args):
    """
    Calculates the kmer distance between contigs and bins that they belong to.
    """
    # Intialize a distance dictionary
    distances = {}
    
    # Find contigs that matches multiple bins based on methylation pattern
    contigs_in_bins_matching_multiple_bins = belonging_score[(belonging_score["contamination_group"] == "multiple_contamination") & (belonging_score["bin_id"] != "unbinned")]["contig"].unique()
    
    for contig in contigs_in_bins_matching_multiple_bins:
        # Initialize a dictionary for this contig's distances
        contig_distances = {}
        
        # Find matching bins
        matching_bins = belonging_score[belonging_score["contig"] == contig]["bin"].unique()

        # Generate kmer_df
        kmers = generate_kmers(args.kmer_window_size)
        kmer_df = pd.DataFrame(index=kmers)
        
        # count kmer for each matching bin
        for bin in matching_bins:
            # Initialize a Series to hold summed k-mer counts for the bin
            bin_kmer_counts = pd.Series(0, index=kmers)
            
            # Get contigs in bin
            contigs_in_bin = contig_bins[contig_bins["bin"] == bin]["contig"].unique()
            
            # remove contig of interest from contigs_in_bin
            contigs_in_bin = contigs_in_bin[contigs_in_bin != contig]
            
            # calculate kmer frequency for each contig in bin from the assembly_file dictionary
            for contig_in_bin in contigs_in_bin:
                # Get contig sequence
                contig_sequence = assembly_file.get(contig_in_bin, "")
                # Count kmer for contig
                contig_kmer_counts = count_kmers(contig_sequence, args.kmer_window_size)
                # Add contig_kmer_counts to bin_kmer_counts
                bin_kmer_counts += pd.Series(contig_kmer_counts).reindex(bin_kmer_counts.index, fill_value=0)
            
            # Add summed k-mer counts for the bin to the DataFrame
            kmer_df[bin] = bin_kmer_counts
            
        # Count k-mer for contig of interest
        contig_sequence = assembly_file.get(contig, "")
        contig_kmer_counts = count_kmers(contig_sequence, args.kmer_window_size)
        kmer_df[contig] = pd.Series(contig_kmer_counts).reindex(kmer_df.index, fill_value=0)
        
        # Convert to kmer frequency
        total_kmer_counts = kmer_df.sum(axis=0)
        kmer_df = kmer_df.div(total_kmer_counts, axis=1)
        
        # Calculate distance between contig column and each bin column
        # Iterate over columns in the DataFrame
        for col in kmer_df.columns:
            if col != contig:  # Skip the contig column itself
                # Calculate the distance between the contig column and the current bin column
                distance = euclidean(kmer_df[contig], kmer_df[col])
                # Store the distance in the dictionary
                contig_distances[col] = distance
        
        # After processing all bins, store the contig's distances in the main dictionary
        distances[contig] = contig_distances
        
    return distances



def generate_kmers(k, alphabet='ATCG'):
    """Generate all possible k-mers from the given alphabet."""
    return [''.join(p) for p in itertools.product(alphabet, repeat=k)]

def count_kmers(sequence, k):
    """Count k-mers in a given sequence."""
    kmer_counts = {kmer: 0 for kmer in generate_kmers(k)}
    for i in range(len(sequence) - k + 1):
        kmer = str(sequence[i:i+k])
        if kmer in kmer_counts:  # This check ensures we only count valid k-mers
            kmer_counts[kmer] += 1
    return kmer_counts

