import argparse

def  create_parser():
    parser = argparse.ArgumentParser(description="Motif identification from Nanopore methylation pileup")
    
    # Argument for assembly file
    parser.add_argument("assembly", type=str, help="Path to the assembly file.")
    
    # Argument for pileup file
    parser.add_argument("pileup", type=str, help="Path to the modkit pileup file.")

    # Argument for output file
    parser.add_argument("output", type=str, help="Path to the output file.")

    # Argument for motif length
    parser.add_argument("-k", "--max_motif_length", type=int, default=30, help="Length of the motif to be searched for. Default is 30.")

    #  Argument for minimum read methylation fraction
    parser.add_argument("-f", "--min_fraction", type=float, default=0.8, help="Minimum fraction of reads that must be methylated for a position to be considered modified. Default is 0.8.")

    # min coverage
    parser.add_argument("-c", "--min_coverage", type=int, default=10, help="Minimum coverage for a position to be considered modified. Default is 10.")

    # min KL-divergence
    parser.add_argument("--min_kl_divergence", type=float, default=0.1, help="Minimum KL-divergence for growing motif. Default is 0.1.")

    # min_cdf_score
    parser.add_argument("--min_cdf_score", type=float, default=0.8, help="Minimum CDF score for motif to be considered valid. Default is 0.8.")

    # cdf_position
    parser.add_argument("--cdf_position", type=float, default=0.55, help="Position in the CDF to be used as a threshold for motif identification. Default is 0.55.")

    

    return parser
    