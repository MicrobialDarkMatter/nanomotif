import argparse

def  create_parser():
    parser = argparse.ArgumentParser(description="Motif identification from Nanopore methylation pileup")
    
    # Argument for assembly file
    parser.add_argument("assembly", type=str, help="Path to the assembly file.")
    
    # Argument for pileup file
    parser.add_argument("pileup", type=str, help="Path to the modkit pileup file.")

    # Argument for output file
    parser.add_argument("output", type=str, help="Path to the output file.")

    # Argument for number of threads
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads to use. Default is 1.")

    # Argument for motif length
    parser.add_argument("-k", "--max_motif_length", type=int, default=30, help="Length of the motif to be searched for. Default is 30.")

    #  Argument for minimum read methylation fraction
    parser.add_argument("-f", "--min_fraction", type=float, default=0.8, help="Minimum fraction of reads that must be methylated for a position to be considered modified. Default is 0.8.")

    # min coverage
    parser.add_argument("-c", "--min_coverage", type=int, default=10, help="Minimum coverage for a position to be considered modified. Default is 10.")

    # min KL-divergence
    parser.add_argument("--min_kl_divergence", type=float, default=0.2, help="Minimum KL-divergence for growing motif. Default is 0.1.")

    # Verbosity
    parser.add_argument("-v", "--verbose", action="store_true", help="Increase output verbosity. (set logger to debug level)")

    return parser
    