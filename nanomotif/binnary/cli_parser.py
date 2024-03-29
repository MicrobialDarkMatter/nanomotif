import argparse


def add_common_arguments(subparser):
    """Function to add common arguments to a subparser."""
    subparser.add_argument(
        "--motifs_scored", type=str, help="Path to motifs-scored.tsv from nanomotif", required=True
    )
    subparser.add_argument("--bin_motifs", type=str, help="Path to bin-motifs.tsv file", required=True)
    subparser.add_argument(
        "--contig_bins", type=str, help="Path to bins.tsv file for contig bins", required=True
    )
    subparser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads to use for multiprocessing")

    subparser.add_argument(
        "--mean_methylation_cutoff",
        type=float,
        default=0.25,
        help="Cutoff value for considering a motif as methylated",
    )
    subparser.add_argument(
        "--n_motif_contig_cutoff",
        type=int,
        default=10,
        help="Number of motifs that needs to be observed in a contig before it is considered valid for scoring",
    )
    subparser.add_argument(
        "--n_motif_bin_cutoff",
        type=int,
        default=500,
        help="Number of motifs that needs to be observed in a bin to be considered valid for scoring",
    )
    
    subparser.add_argument(
        "--ambiguous_motif_percentage_cutoff",
        type=float,
        default=0.40,
        help="Percentage of ambiguous motifs defined as mean methylation between 0.05 and 0.40 in a bin. Motifs with an ambiguous methylation percentage of more than this value are removed from scoring. Default is 0.40",
    )
    subparser.add_argument("--out", type=str, help="Path to output directory", required=True)


def get_parser():
    # Create the top-level parser
    parser = argparse.ArgumentParser(description="DNA Methylation Pattern Analysis")
    subparsers = parser.add_subparsers(dest="command", help="Sub-command help")

    # Create subparsers for each command
    commands = ["detect_contamination", "include_contigs"] #"bin_unbinned", "bin_contamination"
    for command in commands:
        subparser = subparsers.add_parser(command, help=f"{command} help")
        add_common_arguments(subparser)  # Add common arguments to each subparser
        
        if command == "include_contigs":
             # Create a mutually exclusive group within the include_contigs subparser
            group = subparser.add_mutually_exclusive_group(required=True)
            
            # Option for providing an existing file
            group.add_argument(
                "--contamination_file",
                type=str,
                help="Path to an existing contamination file to include in the analysis"
            )

            # Option to indicate that the detect_contamination workflow should be run
            group.add_argument(
                "--run_detect_contamination",
                action='store_true',
                help="Indicate that the detect_contamination workflow should be run first"
            )
            
            subparser.add_argument(
                "--write_bins",
                action='store_true',
                help="If specified, new bins will be written to a bins folder. Requires --assembly_file to be specified.",
            )
            
            subparser.add_argument(
                "--assembly_file",
                type=str,
                help="Path to assembly.fasta file"
            )
            
            subparser.add_argument(
                "--min_motif_comparisons",
                type=int,
                default=5,
                help="Minimum number of non-NA motif comparisons required to include a contig in the analysis. Default is 5",
            )
            
    return parser
