import argparse
from nanomotif._version import __version__
import os


def  create_parser():
    formatter = lambda prog: argparse.HelpFormatter(prog,max_help_position=28)
    parser = argparse.ArgumentParser(description="Motif identification and utilisation commands", formatter_class=formatter)
    parser.add_argument("--version", action="version", version="%(prog)s {}".format(__version__))
    subparsers = parser.add_subparsers(help="-- Command descriptions --", dest="command")

    ###########################################################################
    # Shared sub-command arguments
    parser_positional = argparse.ArgumentParser(add_help=False)
    parser_positional.add_argument("assembly", type=str, help="path to the assembly file.")
    parser_positional.add_argument("pileup", type=str, help="path to the modkit pileup file.")

    parser_optional = argparse.ArgumentParser(add_help=False)
    parser_optional.add_argument("--out", type=str, help="path to the output folder", default="nanomotif")
    parser_optional.add_argument("-t", "--threads", type=int, default=1, help="number of threads to use. Default is 1")
    parser_optional.add_argument("-v", "--verbose", action="store_true", help="increase output verbosity. (set logger to debug level)")
    parser_optional.add_argument("--seed", type=int, default=1, help="seed for random number generator. Default: %(default)s")
    parser_optional.add_argument("--threshold_methylation_general", type=float, default=0.70, help="minimum fraction of reads that must be methylated at a position for the position to be methylated. These position are used for counting number of methylated position of a motif. Default: %(default)s")

    ###########################################################################
    # Find Motifs
    parser_shared_find_motifs = argparse.ArgumentParser(add_help=False)
    parser_shared_find_motifs.add_argument("--search_frame_size", type=int, default=40, help="length of the sequnces sampled around confident methylatyion sites. Default: %(default)s")
    parser_shared_find_motifs.add_argument("--threshold_methylation_confident", type=float, default=0.80, help="minimum fraction of reads that must be methylated at a position for the position to be considered confiently methylated. These position are used to search for candidate motifs. Default: %(default)s")
    parser_shared_find_motifs.add_argument("--threshold_valid_coverage", type=int, default=5, help="minimum valid base coverage for a position to be considered. Default: %(default)s")
    parser_shared_find_motifs.add_argument("--minimum_kl_divergence", type=float, default=0.05, help="minimum KL-divergence for a position to considered for expansion in  motif search. Higher value means less exhaustive, but faster search. Default: %(default)s")
    parser_shared_find_motifs.add_argument("--min_motifs_contig", type=int, default=20, help="minimum number of times a motif has to have been oberserved in a contig. Default: %(default)s")
    parser_find_motifs = subparsers.add_parser(
        'find_motifs', 
        parents=[parser_positional, parser_optional, parser_shared_find_motifs], 
        help="Finds motifs directly on contig level in provided assembly"
    )

    ###########################################################################
    # Score motifs
    parser_score_motifs = subparsers.add_parser(
        'score_motifs', 
        parents=[parser_positional, parser_optional],
        help="Find motifs indirectly in contigs by scoring with motifs found in other contigs"
    )
    parser_score_motifs.add_argument("motifs", type=str, help="path to the motifs file.")
    
    ###########################################################################
    # Bin consensus
    parser_shared_bin_consensus = argparse.ArgumentParser(add_help=False, conflict_handler="resolve")
    parser_shared_bin_consensus.add_argument("--min_motifs_bin", type=int, default=50, help="minimum number of times a motif has to have been oberserved in a bin. Default: %(default)s")
    parser_bin_consensus = subparsers.add_parser(
        'bin_consensus', 
        parents=[parser_positional, parser_optional, parser_shared_bin_consensus],
        help="Indentifies highly methylated motifs in bins"
    )
    parser_bin_consensus.add_argument("motifs", type=str, help="path to the motifs file.")
    parser_bin_consensus.add_argument("bins", type=str, help="tsv file specifying which bin contigs belong.")
    parser_bin_consensus.add_argument("motifs_scored", metavar="motifs-scored", type=str, help="path to the motif-scored file.")

    ###########################################################################
    # Complete workflow
    parser_complete_workflow = subparsers.add_parser('motif_discovery', help='Runs find_motifs, score_motifs and bin_consensus', parents=[parser_positional, parser_optional, parser_shared_find_motifs, parser_shared_bin_consensus], conflict_handler="resolve")
    parser_complete_workflow.add_argument("bins", type=str, help="tsv file specifying which bin contigs belong.")

    ###########################################################################
    # Bin contamination and inclusion
    parser_binnary_shared = argparse.ArgumentParser(description="Contamination DNA Methylation Pattern", add_help=False)
    """Function to add common arguments to a subparser."""
    parser_binnary_shared.add_argument(
        "--motifs_scored", type=str, help="Path to motifs-scored.tsv from nanomotif", required=True
    )
    parser_binnary_shared.add_argument("--bin_motifs", type=str, help="Path to bin-motifs.tsv file", required=True)
    parser_binnary_shared.add_argument(
        "--contig_bins", type=str, help="Path to bins.tsv file for contig bins", required=True
    )
    parser_binnary_shared.add_argument("-t", "--threads", type=int, default=1, help="Number of threads to use for multiprocessing")
    parser_binnary_shared.add_argument(
        "--mean_methylation_cutoff",
        type=float,
        default=0.25,
        help="Cutoff value for considering a motif as methylated",
    )
    parser_binnary_shared.add_argument(
        "--n_motif_contig_cutoff",
        type=int,
        default=10,
        help="Number of motifs that needs to be observed in a contig before it is considered valid for scoring",
    )
    parser_binnary_shared.add_argument(
        "--n_motif_bin_cutoff",
        type=int,
        default=500,
        help="Number of motifs that needs to be observed in a bin to be considered valid for scoring",
    )
    
    parser_binnary_shared.add_argument(
        "--ambiguous_motif_percentage_cutoff",
        type=float,
        default=0.40,
        help="Percentage of ambiguous motifs defined as mean methylation between 0.05 and 0.40 in a bin. Motifs with an ambiguous methylation percentage of more than this value are removed from scoring. Default is 0.40",
    )
    parser_binnary_shared.add_argument(
        "--write_bins",
        action='store_true',
        help="If specified, new bins will be written to a bins folder. Requires --assembly_file to be specified.",
    )
    parser_binnary_shared.add_argument(
        "--assembly_file",
        type=str,
        help="Path to assembly.fasta file"
    )
    parser_binnary_shared.add_argument(
        "--save_scores",
        action='store_true',
        help="If specified, the scores for each comparison will be saved to a scores folder in the output directory"
    )
    
    parser_binnary_shared.add_argument("--out", type=str, help="Path to output directory", required=True, default="nanomotif")
    
    # Binnary contamination
    parser_contamination = subparsers.add_parser(
        'detect_contamination', 
        help="Detect contamination in bins",
        parents=[parser_binnary_shared]
    )  
    
    parser_contamination.add_argument(
        '--contamination_file',
        type=str,
        help="Path to an existing contamination file if bins should be outputtet as a post-processing step"
    )
    
    # Binnary inclusion
    parser_inclusion = subparsers.add_parser(
        'include_contigs', 
        help="Include contigs in bins",
        parents=[parser_binnary_shared]
    )
    
    group = parser_inclusion.add_mutually_exclusive_group(required=False)
    group.add_argument(
        "--contamination_file",
        type=str,
        help="Path to an existing contamination file to include in the analysis"
    )
    group.add_argument(
        "--run_detect_contamination",
        action='store_true',
        help="Indicate that the detect_contamination workflow should be run first"
    )
    parser_inclusion.add_argument(
        "--min_motif_comparisons",
        type=int,
        default=5,
        help="Minimum number of non-NA motif comparisons required to include a contig in the analysis. Default is 5",
    )
    
    ###########################################################################
    # MTase-linker
    parser_mtase_linker = subparsers.add_parser(
        'MTase-linker', 
        help="Commands related to MTase-linker"
    )

    mtase_linker_subparsers = parser_mtase_linker.add_subparsers(
        title="MTase-linker commands", 
        dest="mtase_linker_command"
    )

    parser_mtase_linker_run = mtase_linker_subparsers.add_parser(
        'run', 
        help="Run the MTase-linker workflow"
    )
    parser_mtase_linker_run.add_argument("-t", "--threads", type=int, default=1, help="Number of threads to use. Default is 1")
    parser_mtase_linker_run.add_argument("--forceall", type=bool, default=False, help="Flag for snakemake. Forcerun workflow regardless of already created output (default = False)")
    parser_mtase_linker_run.add_argument("--dryrun", type=bool, default=False, help="Flag for snakemake. Dry-run the workflow. Default is False")
    parser_mtase_linker_run.add_argument("--assembly", type=str, required=True, help="Path to assembly file.")
    parser_mtase_linker_run.add_argument("--contig_bin", type=str, required=True, help="tsv file specifying which bin contigs belong.")
    parser_mtase_linker_run.add_argument("--bin_motifs", type=str, required=True, help="bin-motifs.tsv output from nanomotif.")
    parser_mtase_linker_run.add_argument("-d", "--dependency_dir", type=str, default=os.path.join(os.getcwd(), "ML_dependencies"), help="Path to the ML_dependencies directory created during installation of the MTase-linker module. Default is cwd/ML_dependencies")
    parser_mtase_linker_run.add_argument("-o", "--out", type=str, default=os.path.join(os.getcwd(), "mtase_linker"), help="Path to output directory. Default is cwd")
    parser_mtase_linker_run.add_argument("--identity", type=str, default=80, help="Identity threshold for motif prediction. Default is 80")
    parser_mtase_linker_run.add_argument("--qcovs", type=str, default=80, help="Query coverage for motif prediction. Default is 80")

    parser_mtase_linker_install = mtase_linker_subparsers.add_parser(
        'install', 
        help="Install additional dependencies for MTase-linker"
    )
    parser_mtase_linker_install.add_argument("-d", "--dependency_dir", type=str, default=os.getcwd(), help="Path to installation location of dependencies. A folder named ML_dependencies will be generated. Default is cwd.")
    


    ###########################################################################
    # Check installation
    parser_check_installation = subparsers.add_parser('check_installation', parents=[parser_optional, parser_shared_find_motifs, parser_shared_bin_consensus], add_help=False, help="Performs small test run to verify that the installation is correct.")
    
    
    return parser
    