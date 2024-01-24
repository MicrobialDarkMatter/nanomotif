import argparse
from nanomotif._version import __version__


def  create_parser():
    parser = argparse.ArgumentParser(description="Motif identification and utilisation commands")
    parser.add_argument("--version", action="version", version="%(prog)s {}".format(__version__))
    subparsers = parser.add_subparsers(help="Command descriptions", dest="command")

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
    parser_optional.add_argument("--threshold_methylation_general", type=float, default=0.6, help="minimum fraction of reads that must be methylated at a position for the position to be methylated. These position are used for counting number of methylated position of a motif. Default: %(default)s")

    ###########################################################################
    # Find Motifs
    parser_shared_find_motifs = argparse.ArgumentParser(add_help=False)
    parser_shared_find_motifs.add_argument("--search_frame_size", type=int, default=40, help="length of the sequnces sampled around confident methylatyion sites. Default: %(default)s")
    parser_shared_find_motifs.add_argument("--threshold_methylation_confident", type=float, default=0.8, help="minimum fraction of reads that must be methylated at a position for the position to be considered confiently methylated. These position are used to search for candidate motifs. Default: %(default)s")
    parser_shared_find_motifs.add_argument("--threshold_valid_coverage", type=int, default=5, help="minimum valid base coverage for a position to be considered. Default: %(default)s")
    parser_shared_find_motifs.add_argument("--minimum_kl_divergence", type=float, default=0.2, help="minimum KL-divergence for a position to considered for expansion in  motif search. Higher value means less exhaustive, but faster search. Default: %(default)s")
    parser_find_motifs = subparsers.add_parser(
        'find-motifs', 
        parents=[parser_positional, parser_optional, parser_shared_find_motifs], 
        help="identifies motifs in contigs"
    )

    ###########################################################################
    # Score motifs
    parser_score_motifs = subparsers.add_parser(
        'score-motifs', 
        parents=[parser_positional, parser_optional],
        help="generate feature complete output"
    )
    parser_score_motifs.add_argument("motifs", type=str, help="path to the motifs file.")
    
    ###########################################################################
    # Bin consensus
    parser_shared_bin_consensus = argparse.ArgumentParser(add_help=False, conflict_handler="resolve")
    parser_shared_bin_consensus.add_argument("bins", type=str, help="tsv file specifying which bin contigs belong.")

    parser_bin_consensus = subparsers.add_parser(
        'bin-consensus', 
        parents=[parser_positional, parser_optional, parser_shared_bin_consensus],
        help="generate consensus set of motif for each bin"
    )
    parser_bin_consensus.add_argument("motifs", type=str, help="path to the motifs file.")
    parser_bin_consensus.add_argument("motifs_scored", metavar="motifs-scored", type=str, help="path to the motif-scored file.")

    ###########################################################################
    # Complete workflow
    parser_complete_workflow = subparsers.add_parser('complete-workflow', help='run find-motifs, score-motifs and bin-consensus', parents=[parser_positional, parser_optional, parser_shared_find_motifs, parser_shared_bin_consensus], conflict_handler="resolve")

    ###########################################################################
    # Check installation
    parser_check_installation = subparsers.add_parser('check-installation', parents=[parser_optional, parser_shared_find_motifs], conflict_handler="resolve")
    parser_check_installation.add_argument("--out", type=str, help="path to the output folder. Default: %(default)s", default="nanomotif-install-check")
    return parser
    