import argparse
from nanomotif._version import __version__


def  create_parser():
    parser = argparse.ArgumentParser(description="Motif identification and utilisation commands")
    parser.add_argument("--version", action="version", version="%(prog)s {}".format(__version__))

    subparsers = parser.add_subparsers(help="sub-command help", dest="command")

    ###########################################################################
    # Find Motifs
    parser_find_motifs = subparsers.add_parser('find-motifs', help='Find motifs in assembly using modkit pileup')

    parser_find_motifs.add_argument("assembly", type=str, help="Path to the assembly file.")
    parser_find_motifs.add_argument("pileup", type=str, help="Path to the modkit pileup file.")
    parser_find_motifs.add_argument("output", type=str, help="Path to the output file.")

    parser_find_motifs.add_argument("--run_score_motifs", action="store_true", help="Run score-motifs after find-motifs.")

    parser_find_motifs.add_argument("-t", "--threads", type=int, default=1, help="Number of threads to use. Default is 1.")
    parser_find_motifs.add_argument("-v", "--verbose", action="store_true", help="Increase output verbosity. (set logger to debug level)")

    parser_find_motifs.add_argument("--search_frame_size", type=int, default=40, help="Length of the sequnces sampled around methylatyion sites. Default: %(default)s")
    parser_find_motifs.add_argument("--threshold_methylation_confident", type=float, default=0.8, help="Minimum fraction of reads that must be methylated to be considered in motif search. Default: %(default)s")
    parser_find_motifs.add_argument("--threshold_methylation_general", type=float, default=0.5, help="Minimum fraction of reads that must be methylated for a position to be considered modified. Default: %(default)s")
    parser_find_motifs.add_argument("--threshold_valid_coverage", type=int, default=5, help="Minimum valid base coverage for a position to be considered. Default: %(default)s")
    parser_find_motifs.add_argument("--minimum_kl_divergence", type=float, default=0.2, help="Minimum KL-divergence for a position to considered for expansion in  motif search. Default: %(default)s")

    ###########################################################################
    # Score motifs
    parser_score_motifs = subparsers.add_parser('score-motifs', help='Find degree of methylation of all identified motifs for all sequences in the assembly')

    parser_score_motifs.add_argument("assembly", type=str, help="Path to the assembly file.")
    parser_score_motifs.add_argument("pileup", type=str, help="Path to the modkit pileup file.")
    parser_score_motifs.add_argument("motifs", type=str, help="Path to the motifs file.")
    parser_score_motifs.add_argument("output", type=str, help="Path to the output file.")

    parser_score_motifs.add_argument("-t", "--threads", type=int, default=1, help="Number of threads to use. Default: %(default)s")
    parser_score_motifs.add_argument("-v", "--verbose", action="store_true", help="Increase output verbosity")
    parser_score_motifs.add_argument("--threshold_methylation_general", type=float, default=0.5, help="Minimum fraction of reads that must be methylated for a position to be considered modified. Default: %(default)s")

    ###########################################################################
    # Clean bins using motifs
    parser_clean_bins = subparsers.add_parser('clean-bins', help='Clean bins using methylation motifs [NOT IMPLEMENTED]')

    ###########################################################################
    # place MGEs
    parser_associate_mges = subparsers.add_parser('associate-mges', help='Associate plasmids and viruses with bin of similair motif profile [NOT IMPLEMENTED]')

    return parser
    