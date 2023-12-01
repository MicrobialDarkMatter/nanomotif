import argparse
from nanomotif._version import __version__


def  create_parser():
    parser = argparse.ArgumentParser(description="Motif identification and utilisation commands")
    parser.add_argument("--version", action="version", version="%(prog)s {}".format(__version__))

    subparsers = parser.add_subparsers(help="sub-command help", dest="command")

    ###########################################################################
    # Complete workflow
    parser_complete_workflow = subparsers.add_parser('complete-workflow', help='Run the complete workflow')

    parser_complete_workflow.add_argument("assembly", type=str, help="Assembly file.")
    parser_complete_workflow.add_argument("pileup", type=str, help="Modkit pileup file.")
    parser_complete_workflow.add_argument("bins", type=str, help="File specifying to which bin contigs belong. (tsv file with no header and columns: contig, bin)")
    parser_complete_workflow.add_argument("output", type=str, help="Output directory.")

    parser_complete_workflow.add_argument("-t", "--threads", type=int, default=1, help="Number of threads to use. Default is 1.")
    parser_complete_workflow.add_argument("-v", "--verbose", action="store_true", help="Increase output verbosity. (set logger to debug level)")

    parser_complete_workflow.add_argument("--search_frame_size", type=int, default=40, help="Length of the sequnces sampled around methylatyion sites. Default: %(default)s")
    parser_complete_workflow.add_argument("--threshold_methylation_confident", type=float, default=0.8, help="Minimum fraction of reads that must be methylated to be considered in motif search. Default: %(default)s")
    parser_complete_workflow.add_argument("--threshold_methylation_general", type=float, default=0.5, help="Minimum fraction of reads that must be methylated for a position to be considered modified. Default: %(default)s")
    parser_complete_workflow.add_argument("--threshold_valid_coverage", type=int, default=5, help="Minimum valid base coverage for a position to be considered. Default: %(default)s")
    parser_complete_workflow.add_argument("--minimum_kl_divergence", type=float, default=0.2, help="Minimum KL-divergence for a position to considered for expansion in  motif search. Default: %(default)s")



    ###########################################################################
    # Find Motifs
    parser_find_motifs = subparsers.add_parser('find-motifs', help='Find motifs in contigs of the assembly')

    parser_find_motifs.add_argument("assembly", type=str, help="Path to the assembly file.")
    parser_find_motifs.add_argument("pileup", type=str, help="Path to the modkit pileup file.")
    parser_find_motifs.add_argument("output", type=str, help="Path to the output file.")


    parser_find_motifs.add_argument("-t", "--threads", type=int, default=1, help="Number of threads to use. Default is 1.")
    parser_find_motifs.add_argument("-v", "--verbose", action="store_true", help="Increase output verbosity. (set logger to debug level)")

    parser_find_motifs.add_argument("--search_frame_size", type=int, default=40, 
                                    help="Length of the sequnces sampled around confident methylatyion sites. Default: %(default)s")
    parser_find_motifs.add_argument("--threshold_methylation_confident", type=float, default=0.8, 
                                    help="Minimum fraction of reads that must be methylated at a position for the position to be considered confiently methylated. These position are used to search for candidate motifs. Default: %(default)s")
    parser_find_motifs.add_argument("--threshold_methylation_general", type=float, default=0.6, 
                                    help="Minimum fraction of reads that must be methylated at a position for the position to be methylated. These position are used for counting number of methylated position of a motif. Default: %(default)s")
    parser_find_motifs.add_argument("--threshold_valid_coverage", type=int, default=5, 
                                    help="Minimum valid base coverage for a position to be considered. Default: %(default)s")
    parser_find_motifs.add_argument("--minimum_kl_divergence", type=float, default=0.2, 
                                    help="Minimum KL-divergence for a position to considered for expansion in  motif search. Higher value means less exhaustive, but faster search. Default: %(default)s")

    ###########################################################################
    # Score motifs
    parser_score_motifs = subparsers.add_parser('score-motifs', help='Find degree of methylation of all identified motifs for all contigs in the assembly')

    parser_score_motifs.add_argument("assembly", type=str, help="Path to the assembly file.")
    parser_score_motifs.add_argument("pileup", type=str, help="Path to the modkit pileup file.")
    parser_score_motifs.add_argument("motifs", type=str, help="Path to the motifs file.")
    parser_score_motifs.add_argument("output", type=str, help="Path to the output file.")

    parser_score_motifs.add_argument("-t", "--threads", type=int, default=1, help="Number of threads to use. Default: %(default)s")
    parser_score_motifs.add_argument("-v", "--verbose", action="store_true", help="Increase output verbosity")
    parser_score_motifs.add_argument("--threshold_methylation_general", type=float, default=0.5, help="Minimum fraction of reads that must be methylated for a position to be considered modified. Default: %(default)s")
    ###########################################################################
    # Bin consensus
    parser_bin_consensus = subparsers.add_parser('bin-consensus', help='Find consensus motif in each bin')

    parser_bin_consensus.add_argument("bins", type=str, help="TSV file specifying which bin contigs belong.")
    parser_bin_consensus.add_argument("assembly", type=str, help="Path to the assembly file.")
    parser_bin_consensus.add_argument("pileup", type=str, help="Path to the modkit pileup file.")
    parser_bin_consensus.add_argument("motifs", type=str, help="Path to the motifs file.")
    parser_bin_consensus.add_argument("motifs_scored", metavar="motifs-scored", type=str, help="Path to the motif-scored file.")
    parser_bin_consensus.add_argument("output", type=str, help="Path to the output file.")

    parser_bin_consensus.add_argument("-t", "--threads", type=int, default=1, help="Number of threads to use. Default: %(default)s")
    parser_bin_consensus.add_argument("-v", "--verbose", action="store_true", help="Increase output verbosity")
    parser_bin_consensus.add_argument("--threshold_methylation_general", type=float, default=0.6, help="Minimum fraction of reads that must be methylated for a position to be considered modified. Default: %(default)s")

#    ###########################################################################
#    # Clean bins using motifs
#    parser_clean_bins = subparsers.add_parser('clean-bins', help='Clean bins using methylation motifs [NOT IMPLEMENTED]')
#
#    ###########################################################################
#    # place MGEs
#    parser_associate_mges = subparsers.add_parser('associate-mges', help='Associate plasmids and viruses with bin of similair motif profile [NOT IMPLEMENTED]')

    return parser
    