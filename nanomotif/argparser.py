import argparse
from nanomotif._version import __version__
import os


def create_parser():
    formatter = lambda prog: argparse.HelpFormatter(prog, max_help_position=28)
    parser = argparse.ArgumentParser(description="Motif identification and utilisation commands", formatter_class=formatter)
    parser.add_argument("--version", action="version", version="%(prog)s {}".format(__version__))
    subparsers = parser.add_subparsers(help="-- Command descriptions --", dest="command", title = "commands",
                                       metavar='{motif_discovery, MTase-linker, check_installation}')

    def add_general_arguments(parser):
        """
        Function to add general arguments to a parser.
        """
        group = parser.add_argument_group("general arguments")
        group.add_argument(
            "-t", "--threads", type=int, default=1, 
            help="Number of threads to use. Default is 1"
        )

        group.add_argument(
            "-v", "--verbose", action="store_true", 
            help="Increase output verbosity. (set logger to debug level)"
        )

        group.add_argument(
            "--seed", type=int, default=1, 
            help="Seed for random number generator. Default: %(default)s"
        )

        group.add_argument(
            "-h", "--help", action="help", 
            help="show this help message and exit"
        )

    def add_contig_bin_arguments(parser):
        """
        Function to add arguments for contig bin association.
        """
        group_major = parser.add_argument_group("contig bin arguments, use one of:")
        group = group_major.add_mutually_exclusive_group(required=True)


        group.add_argument(
            "-c", "--contig_bin", type=str, 
            help="Two-column TSV file (contig, bin) without header: the bin of each contig."
        )

        group.add_argument(
            "-f", "--files", nargs='+', 
            help="List of bin FASTA files with contig names as headers."
        )

        group.add_argument(
            "-d", "--directory", 
            help="Directory containing bin FASTA files with contig names as headers."
        )

        group_major.add_argument(
            "--extension", type=str, default=".fasta",
            help="File extension of the bin FASTA files if using -d (DIRECTORY) argument. Default is '.fasta'."
        )


    
    def add_methylation_threshold_arguments(parser):
        """
        Function to add arguments for methylation threshold.
        """
        parser.add_argument(
            "--methylation_threshold_low", type=float, default=0.30, 
            help="A position is considered unmethylated if its fraction of modified reads is at most this. Default: %(default)s"
        )

        parser.add_argument(
            "--methylation_threshold_high", type=float, default=0.70, 
            help="A position is considered methylated if its fraction of modified reads is at least this. Default: %(default)s"
        )

    def add_motif_search_arguments(parser):
        parser.add_argument(
            "--search_frame_size", type=int, default=40, 
            help="Length of the sequences sampled around confidently methylated positions. Default: %(default)s"
        )
        
        parser.add_argument(
            "--minimum_kl_divergence", type=float, default=0.05, 
            help="Minimum KL divergence for a position to be considered for expansion in the motif search. A higher value means a less exhaustive, but faster search. Default: %(default)s"
        )
        
        parser.add_argument(
            "--min_motif_score", type=float, default=1.5, 
            help="Minimum score for a motif to be kept after identification. Default: %(default)s"
        )   


    ###########################################################################
    # Find motifs on bin level
    parser_find_motifs_bin = subparsers.add_parser(
        'motif_discovery',
        help="Finds motifs directly on bin level in provided assembly",
        add_help=False
    )


    parser_find_motifs_bin.add_argument("assembly", type=str, 
        help="Path to the assembly (FASTA).")
    
    parser_find_motifs_bin.add_argument("pileup", type=str, 
        help="Path to the modkit pileup (plain or bgzipped).") 
    
    add_contig_bin_arguments(parser_find_motifs_bin)

    parser_find_motifs_bin_options = parser_find_motifs_bin.add_argument_group("Options")
    parser_find_motifs_bin_options.add_argument(
        "--out", type=str, 
        help="Path to the output folder. Default: %(default)s", default="nanomotif"
    )
    add_methylation_threshold_arguments(parser_find_motifs_bin_options)

    add_motif_search_arguments(parser_find_motifs_bin_options)
    

    parser_find_motifs_bin_options.add_argument(
        "--threshold_valid_coverage", type=int, default=5, 
        help="Minimum number of valid reads (Nvalid_cov) for a position to be used. Default: %(default)s"
    )

    parser_find_motifs_bin_options.add_argument(
        "--pileup_cache", type=str, default=None, metavar="FILE",
        help="Keep the filtered pileup in this file (one compressed Arrow file) and reuse it in later runs "
             "with the same pileup, contig bins, --threshold_valid_coverage and --methylation_threshold_high, "
             "skipping the pileup pass. Off by default (the filtered pileup is temporary)."
    )

    parser_find_motifs_bin_options.add_argument(
        "--min_motifs_bin", type=int, default=50, 
        help="Minimum number of scored motif occurrences (n_mod + n_nomod) in a bin for a motif to be reported. Default: %(default)s"
    )

    add_general_arguments(parser_find_motifs_bin)


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
    parser_mtase_linker_run.add_argument("--minimum_motif_methylation", type=float, default=0.5, help="Minimum fraction of motif occurrences in the bin that must be methylated for the motif to be considered. Default is 0.5")

    parser_mtase_linker_install = mtase_linker_subparsers.add_parser(
        'install', 
        help="Install additional dependencies for MTase-linker"
    )
    parser_mtase_linker_install.add_argument("-d", "--dependency_dir", type=str, default=os.getcwd(), help="Path to installation location of dependencies. A folder named ML_dependencies will be generated. Default is cwd.")
    


    ###########################################################################
    # Check installation
    parser_check_installation = subparsers.add_parser('check_installation', help="Performs small test run to verify that the installation is correct.")
    
    
    return parser
    

class CustomHelpFormatter(argparse.HelpFormatter):
    def add_arguments(self, actions):
        # Group by container to maintain ordering
        action_groups = {}
        for action in actions:
            container = getattr(action, 'container', None)
            if container not in action_groups:
                action_groups[container] = []
            action_groups[container].append(action)

        for container, group_actions in action_groups.items():
            super().add_arguments(group_actions)
