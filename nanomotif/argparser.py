import argparse
from nanomotif._version import __version__
import os


def create_parser():
    formatter = lambda prog: argparse.HelpFormatter(prog, max_help_position=28)
    parser = argparse.ArgumentParser(description="Motif identification and utilisation commands", formatter_class=formatter)
    parser.add_argument("--version", action="version", version="%(prog)s {}".format(__version__))
    subparsers = parser.add_subparsers(help="-- Command descriptions --", dest="command")

    def add_general_arguments(parser):
        """
        Function to add general arguments to a parser.
        """
        group = parser.add_argument_group("General Arguments")
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
        group_major = parser.add_argument_group("Contig Bin Arguments")
        group = group_major.add_mutually_exclusive_group(required=True)


        group.add_argument(
            "-c", "--contig_bin", type=str, 
            help="TSV file specifying which bin contigs belong."
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
            help="File extension of the bin FASTA files. Default is '.fasta'."
        )


    
    def add_methylation_threshold_arguments(parser):
        """
        Function to add arguments for methylation threshold.
        """
        parser.add_argument(
            "--methylation_threshold_low", type=float, default=0.30, 
            help="A position is considered non-methylated if fraction of methylation is below this threshold. Default: %(default)s"
        )

        parser.add_argument(
            "--methylation_threshold_high", type=float, default=0.70, 
            help="A position is considered methylated if fraction of methylated reads is above this threshold. Default: %(default)s"
        )

    def add_motif_search_arguments(parser):
        parser.add_argument(
            "--search_frame_size", type=int, default=40, 
            help="length of the sequnces sampled around confident methylation sites. Default: %(default)s"
        )
        
        parser.add_argument(
            "--minimum_kl_divergence", type=float, default=0.05, 
            help="Minimum KL-divergence for a position to considered for expansion in  motif search. Higher value means less exhaustive, but faster search. Default: %(default)s"
        )
        
        parser.add_argument(
            "--min_motif_score", type=float, default=0.2, 
            help="Minimum score for a motif to be kept after identification. Default: %(default)s"
        )   

    ###########################################################################
    # Find Motifs
    parser_find_motifs = subparsers.add_parser(
        'find_motifs', 
        help="Finds motifs directly on contig level in provided assembly",
        add_help=False
    )

    parser_find_motifs_options = parser_find_motifs.add_argument_group("Options") 

    parser_find_motifs.add_argument(
        "assembly", type=str, 
        help="path to the assembly file."
    )

    parser_find_motifs.add_argument(
        "pileup", type=str, 
        help="path to the modkit pileup file."
    )   
    
    parser_find_motifs_options.add_argument(
        "--out", type=str, 
        help="path to the output folder", default="nanomotif"
    )
    
    add_methylation_threshold_arguments(parser_find_motifs_options)

    add_motif_search_arguments(parser_find_motifs_options)

    parser_find_motifs_options.add_argument(
        "--threshold_valid_coverage", type=int, default=5, 
        help="minimum valid base coverage for a position to be considered. Default: %(default)s"
    )
    
    parser_find_motifs_options.add_argument(
        "--min_motifs_contig", type=int, default=20, 
        help="minimum number of times a motif has to have been oberserved in a contig. Default: %(default)s"
    )
    
    parser_find_motifs_options.add_argument(
        "--read_level_methylation", action="store_true", 
        help="If specified, methylation is calculated on read level instead of contig level. This is slower but produces more stable motifs."
    )
    
    add_general_arguments(parser_find_motifs)

    ###########################################################################
    # Find motifs on bin level
    parser_find_motifs_bin = subparsers.add_parser(
        'motif_discovery',
        help="Finds motifs directly on bin level in provided assembly",
        add_help=False
    )


    parser_find_motifs_bin.add_argument("assembly", type=str, 
        help="path to the assembly file.")
    
    parser_find_motifs_bin.add_argument("pileup", type=str, 
        help="path to the modkit pileup file.") 
    
    add_contig_bin_arguments(parser_find_motifs_bin)

    parser_find_motifs_bin_options = parser_find_motifs_bin.add_argument_group("Options")
    parser_find_motifs_bin_options.add_argument(
        "--out", type=str, 
        help="path to the output folder", default="nanomotif"
    )
    add_methylation_threshold_arguments(parser_find_motifs_bin_options)

    add_motif_search_arguments(parser_find_motifs_bin_options)
    

    parser_find_motifs_bin_options.add_argument(
        "--threshold_valid_coverage", type=int, default=5, 
        help="Minimum valid base coverage (Nvalid_cov) for a position to be considered. Default: %(default)s"
    )

    parser_find_motifs_bin_options.add_argument(
        "--min_motifs_bin", type=int, default=20, 
        help="Minimum number of motif observations in a bin. Default: %(default)s"
    )

    add_general_arguments(parser_find_motifs_bin)


    ###########################################################################
    # Bin contamination and inclusion   
    parser_binnary_shared = argparse.ArgumentParser(description="Contamination DNA Methylation Pattern", add_help=False)
    """Function to add common arguments to a subparser."""
    
    parser_binnary_shared_mandatory = parser_binnary_shared.add_argument_group("Mandatory Arguments")
    
    
    parser_binnary_shared_mandatory.add_argument(
        "--pileup", type=str, help="Path to pileup.bed", required=True
    )
    parser_binnary_shared_mandatory.add_argument(
        "--assembly", type=str, help="Path to assembly file [fasta format required]", required=True
    )
    parser_binnary_shared_mandatory.add_argument("--bin_motifs", type=str, help="Path to bin-motifs.tsv file", required=True)
    parser_binnary_shared_mandatory.add_argument(
        "--contig_bins", type=str, help="Path to bins.tsv file for contig bins", required=True
    )
    parser_binnary_shared.add_argument("-t", "--threads", type=int, default=1, help="Number of threads to use for multiprocessing")

    parser_binnary_shared.add_argument(
        "--min_valid_read_coverage",
        type=int,
        default=3,
        help="Minimum read coverage for calculating methylation [used with methylation_util executable]",
    )
    parser_binnary_shared.add_argument(
        "--methylation_threshold",
        type=int,
        default=24,
        help="Filtering criteria for trusting contig methylation. It is the product of mean_read_coverage and N_motif_observation. Higher value means stricter criteria. [default: 24]",
    )
    
    parser_binnary_shared.add_argument(
        "--num_consensus",
        type=int,
        default=4,
        help="Number of models that has to agree for classifying as contaminant",
    )
    parser_binnary_shared.add_argument(
        "--force",
        action='store_true',
        help="Force override of motifs-scored-read-methylation.tsv. If not set existing file will be used.",
    )
    
    parser_binnary_shared.add_argument(
        "--write_bins",
        action='store_true',
        help="If specified, new bins will be written to a bins folder. Requires --assembly_file to be specified.",
    )
    
    parser_binnary_shared_mandatory.add_argument("--out", type=str, help="Path to output directory", required=True, default="nanomotif")
    
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
    parser_inclusion.add_argument(
        "--mean_model_confidence",
        type = float,
        help = "Mean probability between models for including contig. Contigs above this value will be included. [default: 0.8]",
        default = 0.8
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