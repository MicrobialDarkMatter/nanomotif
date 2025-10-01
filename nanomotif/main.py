import nanomotif as nm
from nanomotif.logger import configure_logger
from nanomotif._version import __version__
import logging as log
import os
import sys
import shutil
from pathlib import Path
import json
import time
import numpy as np
import random
import warnings
import subprocess


def shared_setup(args, working_dir):
    # Check if output directory exists
    if not os.path.exists(args.out):
        os.makedirs(args.out)
    else:
        log.warning(f"Output directory {args.out} already exists")

    # Set up logging
    LOG_DIR = working_dir + "/logs"
    Path(LOG_DIR).mkdir(parents=True, exist_ok=True)
    configure_logger(LOG_DIR + f"/{args.command}.main.log", args.verbose, stdout=True)
    warnings.filterwarnings("ignore")


    # Log arguments
    log.info(f"nanomotif version: {__version__}")
    log.debug("Arguments:")
    for arg in vars(args):
        log.debug(f"{arg}: {getattr(args, arg)}")
    
    # Save settings
    with open(working_dir + f"/args.{args.command}.json", "w") as f:
        json.dump(vars(args), f, indent=2)
    
    # Set seeds
    nm.seed.set_seed(args.seed)


def find_motifs_bin(args, pl, min_mods_pr_contig = 50, min_mod_frequency = 10000):
    """
    Nanomotif motif finder module

    Args:
        args (argparse.Namespace): Arguments
        assembly (nanomotif.Assembly): Assembly data
    
    Returns:
        pandas.DataFrame: Motif data
    """
    import polars as pl

    log.info("Starting nanomotif motif finder")
    bin_contig = nm.fasta.generate_contig_bin(args)
    if bin_contig is None or len(bin_contig) == 0:
        log.error("No bin contig mapping found")
        return

    log.info("Loading assembly")
    assembly = nm.fasta.load_fasta_fastx(args.assembly)

    log.info("Identifying motifs")
    config = nm.find_motifs_bin.ProcessorConfig(
        assembly = assembly,
        pileup_path = args.pileup,
        bin_contig = bin_contig,
        threads = args.threads,
        search_frame_size = args.search_frame_size,
        methylation_threshold_low = args.methylation_threshold_low,
        methylation_threshold_high = args.methylation_threshold_high,
        minimum_kl_divergence = args.minimum_kl_divergence,
        score_threshold = args.min_motif_score,
        verbose = args.verbose,
        log_dir = args.out + "/logs",
        seed = args.seed,
        output_dir = args.out
    )
    motif_discovery_parralel = nm.find_motifs_bin.BinMotifProcessorBuilder(config).build()
    motifs = motif_discovery_parralel.run()
    if motifs is None or len(motifs) == 0:
        log.info("No motifs were identified")
        return None
    else:
        log.debug(f"Identified {len(motifs)} motifs in {len(motifs['reference'].unique())} bins")

    motifs = motifs.filter(pl.col("n_mod") + pl.col("n_nomod") >= args.min_motifs_bin)
    log.info("Writing motifs")
    motifs.write_motif_formatted(args.out + "/bin-motifs.tsv")
    log.info("Done finding motifs")
    log.info(f"Identified {len(motifs)} motifs in {len(motifs['reference'].unique())} bins")
    return motifs


####################################################################################################
# Binnary - contamination and inclusion
from nanomotif.binnary import data_processing, detect_contamination, include_contigs
from nanomotif.binnary.logging import set_logger_config
from epymetheus.epymetheus import methylation_pattern, MethylationOutput


def binnary(args, pl):
    """
    binnary entry point for the DNA Methylation Pattern Analysis tool.
    Orchestrates the workflow of the tool based on the provided arguments.
    """
    
    print("Starting Binnary ", args.command, " analysis...")

    
    # Settting up the logger
    log.info("Starting Binnary analysis...")
    
    # Step 1: Load and preprocess data
    # These functions would be defined in your data_processing module
    (
        bin_motifs,
        contig_bins,
    ) = data_processing.load_data(args)

    

    # Step 2: create motifs_scored_in_bins and bin_consensus_motifs_filtered
    motifs_in_bin_consensus = bin_motifs\
        .with_columns(
            (pl.col("motif") + "_" + pl.col("mod_type") + "_" + pl.col("mod_position").cast(pl.Utf8)).alias("motif_mod")            
        )\
        .get_column("motif_mod").unique()

    if args.methylation_output_type == "median":
        output_type = MethylationOutput.Median
    else:
        output_type = MethylationOutput.WeightedMean

    contig_methylation_file = f"motifs-scored-read-methylation_{args.methylation_output_type}.tsv"
    if not args.force:
        log.info(f"Check if {contig_methylation_file} exists")

    if os.path.isfile(os.path.join(args.out,contig_methylation_file)) and not args.force:
        log.info("motifs-scored-read-methylation.tsv exists. Using existing file! Use --force to override this.")

    elif not os.path.isfile(os.path.join(args.out, contig_methylation_file)) or args.force:
        log.info(f"Running epymetheus to create {contig_methylation_file}")
        # Create motifs-scored-read-methylation
        return_code = methylation_pattern(
            pileup = args.pileup,
            assembly = args.assembly,
            motifs = motifs_in_bin_consensus,
            threads = args.threads,
            min_valid_read_coverage = args.min_valid_read_coverage,
            batch_size=1000,
            min_valid_cov_to_diff_fraction=0.8,
            output = os.path.join(args.out,contig_methylation_file),
            allow_assembly_pileup_mismatch=False,
            output_type = output_type
        )

        if return_code != 0:
            log.error("Error running epymetheus")


    log.info("Loading assembly file...")
    assembly = data_processing.read_fasta(args.assembly)
    contig_lengths = data_processing.find_contig_lengths(assembly) 

    # Load motifs-scored-read-methylation.tsv
    contig_methylation = pl.read_csv(
        os.path.join(args.out, contig_methylation_file), separator="\t", has_header = True, schema = {
            'contig': pl.String(), 'motif': pl.String(), 'mod_type': pl.String(), 'mod_position': pl.Int8(), 'methylation_value': pl.Float64(),  'mean_read_cov': pl.Float64(), 'n_motif_obs': pl.Int32(),
        }
    )

    contig_methylation = contig_methylation\
        .filter((pl.col("n_motif_obs").cast(pl.Float64) * pl.col("mean_read_cov")) >= args.methylation_threshold)

    # Setting up the contamination analysis
    if (args.command == "detect_contamination" and not args.contamination_file) or (args.command == "include_contigs" and args.run_detect_contamination):
    
        contig_methylation_cont = data_processing.add_bin(
            contig_methylation,
            contig_bins
        )

        contamination = detect_contamination.detect_contamination(
            contig_methylation_cont, contig_lengths, args.num_consensus, args.threads
        )


        data_processing.generate_output(contamination.to_pandas(), args.out, "bin_contamination.tsv")

    if args.command == "detect_contamination":
        # Create a decontaminated contig_bin file
        new_contig_bins = data_processing.create_contig_bin_file(
            contig_bins=contig_bins.to_pandas(), 
            contamination=contamination.to_pandas(),
            include=None
        )
        data_processing.generate_output(new_contig_bins, args.out, "decontaminated_contig_bin.tsv")

    if args.command == "include_contigs":
        # User provided contamination file
        if args.contamination_file:
            log.info("Loading contamination file...")
            contamination = data_processing.load_contamination_file(args.contamination_file)
        # If run_detect_contamination is false and contamination file is not provided, then set contamination to None
        if not args.run_detect_contamination and not args.contamination_file:
            contamination = pl.DataFrame(
                {
                    "contig": []
                }
            )
        
        # move contigs in bin_contamination to unbinned
        contamination_contigs = contamination.get_column("contig")

        log.info("Removing contaminants from bins")
        contig_bins = contig_bins.filter(~pl.col("contig").is_in(contamination_contigs))
        
        contig_methylation_inc = data_processing.add_bin(
            contig_methylation,
            contig_bins
        )

        # Run the include_contigs analysis    
        include_contigs_df = include_contigs.include_contigs(
            contig_methylation_inc,
            contig_lengths,
            mean_probability = args.mean_model_confidence
        )
        
        # Save the include_contigs_df results
        data_processing.generate_output(include_contigs_df.to_pandas(), args.out, "include_contigs.tsv")
        
        # Create a new contig_bin file
        uniquely_assigned_contigs = include_contigs_df\
            .filter(pl.col("confidence") == "high_confidence")\
            .drop("bin")\
            .rename({"assigned_bin": "bin"})\
            .select(["contig", "bin"])\
            .unique()

        new_contig_bins = data_processing.create_contig_bin_file(
            contig_bins=contig_bins.to_pandas(), 
            contamination= contamination.to_pandas(),
            include=uniquely_assigned_contigs.to_pandas()
        )
        data_processing.generate_output(new_contig_bins, args.out, "new_contig_bin.tsv")
        
    if args.write_bins:
        log.info("Write bins flag is set. Writing bins to file...")
        
        bin_dir = os.path.join(args.out, args.command + "_bins")
        
        data_processing.write_bins_from_contigs(new_contig_bins, assembly, bin_dir)
    
    log.info(f"Analysis Completed. Results are saved to: {args.out}")
    print("Analysis Completed. Results are saved to:", args.out)

####################################################################################################
# MTaser-linker

from nanomotif.mtase_linker.dependencies import snakemake_create_environments, get_models, defensefinder_update, check_installation_MTase_linker
from nanomotif.mtase_linker.command import run_MTase_linker

def mtase_linker(args):
    if args.mtase_linker_command == "install":
        snakemake_create_environments(args)
        get_models(args)
        defensefinder_update(args)
        check_installation_MTase_linker(args)
    elif args.mtase_linker_command == "run":
        run_MTase_linker(args)
    else:
        print(f"Unknown MTase-linker command: {args.mtase_linker_command}")
       



def main():
    # Parse arguments
    parser = nm.argparser.create_parser()
    args = parser.parse_args()

    try: 
        os.environ["POLARS_MAX_THREADS"] = str(args.threads)
        os.environ.setdefault("RAYON_NUM_THREADS", "1")
    except:
        os.environ["POLARS_MAX_THREADS"] = "1"
        os.environ.setdefault("RAYON_NUM_THREADS", "1")
    import polars as pl
    if args.command in ["detect_contamination", "include_contigs", "MTase-linker"]:
        args.verbose = False
        args.seed = 1
    
    if args.command == "motif_discovery":
        shared_setup(args, args.out)
        find_motifs_bin(args, pl)

    elif args.command in ["detect_contamination", "include_contigs"]:
        shared_setup(args, args.out)
        binnary(args, pl)

    elif args.command == "MTase-linker":
        mtase_linker(args)

    elif args.command == "check_installation":
        fasta = nm.datasets.geobacillus_plasmids_assembly_path()
        pileup = nm.datasets.geobacillus_plasmids_pileup_path()
        contig_bin = nm.datasets.geobacillus_plasmids_bin_path()
        outdir = "nanomotif_check_installation_output"
        cmd = [
            "nanomotif", "motif_discovery",
            "-t", "1",
            fasta,
            pileup,
            "-c", contig_bin,
            "--out", outdir
        ]
        subprocess.run(cmd)
        shutil.rmtree(outdir, ignore_errors=True)

    else:
        parser.print_help()
        exit()
if __name__ == "__main__":
    main()
