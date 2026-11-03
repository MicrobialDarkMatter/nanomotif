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
import multiprocessing as mp


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
    motif_discovery_parralel = nm.find_motifs_bin.BinMotifProcessorBuilder(config)
    if not motif_discovery_parralel:
        log.error("No motifs were identified")
        return None
    motif_discovery_parralel = motif_discovery_parralel.build()
    motifs = motif_discovery_parralel.run()
    if motifs is None or len(motifs) == 0:
        log.info("No motifs were identified")
        return None
    else:
        log.debug(f"Identified {len(motifs)} motifs in {len(motifs['reference'].unique())} bins")

    motifs = motifs.filter(pl.col("n_mod") + pl.col("n_nomod") >= args.min_motifs_bin)
    if len(motifs) == 0:
        log.info("Motif frequency of all motifs too low")
        return None
    log.info("Writing motifs")
    motifs.write_motif_formatted(args.out + "/bin-motifs.tsv")
    log.info("Done finding motifs")
    log.info(f"Identified {len(motifs)} motifs in {len(motifs['reference'].unique())} bins")
    return motifs


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
    if args.command == "MTase-linker":
        args.verbose = False
        args.seed = 1
    
    if args.command == "motif_discovery":
        shared_setup(args, args.out)
        result = find_motifs_bin(args, pl)
        if result is None:
            # Write empty motif file
            empty_motif_path = os.path.join(args.out, "bin-motifs.tsv")
            with open(empty_motif_path, 'w') as f:
                f.write("reference\tmotif\tmod_position\tmod_type\tn_mod\tn_nomod\tmotif_type\tmotif_complement\tmod_position_complement\tn_mod_complement\tn_nomod_complement\n")

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
    mp.set_start_method("spawn", force=True)
    main()
