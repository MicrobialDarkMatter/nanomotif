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




def find_motifs(args, pl,  pileup = None, assembly = None, min_mods_pr_contig = 50, min_mod_frequency = 10000):
    """
    Nanomotif motif finder module

    Args:
        args (argparse.Namespace): Arguments
        pileup (pandas.DataFrame): Pileup data
        assembly (nanomotif.Assembly): Assembly data
    
    Returns:
        pandas.DataFrame: Motif data
    """
    import polars as pl

    log.info("Starting nanomotif motif finder")

    # Assembly
    if assembly is None:
        log.info("Loading assembly")
        assembly = nm.fasta.load_fasta(args.assembly)

    # Pileup 
    if pileup is None:
        log.info("Loading pileup")
        if not os.path.exists(args.pileup):
            log.error(f"File {args.pileup} does not exist")
            return None
        pileup = nm.load_pileup(args.pileup, min_coverage = args.threshold_valid_coverage, min_fraction = 0)


    
    pileup = pileup.pileup.with_columns([
        (pl.col("contig") + "_" + pl.col("mod_type")).alias("contig_mod")
    ])
    contig_mods_to_keep, contig_mods_to_remove = nm.dataload.extract_contig_mods_with_sufficient_information(pileup.filter(pl.col("fraction_mod") > args.methylation_threshold_high), assembly, min_mods_pr_contig, min_mod_frequency)
    if len(contig_mods_to_keep) == 0:
        log.info("No contigs with sufficient information")
        return None

    log.debug(f"Filtering pileup to keep contigs with more than {min_mods_pr_contig} mods and mod frequency of 1 pr. {min_mod_frequency}")
    pileup = pileup.filter(pl.col("contig_mod").is_in(contig_mods_to_keep))

    # Writing temperary files
    os.makedirs(args.out + "/temp/", exist_ok=True)
    pileup.write_csv(args.out + "/temp/filtered_pileup.tsv", separator="\t")
    with open(args.out + '/temp/contig_mod_combinations_not_processed.tsv', 'w') as f:
        for line in contig_mods_to_remove:
            f.write(f"{line}\n")
    with open(args.out + '/temp/contig_mod_combinations_processed.tsv', 'w') as f:
        for line in contig_mods_to_keep:
            f.write(f"{line}\n")
    log.info("Identifying motifs")
    motifs = nm.find_motifs.process_contig_sample_parallel(
            assembly, pileup, 
            threads = args.threads,
            search_frame_size = args.search_frame_size,

            minimum_kl_divergence = args.minimum_kl_divergence,
            verbose = args.verbose,
            log_dir = args.out + "/logs",
            seed = args.seed
        )
    motifs = pl.DataFrame(motifs)
    if motifs is None or len(motifs) == 0:
        log.info("No motifs were identified")
        return

    log.info("Writing motifs")
    def format_motif_df(df):
        if "model" in df.columns:
            n_mod = [x._alpha for x in df["model"]]
            n_nomod = [x._beta for x in df["model"]]
        motif_iupac = [nm.seq.regex_to_iupac(x) for x in df["motif"]]
        motif_type = [nm.utils.motif_type(x) for x in motif_iupac]

        df_out = df.with_columns([
            pl.Series("motif", motif_iupac),
            pl.Series("motif_type", motif_type)
        ])
        if "model" in df.columns:
            df_out = df_out.with_columns([
                pl.Series("n_mod", n_mod),
                pl.Series("n_nomod", n_nomod)
            ])
        try:
            df_out = df_out.select([
                "contig", "motif", "mod_position", "mod_type", "n_mod", "n_nomod", "motif_type",
                "motif_complement", "mod_position_complement", "n_mod_complement", "n_nomod_complement"
            ])
        except:
            df_out = df_out.select([
                "contig", "motif", "mod_position", "mod_type", "n_mod", "n_nomod", "motif_type",
            ])
        # Subtract prior alpha and beta from n_mod and n_nomod
        if "model" in df.columns:
            df_out = df_out.with_columns([
                pl.col("n_mod") - nm.model.DEFAULT_PRIOR_BETA,
                pl.col("n_nomod") - nm.model.DEFAULT_PRIOR_ALPHA
            ])
        return df_out
    def save_motif_df(df, name):
        df = format_motif_df(df)
        df = df.sort(["contig", "mod_type", "motif"])
        df.write_csv(args.out + "/" + name + ".tsv", separator="\t")
    os.makedirs(args.out + "/precleanup-motifs/", exist_ok=True)
    motifs.drop("model").write_csv(args.out + "/precleanup-motifs/motifs-raw-unformatted.tsv", separator="\t")
    save_motif_df(motifs, "precleanup-motifs/motifs-raw")

    log.info("Postprocessing motifs")
    motifs_file_name = "precleanup-motifs/motifs"

    log.info(" - Writing motifs")
    motifs = motifs.filter(pl.col("score") > args.min_motif_score)
    if len(motifs) == 0:
        log.info("No motifs found")
        return
    
    motifs_file_name = motifs_file_name + "-score"
    save_motif_df(motifs, motifs_file_name)

    log.info(" - Removing sub motifs")
    motifs = nm.postprocess.remove_sub_motifs(motifs)
    if len(motifs) == 0:
        log.info("No motifs found")
        return
    motifs_file_name = motifs_file_name +   "-sub"
    save_motif_df(motifs, motifs_file_name)

    log.info(" - Removing noisy motifs")
    motifs = nm.postprocess.remove_noisy_motifs(motifs)
    if len(motifs) == 0:
        log.info("No motifs found")
        return
    motifs_file_name = motifs_file_name +   "-noise"
    save_motif_df(motifs, motifs_file_name)

    log.info(" - Merging motifs")
    motifs = nm.postprocess.merge_motifs_in_df(motifs, pileup, assembly)
    if len(motifs) == 0:
        log.info("No motifs found")
        return
    motifs_file_name = motifs_file_name +   "-merge"
    save_motif_df(motifs, motifs_file_name)

    log.info(" - Joining motif complements")
    motifs = nm.postprocess.join_motif_complements(motifs)
    if len(motifs) == 0:
        log.info("No motifs found")
        return
    motifs_file_name = motifs_file_name +   "-complement"
    save_motif_df(motifs, motifs_file_name)

    log.info(" - Removing motifs observed less than min count")
    motifs = motifs.filter((pl.col("n_mod") + pl.col("n_nomod")) > args.min_motifs_contig)
    if len(motifs) == 0:
        log.info("No motifs found")
        return
    motifs_file_name = motifs_file_name +   "-mincount"
    save_motif_df(motifs, motifs_file_name)



    save_motif_df(motifs, "motifs")

    log.info("Done finding motifs")
    return format_motif_df(motifs)

def find_motifs_bin(args, pl,  pileup = None, assembly = None, min_mods_pr_contig = 50, min_mod_frequency = 10000):
    """
    Nanomotif motif finder module

    Args:
        args (argparse.Namespace): Arguments
        pileup (pandas.DataFrame): Pileup data
        assembly (nanomotif.Assembly): Assembly data
    
    Returns:
        pandas.DataFrame: Motif data
    """
    import polars as pl

    log.info("Starting nanomotif motif finder")
    # Bin contig relationsship
    bin_contig = nm.fasta.generate_contig_bin(args)

    # Assembly
    if assembly is None:
        log.info("Loading assembly")
        assembly = nm.fasta.load_fasta(args.assembly)

    # Pileup 
    if pileup is None:
        log.info("Loading pileup")
        if not os.path.exists(args.pileup):
            log.error(f"File {args.pileup} does not exist")
            return None
        pileup = nm.load_pileup(args.pileup, min_coverage = args.threshold_valid_coverage, min_fraction = 0)


    
    pileup = pileup.pileup.with_columns([
        (pl.col("contig") + "_" + pl.col("mod_type")).alias("contig_mod")
    ])
    contig_mods_to_keep, contig_mods_to_remove = nm.dataload.extract_contig_mods_with_sufficient_information(pileup.filter(pl.col("fraction_mod") > args.methylation_threshold_high), assembly, min_mods_pr_contig, min_mod_frequency)
    if len(contig_mods_to_keep) == 0:
        log.info("No contigs with sufficient information")
        return None

    log.debug(f"Filtering pileup to keep contigs with more than {min_mods_pr_contig} mods and mod frequency of 1 pr. {min_mod_frequency}")
    pileup = pileup.filter(pl.col("contig_mod").is_in(contig_mods_to_keep))

    # Writing temperary files
    os.makedirs(args.out + "/temp/", exist_ok=True)
    pileup.write_csv(args.out + "/temp/filtered_pileup.tsv", separator="\t")
    with open(args.out + '/temp/contig_mod_combinations_not_processed.tsv', 'w') as f:
        for line in contig_mods_to_remove:
            f.write(f"{line}\n")
    with open(args.out + '/temp/contig_mod_combinations_processed.tsv', 'w') as f:
        for line in contig_mods_to_keep:
            f.write(f"{line}\n")

    log.info("Identifying motifs")
    motifs = nm.find_motifs_bin.process_binned_sample_parallel(
            assembly, pileup, bin_contig,
            threads = args.threads,
            search_frame_size = args.search_frame_size,
            methylation_threshold_high = args.methylation_threshold_high,
            methylation_threshold_low = args.methylation_threshold_low,
            minimum_kl_divergence = args.minimum_kl_divergence,
            verbose = args.verbose,
            log_dir = args.out + "/logs",
            seed = args.seed
        )
    motifs = pl.DataFrame(motifs)
    if motifs is None or len(motifs) == 0:
        log.info("No motifs were identified")
        return

    log.info("Writing motifs")
    def format_motif_df(df):
        if "model" in df.columns:
            n_mod = [x._alpha for x in df["model"]]
            n_nomod = [x._beta for x in df["model"]]
        motif_iupac = [nm.seq.regex_to_iupac(x) for x in df["motif"]]
        motif_type = [nm.utils.motif_type(x) for x in motif_iupac]

        df_out = df.with_columns([
            pl.Series("motif", motif_iupac),
            pl.Series("motif_type", motif_type)
        ])
        if "model" in df.columns:
            df_out = df_out.with_columns([
                pl.Series("n_mod", n_mod),
                pl.Series("n_nomod", n_nomod)
            ])
        try:
            df_out = df_out.select([
                "bin", "motif", "mod_position", "mod_type", "n_mod", "n_nomod", "motif_type",
                "motif_complement", "mod_position_complement", "n_mod_complement", "n_nomod_complement"
            ])
        except:
            df_out = df_out.select([
                "bin", "motif", "mod_position", "mod_type", "n_mod", "n_nomod", "motif_type",
            ])
        # Subtract prior alpha and beta from n_mod and n_nomod
        if "model" in df.columns:
            df_out = df_out.with_columns([
                pl.col("n_mod") - nm.model.DEFAULT_PRIOR_BETA,
                pl.col("n_nomod") - nm.model.DEFAULT_PRIOR_ALPHA
            ])
        return df_out
    def save_motif_df(df, name):
        df = format_motif_df(df)
        df = df.sort(["bin", "mod_type", "motif"])
        df.write_csv(args.out + "/" + name + ".tsv", separator="\t")
    os.makedirs(args.out + "/precleanup-motifs/", exist_ok=True)
    motifs.drop("model").write_csv(args.out + "/precleanup-motifs/motifs-raw-unformatted.tsv", separator="\t")
    save_motif_df(motifs, "precleanup-motifs/motifs-raw")

    log.info("Postprocessing motifs")
    motifs_file_name = "precleanup-motifs/motifs"

    log.info(" - Writing motifs")
    motifs = motifs.filter(pl.col("score") > args.min_motif_score)
    if len(motifs) == 0:
        log.info("No motifs found")
        return
    
    motifs_file_name = motifs_file_name + "-score"
    save_motif_df(motifs, motifs_file_name)


    log.info(" - Removing noisy motifs")
    motifs = nm.postprocess.remove_noisy_motifs(motifs)
    if len(motifs) == 0:
        log.info("No motifs found")
        return
    motifs_file_name = motifs_file_name +   "-noise"
    save_motif_df(motifs, motifs_file_name)

    log.info(" - Merging motifs")
    motifs = nm.find_motifs_bin.merge_motifs_in_df(motifs, pileup, assembly, bin_contig)
    if len(motifs) == 0:
        log.info("No motifs found")
        return
    motifs_file_name = motifs_file_name +   "-merge"
    save_motif_df(motifs, motifs_file_name)

    log.info(" - Removing sub motifs")
    motifs = motifs.rename({"bin":"contig"})
    motifs = nm.postprocess.remove_sub_motifs(motifs)
    if len(motifs) == 0:
        log.info("No motifs found")
        return
    motifs_file_name = motifs_file_name +   "-sub"
    motifs = motifs.rename({"contig":"bin"})
    save_motif_df(motifs, motifs_file_name)

    log.info(" - Joining motif complements")
    motifs = motifs.rename({"bin":"contig"})
    motifs = nm.postprocess.join_motif_complements(motifs)
    if len(motifs) == 0:
        log.info("No motifs found")
        return
    motifs_file_name = motifs_file_name +   "-complement"
    motifs = motifs.rename({"contig":"bin"})
    save_motif_df(motifs, motifs_file_name)

    log.info(" - Removing motifs observed less than min count")
    motifs = motifs.filter((pl.col("n_mod") + pl.col("n_nomod")) > args.min_motifs_bin)
    if len(motifs) == 0:
        log.info("No motifs found")
        return
    motifs_file_name = motifs_file_name +   "-mincount"
    save_motif_df(motifs, motifs_file_name)



    save_motif_df(motifs, "bin-motifs")

    log.info("Done finding motifs")
    return format_motif_df(motifs)


def motif_discovery_legacy(args, pl):
    # Check if all required files exist
    if not os.path.exists(args.pileup):
        log.error(f"File {args.pileup} does not exist")
        return
    if not os.path.exists(args.assembly):
        log.error(f"File {args.assembly} does not exist")
        return
    if not os.path.exists(args.bins):
        log.error(f"File {args.bins} does not exist")
        return
    empty_find_motifs = pl.DataFrame({
        "contig": [],
        "motif": [],
        "mod_position": [],
        "mod_type": [],
        "n_mod": [],
        "n_nomod": [],
        "motif_type": [],
        "motif_complement": [],
        "mod_position_complement": [],
        "n_mod_complement": [],
        "n_nomod_complement": []
    })
    empty_bin_motifs = pl.DataFrame({
        "bin": [],
        "mod_type": [],
        "motif": [],
        "mod_position": [],
        "n_mod_bin": [],
        "n_nomod_bin": [],
        "motif_type": [],
        "motif_complement": [],
        "mod_position_complement": [],
        "n_mod_complement": [],
        "n_nomod_complement": []
    })
    
    # Check if output directory exists
    log.info("Loading required files")
    log.debug("Loading pileup")
    if args.read_level_methylation:
        pileup = nm.load_pileup(args.pileup,min_coverage = args.threshold_valid_coverage, min_fraction = 0)
    else:
        pileup = nm.load_pileup(args.pileup,min_coverage = args.threshold_valid_coverage, min_fraction = args.threshold_methylation_general)
        log.debug("Loading assembly")
    assembly = nm.fasta.load_fasta(args.assembly)

    # Find motifs
    log.info("Finding motifs")
    motifs = find_motifs(args, pl, pileup=pileup, assembly=assembly)
    if motifs is None:
        empty_find_motifs.write_csv(args.out + "/motifs.tsv", separator="\t")
        empty_bin_motifs.write_csv(args.out + "/bin-motifs.tsv", separator="\t")
        log.info("Stopping workflow")
        return

    # Score all motifs
    log.info("Scoring motifs")
    scored_all = score_motifs(args, pl, pileup=pileup, assembly=assembly, motifs=motifs)

    # Bin consensus
    log.info("Finding bin consensus motifs")
    bin_consensus(args, pl, pileup=pileup, assembly=assembly, motifs=motifs, motifs_scored=scored_all)

    log.info("Done")

def check_install(args, pl):
    
    # Check if output directory exists
    log.info("Loading required files")
    args.out = "nanomotif_install_check"
    args.save_motif_positions = False
    args.pileup = nm.datasets.geobacillus_plasmids_pileup_path()

    pileup = nm.datasets.geobacillus_plasmids_pileup()
    assembly = nm.datasets.geobacillus_plasmids_assembly()

    # Find motifs
    log.info("Finding motifs")
    motifs = find_motifs(args, pl, pileup=pileup, assembly=assembly)

    # Score all motifs
    log.info("Scoring motifs")
    scored_all = score_motifs(args, pl, pileup=pileup, assembly=assembly, motifs=motifs)

    # Bin consensus
    log.info("Finding bin consensus motifs")
    args.bins = nm.datasets.geobacillus_plasmids_bin_path()
    bin_consensus(args, pl, pileup=pileup, assembly=assembly, motifs=motifs, motifs_scored=scored_all)

    # Find motifs bin
    log.info("Finding bin motifs")
    find_motifs_bin(args, pl, pileup=pileup, assembly=assembly)
    
    log.info("Done")
    for _ in range(3):
        try:
            shutil.rmtree(args.out)
            break
        except OSError as e:
            time.sleep(0.5) 



####################################################################################################
# Binnary - contamination and inclusion
from nanomotif.binnary import data_processing, detect_contamination, include_contigs
from nanomotif.binnary.logging import set_logger_config
from pymethylation_utils.utils import run_epimetheus


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

    contig_methylation_file = "motifs-scored-read-methylation.tsv"
    if not args.force:
        log.info(f"Check if {contig_methylation_file} exists")

    if os.path.isfile(os.path.join(args.out,contig_methylation_file)) and not args.force:
        log.info("motifs-scored-read-methylation.tsv exists. Using existing file! Use --force to override this.")
    elif not os.path.isfile(os.path.join(args.out, contig_methylation_file)) or args.force:
        log.info(f"Running epimetheus to create {contig_methylation_file}")
        # Create motifs-scored-read-methylation
        return_code = run_epimetheus(
            pileup = args.pileup,
            assembly = args.assembly,
            motifs = motifs_in_bin_consensus,
            threads = args.threads,
            min_valid_read_coverage = args.min_valid_read_coverage,
            output = os.path.join(args.out,contig_methylation_file)
        )

        if return_code != 0:
            log.error("Error running epimetheus")


    log.info("Loading assembly file...")
    assembly = data_processing.read_fasta(args.assembly)
    contig_lengths = data_processing.find_contig_lengths(assembly) 

    # Load motifs-scored-read-methylation.tsv
    contig_methylation = pl.read_csv(
        os.path.join(args.out, contig_methylation_file), separator="\t", has_header = True, schema = {
            'contig': pl.String(), 'motif': pl.String(), 'mod_type': pl.String(), 'mod_position': pl.Int8(), 'median': pl.Float64(),  'mean_read_cov': pl.Float64(), 'N_motif_obs': pl.Int32(), 'motif_occurences_total': pl.Int32(),
        }
    )

    contig_methylation = contig_methylation\
        .filter((pl.col("N_motif_obs").cast(pl.Float64) * pl.col("mean_read_cov")) >= args.methylation_threshold)

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
    except:
        os.environ["POLARS_MAX_THREADS"] = "1"
    import polars as pl
    if args.command in ["detect_contamination", "include_contigs", "MTase-linker"]:
        args.verbose = False
        args.seed = 1
    
    if args.command == "find_motifs":
        shared_setup(args, args.out)
        find_motifs(args, pl)
    elif args.command == "motif_discovery":
        shared_setup(args, args.out)
        find_motifs_bin(args, pl)
    

    elif args.command in ["detect_contamination", "include_contigs"]:
        shared_setup(args, args.out)
        binnary(args, pl)

    elif args.command == "MTase-linker":
        mtase_linker(args)

    elif args.command == "check_installation":
        
        outdir = "tests/cli_test_motif_discovery"
        cmd = [
            "nanomotif", "find_motifs",
            "-t", "1",
            "nanomotif/datasets/geobacillus-plasmids.assembly.fasta",
            "nanomotif/datasets/geobacillus-plasmids.pileup.bed",
            # "-c", "nanomotif/datasets/geobacillus-contig-bin.tsv",
            "--out", outdir
        ]
        subprocess.run(cmd)
        shutil.rmtree(outdir, ignore_errors=True)

    else:
        parser.print_help()
        exit()
if __name__ == "__main__":
    main()
