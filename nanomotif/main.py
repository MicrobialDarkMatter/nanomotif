import nanomotif as nm
from nanomotif.logger import configure_logger
import logging as log
import os
import sys
import shutil
from pathlib import Path
import json
os.environ["POLARS_MAX_THREADS"] = "1"
import polars as pl
from polars import col
import numpy as np
import random
import warnings
from nanomotif._version import __version__


def shared_setup(args, working_dir):
    warnings.filterwarnings("ignore")
    # Check if output directory exists
    if not os.path.exists(args.out):
        os.makedirs(args.out)
    else:
        log.warn(f"Output directory {args.out} already exists")

    # Set up logging
    LOG_DIR = working_dir + "/logs"
    Path(LOG_DIR).mkdir(parents=True, exist_ok=True)
    configure_logger(LOG_DIR + f"/{args.command}.main.log", args.verbose)

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




def find_motifs(args, pileup = None, assembly = None):
    # Run nanomotif
    log.info("Starting nanomotif motif finder")
    if pileup is None:
        log.info("Loading pileup")
        pileup = nm.load_pileup(args.pileup, threads = args.threads, min_fraction = args.threshold_methylation_general)
    if assembly is None:
        log.info("Loading assembly")
        assembly = nm.load_assembly(args.assembly)
    assm_lengths = pl.DataFrame({
        "contig": list(assembly.assembly.keys()),
        "length": [len(contig) for contig in assembly.assembly.values()]
    })
    log.info("Filtering pileup")
    # Filter pileup to contigs with mods, minimum 1 mod per 10kb
    contigs_with_mods = pileup.pileup \
        .groupby(["contig", "mod_type"]) \
        .agg(pl.count()) \
        .join(assm_lengths, on = "contig") \
        .filter(pl.col("count") > pl.col("length")/10000) \
        .get_column("contig").unique().to_list()

    total_contigs = len(assembly.assembly.keys())

    contigs_to_process = [contig for contig in assembly.assembly.keys() if contig in contigs_with_mods]
    pileup = pileup.pileup.filter(pl.col("contig").is_in(contigs_to_process))
    remaining_contigs = pileup.get_column("contig").unique().to_list()

    os.makedirs(args.out + "/temp/", exist_ok=True)
    pileup.write_csv(args.out + "/temp/filtered_pileup.tsv", separator="\t")
    log.info(f"Processing {len(remaining_contigs)} of {total_contigs} contigs")
    log.info("Identifying motifs")
    motifs = nm.evaluate.process_sample_parallel(
            assembly, pileup, 
            threads = args.threads,
            search_frame_size = args.search_frame_size,
            threshold_methylation_confident = args.threshold_methylation_confident,
            threshold_methylation_general = args.threshold_methylation_general,
            threshold_valid_coverage = args.threshold_valid_coverage,
            minimum_kl_divergence = args.minimum_kl_divergence,
            verbose = args.verbose,
            log_dir = args.out + "/logs",
            seed = args.seed
        )
    motifs = pl.DataFrame(motifs)
    if motifs is None:
        log.info("No motifs found")
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
    save_motif_df(motifs, "precleanup-motifs/motifs-raw")

    log.info("Postprocessing motifs")
    motifs_file_name = "precleanup-motifs/motifs"

    log.info(" - Writing motifs")
    motifs = motifs.filter(col("score") > 0.1)
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

def score_motifs(args, pileup = None, assembly = None, motifs = None):
    log.info("Starting nanomotif motif scorer")
    if pileup is None:
        log.info("Loading pileup")
        pileup =  nm.load_pileup(args.pileup, threads = args.threads, min_fraction = args.threshold_methylation_general)
    if assembly is None:
        log.info("Loading assembly")
        assembly = nm.load_assembly(args.assembly)
    if motifs is None:
        log.info("Loading motifs")
        motifs = pl.read_csv(args.motifs, separator="\t")
    
    # Ensure motif are iupac
    motifs.with_columns([
        pl.col("motif").map_elements(lambda x: nm.seq.regex_to_iupac(x)).alias("motif")
    ])

    
    if any([item for item in motifs.columns if "_complement" in item]):
        log.debug("Reformatting motifs to include scoring of complementary motifs")
        motifs_fwd = motifs[[item for item in motifs.columns if "complement" not in item]]
        motifs_rev = motifs[["contig", "motif_type", "mod_type"] + [item for item in motifs.columns if "complement" in item]] \
            .select(pl.all().name.map(lambda col_name: col_name.replace('_complement', ''))) \
            .select(motifs_fwd.columns) \
            .filter(pl.col("motif").is_not_null())
        
        motifs = pl.concat([motifs_fwd, motifs_rev]).unique(["motif", "contig", "mod_type", "mod_position"])
    
    # Convert to regex for motif scoring
    motifs.with_columns([
        pl.col("motif").map_elements(lambda x: nm.seq.iupac_to_regex(x)).alias("motif")
    ])

    log.info("Scoring motifs")
    scored_all = nm.scoremotifs.score_sample_parallel(
        assembly, pileup.pileup, motifs,
        threads = args.threads,
        threshold_methylation_general = args.threshold_methylation_general,
        threshold_valid_coverage = 1,
        verbose = args.verbose,
        log_dir = args.out + "/logs",
        seed = args.seed
        )
    scored_all = nm.postprocess.join_motif_complements(scored_all)

    scored_all = scored_all.unique([
            "motif", "contig", "mod_type", "mod_position"
        ])
    
    scored_all = scored_all.sort(["contig", "mod_type", "motif"])
    scored_all.write_csv(args.out + "/motifs-scored.tsv", separator="\t")
    return scored_all

def bin_consensus(args, pileup = None, assembly = None, motifs = None, motifs_scored = None):
    bins = pl.read_csv(args.bins, separator="\t", has_header=False) \
        .rename({"column_1":"contig", "column_2":"bin"})
    if motifs is None:
        log.info("Loading motifs-scored")
        motifs = pl.read_csv(args.motifs, separator="\t")
    if motifs_scored is None:
        log.info("Loading motifs")
        motifs_scored = pl.read_csv(args.motifs_scored, separator="\t")
    if pileup is None:
        log.info("Loading pileup")
        pileup =  nm.load_pileup(args.pileup, threads = args.threads, min_fraction = args.threshold_methylation_general)
    if assembly is None:
        log.info("Loading assembly")
        assembly = nm.load_assembly(args.assembly)

    
    if any([item for item in motifs_scored.columns if "_complement" in item]):
        motifs_fwd = motifs_scored[[item for item in motifs_scored.columns if "complement" not in item]]
        motifs_rev = motifs_scored[["contig", "motif_type", "mod_type"] + [item for item in motifs_scored.columns if "complement" in item]] \
            .select(pl.all().name.map(lambda col_name: col_name.replace('_complement', ''))) \
            .select(motifs_fwd.columns) \
            .filter(pl.col("motif").is_not_null())
        motifs_scored = pl.concat([motifs_fwd, motifs_rev]).unique(["motif", "contig", "mod_type", "mod_position"])
    
    if any([item for item in motifs.columns if "_complement" in item]):
        motifs_fwd = motifs[[item for item in motifs.columns if "complement" not in item]]
        motifs_rev = motifs[["contig", "motif_type", "mod_type"] + [item for item in motifs.columns if "complement" in item]] \
            .select(pl.all().name.map(lambda col_name: col_name.replace('_complement', ''))) \
            .select(motifs_fwd.columns) \
            .filter(pl.col("motif").is_not_null())
        motifs = pl.concat([motifs_fwd, motifs_rev]).unique(["motif", "contig", "mod_type", "mod_position"])


    output = nm.bin_consensus.within_bin_motifs_consensus(pileup.pileup, assembly, motifs, motifs_scored, bins)
    output = nm.bin_consensus.merge_bin_motifs(output, bins, pileup, assembly)

    output = output.rename({"bin":"contig", "n_mod_bin":"n_mod", "n_nomod_bin":"n_nomod"})

    output = nm.postprocess.join_motif_complements(output)

    output = output.rename({"contig":"bin", "n_mod":"n_mod_bin", "n_nomod":"n_nomod_bin"})
    output = output.sort(["bin", "mod_type", "motif"])

    output = output.filter(pl.col("n_mod_bin") + pl.col("n_nomod_bin") > args.min_motifs_bin)
    output.write_csv(args.out + "/bin-motifs.tsv", separator="\t")

def motif_discovery(args):
    # Check if output directory exists
    log.info("Loading required files")
    pileup = nm.load_pileup(args.pileup, threads = args.threads, min_fraction = args.threshold_methylation_general)
    assembly = nm.load_assembly(args.assembly)

    # Find motifs
    log.info("Finding motifs")
    motifs = find_motifs(args, pileup=pileup, assembly=assembly)
    if motifs is None:
        log.info("Stopping workflow")
        return

    # Score all motifs
    log.info("Scoring motifs")
    scored_all = score_motifs(args, pileup=pileup, assembly=assembly, motifs=motifs)

    # Bin consensus
    log.info("Finding bin consensus motifs")
    bin_consensus(args, pileup=pileup, assembly=assembly, motifs=motifs, motifs_scored=scored_all)

    log.info("Done")

def check_install(args):
    
    # Check if output directory exists
    log.info("Loading required files")
    args.out = "nanomotif_install_check"

    pileup = nm.datasets.geobacillus_plasmids_pileup()
    assembly = nm.datasets.geobacillus_plasmids_assembly()

    # Find motifs
    log.info("Finding motifs")
    motifs = find_motifs(args, pileup=pileup, assembly=assembly)

    # Score all motifs
    log.info("Scoring motifs")
    scored_all = score_motifs(args, pileup=pileup, assembly=assembly, motifs=motifs)

    # Bin consensus
    log.info("Finding bin consensus motifs")
    args.bins = nm.datasets.geobacillus_plasmids_bin_path()
    bin_consensus(args, pileup=pileup, assembly=assembly, motifs=motifs, motifs_scored=scored_all)
    
    log.info("Done")
    shutil.rmtree(args.out)


####################################################################################################
# Binnary - contamination and inclusion
from nanomotif.binnary import data_processing, detect_contamination, include_contigs
from nanomotif.binnary.logging import set_logger_config



def binnary(args):
    """
    binnary entry point for the DNA Methylation Pattern Analysis tool.
    Orchestrates the workflow of the tool based on the provided arguments.
    """
    
    # Conditional check for --write_bins and --assembly_file
    if args.write_bins and not args.assembly_file:
        print("Error: --assembly_file must be specified when --write_bins is used.")
        sys.exit(1)
    elif not args.write_bins and args.assembly_file:
        print("Error: --assembly_file can only be used when --write_bins is specified.")
        sys.exit(1)
    
    print("Starting Binnary ", args.command, " analysis...")

    
    # Settting up the logger
    # set_logger_config(args)
    # logger = log.getLogger(__name__)
    log.info("Starting Binnary analysis...")
    
    # Set number of threads for polars
    os.environ['POLARS_MAX_THREADS'] = str(args.threads)
    
    POLAR_THREADS = pl.threadpool_size()
    log.info(f"Polars is using {POLAR_THREADS} threads.")
    
    # Step 1: Load and preprocess data
    # These functions would be defined in your data_processing module
    (
        motifs_scored,
        bin_motifs,
        contig_bins,
    ) = data_processing.load_data(args)

    

    # Step 2: create motifs_scored_in_bins and bin_motif_binary
    bin_motif_binary = data_processing.prepare_bin_consensus(bin_motifs, args)
    
    motifs_in_bin_consensus = bin_motif_binary.select("motif_mod").unique()["motif_mod"]
    
    motifs_scored_in_bins = data_processing.prepare_motifs_scored_in_bins(
        motifs_scored,
        motifs_in_bin_consensus,
        contig_bins
    )
    
    
    # Create the bin_consensus dataframe for scoring
    log.info("Creating bin_consensus dataframe for scoring...")
    bin_motifs_from_motifs_scored_in_bins = data_processing.construct_bin_consensus_from_motifs_scored_in_bins(
        motifs_scored_in_bins,
        args
    )

    # Setting up the contamination analysis
    if (args.command == "detect_contamination" and not args.contamination_file) or (args.command == "include_contigs" and args.run_detect_contamination):
        contamination = detect_contamination.detect_contamination(
            motifs_scored_in_bins, bin_motifs_from_motifs_scored_in_bins, args
        )
        data_processing.generate_output(contamination.to_pandas(), args.out, "bin_contamination.tsv")
    elif args.contamination_file:
        log.info("Loading contamination file...")
        contamination = data_processing.load_contamination_file(args.contamination_file)
        
    if args.command == "detect_contamination":
        # Create a decontaminated contig_bin file
        new_contig_bins = data_processing.create_contig_bin_file(
            contig_bins=contig_bins.to_pandas(), 
            contamination=contamination.to_pandas(),
            include=None
        )
        data_processing.generate_output(new_contig_bins, args.out, "decontaminated_contig_bin.tsv")

    if args.command == "include_contigs":
        # # User provided contamination file
        # if args.contamination_file:
        #     log.info("Loading contamination file...")
        #     contamination = data_processing.load_contamination_file(args.contamination_file)
        
        ## If run_detect_contamination is false and contamination file is not provided, then set contamination to None
        if not args.run_detect_contamination and not args.contamination_file:
            contamination = pl.DataFrame(
                {
                    "bin": [],
                    "bin_contig_compare": [],
                    "binary_methylation_missmatch_score": [],
                    "non_na_comparisons": [],
                    "contig": []
                }
            )
        
        
        
        # Run the include_contigs analysis    
        include_contigs_df = include_contigs.include_contigs(
            motifs_scored_in_bins, bin_motifs_from_motifs_scored_in_bins, contamination, args
        )
        
        # Save the include_contigs_df results
        data_processing.generate_output(include_contigs_df.to_pandas(), args.out, "include_contigs.tsv")
        
        # Create a new contig_bin file
        
        
        new_contig_bins = data_processing.create_contig_bin_file(
            contig_bins=contig_bins.to_pandas(), 
            contamination= contamination.to_pandas(),
            include=include_contigs_df.to_pandas()
        )
        data_processing.generate_output(new_contig_bins, args.out, "new_contig_bin.tsv")
        
    if args.write_bins:
        log.info("Write bins flag is set. Writing bins to file...")
        log.info("Loading assembly file...")
        assembly = data_processing.read_fasta(args.assembly_file)
        
        bin_dir = os.path.join(args.out, args.command + "_bins")
        
        data_processing.write_bins_from_contigs(new_contig_bins, assembly, bin_dir)
    
    log.info(f"Analysis Completed. Results are saved to: {args.out}")
    print("Analysis Completed. Results are saved to:", args.out)


from nanomotif.mtase_linker.dependencies import snakemake_create_environments, get_models, defensefinder_update
from nanomotif.mtase_linker.command import run_MTase_linker

def mtase_linker(args):
    if args.mtase_linker_command == "install":
        snakemake_create_environments(args)
        get_models(args)
        defensefinder_update(args)
    elif args.mtase_linker_command == "run":
        run_MTase_linker(args)
    else:
        print(f"Unknown MTase-linker command: {args.mtase_linker_command}")
       



def main():
    # Parse arguments
    parser = nm.argparser.create_parser()
    args = parser.parse_args()
    
    if args.command in ["detect_contamination", "include_contigs", "MTase-linker"]:
        args.verbose = False
        args.seed = 1
    
    if args.command == "find_motifs":
        shared_setup(args, args.out)
        find_motifs(args)
    elif args.command == "score_motifs":
        shared_setup(args, args.out)
        score_motifs(args)
    elif args.command == "bin_consensus":
        shared_setup(args, args.out)
        bin_consensus(args)
    elif args.command == "motif_discovery":
        shared_setup(args, args.out)
        motif_discovery(args)

    elif args.command in ["detect_contamination", "include_contigs"]:
        shared_setup(args, args.out)
        binnary(args)

    elif args.command == "MTase-linker":
        mtase_linker(args)

    elif args.command == "check_installation":
        args.out = "nanomotif_install_check"
        shared_setup(args, args.out)
        check_install(args)

    else:
        parser.print_help()
        exit()
if __name__ == "__main__":
    main()
