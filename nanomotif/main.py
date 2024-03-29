import nanomotif as nm
from nanomotif.logger import configure_logger
import logging as log
import os
import shutil
from pathlib import Path
import json
os.environ["POLARS_MAX_THREADS"] = "1"
import polars as pl
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

    if motifs is None:
        log.info("No motifs found")
        return

    log.info("Writing motifs")

    def format_motif_df(df):
        if "model" in df.columns:
            df = df.with_columns(
                pl.col("model").apply(lambda x: x._alpha).alias("n_mod"),
                pl.col("model").apply(lambda x: x._beta).alias("n_nomod")
            )
        df_out = df.with_columns([
            pl.col("motif").apply(lambda x: nm.seq.regex_to_iupac(x)).alias("motif")
        ]).with_columns([
            pl.col("motif").apply(lambda x: nm.utils.motif_type(x)).alias("motif_type")
        ]).unique(["motif", "contig", "mod_type", "mod_position"])
        try:
            df_out = df_out.select([
                "contig", "motif", "mod_position", "mod_type", "n_mod", "n_nomod", "motif_type",
                "motif_complement", "mod_position_complement", "n_mod_complement", "n_nomod_complement"
            ])
        except:
            df_out = df_out.select([
                "contig", "motif", "mod_position", "mod_type", "n_mod", "n_nomod", "motif_type",
            ])
        return df_out
    def save_motif_df(df, name):
        df = format_motif_df(df)
        df = df.sort(["contig", "mod_type", "motif"])
        df.write_csv(args.out + "/" + name + ".tsv", separator="\t")
    os.makedirs(args.out + "/precleanup-motifs/", exist_ok=True)
    save_motif_df(motifs, "precleanup-motifs/motifs-raw")

    log.info("Postprocessing motifs")

    log.info(" - Writing motifs")
    motifs = motifs.filter(pl.col("score") > 0.1)
    if len(motifs) == 0:
        log.info("No motifs found")
        return
    save_motif_df(motifs, "precleanup-motifs/motifs-score")

    log.info(" - Removing sub motifs")
    motifs = nm.postprocess.remove_sub_motifs(motifs)
    if len(motifs) == 0:
        log.info("No motifs found")
        return
    save_motif_df(motifs, "precleanup-motifs/motifs-score-sub")

    log.info(" - Removing noisy motifs")
    motifs = nm.postprocess.remove_noisy_motifs(motifs)
    if len(motifs) == 0:
        log.info("No motifs found")
        return
    save_motif_df(motifs, "precleanup-motifs/motifs-score-sub-noise")

    log.info(" - Merging motifs")
    motifs = nm.postprocess.merge_motifs_in_df(motifs, pileup, assembly)
    if len(motifs) == 0:
        log.info("No motifs found")
        return
    save_motif_df(motifs, "precleanup-motifs/motifs-score-sub-noise-merge")

    log.info(" - Joining motif complements")
    motifs = nm.postprocess.join_motif_complements(motifs)
    save_motif_df(motifs, "precleanup-motifs/motifs-score-sub-noise-merge-complement")

    log.info(" - Removing motifs observed less than min count")
    motifs = motifs.filter(pl.col("n_mod") + pl.col("n_nomod") > 50)
    if len(motifs) == 0:
        log.info("No motifs found")
        return
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
    output = output.rename({"bin":"contig", "n_mod_bin":"n_mod", "n_nomod_bin":"n_nomod"})

    output = nm.postprocess.join_motif_complements(output)

    output = output.rename({"contig":"bin", "n_mod":"n_mod_bin", "n_nomod":"n_nomod_bin"})
    output = output.sort(["bin", "mod_type", "motif"])
    output.write_csv(args.out + "/bin-motifs.tsv", separator="\t")

def metagenomic_workflow(args):
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

def main():
    # Parse arguments
    parser = nm.argparser.create_parser()
    args = parser.parse_args()
    shared_setup(args, args.out)
    if args.command == "find-motifs":
        find_motifs(args)
    elif args.command == "score-motifs":
        score_motifs(args)
    elif args.command == "bin-consensus":
        bin_consensus(args)
    elif args.command == "complete-workflow":
        metagenomic_workflow(args)
    elif args.command == "check-installation":
        check_install(args)
    #    pass
    #elif args.command == "associate-mges":
    #    pass
    else:
        parser.print_help()
        exit()
if __name__ == "__main__":
    main()
