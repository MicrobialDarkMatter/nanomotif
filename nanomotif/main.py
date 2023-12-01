import nanomotif as nm
from nanomotif.logger import configure_logger
import logging as log
import os
from pathlib import Path
import json
os.environ["POLARS_MAX_THREADS"] = "1"
import polars as pl
import warnings

def shared_setup(args, working_dir):
    warnings.filterwarnings("ignore")

    # Set up logging
    LOG_DIR = working_dir + "/logs"
    Path(LOG_DIR).mkdir(parents=True, exist_ok=True)
    configure_logger(LOG_DIR + f"/{args.command}.main.log", args.verbose)

    # Log arguments
    log.debug("Arguments:")
    for arg in vars(args):
        log.debug(f"{arg}: {getattr(args, arg)}")
    
    # Save settings
    with open(working_dir + f"/args.{args.command}.json", "w") as f:
        json.dump(vars(args), f, indent=2)
    
    # Check if output directory exists
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    else:
        log.error(f"Output directory {args.output} already exists")
        return





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
            log_dir = args.output + "/logs"
        )

    if motifs is None:
        log.info("No motifs found")
        return

    log.info("Writing motifs")
    def save_motif_df(df, name):
        df.with_columns(
            pl.col("model").apply(lambda x: x._alpha).alias("alpha"),
            pl.col("model").apply(lambda x: x._beta).alias("beta")
        ).drop("model").with_columns([
            pl.col("motif").apply(lambda x: nm.seq.regex_to_iupac(x)).alias("motif_iupac")
        ]) \
        .with_columns([
            pl.col("motif").apply(lambda x: nm.utils.motif_type(x)).alias("motif_type")
        ]) \
        .unique(["motif_iupac", "contig", "mod_type", "mod_position"]) \
        .drop("sequence", "score") \
        .write_csv(args.output + "/" + name + ".tsv", separator="\t")
    os.makedirs(args.output + "/precleanup-motifs/", exist_ok=True)
    save_motif_df(motifs, "precleanup-motifs/motifs-raw")

    log.info("Postprocessing motifs")

    log.info(" - Writing motifs")
    motifs = motifs.filter(pl.col("score") > 0.1)
    if len(motifs) == 0:
        log.info("No motifs found")
        return
    save_motif_df(motifs, "precleanup-motifs/motifs-score-filtered")

    log.info(" - Removing sub motifs")
    motifs = nm.postprocess.remove_sub_motifs(motifs)
    if len(motifs) == 0:
        log.info("No motifs found")
        return
    save_motif_df(motifs, "precleanup-motifs/motifs-score-sub-filtered")

    log.info(" - Removing noisy motifs")
    motifs = nm.postprocess.remove_noisy_motifs(motifs)
    if len(motifs) == 0:
        log.info("No motifs found")
        return
    save_motif_df(motifs, "precleanup-motifs/motifs-score-sub-noise-filtered")

    log.info(" - Merging motifs")
    motifs = nm.postprocess.merge_motifs_in_df(motifs, pileup, assembly)
    if len(motifs) == 0:
        log.info("No motifs found")
        return
    save_motif_df(motifs, "motifs")

    log.info("Done finding motifs")
    return motifs.drop("model")

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
    log.info("Scoring motifs")
    scored_all = nm.scoremotifs.score_sample_parallel(
        assembly, pileup.pileup, motifs,
        threads = args.threads,
        threshold_methylation_general = args.threshold_methylation_general,
        threshold_valid_coverage = 1,
        verbose = args.verbose,
        log_dir = args.output + "/logs"
        )
    scored_all = scored_all \
        .with_columns(
            pl.col("motif").apply(lambda x: nm.seq.regex_to_iupac(x)).alias("motif_iupac")
        ) \
        .with_columns([
            pl.col("motif").apply(lambda x: nm.utils.motif_type(x)).alias("motif_type")
        ]).unique(["motif_iupac", "contig", "mod_type", "mod_position"])
    scored_all.drop("model").write_csv(args.output + "/motifs-scored.tsv", separator="\t")
    return scored_all.drop("model")

def bin_consensus(args, pileup = None, assembly = None, motifs = None, motifs_scored = None):
    bins = pl.read_csv(args.bins, separator="\t", has_header=False) \
        .rename({"column_1":"contig", "column_2":"bin"})
    log.info("Loading motifs-scored")
    if motifs is None:
        motifs = pl.read_csv(args.motifs, separator="\t")
    log.info("Loading motifs")
    if motifs_scored is None:
        motifs_scored = pl.read_csv(args.motifs_scored, separator="\t")
    log.info("Loading pileup")
    if pileup is None:
        pileup =  nm.load_pileup(args.pileup, threads = args.threads, min_fraction = args.threshold_methylation_general)
    log.info("Loading assembly")
    if assembly is None:
        assembly = nm.load_assembly(args.assembly)
    motifs_scored = motifs_scored.with_columns([
        pl.col("motif").apply(lambda x: nm.seq.regex_to_iupac(x)).alias("motif_iupac")
    ])
    motifs = motifs.with_columns([
        pl.col("motif").apply(lambda x: nm.seq.regex_to_iupac(x)).alias("motif_iupac")
    ])

    output = nm.bin_consensus.within_bin_motifs_consensus(pileup.pileup, assembly, motifs, motifs_scored, bins)
    output.write_csv(args.output + "/bin-motifs.tsv", separator="\t")

def metagenomic_workflow(args):
    # Check if output directory exists
    pileup = nm.load_pileup(args.pileup, threads = args.threads, min_fraction = args.threshold_methylation_general)
    assembly = nm.load_assembly(args.assembly)

    # Find motifs
    motifs = find_motifs(args, pileup=pileup, assembly=assembly)

    # Score all motifs
    scored_all = score_motifs(args, pileup=pileup, assembly=assembly, motifs=motifs)

    # Bin consensus
    bin_consensus(args, pileup=pileup, assembly=assembly, motifs=motifs, motifs_scored=scored_all)

    log.info("Completely done")


def main():
    # Parse arguments
    parser = nm.argparser.create_parser()
    args = parser.parse_args()
    shared_setup(args, args.output)
    if args.command == "find-motifs":
        find_motifs(args)
    elif args.command == "score-motifs":
        score_motifs(args)
    elif args.command == "bin-consensus":
        bin_consensus(args)
    elif args.command == "complete-workflow":
        metagenomic_workflow(args)
    #elif args.command == "clean-bins":
    #    pass
    #elif args.command == "associate-mges":
    #    pass
    else:
        parser.print_help()
        exit()
if __name__ == "__main__":
    main()
