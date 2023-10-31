import nanomotif as nm
import logging as log
import polars as pl
import os
import json

def main():
    # Parse arguments
    parser = nm.argparser.create_parser()
    args = parser.parse_args()

    # Check if output directory exists
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    else:
        log.error(f"Output directory {args.output} already exists")
        return
    
    # Save settings
    with open(args.output + "/args.json", "w") as f:
        json.dump(vars(args), f, indent=2)

    # Set up logging
    log.basicConfig(
        filename=args.output + "/nanomotif.log",
        level=log.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%d-%b-%y %H:%M:%S",
    )
    log.getLogger().addHandler(log.StreamHandler())
    if args.verbose:
        log.getLogger().setLevel(log.DEBUG)


    # Log arguments
    log.info("Arguments:")
    for arg in vars(args):
        log.debug(f"{arg}: {getattr(args, arg)}")

    # Run nanomotif
    log.info("Starting nanomotif")
    log.info("Loading assembly")
    assembly = nm.load_assembly(args.assembly)
    assm_lengths = pl.DataFrame({
        "contig": list(assembly.assembly.keys()),
        "length": [len(contig) for contig in assembly.assembly.values()]
    })
    log.info("Loading pileup")
    pileup = nm.load_pileup(args.pileup)

    # Filter pileup to contigs with mods
    contigs_with_mods = pileup.pileup \
        .filter(pl.col("fraction_mod") > args.min_fraction) \
        .groupby(["contig", "mod_type"]) \
        .agg(pl.count()) \
        .join(assm_lengths, on = "contig") \
        .filter(pl.col("count") > pl.col("length")/10000) \
        .get_column("contig").unique().to_list()

    total_contigs = len(assembly.assembly.keys())
    log.info(f"Processing {len(contigs_with_mods)} of {total_contigs} contigs")

    contigs_to_process = [contig for contig in assembly.assembly.keys() if contig in contigs_with_mods]
    pileup = pileup.pileup.filter(pl.col("contig").is_in(contigs_to_process))

    log.info("Identifying motifs")
    motifs = nm.evaluate.process_sample(
        assembly,
        pileup,
        max_candidate_size = args.max_motif_length,
        min_read_methylation_fraction = args.min_fraction,
        min_valid_coverage = args.min_coverage,
        min_kl_divergence = args.min_kl_divergence
    )
    if motifs is None:
        log.info("No motifs found")
        return

    log.info("Writing motifs")
<<<<<<< HEAD
    def save_motif_df(df, name):
        df.with_columns(
            pl.col("model").apply(lambda x: x._alpha).alias("alpha"),
            pl.col("model").apply(lambda x: x._beta).alias("beta")
        ).drop("model").write_csv(args.output + "/" + name + ".tsv", separator="\t")
    save_motif_df(motifs, "motifs-raw")

    log.info("Postprocessing motifs")

    log.info(" - Writing motifs")
    motifs = motifs.filter(pl.col("score") > 0.15)
    save_motif_df(motifs, "motifs-score-filtered")

    log.info(" - Removing sub motifs")
    motifs = nm.postprocess.remove_sub_motifs(motifs)
    save_motif_df(motifs, "motifs-score-sub-filtered")

    log.info(" - Removing noisy motifs")
    motifs = nm.postprocess.remove_noisy_motifs(motifs)
    save_motif_df(motifs, "motifs-score-sub-noise-filtered")

    log.info(" - Merging motifs")
    motifs = nm.postprocess.merge_motifs_in_df(motifs, pileup, assembly)
    save_motif_df(motifs, "motifs")
=======
    motifs.with_columns(
        pl.col("model").apply(lambda x: x._alpha).alias("alpha"),
        pl.col("model").apply(lambda x: x._beta).alias("beta")
    ).drop("model").write_csv(args.output + "/motifs-raw.tsv", separator="\t")

    log.info("Postprocessing motifs")
    motifs = motifs.filter(pl.col("score") > 0.15)
    motifs = nm.postprocess.remove_sub_motifs(motifs)
    motifs = nm.postprocess.remove_noisy_motifs(motifs)
    motifs.with_columns(
        pl.col("model").apply(lambda x: x._alpha).alias("alpha"),
        pl.col("model").apply(lambda x: x._beta).alias("beta")
    ).drop("model").write_csv(args.output + "/motifs.tsv", separator="\t")
>>>>>>> 756ec0aaa02f67bccb18ad0f571059b2960f54b5

if __name__ == "__main__":
    main()