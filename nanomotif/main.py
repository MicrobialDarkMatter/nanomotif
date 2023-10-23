import nanomotif as nm
import logging as log
import polars as pl
def main():
    parser = nm.argparser.create_parser()
    args = parser.parse_args()

    log.info("Loading assembly")
    assembly = nm.load_assembly(args.assembly)
    assm_lengths = pl.DataFrame({
        "contig": list(assembly.assembly.keys()),
        "length": [len(contig) for contig in assembly.assembly.values()]
    })
    log.info("Loading pileup")
    pileup = nm.load_pileup(args.pileup)

    contigs_with_mods = pileup.pileup \
        .filter(pl.col("fraction_mod") > args.min_fraction) \
        .groupby(["contig", "mod_type"]) \
        .agg(pl.count()) \
        .join(assm_lengths, on = "contig") \
        .filter(pl.col("count") > pl.col("length")/10000) \
        .get_column("contig").unique().to_list()
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
    motifs.with_columns(
        pl.col("model").apply(lambda x: x._alpha).alias("alpha"),
        pl.col("model").apply(lambda x: x._beta).alias("beta")
    ).drop("model").write_csv(args.output, separator="\t")

if __name__ == "__main__":
    main()