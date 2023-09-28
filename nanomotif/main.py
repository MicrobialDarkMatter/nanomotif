import nanomotif as nm
import logging as log
import polars as pl
def main():
    parser = nm.argparser.create_parser()
    args = parser.parse_args()

    log.info("Loading assembly")
    assembly = nm.load_assembly(args.assembly)
    log.info("Loading pileup")
    pileup = nm.load_pileup(args.pileup)

    contigs_with_mods = pileup.pileup.groupby("contig") \
        .agg(pl.count()) \
        .filter(pl.col("count") > 50) \
        .get_column("contig").to_list()
    contigs_to_process = [contig for contig in assembly.assembly.keys() if contig in contigs_with_mods]
    pileup.pileup = pileup.pileup.filter(pl.col("contig").is_in(contigs_to_process))

    log.info("Identifying motifs")
    motifs = nm.evaluate.process_sample(
        assembly,
        pileup.pileup,
        args.max_motif_length,
        args.min_fraction,
        args.min_coverage,
        args.min_kl_divergence,
        args.min_cdf_score,
        args.cdf_position
    )

    log.info("Writing motifs")
    motifs.with_columns(
        pl.col("model").apply(lambda x: x._alpha).alias("alpha"),
        pl.col("model").apply(lambda x: x._beta).alias("beta")
    ).drop("model").write_csv(args.output, separator="\t")

if __name__ == "__main__":
    main()