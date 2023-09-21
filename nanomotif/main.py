import nanomotif as nm
import polars as pl
def main():
    parser = nm.argparser.create_parser()
    args = parser.parse_args()

    assembly = nm.load_assembly(args.assembly)
    pileup = nm.load_pileup(args.pileup)

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
    motifs.with_columns(
        pl.col("model").apply(lambda x: x._alpha).alias("alpha"),
        pl.col("model").apply(lambda x: x._beta).alias("beta")
    ).drop("model").write_csv(args.output, separator="\t")

if __name__ == "__main__":
    main()