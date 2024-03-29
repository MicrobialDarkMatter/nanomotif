import polars as pl
import nanomotif as nm


def within_bin_motifs_consensus(pileup, assembly, motifs, motifs_scored, bins):
    motifs = motifs.with_columns([
        pl.lit(True).alias("directly_detected")
    ])
    motifs_scored = motifs_scored \
    .join(bins, on="contig", how = "outer") \
    .with_columns(
        pl.when(pl.col("bin").is_null()).then(False).otherwise(pl.lit(True)).alias("binned")
    ).join(
        motifs.drop("sequence",	"score", "alpha", "beta"), on=["contig", "motif", "mod_type", "mod_position"], how="left"
    ).with_columns(
        pl.when(pl.col("directly_detected").is_null()).then(pl.lit(False)).otherwise(pl.col("directly_detected")).alias("directly_detected")
    )

    # keep directly detected motifs
    bin_consensus = motifs_scored.groupby("bin", "motif", "mod_position", "mod_type") \
        .apply(lambda group: group.filter(
            pl.col("directly_detected").any()
        ))
    # Merge motifs
    merged_motifs = merge_motifs_in_df(bin_consensus.select(["contig","mod_type", "alpha", "beta","motif","mod_position",]), pileup, assembly)
    # Remove submotifs
    merged_motifs = nm.postprocess.remove_sub_motifs(merged_motifs)

    # Join merged motifs with unmerged motifs
    merged_motifs = merged_motifs \
        .with_columns(
            pl.col("motif").apply(lambda x: nm.seq.regex_to_iupac(x)).alias("motif")
        ) \
        .join(bin_consensus, on=["contig", "motif", "mod_position", "mod_type"], how="left")


    bin_consensus = merged_motifs.unique(["bin", "motif", "mod_position", "mod_type", "contig"]) \
        .with_columns(
            (pl.when(pl.col("mean") > 0.5).then(True).otherwise(False)).alias("count_mean_threshold")
        )
    bin_consensus = bin_consensus \
        .groupby("bin", "motif", "mod_position", "mod_type") \
        .agg(
            pl.col("alpha").sum().alias("alpha_sum"),
            pl.col("beta").sum().alias("beta_sum"),
            pl.col("directly_detected").sum().alias("count_directly_detected"),
            pl.col("count_mean_threshold").sum().alias("count_mean_threshold"),
            pl.col("contig").count().alias("contig_count")
        ).with_columns( 
            (pl.col("alpha_sum") / (pl.col("alpha_sum") + pl.col("beta_sum"))).alias("mean_sum"),
            (pl.col("count_mean_threshold") / pl.col("contig_count").cast(pl.Float32)).alias("fraction_mean_threshold"),
        )
    bin_consensus = bin_consensus \
        .with_columns(
            pl.col("motif").apply(lambda x: nm.utils.motif_type(x)).alias("motif_type")
        )
        

    output = bin_consensus.select(["bin", "motif", "mod_position", "mod_type", "contig_count", "count_mean_threshold", "motif_type", "alpha_sum", "beta_sum", "mean_sum"]) \
        .rename({"alpha_sum":"methylated_count", "beta_sum":"non_methylated_count", "count_mean_threshold":"contigs_with_motif", "mean_sum":"mean_methylation_bin"}) \
        .filter(pl.col("contigs_with_motif") > 0) \
        .filter(pl.col("motif_type") != "ambiguous") \
        .filter(pl.col("methylated_count") > 10) \
        .with_columns(
            pl.when(pl.col("bin").is_null()).then("unbinned").otherwise(pl.col("bin")).alias("bin")
        )
    return output



def merge_motifs_in_df(motif_df, pileup, assembly):
    new_df = []
    for (contig, mod_type), df in motif_df.groupby("contig", "mod_type"):
        # Get list of motifs
        motif_seq = df["motif"].to_list()
        motif_pos = df["mod_position"].to_list()
        motifs = [nm.candidate.Motif(seq, pos) for seq, pos in zip(motif_seq, motif_pos)]

        # Merge motifs
        merged_motifs, pre_merge_motifs = nm.candidate.merge_motifs(motifs)
        
        # Create a new dataframe with the non merged motifs
        if  len(pre_merge_motifs) == 0:
            # No motifs were merged
            new_df.append(df)
            continue
        single_df = df.filter(pl.col("motif").is_in(pre_merge_motifs).not_())
        new_df.append(single_df)
        merged_df = []
        for motif in merged_motifs:
            merged_model = nm.evaluate.motif_model_contig(
                pileup.filter((pl.col("contig") == contig) & (pl.col("mod_type") == mod_type)), 
                assembly.assembly[contig].sequence,
                motif
            )
            merged_df.append(pl.DataFrame({
                "contig": contig,
                "mod_type": mod_type,
                "alpha": merged_model._alpha,
                "beta": merged_model._beta,
                "motif": motif.string,
                "mod_position": motif.mod_position
            }))
        new_df.append(pl.concat(merged_df))
    new_df = pl.concat(new_df)

    return new_df

