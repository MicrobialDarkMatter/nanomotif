import polars as pl
from polars import col
import nanomotif as nm
import logging as log
import numpy as np

def within_bin_motifs_consensus(pileup, assembly, motifs, motifs_scored, bins):
    log.debug("Starting within bin motif consensus")
    motifs = motifs.with_columns([
        pl.lit(True).alias("directly_detected")
    ]).with_columns([
        pl.col("motif") \
            .map_elements(lambda x: nm.seq.regex_to_iupac(x)) \
            .map_elements(lambda x: nm.seq.iupac_to_regex(x)) \
            .alias("motif")
    ])

    log.debug("Joining motifs with bins")
    motifs_scored = motifs_scored.with_columns([
        pl.col("motif") \
            .map_elements(lambda x: nm.seq.regex_to_iupac(x)) \
            .map_elements(lambda x: nm.seq.iupac_to_regex(x)).alias("motif")
    ]) \
    .join(bins, on="contig", how = "left") \
    .filter(pl.col("bin").is_not_null()) \
    .join(
        motifs.select("contig", "motif", "mod_type", "mod_position", "directly_detected"), 
        on=["contig", "motif", "mod_type", "mod_position"], 
        how="left"
    ).with_columns(
        pl.when(pl.col("directly_detected").is_null()).then(pl.lit(False)).otherwise(pl.col("directly_detected")).alias("directly_detected")
    )

    log.debug("Calculating motif scores")
    # keep directly detected motifs
    bin_consensus = motifs_scored.groupby("bin", "motif", "mod_position", "mod_type") \
        .apply(lambda group: group.filter(pl.col("directly_detected").any()))
    log.debug("Mergin similar motifs within bin")
    # Merge motifs
    bin_consensus = bin_consensus.with_columns(
        pl.lit(-1).alias("score")
    )
    merged_motifs = nm.postprocess.merge_motifs_in_df(bin_consensus.select(["contig","mod_type", "n_mod", "n_nomod","motif","mod_position"]), pileup, assembly)
    log.debug("Removing submotifs")
    bin_consensus = bin_consensus.drop("score")
    # Remove submotifs
    merged_motifs = nm.postprocess.remove_sub_motifs(merged_motifs)
    
    log.debug("Calculating motif scores")
    # Join merged motifs with unmerged motifs
    bin_consensus = merged_motifs \
        .join(
            bin_consensus, on=["contig", "motif", "mod_position", "mod_type"], 
            how="left"
        )
    bin_consensus = bin_consensus.with_columns(
            pl.col("motif").apply(lambda x: nm.seq.regex_to_iupac(x)).alias("motif"),
            (pl.col("n_mod")  / (pl.col("n_mod") + pl.col("n_nomod"))).alias("mean")
        )
    bin_consensus=bin_consensus.unique([
            "bin", "motif", "mod_position", "mod_type", "contig"
        ])
    bin_consensus = bin_consensus.with_columns(
            (pl.when(pl.col("mean") > 0.5).then(True).otherwise(False)).alias("n_contigs_above_mean")
        )
    bin_consensus = bin_consensus.groupby("bin", "motif", "mod_position", "mod_type") \
        .agg(
            pl.col("n_mod").sum().alias("n_mod_bin"),
            pl.col("n_nomod").sum().alias("n_nomod_bin"),
            pl.col("n_contigs_above_mean").sum().alias("n_contigs_above_mean"),
            pl.col("contig").count().alias("contig_count")
        ).with_columns( 
            pl.col("motif").apply(lambda x: nm.utils.motif_type(x)).alias("motif_type"),
            (pl.col("n_mod_bin") / (pl.col("n_mod_bin") + pl.col("n_nomod_bin"))).alias("mean_sum")
        )
        
    log.debug("Filtering motifs")
    output = bin_consensus.select(["bin", "motif", "mod_position", "mod_type", "contig_count", "n_contigs_above_mean", "motif_type", "n_mod_bin", "n_nomod_bin", "mean_sum"]) \
        .rename({"n_contigs_above_mean":"contigs_with_motif", "mean_sum":"mean_methylation_bin"}) \
        .filter(pl.col("contigs_with_motif") > 0) \
        .filter(pl.col("n_mod_bin") > 50)
    return output



def merge_motifs_in_df(motif_df, pileup, assembly, mean_shift_threshold = -0.25):
    new_df = []
    for (contig, mod_type), df in motif_df.groupby("contig", "mod_type"):
        # Get list of motifs
        motif_seq = df["motif"].to_list()
        motif_pos = df["mod_position"].to_list()
        motifs = [nm.candidate.Motif(seq, pos) for seq, pos in zip(motif_seq, motif_pos)]

        # Merge motifs
        merged_motifs = nm.candidate.merge_motifs(motifs)
        all_merged_motifs = []
        all_premerge_motifs = []
        # Check mean shift of premerge motifs to merged motif is high enough
        for cluster, motifs in merged_motifs.items():
            merged_motif = motifs[0]
            premerge_motifs = motifs[1]
            merge_mean = nm.evaluate.motif_model_contig(
                pileup.filter((col("contig") == contig) & (col("mod_type") == mod_type)), 
                assembly.assembly[contig].sequence,
                merged_motif
            ).mean()
            pre_merge_means = []
            for pre_merge_motif in premerge_motifs:
                pre_merge_means.append(nm.evaluate.motif_model_contig(
                    pileup.filter((col("contig") == contig) & (col("mod_type") == mod_type)), 
                    assembly.assembly[contig].sequence,
                    pre_merge_motif
                ).mean())
            
            pre_merge_mean = sum(np.array(pre_merge_means)) / len(pre_merge_means)
            mean_shift = merge_mean - pre_merge_mean
            if mean_shift < mean_shift_threshold:
                log.info(f"Mean shift of merged motif {merged_motif} is {mean_shift}, keeping original motifs")
            
            else:
                log.info(f"Mean shift of merged motif {merged_motif} is {mean_shift}")
                all_merged_motifs.append(merged_motif)
                all_premerge_motifs.extend(premerge_motifs)

        # Create a new dataframe with the non merged motifs
        if  len(all_premerge_motifs) == 0:
            # No motifs were merged
            new_df.append(df)
            continue
        single_df = df.filter(col("motif").is_in(all_premerge_motifs).not_())
        new_df.append(single_df)
        merged_df = []
        for motif in all_merged_motifs:
            merged_model = nm.evaluate.motif_model_contig(
                pileup.filter((col("contig") == contig) & (col("mod_type") == mod_type)), 
                assembly.assembly[contig].sequence,
                motif
            )
            merged_df.append(pl.DataFrame({
                "contig": contig,
                "mod_type": mod_type,
                "n_mod": merged_model._alpha,
                "n_nomod": merged_model._beta - 1,
                "motif": nm.candidate.regex_to_iupac(motif.string),
                "mod_position": motif.mod_position
            }))
        new_df.append(pl.concat(merged_df))
    new_df = pl.concat(new_df)
    return new_df
