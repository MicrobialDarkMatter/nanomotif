import polars as pl
from polars import col
import nanomotif as nm
import logging as log
import numpy as np

def within_bin_motifs_consensus(pileup, assembly, motifs, motifs_scored, bins, minimum_contig_motif_methylation=0.25, minimum_bin_methylation=0.75):
    log.debug("Starting within bin motif consensus")

    motifs = motifs.with_columns(pl.lit(True).alias("directly_detected"))
    def convert_motifs_to_regex(motifs):
        return motifs.with_columns(pl.col("motif").map_elements(lambda x: nm.seq.regex_to_iupac(x)).map_elements(lambda x: nm.seq.iupac_to_regex(x)))
    motifs = convert_motifs_to_regex(motifs)
    motifs_scored = convert_motifs_to_regex(motifs_scored)

    motifs_scored = motifs_scored \
        .join(bins, on="contig", how = "left") \
        .filter(pl.col("bin").is_not_null())
    motifs_scored = motifs_scored.join(
        motifs.select("contig", "motif", "mod_type", "mod_position", "directly_detected"), on=["contig", "motif", "mod_type", "mod_position"], how="left"
    ).with_columns(
        pl.col("directly_detected").fill_null(False)
    ).groupby("bin", "motif", "mod_position", "mod_type") \
            .apply(lambda group: group.filter(pl.col("directly_detected").any()))
    
    methylated_contig_mass = motifs_scored.filter(
            ((pl.col("n_mod") / (pl.col("n_mod") + pl.col("n_nomod"))) > minimum_contig_motif_methylation)
        ).groupby(["bin", "motif", "mod_type", "mod_position"]).agg([
            (pl.col("n_mod") + pl.col("n_nomod")).sum().alias("n_motif_meth_contigs")
        ])
    all_contig_mass = motifs_scored.groupby(["bin", "motif", "mod_type", "mod_position"]).agg([
            (pl.col("n_mod") + pl.col("n_nomod")).sum().alias("n_motif_contigs"),
        ])
    result = methylated_contig_mass.join(all_contig_mass, on=["bin", "motif", "mod_type", "mod_position"]) \
        .with_columns(
            (pl.col("n_motif_meth_contigs") / pl.col("n_motif_contigs")).alias("methylated_proportion")
        )
    motifs_scored_filt = motifs_scored.join(result, on=["bin", "motif", "mod_type", "mod_position"]).filter(pl.col("methylated_proportion") > minimum_bin_methylation)

    
    log.debug("Removing submotifs")
    # Remove submotifs
    contig_motifs = nm.postprocess.remove_sub_motifs(motifs_scored_filt)

    contig_motifs = contig_motifs.with_columns(
            pl.col("motif").apply(lambda x: nm.seq.regex_to_iupac(x)).alias("motif"),
            (pl.col("n_mod")  / (pl.col("n_mod") + pl.col("n_nomod"))).alias("mean")
        )
    bin_motifs = contig_motifs.groupby("bin", "motif", "mod_position", "mod_type") \
        .agg(
            pl.col("n_mod").sum().alias("n_mod_bin"),
            pl.col("n_nomod").sum().alias("n_nomod_bin"),
            pl.col("contig").count().alias("contig_count")
        ).with_columns( 
            pl.col("motif").apply(lambda x: nm.utils.motif_type(x)).alias("motif_type"),
            (pl.col("n_mod_bin") / (pl.col("n_mod_bin") + pl.col("n_nomod_bin"))).alias("mean_sum")
        )
    return bin_motifs



def merge_motifs_in_df(motif_df, pileup, assembly, mean_shift_threshold = -0.25):
    # DEPRECATED CURRENTLY
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
