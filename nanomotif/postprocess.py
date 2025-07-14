import nanomotif as nm
import polars as pl
import numpy as np
from polars import col
import logging as log

def remove_noisy_motifs(motif_df):
    """
    Remove motifs that have isolated bases
    """
    assert "motif" in motif_df.columns
    assert "mod_position" in motif_df.columns
    assert len(motif_df) > 0
    motif_strings = motif_df.get_column("motif").to_list()
    positions = motif_df.get_column("mod_position").to_list()
    motifs = [nm.motif.Motif(motif_string, pos) for motif_string, pos in zip(motif_strings, positions)]
    clean_motifs = []
    for motif in motifs:
        if not motif.have_isolated_bases(isolation_size = 3):
            clean_motifs.append(motif.string)
    if not clean_motifs:
        return motif_df.clear()
    else: 
        motif_df_clean = motif_df.filter(pl.col("motif").is_in(clean_motifs))
        return motif_df_clean

def remove_child_motifs(motifs):
    parent_motifs = []
    for i, motif in enumerate(motifs):
        parent = True
        for j, other in enumerate(motifs):
            if i == j:
                continue
            if motif.sub_string_of(other):
                parent = False
                break
        if parent:
            parent_motifs.append(motif)
    return parent_motifs

def remove_sub_motifs(motif_df):
    motif_df_clean = []
    for contig, df in motif_df.group_by("contig"):
        motif_strings = df.get_column("motif").to_list()
        positions = df.get_column("mod_position").to_list()
        motifs = [nm.motif.Motif(motif_string, pos) for motif_string, pos in zip(motif_strings, positions)]
        parent_motifs = remove_child_motifs(motifs)
        df_clean = df.filter(col("motif").is_in(parent_motifs))
        motif_df_clean.append(df_clean)
    motif_df_clean = pl.concat(motif_df_clean)
    return motif_df_clean

def merge_motifs_in_df(motif_df, pileup, assembly, mean_shift_threshold = -0.2):
    new_df = []
    for (contig, mod_type), df in motif_df.group_by("contig", "mod_type"):
        # Get list of motifs
        motif_seq = df["motif"].to_list()
        motif_pos = df["mod_position"].to_list()
        motifs = [nm.motif.Motif(seq, pos) for seq, pos in zip(motif_seq, motif_pos)]

        # Merge motifs
        merged_motifs = nm.motif.merge_motifs(motifs)
        all_merged_motifs = []
        all_premerge_motifs = []
        # Check mean shift of premerge motifs to merged motif is high enough
        for cluster, motifs in merged_motifs.items():
            merged_motif = motifs[0]
            premerge_motifs = motifs[1]
            merge_mean = nm.find_motifs.motif_model_contig(
                pileup.filter((col("contig") == contig) & (col("mod_type") == mod_type)), 
                assembly.assembly[contig].sequence,
                merged_motif
            ).mean()
            pre_merge_means = []
            for pre_merge_motif in premerge_motifs:
                pre_merge_means.append(nm.find_motifs.motif_model_contig(
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
            merged_model = nm.find_motifs.motif_model_contig(
                pileup.filter((col("contig") == contig) & (col("mod_type") == mod_type)), 
                assembly.assembly[contig].sequence,
                motif
            )
            merged_df.append(pl.DataFrame({
                "sequence": motif.string,
                "score": -1.,
                "contig": contig,
                "mod_type": mod_type,
                "model": merged_model,
                "motif": motif.string,
                "mod_position": motif.mod_position
            }))
        new_df.append(pl.concat(merged_df))
    new_df = pl.concat(new_df)
    return new_df

def join_motif_complements(motif_df):
    if "model" in motif_df.columns:
        if "model" in motif_df.columns:
            n_mod = [x._alpha for x in motif_df["model"]]
            n_nomod = [x._beta for x in motif_df["model"]]

            motif_df = motif_df.with_columns([
                pl.Series("n_mod", n_mod),
                pl.Series("n_nomod", n_nomod)
            ]).drop(["model"])
    motif_df = motif_df.select([
        "contig", "mod_type", "motif", "mod_position", "n_mod", "n_nomod"
    ])
    motif_df = motif_df.with_columns([
        pl.col("motif").map_elements(lambda x: nm.motif.regex_to_iupac(x)).alias("motif")
    ]).with_columns([
        pl.col("motif").map_elements(lambda x: nm.utils.motif_type(x)).alias("motif_type")
    ]).with_columns([
        pl.col("motif").map_elements(lambda x: nm.motif.reverse_compliment(x)).alias("motif_complement")
    ])
    
    
    motifs_out = motif_df \
        .join(
            motif_df.select([
                "contig", "mod_type", "motif_complement", "mod_position", "n_mod", "n_nomod"
            ]),
            left_on = ["contig", "mod_type", "motif"],
            right_on = ["contig", "mod_type", "motif_complement"],
            how = "left",
            suffix = "_complement"
        ) \
        .filter(
            (pl.col("motif") >= pl.col("motif_complement")) | pl.col("mod_position_complement").is_null()
        ) \
        .with_columns([
            pl.when(pl.col("mod_position_complement").is_not_null()) \
                .then(pl.col("motif_complement")) \
                .otherwise(None) \
                .alias("motif_complement")
        ])
    return motifs_out