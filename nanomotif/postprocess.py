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
        return motif_df
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

def get_motif_parental_relationship(motifs):
    relationships = []
    for i, motif1 in enumerate(motifs):
        for j, motif2 in enumerate(motifs):
            if i != j and motif1.sub_string_of(motif2):
                # Avoid adding duplicates
                if (motif2, motif1) not in relationships:
                    relationships.append((motif2, motif1))
    return relationships


def remove_sub_motifs(motif_df: nm.motif.MotifSearchResult):
    assert "motif" in motif_df.columns
    assert "mod_position" in motif_df.columns
    assert "model" in motif_df.columns
    assert len(motif_df) > 0
    for (contig, mod_type), df in motif_df.group_by("reference", "mod_type"):
        log.debug(f"Processing evaluating sub-motifs for: {contig}, mod_type {mod_type}")
        motif_strings = df.get_column("motif").to_list()
        positions = df.get_column("mod_position").to_list()
        motifs = [nm.motif.Motif(motif_string, pos) for motif_string, pos in zip(motif_strings, positions)]
        parent_motifs = get_motif_parental_relationship(motifs)
        log.debug(f"Found {len(parent_motifs)} sub-motifs motifs from {len(motifs)} total motifs")
        for parent, child in parent_motifs:
            if child.string not in df["motif"].to_list() or child.mod_position not in df["mod_position"].to_list():
                log.warning(f"Child motif {child} not found in dataframe")
                continue
            model_child = df.filter(col("motif").eq(child.string) & col("mod_position").eq(child.mod_position))["model"][0]
            model_parent = df.filter(col("motif").eq(parent.string) & col("mod_position").eq(parent.mod_position))["model"][0]
            descend_score = nm.find_motifs_bin.predictive_evaluation_score(model_child, model_parent)
            if descend_score > 0.5:
                log.info(f"Keeping sub motif {child} as it has a high predictive score of {descend_score} against parent motif {parent}")
                motif_to_discard = parent
            else:
                log.info(f"Removing sub motif {child} as it has a low predictive score of {descend_score} against parent motif {parent}")
                motif_to_discard = child
            motif_df = motif_df.remove(
                col("motif").eq(motif_to_discard.string) & 
                col("mod_position").eq(motif_to_discard.mod_position) &
                col("mod_type").eq(mod_type)
            )
    return motif_df

def merge_motifs_in_df(motif_df: nm.motif.MotifSearchResult, pileup, assembly, mean_shift_threshold = -0.2):
    new_df = []
    for (contig, mod_type), df in motif_df.group_by("reference", "mod_type"):
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

def join_motif_complements(motif_df: nm.motif.MotifSearchResult):
    motif_rev_complement = motif_df.with_columns([
        pl.col("motif_iupac").map_elements(
            lambda x: nm.motif.reverse_compliment(x), 
            return_dtype=pl.Utf8
        ).alias("motif_reverse")
    ])


    motifs_out = motif_df \
        .join(
            motif_rev_complement,
            left_on = ["reference", "mod_type", "motif_iupac"],
            right_on = ["reference", "mod_type", "motif_reverse"],
            how = "left",
            suffix = "_complement"
        ).filter(
            (pl.col("motif_iupac") >= pl.col("motif_complement")) | pl.col("mod_position_complement").is_null()
        ).with_columns([
            pl.when(pl.col("mod_position_complement").is_not_null()) \
                .then(pl.col("motif_complement")) \
                .otherwise(None) \
                .alias("motif_complement")
        ])
    return motifs_out