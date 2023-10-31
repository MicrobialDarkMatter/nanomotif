import nanomotif as nm
import polars as pl
from polars import col
import logging as log

def remove_noisy_motifs(motif_df):
    motif_df_clean = []
    for contig, df in motif_df.groupby("contig"):
        motif_strings = df.get_column("motif").to_list()
        positions = df.get_column("mod_position").to_list()
        motifs = [nm.candidate.Motif(motif_string, pos) for motif_string, pos in zip(motif_strings, positions)]
        clean_motifs = []
        for motif in motifs:
            if not motif.have_isolated_bases():
                clean_motifs.append(motif)
        df_clean = df.filter(col("motif").is_in(clean_motifs))
        motif_df_clean.append(df_clean)
    motif_df_clean = pl.concat(motif_df_clean)
    return motif_df_clean


def remove_sub_motifs(motif_df):
    motif_df_clean = []
    for contig, df in motif_df.groupby("contig"):
        motif_strings = df.get_column("motif").to_list()
        positions = df.get_column("mod_position").to_list()
        motifs = [nm.candidate.Motif(motif_string, pos) for motif_string, pos in zip(motif_strings, positions)]
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
        df_clean = df.filter(col("motif").is_in(parent_motifs))
        motif_df_clean.append(df_clean)
    motif_df_clean = pl.concat(motif_df_clean)
    return motif_df_clean
<<<<<<< HEAD

def merge_motifs_in_df(motif_df, pileup, assembly):
    new_df = []
    for (contig, mod_type), df in motif_df.groupby("contig", "mod_type"):
        # Get list of motifs
        motif_seq = df["motif"].to_list()
        motif_pos = df["mod_position"].to_list()
        motifs = [nm.candidate.Motif(seq, pos) for seq, pos in zip(motif_seq, motif_pos)]

        # Merge motifs
        merged_motifs, single_motifs = nm.candidate.merge_motifs(motifs)
        
        # Create a new dataframe with the non merged motifs
        single_df = df.filter(col("motif").is_in(single_motifs))
        new_df.append(single_df)

        if len(merged_motifs) == 0:
            # No motifs were mergged
            continue
        merged_df = []
        for motif in merged_motifs:
            merged_model = nm.evaluate.motif_model_contig(
                pileup.filter((col("contig") == contig) & (col("mod_type") == mod_type)), 
                assembly.assembly[contig].sequence,
                motif
            )
            merged_df.append(pl.DataFrame({
                "sequence": motif.string,
                "model": merged_model,
                "score": -1.,
                "contig": contig,
                "mod_type": mod_type,
                "motif": motif.string,
                "mod_position": motif.mod_position
            }))
        new_df.append(pl.concat(merged_df))
    new_df = pl.concat(new_df)
    return new_df
=======
>>>>>>> 756ec0aaa02f67bccb18ad0f571059b2960f54b5
