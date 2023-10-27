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
