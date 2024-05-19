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

    
    #log.debug("Removing submotifs")
    # Remove submotifs
    #contig_motifs = nm.postprocess.remove_sub_motifs(motifs_scored_filt)

    contig_motifs = motifs_scored_filt.with_columns(
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
            (pl.col("n_mod_bin") / (pl.col("n_mod_bin") + pl.col("n_nomod_bin"))).alias("mean_methylation")
        )
    return bin_motifs

def merge_bin_motifs(bin_motifs, bins, pileup, assembly):
    assert list(bin_motifs.schema.keys()) == ['bin', 'motif', 'mod_position', 'mod_type', 'n_mod_bin', 'n_nomod_bin', 'contig_count', 'motif_type', 'mean_methylation']

    for (bin, mod_type), df in bin_motifs.groupby("bin", "mod_type"):
        contig_count = df.get_column("contig_count").max()
        
        # Get list of motifs
        motif_seq = df["motif"].to_list()
        motif_pos = df["mod_position"].to_list()
        motifs = [nm.candidate.Motif(nm.candidate.iupac_to_regex(seq), pos) for seq, pos in zip(motif_seq, motif_pos)]
        
        merged_motifs = nm.candidate.merge_motifs(motifs)
        for cluster, motifs in merged_motifs.items():
            merged_motif = motifs[0]
            premerge_motifs = motifs[1]

            premerge_motifs_iupac = [motif.iupac() for motif in premerge_motifs]
            previous_motif_mean_max = bin_motifs.filter(pl.col("motif").is_in(premerge_motifs_iupac)).get_column("mean_methylation").max()

            merge_motif_n_mod = 0
            merge_motif_n_nomod = 0
            for contig in bins.filter(pl.col("bin") == bin).get_column("contig").unique():
                merge_motif_model = nm.evaluate.motif_model_contig(
                            pileup.pileup.filter((pl.col("contig") == contig) & (pl.col("mod_type") == mod_type)), 
                            assembly.assembly[contig].sequence,
                            merged_motif
                        )
                
                merge_motif_n_mod += merge_motif_model._alpha
                merge_motif_n_nomod += merge_motif_model._beta
            merge_motif_mean = merge_motif_n_mod / (merge_motif_n_mod + merge_motif_n_nomod)

            if merge_motif_mean - previous_motif_mean_max > -0.1:
                bin_motifs = bin_motifs.filter(pl.col("motif").is_in(premerge_motifs_iupac).not_())
                bin_motifs = pl.concat(
                    [bin_motifs,
                    pl.DataFrame(
                        {
                            "bin": bin,
                            "motif": merged_motif.iupac(),
                            "mod_position": merged_motif.mod_position,
                            "mod_type": mod_type,
                            "n_mod_bin": merge_motif_n_mod,
                            "n_nomod_bin": merge_motif_n_nomod,
                            "contig_count": contig_count,
                            "motif_type": nm.utils.motif_type(merged_motif.iupac()),
                            "mean_methylation": merge_motif_mean
                        },
                        schema=bin_motifs.schema
                    )]
                )
    return bin_motifs