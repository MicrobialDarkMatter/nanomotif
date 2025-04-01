#load packages
import pandas as pd
import os

# Defining config file location
#configfile: "config/mtase_linker_config.yaml"

# Loading unique bin names
#bin_ids = glob_wildcards(os.path.join(config["BINSDIR"], "{binname}.fa")).binname

# Defining output directory
OUTPUTDIR = config["OUTPUTDIRECTORY"]
DEPENDENCYDIR = config["DEPENDENCY_PATH"]
ASSEMBLY_PATH = config["ASSEMBLY"]
BASENAME = os.path.splitext(os.path.basename(ASSEMBLY_PATH))[0]

#################################################################################
                        #MTase-Linker Pipeline
#################################################################################


rule all:
    input:
        os.path.join(OUTPUTDIR, "mtase_assignment_table.tsv")



# Running prodigal to extract protein sequences

rule prodigal:
     input:
          ASSEMBLY_PATH
     output:
          gene_coords = os.path.join(OUTPUTDIR, "prodigal", f"{BASENAME}.gff"),
          protein_translations = os.path.join(OUTPUTDIR, "prodigal", f"{BASENAME}.faa")
     threads: config["THREADS"]
     conda: 
          "envs/prodigal-2.6.3.yaml"
     shell:
          """
          prodigal -i {input} -o {output.gene_coords} -a {output.protein_translations} -f gff
          """


# Removing asterisk from AA sequence ends in prodigal sequence file
rule process_prodigal:
    input:
        os.path.join(OUTPUTDIR, "prodigal", f"{BASENAME}.faa")
    output:
        os.path.join(OUTPUTDIR, "prodigal", f"{BASENAME}_processed.faa")
    shell:
        """
        sed '/^[^>]/ s/\*$//' {input} > {output}
        """


# Defensefinder to extract MTase genes
rule defenseFinder:
    input:
        faa = os.path.join(OUTPUTDIR, "prodigal", f"{BASENAME}_processed.faa"),
        #dummy = os.path.join(OUTPUTDIR,"defense-finder_update.done")
    output: 
        hmmer_file = os.path.join(OUTPUTDIR, "defensefinder", f"{BASENAME}_processed_defense_finder_hmmer.tsv"),
        system_file = os.path.join(OUTPUTDIR, "defensefinder", f"{BASENAME}_processed_defense_finder_systems.tsv")
    params:
        modelsdir = os.path.join(DEPENDENCYDIR, "df_models"),
        outputdir = os.path.join(OUTPUTDIR, "defensefinder")
    threads: config["THREADS"]
    conda:
        "envs/defensefinder-2.0.0.yaml"
    shell:
        """
        defense-finder run --db-type gembase --workers {threads} --models-dir {params.modelsdir} --out-dir {params.outputdir} {input.faa}
        """


# Extracting MTase protein ids and AA sequences with defensefinder output
rule extract_MTase_protein_seqs:
    input:
        hmmer_file = os.path.join(OUTPUTDIR, "defensefinder", f"{BASENAME}_processed_defense_finder_hmmer.tsv"),
        system_file = os.path.join(OUTPUTDIR, "defensefinder", f"{BASENAME}_processed_defense_finder_systems.tsv"),
        prodigal_AAs = os.path.join(OUTPUTDIR, "prodigal", f"{BASENAME}_processed.faa")
    output:
        DF_MTase_table = os.path.join(OUTPUTDIR, "defensefinder", f"{BASENAME}_processed_defense_finder_mtase.tsv"),
        DF_MTase_AA = os.path.join(OUTPUTDIR, "defensefinder", f"{BASENAME}_processed_defense_finder_mtase.faa"),
    conda:
        "envs/python_env-3.12.0.yaml",
    script:
        "src/extract_MTase_genes.py"



# Blasting MTase sequences against REbase
rule blastp_REbase:
    input: 
        query = os.path.join(OUTPUTDIR, "defensefinder", f"{BASENAME}_processed_defense_finder_mtase.faa")
    output:
        blast = temp(os.path.join(OUTPUTDIR, "blastp", f"{BASENAME}_rebase_mtase_alignment_tmp.tsv")),
        blast_header = os.path.join(OUTPUTDIR, "blastp", f"{BASENAME}_rebase_mtase_alignment.tsv")
    params: 
        max_target_seqs=100,
        db = os.path.join(DEPENDENCYDIR, "REbase_RecSeq_all", "REbase_RecSeq_all")
    conda:
        "envs/blast-2.14.1.yaml"
    threads: config["THREADS"]
    shell:
        """
        blastp -query {input.query} -db {params.db} -out {output.blast} \
        -outfmt "6 std qcovs salltitles" \
        -max_target_seqs {params.max_target_seqs} \
        -num_threads {threads}

        if [ -s {output.blast} ]; then
            echo -e "qseqid\\tsseqid\\tpident\\tlength\\tmismatch\\tgapopen\\tqstart\\tqend\\tsstart\\tsend\\tevalue\\tbitscore\\tqcovs\\tsalltitles" \
            | cat - {output.blast} > {output.blast_header}
        else
            echo -e "qseqid\\tsseqid\\tpident\\tlength\\tmismatch\\tgapopen\\tqstart\\tqend\\tsstart\\tsend\\tevalue\\tbitscore\\tqcovs\\tsalltitles" > {output.blast_header}
        fi
        """



# Retrieving significant BLASTP hits for motif guess
rule blastp_sign_hits:
    input:
        os.path.join(OUTPUTDIR, "blastp", f"{BASENAME}_rebase_mtase_alignment.tsv")
    output:
        os.path.join(OUTPUTDIR, "blastp", f"{BASENAME}_rebase_mtase_sign_alignment.tsv")
    params: 
        identity = config["IDENTITY"],
        qcovs = config["QCOVS"]
    conda:
        "envs/python_env-3.12.0.yaml"
    script:
        "src/blastp_processing.py"




# Mod type predictions
rule run_pfam_hmm:
    input:
        DF_MTase_AA = os.path.join(OUTPUTDIR, "defensefinder", f"{BASENAME}_processed_defense_finder_mtase.faa")
    output:
        hmmsearch = os.path.join(OUTPUTDIR, "pfam_hmm_hits", f"{BASENAME}_raw_hmm_hits_mtase_aaseqs.hits"),
        hits = os.path.join(OUTPUTDIR, "pfam_hmm_hits", f"{BASENAME}_hmm_hits_mtase_aaseqs.tsv")
    params: 
        hmm_profiles = os.path.join(DEPENDENCYDIR, "PFAM_MTase_profiles.hmm")
    threads: config["THREADS"]
    conda: 
        "envs/defensefinder-2.0.0.yaml"
    shell:
        """
        if [ -s {input.DF_MTase_AA} ]; then
            hmmsearch --cut_ga --cpu {threads} --tblout {output.hmmsearch} {params.hmm_profiles} {input.DF_MTase_AA}

            if [ -s {output.hmmsearch} ] && grep -q -v '^#' {output.hmmsearch}; then
                grep -v '^#' {output.hmmsearch} | awk -v OFS='\\t' '{{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19}}' > {output.hits}
            else
                echo "No matches found or only header lines present in hmmsearch output." >&2
                touch {output.hits}
            fi
        else
            touch {output.hmmsearch} {output.hits}
        fi
        """



# Converting pfam hits to mod type predictions:
rule mod_type_prediction_table:
    input:
        os.path.join(OUTPUTDIR, "pfam_hmm_hits", f"{BASENAME}_hmm_hits_mtase_aaseqs.tsv")
    output:
        os.path.join(OUTPUTDIR, "pfam_hmm_hits", f"{BASENAME}_gene_id_mod_table.tsv")
    conda:
        "envs/python_env-3.12.0.yaml"
    script:
        "src/mod_typing.py"



# Generating final output
rule motif_assignment:
    input:
        mod_table_path = os.path.join(OUTPUTDIR, "pfam_hmm_hits", f"{BASENAME}_gene_id_mod_table.tsv"),
        DefenseFinder_path = os.path.join(OUTPUTDIR, "defensefinder", f"{BASENAME}_processed_defense_finder_mtase.tsv"),
        BLASTP_path = os.path.join(OUTPUTDIR, "blastp", f"{BASENAME}_rebase_mtase_sign_alignment.tsv"),
        nanomotif_table_path = config["NANOMOTIF"],
        contig_bin = config["CONTIG_BIN"]
    output:
        MTase_assignment_table = os.path.join(OUTPUTDIR, "mtase_assignment_table.tsv"),
        nanomotif_assignment_table = os.path.join(OUTPUTDIR, "nanomotif_assignment_table.tsv")
    params:
        MINIMUM_METHYLATION = config["MINIMUM_METHYLATION"],
    conda:
        "envs/python_env-3.12.0.yaml"
    script:
        "src/motif_assignment.py"

