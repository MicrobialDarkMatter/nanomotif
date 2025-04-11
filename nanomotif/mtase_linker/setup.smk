##### Setup rules #####

rule all:
    input:
        "envs/prodigal-2.6.3.yaml",
        "envs/defensefinder-2.0.0.yaml",
        "envs/python_env-3.12.0.yaml",
        "envs/blast-2.14.1.yaml"

# Rule for specifying prodigal environment
rule prodigal:
    output:
        touch("envs/prodigal-2.6.3.yaml")
    conda:
        "envs/prodigal-2.6.3.yaml"
    shell:
        """
        touch {output}
        """

# Rule for specifying defensefinder environment
rule defensefinder:
    output:
        touch("envs/defensefinder-2.0.0.yaml")
    conda:
        "envs/defensefinder-2.0.0.yaml"
    shell:
        """
        touch {output}
        """

# Rule for specifying python_env environment
rule python_env:
    output:
        touch("envs/python_env-3.12.0.yaml")
    conda:
        "envs/python_env-3.12.0.yaml"
    shell:
        """
        touch {output}
        """

# Rule for specifying blast-2.14.1.yaml environment
rule blastp:
    output:
        touch("envs/blast-2.14.1.yaml")
    conda:
        "envs/blast-2.14.1.yaml"
    shell:
        """
        touch {output}
        """

# Defensefinder to extract MTase genes

rule defense_finder_update:
    conda:
        "envs/defensefinder-2.0.0.yaml"
    output:
        dummy=touch(os.path.join(config["DEPENDENCYDIR"], "df_models", "defense-finder_update.done"))
    params: 
        models_dir= os.path.join(config["DEPENDENCYDIR"], "df_models")
    shell:
        """
        defense-finder update --models-dir {params.models_dir}
        """

rule BLASTP_dbget:
    output:
        os.path.join(config["DEPENDENCYDIR"], "REbase_RecSeq_all.fa")
    params: 
        url1="https://raw.githubusercontent.com/JSBoejer/MTase-Linker_models/main/REbase_all_protein_seqs_RecSeqs_anno_format_08_01_2024.fa.gz",
        url2="https://raw.githubusercontent.com/JSBoejer/MTase-Linker_models/main/PFAM_MTase_profiles.hmm",
        temp_output1=os.path.join(config["DEPENDENCYDIR"], "REbase_all_protein_seqs_RecSeqs_anno_format_08_01_2024.fa.gz"),
        temp_output2=os.path.join(config["DEPENDENCYDIR"], "PFAM_MTase_profiles.hmm")
    shell:
        """
        wget -O {params.temp_output1} {params.url1} && gunzip -c {params.temp_output1} > {output}
        wget -O {params.temp_output2} {params.url2}
        """

# Make BLASTP database
rule make_REbase_db:
    input:
      os.path.join(config["DEPENDENCYDIR"], "REbase_RecSeq_all.fa")
    output:
        os.path.join(config["DEPENDENCYDIR"], "REbase_RecSeq_all", "REbase_RecSeq_all.pdb")
    params:
        os.path.join(config["DEPENDENCYDIR"], "REbase_RecSeq_all", "REbase_RecSeq_all")
    conda:
        "envs/blast-2.14.1.yaml"
    shell:
        """
      makeblastdb -in {input} -dbtype prot -title REbase_RecSeq_all -out {params}
        """