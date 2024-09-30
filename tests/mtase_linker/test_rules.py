import subprocess
import os
import filecmp


def run_snakemake():
        cmd = [
            "snakemake", 
            "--snakefile", "nanomotif/mtase_linker/MTase_linker.smk",
            "--cores", "1",
            "--config", 
            "THREADS=10",
            "ASSEMBLY=nanomotif/datasets/mtase_linker_check_installation/assembly.polished.fa", 
            "CONTIG_BIN=nanomotif/datasets/mtase_linker_check_installation/contig_bin.tsv", 
            "OUTPUTDIRECTORY=datasets/mtase_linker_testdata/results", 
            "DEPENDENCY_PATH=../restriction-modification-annotation/data/ML_dependencies", "IDENTITY=80",
            "QCOVS=80",
            "NANOMOTIF=nanomotif/datasets/mtase_linker_check_installation/bin-motifs.tsv",            
            "--use-conda",
            "--conda-prefix", 
            "../restriction-modification-annotation/data/ML_dependencies/ML_envs",
            "--forceall"
                
            ]
        subprocess.run(cmd, cwd=os.path.join(os.getcwd()), check=True)
        
def compare_files(file1, file2):
    """Compare two files and return True if they are the same, False otherwise."""
    return filecmp.cmp(file1, file2, shallow=False)




def test_snakemakepipeline():


    run_snakemake()

    #### PRODIGAL ####
    # Result files
    prodigal_result_file1 = "datasets/mtase_linker_testdata/results/prodigal/assembly.polished.faa"
    prodigal_result_file2 = "datasets/mtase_linker_testdata/results/prodigal/assembly.polished.gff"

    # Expected files
    prodigal_expected_file1 = "datasets/mtase_linker_testdata/expected_results/prodigal/assembly.polished.faa"
    prodigal_expected_file2 = "datasets/mtase_linker_testdata/expected_results/prodigal/assembly.polished.gff"

    # Prodigal comparison
    assert compare_files(prodigal_result_file1, prodigal_expected_file1), "Test failed: Prodigal output file1 do not match."
    assert compare_files(prodigal_result_file2, prodigal_expected_file2), "Test failed: Prodigal output file2 do not match."

    #### PROCESS PRODIGAL ####

    # Result files
    process_prodigal_result_file = "datasets/mtase_linker_testdata/results/prodigal/assembly.polished_processed.faa"
    
    # Expected files
    process_prodigal_expected_file = "datasets/mtase_linker_testdata/expected_results/prodigal/assembly.polished_processed.faa"

    # Process prodigal comparison
    assert compare_files(process_prodigal_result_file, process_prodigal_expected_file), "Test failed: Process prodigal output file do not match."

    #### DEFENSE FINDER ####

    # Result files
    defense_finder_result_file1 = "datasets/mtase_linker_testdata/results/defensefinder/assembly.polished_processed_defense_finder_hmmer.tsv"
    defense_finder_result_file2 = "datasets/mtase_linker_testdata/results/defensefinder/assembly.polished_processed_defense_finder_systems.tsv"

    # Expected files
    defense_finder_expected_file1 = "datasets/mtase_linker_testdata/expected_results/defensefinder/assembly.polished_processed_defense_finder_hmmer.tsv"
    defense_finder_expected_file2 = "datasets/mtase_linker_testdata/expected_results/defensefinder/assembly.polished_processed_defense_finder_systems.tsv"

    # Defense finder comparison
    assert compare_files(defense_finder_result_file1, defense_finder_expected_file1), "Test failed: Defense finder output file1 do not match."
    assert compare_files(defense_finder_result_file2, defense_finder_expected_file2), "Test failed: Defense finder output file2 do not match."

    #### extract_mtase ####

    # Result files
    extract_mtase_result_file1 = "datasets/mtase_linker_testdata/results/defensefinder/assembly.polished_processed_defense_finder_mtase.faa"
    extract_mtase_result_file2 = "datasets/mtase_linker_testdata/results/defensefinder/assembly.polished_processed_defense_finder_mtase.tsv"

    # Expected files
    extract_mtase_expected_file1 = "datasets/mtase_linker_testdata/expected_results/defensefinder/assembly.polished_processed_defense_finder_mtase.faa"
    extract_mtase_expected_file2 = "datasets/mtase_linker_testdata/expected_results/defensefinder/assembly.polished_processed_defense_finder_mtase.tsv"

    # Extract mtase comparison
    assert compare_files(extract_mtase_result_file1, extract_mtase_expected_file1), "Test failed: extract_MTase_protein_seqs output file1 do not match."
    assert compare_files(extract_mtase_result_file2, extract_mtase_expected_file2), "Test failed: extract_MTase_protein_seqs output file2 do not match."

    #### blastp Rebase ####

    # Result files
    blastp_rebase_result_file = "datasets/mtase_linker_testdata/results/blastp/assembly.polished_rebase_mtase_alignment.tsv"

    # Expected files
    blastp_rebase_expected_file = "datasets/mtase_linker_testdata/expected_results/blastp/assembly.polished_rebase_mtase_alignment.tsv"

    # Blastp Rebase comparison
    assert compare_files(blastp_rebase_result_file, blastp_rebase_expected_file), "Test failed: blastp_REbase output file do not match."

    #### blastp Sign ####

    # Result files
    blastp_sign_result_file = "datasets/mtase_linker_testdata/results/blastp/assembly.polished_rebase_mtase_sign_alignment.tsv"

    # Expected files
    blastp_sign_expected_file = "datasets/mtase_linker_testdata/expected_results/blastp/assembly.polished_rebase_mtase_sign_alignment.tsv"

    # Blastp Sign comparison
    assert compare_files(blastp_sign_result_file, blastp_sign_expected_file), "Test failed: blastp_Sign output file do not match."

    #### run_pfam_hmm ####

    # Result files
    run_pfam_hmm_result_file1 = "datasets/mtase_linker_testdata/results/pfam_hmm_hits/assembly.polished_raw_hmm_hits_mtase_aaseqs.hits"
    run_pfam_hmm_result_file2 = "datasets/mtase_linker_testdata/results/pfam_hmm_hits/assembly.polished_hmm_hits_mtase_aaseqs.tsv"

    # Expected files
    run_pfam_hmm_expected_file1 = "datasets/mtase_linker_testdata/expected_results/pfam_hmm_hits/assembly.polished_raw_hmm_hits_mtase_aaseqs.hits"
    run_pfam_hmm_expected_file2 = "datasets/mtase_linker_testdata/expected_results/pfam_hmm_hits/assembly.polished_hmm_hits_mtase_aaseqs.tsv"

    # Run pfam hmm comparison
    #assert compare_files(run_pfam_hmm_result_file1, run_pfam_hmm_expected_file1), "Test failed: run_pfam_hmm output file1 do not match." #Comparison fails due to timestamps
    assert compare_files(run_pfam_hmm_result_file2, run_pfam_hmm_expected_file2), "Test failed: run_pfam_hmm output file2 do not match."

    #### mod_type_prediction ####

    # Result files
    mod_type_prediction_result_file = "datasets/mtase_linker_testdata/results/pfam_hmm_hits/assembly.polished_gene_id_mod_table.tsv"

    # Expected files
    mod_type_prediction_expected_file = "datasets/mtase_linker_testdata/expected_results/pfam_hmm_hits/assembly.polished_gene_id_mod_table.tsv"

    # Mod type prediction comparison
    assert compare_files(mod_type_prediction_result_file, mod_type_prediction_expected_file), "Test failed: mod_type_prediction output file do not match."

    #### motif_assignment ####

    # Result files
    motif_assignment_result_file1 = "datasets/mtase_linker_testdata/results/mtase_assignment_table.tsv"
    motif_assignment_result_file2 = "datasets/mtase_linker_testdata/results/nanomotif_assignment_table.tsv"

    # Expected files
    motif_assignment_expected_file1 = "datasets/mtase_linker_testdata/expected_results/mtase_assignment_table.tsv"
    motif_assignment_expected_file2 = "datasets/mtase_linker_testdata/expected_results/nanomotif_assignment_table.tsv"

    # Motif assignment comparison
    assert compare_files(motif_assignment_result_file1, motif_assignment_expected_file1), "Test failed: motif_assignment output file1 do not match."
    assert compare_files(motif_assignment_result_file2, motif_assignment_expected_file2), "Test failed: motif_assignment output file2 do not match."

if __name__ == '__main__':
    test_snakemakepipeline()
