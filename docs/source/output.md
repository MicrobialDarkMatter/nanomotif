# Output explanation
## Motif discovery

Motif discovery output two primary files: bin-motifs.tsv and motif.tsv. `bin-motifs.tsv` contains the motifs found in each bin, while `motif.tsv` contains motifs found for each contig. Both files follow the format in the table below.

| **Column**                   | **Description**                                                                                                                                                               |
|---------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **contig/bin**                   | Te contig/bin wherein the motif has been discovered |
| **mod_type**              | Modification type of the motif. Follows single letter code (a=6mA, m=5mC, 21839=4mC)                                                         |
| **motif**                 | The sequence of the motif.  [IUPAC](https://www.bioinformatics.org/sms/iupac.html) letters are used for positions matching multiple nucleotides                                        |
| **mod_position**          | The position of the modified base within the motif sequence following in 0-based index.                            |
| **n_mod(_bin)**             | Number of motif position where the fraction of mapped bases methylated is above the generally methylated threshold (default: 0.7)                                                          |
| **n_nomod(_bin)**             | Number of motif position where the fraction of mapped bases methylated is below the generally methylated threshold (default: 0.7)                                                                                             |
| **motif_type**            | Indicates whether the motif sequence is a palindrome, non-palindromic or bipartite. |
| **motif_complement**      | The corresponding reverse complement motif if identified.                                                 |
| **mod_position_complement** | Same as `mod_position`, but for reverse complement |
| **n_mod_complement**      | Same as `n_mod`, but for reverse complement                                                   |
| **n_nomod_complement**    | Same as `n_nomod`, but for reverse complement                                               |



## Bin improvement

The bin improvement module has two primary outputs: `bin_contamination.tsv` and `include_contigs.tsv`.

**Contamination detection**

Running the nanomotif detect_contamination command will generate a `bin_contamination.tsv` specifying the contigs, which is flagged as contamination. If the --write_bins and the --assembly_file flags are specified new de-contaminated bins will be written to a bins folder. Contigs are only reported in the `bin_contamination.tsv` if 4/4 clustering methods (agg, spectral, gmm, hdbscan) assigns the contig to a different `cluster` than `bin_cluster`.


| **Column**     | **Description**         |
|----|------|
| **contig**          | The identifier of a specific contig.                                                                       |
| **bin**             | The name of the bin that this contig currently belongs to.                                    |
| **method**          | The clustering used to assign the contig to a cluster (e.g., `agg` for agglomerative clustering, `spectral` for spectral clustering, `gmm` for Gaussian Mixture Model, `hdbscan` for HDBSCAN clustering). |
| **cluster**         | The cluster ID assigned to the contig by the specified method. Each method may produce different cluster identifiers.                                                                     |
| **bin_cluster**     | The cluster ID of the current bin.    |
| **bin_length**      | The total length (in base pairs) of all current contigs within the bin.                                                                                                                           |
| **n_contigs_bin**   | The total number of contigs that the bin currently contains.                                                                                                                                        |
| **fraction_contigs** | The fraction of contigs in the bin that is clustered within `bin_cluster`. For example, if the bin has 12 contigs and the `bin_cluster` includes 10 of them, the fraction would be 10/12 ≈ 0.83. |
| **fraction_length**  | The fraction of the total bin length covered by the contigs in the `bin_cluster`. Similar to `fraction_contigs` but based on sequence length.                           |


**Contig inclusion**

The include_contigs command assigns unbinned contigs in the assembly file to bins by training three classifiers, random forest, linear discriminant analysis, and k neighbors classifier, on the methylation pattern of the bins. 

| Header          | Explanation                                                                                                              |
|-----------------|--------------------------------------------------------------------------------------------------------------------------|
| **contig**          | The identifier of a specific contig.                                                                       |
| **bin**             | The name of the bin that this contig previously belongs to. "unbinned" indicates it was initially not assigned to a bin. |
| **assigned_bin** | The bin to which the contig is assigned by the classification.                                     |
| **method**       | The method used to predict the bin assignment (e.g., knn, lda, rf). "knn" = k-Nearest Neighbors, "lda" = Linear Discriminant Analysis, "rf" = Random Forest. |
| **prob**         | The probability assigned by the method for placing the contig into the assigned bin.                          |
| **mean_prob**    | The average probability aggregated from all methods probability for the same contig-bin pair.     |
| **confidence**   | The qualitative confidence level of the assignment (e.g., "high_confidence", "medium_confidence", "low_confidence"). We recommend only using "high_confidence" assignments if no other verifications are performed.   |


## MTase-linker

Running the nanomotif MTase-linker run command will generate two primary output files: mtase_assignment_table.tsv and nanomotif_assignment_table.tsv. 


The mtase_assignment_table.tsv file lists all predicted MTase genes in the genome along with their predicted methylation characteristics and whether the module was able to unambiguously assign any detected motifs to the MTase (`linked` = (True/False)).

| **Column**       | **Description**                                                                                       |
|------------------|-------------------------------------------------------------------------------------------------------|
| **bin**          | bin/genome in which the MTase gene is found                                                                         |
| **gene_id** | A unique identifier for each MTase gene, consisting of the ordinal ID of the contig and an ordinal ID of that gene within the contig. The name is a concatenation of contig ID + "_" + gene numer. |
| **contig** | ID of contig on which the MTase gene was found. |
| **mod_type_pred** | Predicted modification type of the MTase product (ac = 6mA/4mC, m = 5mC). |
| **sub_type** | Predicted restriction modification sub-type of MTase product (I, II, IIG, III). |
| **RM_system** | Binary indicator of whether the MTase gene was found within a complete RM-system (TRUE/FALSE). |
| **motif_type_pred** | Prediction of the motif type of the MTase product. Based on the sub_type. (Palindromic, non-palindromic, and bipartite). |
| **REbase_ID** | ID of closets homolog MTase in REbase if sequence identity and query coverage is >= THRESHOLD (Default: THRESHOLD = 80). |
| **motif_pred** | Predicted recognition motif of MTase based REbase annotation. |
| **linked** | Binary indicator if a motif could be unambiguously linked the MTase gene in the bin (TRUE/FALSE). If TRUE see detected_motif for motif. |
| **detected_motif** | DNA sequence of motif, which could be unambiguously linked to the MTase in the bin. |

nanomotif_assignment_table.tsv files includes data from the bin-motifs.tsv of the nanomotif output with two additional columns `linked` and `candidate_genes`. The `linked` variable is a boolean indicator if the motif could be unambiguously linked to a MTase in the bin/genome (TRUE/FALSE). If True the gene_id of the MTase is provided in `candidate_gene`. If False, the `candidate_gene` variable lists feasible candidate facilitators of the modification based on motif type and modification type predictions.

MTase-linker also provides the outputs from the dependency tools, BLASTP, DefenseFinder and Prodigal.

- If you are interested in RM systems the defensefinder output file ".../defensefinder/{assembly_name}_processed_defense_finder_systems.tsv" may be relevant to determine RM-system genes associated with the MTases. 

- The amino acid sequences of the MTase gene products can be found in ".../defensefinder/{assembly_name}_processed_defense_finder_mtase.faa".

- The amino acid sequences for all proteins predicted can be found in .../prodigal/mmlong2_lite_assembly_processed.faa".