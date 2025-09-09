# Output desciption

In general, nanomotif outputs are provided as tab-separated values (TSV) files. Each TSV file contains one or more columns detailing the results, predictions, and underlying data. 

---

## Motif Discovery

**Overview:**  
The motif discovery step identifies DNA methylation motifs within bins. It outputs `bin-motifs.tsv` . This files contain information about the discovered motifs, their modification types, and their degree of methylation.

**Files:**  
- **motifs.tsv**: Contains motifs detected within individual contigs.
- **bin-motifs.tsv**: Contains motifs aggregated and identified at the bin level, representing a consensus from multiple contigs.

**Columns in bin-motifs.tsv:**

| Column                        | Description                                                                                                                                                                              |
|-------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **contig \| bin**             | The identifier of the contig or bin in which the motif was discovered. |
| **mod_type**                  | The type of base modification associated with the motif. Uses single-letter codes (e.g., `a` for 6mA, `m` for 5mC, `21839` for 4mC).                                                   |
| **motif**                     | The detected DNA motif sequence. Uses [IUPAC degenerate nucleotide codes](https://www.bioinformatics.org/sms/iupac.html) to represent positions where multiple nucleotides are possible. |
| **mod_position**              | The zero-based index of the modified base within the motif. For example, if `mod_position` is 0, the modification occurs on the first nucleotide of the motif.                             |
| **n_mod**        | The count of motif occurrences considered "generally methylated" (fraction of mapped bases methylated, ≥0.7 by default).                                                               |
| **n_nomod**    | The count of motif occurrences considered "generally unmethylated" (fraction of mapped bases methylated, <0.2 by default).                                                              |
| **motif_type**                | Classification of the motif as palindromic, non-palindromic, or bipartite. Palindromic motifs read the same forward and backward, while bipartite motifs have a gap separating two distinct parts. |
| **motif_complement**          | The reverse complement of the motif. Only reported if the reverse complement is identified.                                                                                                 |
| **mod_position_complement**   | The zero-based position of the modified base in the motif complement.                                                                                                                    |
| **n_mod_complement**          | The number of generally methylated occurrences of the complement motif.                                                                                                         |
| **n_nomod_complement**        | The number of generally unmethylated occurrences of the complement motif.                                                                                                       |

---

## Bin Improvement

**Overview:**  
The bin improvement module refines genome bins by identifying contigs that may be contaminants. Suspect contigs are recorded in `bin_contamination.tsv`. If requested, cleaned bins (with contaminants removed) can be created.

Additionally, the `include_contigs` command assigns unbinned contigs to existing bins using classification models trained on methylation patterns, resulting in the `include_contigs.tsv` file.

**Files:**  
- **bin_contamination.tsv**: Lists contigs flagged as contamination.
- **include_contigs.tsv**: Lists previously unbinned contigs and their new bin assignments.

### bin_contamination.tsv

Contigs appear in this file if all four clustering methods (agg, spectral, gmm, hdbscan) assign the contig to a cluster different from the bin’s main cluster, suggesting contamination. The each contaminant will have 4 rows, one for each clustering algorithm, along with the cluster stats.
If the `--write_bins` flag is specified new de-contaminated bins will be written to a bins folder.

| Column           | Description                                                                                                                          |
|------------------|--------------------------------------------------------------------------------------------------------------------------------------|
| **contig**        | The identifier of the contig flagged as contamination.                                                                              |
| **bin**           | The name of the bin that the contig currently belongs to.                                                                            |
| **method**        | The clustering method used (agg, spectral, gmm, hdbscan).                                                                            |
| **cluster**       | The cluster ID assigned by the specified method. Methods produce different numeric cluster IDs.                                       |
| **bin_cluster**   | The cluster ID associated with the bin’s original assignment.                                                                         |
| **bin_length**    | The total length (in base pairs) of all contigs currently in the bin.                                                                 |
| **n_contigs_bin** | The total number of contigs currently assigned to the bin.                                                                            |
| **fraction_contigs** | The fraction of the bin’s contigs grouped under `bin_cluster` (e.g., 10/12 ≈ 0.83 if 10 out of 12 contigs are in `bin_cluster`). |
| **fraction_length** | The fraction of the bin’s length represented by contigs in `bin_cluster`. Similar to fraction_contigs but based on sequence length. |

### include_contigs.tsv

Assigning unbinned contigs to bins uses classification methods (Random Forest, Linear Discriminant Analysis, and k-Nearest Neighbors) trained on methylation data.
If all three classifiers assigns an unbinned contig to the same bin, with a join mean probability above 0.80, the contig is assigned to that bin. This is called a `high_confidence` assignment. Besides the aforementioned `high_confidence` assignment, contigs can also be medium and low confidence assigments.
`medium_confidence` assignments are contigs where all three classifiers agree but the join probability is below 0.8. `low_confidence` is when only two classifiers agree. 
`high_confidence` assignments are outputted in a `new_contig_bin.tsv` file. 

| Column          | Description                                                                                                                         |
|-----------------|-------------------------------------------------------------------------------------------------------------------------------------|
| **contig**       | Identifier of the previously unbinned contig.                                                                                       |
| **bin**          | The original bin assignment. "unbinned" indicates no prior assignment.                                                               |
| **assigned_bin** | The bin predicted for this contig by the classifiers.                                                                                |
| **method**       | The classification method (knn, lda, rf). knn = k-Nearest Neighbors, lda = Linear Discriminant Analysis, rf = Random Forest.          |
| **prob**         | The probability assigned by the method for placing the contig in `assigned_bin`.                                                     |
| **mean_prob**    | The average probability from all methods for the contig-bin pair.                                                                    |
| **confidence**   | Qualitative confidence ("high_confidence", "medium_confidence", "low_confidence"). High confidence predictions are more reliable.     |

---

## MTase-Linker

**Overview:**  
The MTase-linker module attempts to link discovered motifs to predicted MTase genes. MTases are part of Restriction-Modification (RM) systems and recognize specific DNA motifs. By integrating motif discovery results with gene prediction and homology searches, MTase-linker outputs files listing MTase genes, their predicted motif specificity, and whether the discovered motifs can be confidently linked to them.

**Files:**  
- **mtase_assignment_table.tsv**: Lists all predicted MTase genes, their likely modification type and motif specificity, and whether a motif is linked.
- **nanomotif_assignment_table.tsv**: Similar to `bin-motifs.tsv`, with added columns indicating whether a motif is linked to a specific MTase.

### mtase_assignment_table.tsv

| Column          | Description                                                                                                                                                  |
|-----------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **bin**         | The bin/genome in which the MTase gene is located.                                                                                                            |
| **gene_id**     | A unique identifier for each MTase gene, combining the contig ID and a gene number (e.g., "contig_12_5").                                                    |
| **contig**      | The contig on which the MTase gene is found.                                                                                                                  |
| **mod_type_pred** | The predicted modification type (e.g., ac for 6mA/4mC, m for 5mC).                                                                                          |
| **sub_type**    | Predicted RM system subtype (I, II, IIG, III).                                                                                                               |
| **RM_system**   | TRUE/FALSE indicating if the MTase is part of a complete RM system.                                                                                           |
| **motif_type_pred** | Predicted motif type that the MTase may recognize (palindromic, non-palindromic, bipartite).                                                           |
| **REbase_ID**   | Closest homolog MTase from REbase, if identity and coverage meet the threshold (≥80%).                                                                       |
| **motif_pred**  | Predicted recognition motif based on REbase annotation.                                                                                                      |
| **linked**      | TRUE/FALSE indicating if a motif could be unambiguously linked to the MTase gene. If TRUE, see `detected_motif`.                                             |
| **detected_motif** | The exact motif sequence linked to this MTase, if `linked` is TRUE.                                                                                        |

### nanomotif_assignment_table.tsv

This file includes data from `bin-motifs.tsv` plus two additional columns: `linked` and `candidate_genes`.  
- **linked (TRUE/FALSE):** Indicates if the motif can be unambiguously attributed to an MTase.  
- **candidate_genes:** Lists the MTase genes linked to the motif if `linked` = TRUE, or other potential candidates if not.

---

## Additional Outputs From Dependencies

MTase-linker uses external tools for gene prediction and annotation:

- **DefenseFinder output:** `.../defensefinder/{assembly_name}_processed_defense_finder_systems.tsv` lists RM-system gene annotations to verify if MTases form part of a complete RM-system.
- **MTase amino acid sequences:** `.../defensefinder/{assembly_name}_processed_defense_finder_mtase.faa` contains predicted MTase amino acid sequences.
- **All predicted protein sequences:** `.../prodigal/{assembly_name}_processed.faa` contains amino acid sequences for all predicted proteins, including MTases.

---
