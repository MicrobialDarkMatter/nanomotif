# Output explanation
## Motif discovery

COMING SOON

## Bin improvement

COMING SOON

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