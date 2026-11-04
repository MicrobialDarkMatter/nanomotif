# Output description

In general, nanomotif outputs are provided as tab-separated values (TSV) files. Each TSV file contains one or more columns detailing the results, predictions, and underlying data. 

---

## Motif Discovery

**Overview:**
Motif discovery writes the motifs found in each bin to `bin-motifs.tsv`, with how often each motif is
methylated.

**Files in the output folder:**
- **bin-motifs.tsv**: the motifs (columns below). If no motifs are found, the file has only the header.
- **precleanup-motifs/**: per bin and modification type, the motifs after each clean-up step (`motifs.tsv`,
  then `-noise`, `-merge`, `-sub` and `-complement`; see [How it works](how_it_works.md#6-clean-up)).
  Useful to see why a motif was dropped.
- **temp/**: per bin, the sequence windows around methylated positions, the background windows, their
  base frequencies and the search graph (`motif_graph_<mod>.gml`), for debugging.
- **logs/**: the main log (`motif_discovery.main.log`, also printed to the terminal) and one log per worker
  process (`find-motifs.<process id>.log`), with each task's peak memory.
- **args.motif_discovery.json**: the arguments of the run.

**Columns in bin-motifs.tsv:**

| Column                        | Description                                                                                                                                                                              |
|-------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **reference**                 | The bin (or contig, or genome) in which the motif was found. |
| **motif**                     | The motif. Uses [IUPAC degenerate nucleotide codes](https://www.bioinformatics.org/sms/iupac.html) for positions where several bases are possible, and `N` for any base. |
| **mod_position**              | The zero-based position of the modified base in the motif. For example, `GATC` with `mod_position` 1 is methylated at the A. |
| **mod_type**                  | The modification type, as modkit's code: `a` for 6mA, `m` for 5mC, `21839` for 4mC. |
| **n_mod**                     | Motif occurrences (both strands) that are methylated: at least `--methylation_threshold_high` (0.7) of their reads called modified. |
| **n_nomod**                   | Motif occurrences that are unmethylated: at most `--methylation_threshold_low` (0.3) of their reads called modified. Occurrences in between, or with fewer reads than `--threshold_valid_coverage`, are not counted. |
| **motif_type**                | `palindrome` (the motif is its own reverse complement), `non-palindrome`, `bipartite` (two parts separated by three or more `N`) or `ambiguous` (two or more stretches of `N`). |
| **motif_complement**          | The reverse complement motif, when it was found too: the other strand of the same recognition site, methylated at its own position. Empty otherwise. |
| **mod_position_complement**   | The zero-based position of the modified base in the complement motif. |
| **n_mod_complement**          | Methylated occurrences of the complement motif. |
| **n_nomod_complement**        | Unmethylated occurrences of the complement motif. |

The methylated fraction of a motif is `n_mod / (n_mod + n_nomod)`.

---

## Bin improvement

The bin contamination (`detect_contamination`) and inclusion (`include_contigs`) commands were removed in
version 1.2.0, with their output files (`bin_contamination.tsv`, `include_contigs.tsv`,
`new_contig_bin.tsv`). For bin refinement from methylation, see
[EpiRefine](https://github.com/SorenHeidelbach/epirefine).

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
