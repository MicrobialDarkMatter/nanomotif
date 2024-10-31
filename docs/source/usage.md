# Usage

## Motif discovery

Motif discovery is meant for identifying motif at contig and bin level. It consist of three commands `find_motifs`, `score_motifs` & `bin_consensus`. We provide a wrapper command that executes these three commands togehter; `motif_discovery`.

We recommend always using `motif_discovery` unless there is a specific reason for using the seperate commands. 

**QUICK START**
```shell
ASSEMBLY="path/to/assembly.fasta"
PILEUP="path/to/pileup.tsv"
CONTIG_BIN="path/to/contig_bin.tsv"
OUT="path/to/output"
nanomotif motif_discovery $ASSEMBLY $PILEUP $CONTIG_BIN --out $OUT
```

```
usage: nanomotif motif_discovery [-h] [--out OUT] [-t THREADS] [-v] [--seed SEED] [--threshold_methylation_general THRESHOLD_METHYLATION_GENERAL]
                                 [--search_frame_size SEARCH_FRAME_SIZE] [--threshold_methylation_confident THRESHOLD_METHYLATION_CONFIDENT]
                                 [--threshold_valid_coverage THRESHOLD_VALID_COVERAGE] [--minimum_kl_divergence MINIMUM_KL_DIVERGENCE]
                                 assembly pileup bins

positional arguments:
  assembly              path to the assembly file.
  pileup                path to the modkit pileup file.
  bins                  tsv file specifying which bin contigs belong.

optional arguments:
  -h, --help            show this help message and exit
  --out OUT             path to the output folder
  -t THREADS, --threads THREADS
                        number of threads to use. Default is 1
  -v, --verbose         increase output verbosity. (set logger to debug level)
  --seed SEED           seed for random number generator. Default: 1
  --threshold_methylation_general THRESHOLD_METHYLATION_GENERAL
                        minimum fraction of reads that must be methylated at a position for the position to be methylated. These position are used for counting
                        number of methylated position of a motif. Default: 0.7
  --search_frame_size SEARCH_FRAME_SIZE
                        length of the sequnces sampled around confident methylatyion sites. Default: 40
  --threshold_methylation_confident THRESHOLD_METHYLATION_CONFIDENT
                        minimum fraction of reads that must be methylated at a position for the position to be considered confiently methylated. These position are
                        used to search for candidate motifs. Default: 0.8
  --threshold_valid_coverage THRESHOLD_VALID_COVERAGE
                        minimum valid base coverage for a position to be considered. Default: 5
  --minimum_kl_divergence MINIMUM_KL_DIVERGENCE
                        minimum KL-divergence for a position to considered for expansion in motif search. Higher value means less exhaustive, but faster search.
                        Default: 0.05
```


## Bin improvement

### Bin contamination
After motif identification it is possible to identify contamination in bins using the `bin-motifs.tsv`, `contig-bin.tsv` and `motif-scored.tsv` files.

**QUICK START**
```shell
MOTIFS_SCORED="path/to/nanomotif/motifs-scored.tsv"
BIN_MOTIFS="path/to/nanomotif/bin-motifs.tsv"
CONTIG_BIN="path/to/contig_bin.tsv"
OUT="path/to/output"
nanomotif detect_contamination --motifs_scored $MOTIFS_SCORED --bin_motifs $BIN_MOTIFS --contig_bins $CONTIG_BIN --out $OUT
```

This will generate a bin_contamination.tsv specifying the contigs, which is flagged as contamination.

If the `--write_bins` and the `--assembly_file` flags are specified new de-contaminated bins will be written to a bins folder.

If `detect_contamination` was run without the `--write_bins` flag, bins can be written as a post processing step if the `--contamination_file` is specified along with the `--write_bins` flag and the `--assembly_file` flag.


```
usage: nanomotif detect_contamination [-h] --motifs_scored MOTIFS_SCORED --bin_motifs BIN_MOTIFS --contig_bins
                                      CONTIG_BINS [-t THREADS] [--mean_methylation_cutoff MEAN_METHYLATION_CUTOFF]
                                      [--n_motif_contig_cutoff N_MOTIF_CONTIG_CUTOFF]
                                      [--n_motif_bin_cutoff N_MOTIF_BIN_CUTOFF]
                                      [--ambiguous_motif_percentage_cutoff AMBIGUOUS_MOTIF_PERCENTAGE_CUTOFF]
                                      [--write_bins] [--assembly_file ASSEMBLY_FILE] [--save_scores] --out OUT
                                      [--contamination_file CONTAMINATION_FILE]

optional arguments:
  -h, --help            show this help message and exit
  --motifs_scored MOTIFS_SCORED
                        Path to motifs-scored.tsv from nanomotif
  --bin_motifs BIN_MOTIFS
                        Path to bin-motifs.tsv file
  --contig_bins CONTIG_BINS
                        Path to bins.tsv file for contig bins
  -t THREADS, --threads THREADS
                        Number of threads to use for multiprocessing
  --mean_methylation_cutoff MEAN_METHYLATION_CUTOFF
                        Cutoff value for considering a motif as methylated
  --n_motif_contig_cutoff N_MOTIF_CONTIG_CUTOFF
                        Number of motifs that needs to be observed in a contig before it is considered valid for scoring
  --n_motif_bin_cutoff N_MOTIF_BIN_CUTOFF
                        Number of motifs that needs to be observed in a bin to be considered valid for scoring
  --ambiguous_motif_percentage_cutoff AMBIGUOUS_MOTIF_PERCENTAGE_CUTOFF
                        Percentage of ambiguous motifs defined as mean methylation between 0.05 and 0.40 in a bin.
                        Motifs with an ambiguous methylation percentage of more than this value are removed from
                        scoring. Default is 0.40
  --write_bins          If specified, new bins will be written to a bins folder. Requires --assembly_file to be
                        specified.
  --assembly_file ASSEMBLY_FILE
                        Path to assembly.fasta file
  --save_scores         If specified, the scores for each comparison will be saved to a scores folder in the output
                        directory
  --out OUT             Path to output directory
  --contamination_file CONTAMINATION_FILE
                        Path to an existing contamination file if bins should be outputtet as a post-processing step

```

### Include unbinned contigs
The `include_contigs` command assigns unbinned contigs in the assembly file to bins by comparing the methylation pattern of the contig to the bin consensus pattern. The contig must have a unique perfect match to the bin consensus pattern to be assigned to a bin. Additionally, the `include_contigs` assigns all the contigs in the `bin_contamination.tsv` file as unbinned. 

**QUICK START**
```shell
MOTIFS_SCORED="path/to/nanomotif/motifs-scored.tsv"
BIN_MOTIFS="path/to/nanomotif/bin-motifs.tsv"
CONTIG_BIN="path/to/contig_bin.tsv"
OUT="path/to/output"
nanomotif include_contigs --motifs_scored $MOTIFS_SCORED --bin_motifs $BIN_MOTIFS --contig_bins $CONTIG_BIN --run_detect_contamination --out $OUT
```

The output is two files: `include_contigs.tsv` and `new_contig_bin.tsv`. The `include_contigs.tsv` file is the contigs that were assigned based on the methylation pattern and the `new_contig_bin.tsv` is the updated contig_bin file.

If decontamination should not be performed, the `include_contigs` can be run without the `--run_detect_contamination` flag or without the `--contamination_file` flag.

#### save_scores
Save scores will create a csv file for each contig, which contains all the information for scoring the contig. Columns are:
- bin: The bin the contig is compared to
- motif_mod: The motif being compared `motif_string + "_" + mod_type + "_" + mod_position`
- n_mod: The number of modified motifs in bin
- n_nomod: The number of unmodified motifs in bin
- n_motifs_bin: The total number of motifs in bin
- mean_methylation: n_mod / n_motifs_bin
- mean_methylation_bin: If the motif is methylated, this number represent the mean methylation degree only for contigs where the mean mean_methylation was above 0.25
- std_methylation_bin: The std deviation for the above calculation
- n_contigs: The number of contigs used in above calculation
- methylation_binary: Indicates if the bin is methylated or not
- contig_bin: The bin where the contig originally was assigned
- bin_compare: `bin + "_" + contig` for where the contig was originally assigned
- mean: The mean methylation degree for the contig
- methylation_mean_theshold: The threshold for assigning the contig as methylated or not. This depends on the variation of the bin methylation for the given motif
- methylation_binary_compare: The binary methylation for the contig

## MTase-linker

This module links methylation motifs to their corresponding MTase and, when present, their entire RM system.

The MTase-Linker module has additional dependencies that are not automatically installed with Nanomotif. Therefore, before using this module, you must manually install these dependencies using the `MTase-linker install` command.
The `MTase-linker` module requires that conda is available on your system.

```
nanomotif MTase-linker install
```

This will create a folder named `ML_dependencies` in your current working directory, containing the required dependencies for the MTase-linker module. You can use the `--dependency_dir` flag to change the installation location of the `ML_dependencies` folder.

The installation requires conda to generate a few environments. 

When the additional dependencies are installed you can run the workflow using `MTase-linker run`

**QUICK START**
```shell
ASSEMBLY="path/to/assembly.fasta"
CONTIG_BIN="path/to/contig_bin.tsv"
BIN_MOTIFS="path/to/nanomotif/bin-motifs.tsv"
ML_DEPENDICIES="path/to/ML_dependencies"
OUT="path/to/output"
nanomotif MTase-linker run --asembly $ASSEMBLY --contig_bin $CONTIG_BIN --bin_motifs $BIN_MOTIFS -d $ML_DEPENDICIES --out $OUT
```

```
usage: nanomotif MTase-linker run [-h] [-t THREADS] [--forceall FORCEALL] [--dryrun DRYRUN] --assembly ASSEMBLY --contig_bin CONTIG_BIN --bin_motifs BIN_MOTIFS [-d DEPENDENCY_DIR] [-o OUT] [--identity IDENTITY] [--qcovs QCOVS]

optional arguments:
  -h, --help            show this help message and exit
  -t THREADS, --threads THREADS
                        Number of threads to use. Default is 1
  --forceall FORCEALL   Flag for snakemake. Forcerun workflow regardless of already created output (default = False)
  --dryrun DRYRUN       Flag for snakemake. Dry-run the workflow. Default is False
  --assembly ASSEMBLY   Path to assembly file.
  --contig_bin CONTIG_BIN
                        tsv file specifying which bin contigs belong.
  --bin_motifs BIN_MOTIFS
                        bin-motifs.tsv output from nanomotif.
  -d DEPENDENCY_DIR, --dependency_dir DEPENDENCY_DIR
                        Path to the ML_dependencies directory created during installation of the MTase-linker module. Default is cwd/ML_dependencies
  -o OUT, --out OUT     Path to output directory. Default is cwd
  --identity IDENTITY   Identity threshold for motif prediction. Default is 80
  --qcovs QCOVS         Query coverage for motif prediction. Default is 80
```
