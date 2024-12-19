# Usage

---

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


---

## Bin improvement

### Bin contamination
After motif identification it is possible to identify contamination in bins using the `bin-motifs.tsv`, `assembly` and `pileup`.

**QUICK START**
```shell
ASSEMBLY="path/to/assembly.fasta"
PILEUP="path/to/pileup.bed"
BIN_MOTIFS="path/to/nanomotif/bin-motifs.tsv"
CONTIG_BIN="path/to/contig_bin.tsv"
OUT="path/to/output"
nanomotif detect_contamination --pileup $PILEUP --assembly $ASSEMBLY --bin_motifs $BIN_MOTIFS --contig_bins $CONTIG_BIN --out $OUT
```

This will generate a bin_contamination.tsv specifying the contigs, which is flagged as contamination.

If the `--write_bins` flag is specified new de-contaminated bins will be written to a bins folder.

```
usage: nanomotif detect_contamination [-h] --pileup PILEUP --assembly ASSEMBLY
                                      --bin_motifs BIN_MOTIFS --contig_bins
                                      CONTIG_BINS [-t THREADS]
                                      [--min_valid_read_coverage MIN_VALID_READ_COVERAGE]
                                      [--methylation_threshold METHYLATION_THRESHOLD]
                                      [--num_consensus NUM_CONSENSUS]
                                      [--force] [--write_bins] --out OUT
                                      [--contamination_file CONTAMINATION_FILE]

options:
  -h, --help            show this help message and exit
  -t THREADS, --threads THREADS
                        Number of threads to use for multiprocessing
  --min_valid_read_coverage MIN_VALID_READ_COVERAGE
                        Minimum read coverage for calculating methylation
                        [used with methylation_util executable]
  --methylation_threshold METHYLATION_THRESHOLD
                        Filtering criteria for trusting contig methylation. It
                        is the product of mean_read_coverage and
                        N_motif_observation. Higher value means stricter
                        criteria. [default: 24]
  --num_consensus NUM_CONSENSUS
                        Number of models that has to agree for classifying as
                        contaminant
  --force               Force override of motifs-scored-read-methylation.tsv.
                        If not set existing file will be used.
  --write_bins          If specified, new bins will be written to a bins
                        folder. Requires --assembly_file to be specified.
  --contamination_file CONTAMINATION_FILE
                        Path to an existing contamination file if bins should
                        be outputtet as a post-processing step

Mandatory Arguments:
  --pileup PILEUP       Path to pileup.bed
  --assembly ASSEMBLY   Path to assembly file [fasta format required]
  --bin_motifs BIN_MOTIFS
                        Path to bin-motifs.tsv file
  --contig_bins CONTIG_BINS
                        Path to bins.tsv file for contig bins
  --out OUT             Path to output directory```
```

The output is a `bin_contamination.tsv` file. The each contaminant will have 4 rows, one for each clustering algorithm, along with the cluster stats.

### Include unbinned contigs
The `include_contigs` command assigns unbinned contigs in the assembly file to bins by training three supervised classifiers, random forest, linear discriminant analysis, and k-neighbors classifier.
In case all three classifiers assigns a unbinned contig to the same bin with a join mean probability above 0.80, the contig is assigned. This is called a `high_confidence` assignment.


**QUICK START**
```shell
ASSEMBLY="path/to/assembly.fasta"
PILEUP="path/to/pileup.bed"
BIN_MOTIFS="path/to/nanomotif/bin-motifs.tsv"
CONTIG_BIN="path/to/contig_bin.tsv"
OUT="path/to/output"
nanomotif include_contigs --pileup $PILEUP --assembly $ASSEMBLY --bin_motifs $BIN_MOTIFS --contig_bins $CONTIG_BIN --run_detect_contamination --out $OUT
```

The output file is a `include_contigs.tsv`, which will show the classifier assignment stats. Besides the aforementioned `high_confidence` assignment there is also a medium and low confidence assigment.
`medium_confidence` assignments are contigs where all three classifiers agree but the join probability is below 0.8. `low_confidence` is when only two classifiers agree.

`high_confidence` assignments are outputted in a `new_contig_bin.tsv` file. 

If decontamination should not be performed, the `include_contigs` can be run without the `--run_detect_contamination` flag or without the `--contamination_file` flag.

> Note: Assigning contigs based purely on methylation patterns can lead to errors as MAGs can share methylation patterns, which is especially problematic for unrecovered MAGs.

```
usage: nanomotif include_contigs [-h] --pileup PILEUP --assembly ASSEMBLY
                                 --bin_motifs BIN_MOTIFS --contig_bins
                                 CONTIG_BINS [-t THREADS]
                                 [--min_valid_read_coverage MIN_VALID_READ_COVERAGE]
                                 [--methylation_threshold METHYLATION_THRESHOLD]
                                 [--num_consensus NUM_CONSENSUS] [--force]
                                 [--write_bins] --out OUT
                                 [--mean_model_confidence MEAN_MODEL_CONFIDENCE]
                                 [--contamination_file CONTAMINATION_FILE | --run_detect_contamination]

options:
  -h, --help            show this help message and exit
  -t THREADS, --threads THREADS
                        Number of threads to use for multiprocessing
  --min_valid_read_coverage MIN_VALID_READ_COVERAGE
                        Minimum read coverage for calculating methylation
                        [used with methylation_util executable]
  --methylation_threshold METHYLATION_THRESHOLD
                        Filtering criteria for trusting contig methylation. It
                        is the product of mean_read_coverage and
                        N_motif_observation. Higher value means stricter
                        criteria. [default: 24]
  --num_consensus NUM_CONSENSUS
                        Number of models that has to agree for classifying as
                        contaminant
  --force               Force override of motifs-scored-read-methylation.tsv.
                        If not set existing file will be used.
  --write_bins          If specified, new bins will be written to a bins
                        folder. Requires --assembly_file to be specified.
  --mean_model_confidence MEAN_MODEL_CONFIDENCE
                        Mean probability between models for including contig.
                        Contigs above this value will be included. [default:
                        0.8]
  --contamination_file CONTAMINATION_FILE
                        Path to an existing contamination file to include in
                        the analysis
  --run_detect_contamination
                        Indicate that the detect_contamination workflow should
                        be run first

Mandatory Arguments:
  --pileup PILEUP       Path to pileup.bed
  --assembly ASSEMBLY   Path to assembly file [fasta format required]
  --bin_motifs BIN_MOTIFS
                        Path to bin-motifs.tsv file
  --contig_bins CONTIG_BINS
                        Path to bins.tsv file for contig bins
  --out OUT             Path to output directory```
```

---

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
