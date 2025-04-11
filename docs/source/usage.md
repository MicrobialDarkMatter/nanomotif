# Usage

---

## Motif discovery

The motif discovery process identifies motifs at both the contig and bin levels. It consists of three underlying commands: `find_motifs`, `score_motifs`, and `bin_consensus`. To simplify usage, a wrapper command, `motif_discovery`, runs all three steps in sequence. We recommend using `motif_discovery` unless you have a specific reason to run the individual commands separately. See [here](https://nanomotif.readthedocs.io/en/latest/output.html#motifs-tsv) for detailed output description.

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
                                 [--threshold_valid_coverage THRESHOLD_VALID_COVERAGE] [--minimum_kl_divergence MINIMUM_KL_DIVERGENCE] [--min_motifs_contig MIN_MOTIFS_CONTIG]
                                 [--read_level_methylation] [--min_motif_score MIN_MOTIF_SCORE] [--min_motifs_bin MIN_MOTIFS_BIN] [--save-motif-positions]
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
                        minimum fraction of reads that must be methylated at a position for the position to be methylated. These position are used for counting number of
                        methylated position of a motif. Default: 0.7
  --search_frame_size SEARCH_FRAME_SIZE
                        length of the sequnces sampled around confident methylatyion sites. Default: 40
  --threshold_methylation_confident THRESHOLD_METHYLATION_CONFIDENT
                        minimum fraction of reads that must be methylated at a position for the position to be considered confiently methylated. These position are used to
                        search for candidate motifs. Default: 0.8
  --threshold_valid_coverage THRESHOLD_VALID_COVERAGE
                        minimum valid base coverage for a position to be considered. Default: 5
  --minimum_kl_divergence MINIMUM_KL_DIVERGENCE
                        minimum KL-divergence for a position to considered for expansion in motif search. Higher value means less exhaustive, but faster search. Default: 0.05
  --min_motifs_contig MIN_MOTIFS_CONTIG
                        minimum number of times a motif has to have been oberserved in a contig. Default: 20
  --read_level_methylation
                        If specified, methylation is calculated on read level instead of contig level. This is slower but produces more stable motifs.
  --min_motif_score MIN_MOTIF_SCORE
                        minimum score for a motif to be kept after identification considered valid. Default: 0.2
  --min_motifs_bin MIN_MOTIFS_BIN
                        minimum number of times a motif has to have been oberserved in a bin. Default: 50
  --save-motif-positions
                        save motif positions in the output folder
```


---

## Bin improvement

### Bin contamination
After motif identification it is possible to identify contamination in bins using the `bin-motifs.tsv`, `assembly` and `pileup`. This will generate a `bin_contamination.tsv` specifying the contigs, which is flagged as contamination. 
The `detect_contamination` command detects putative contamination using  four clustering methods (agg, spectral, gmm, hdbscan), all of which have to flag the contig as a contaminant. 
See [here](https://nanomotif.readthedocs.io/en/latest/output.html#bin-contamination-tsv) for detailed output description.

**QUICK START**
```shell
ASSEMBLY="path/to/assembly.fasta"
PILEUP="path/to/pileup.bed"
BIN_MOTIFS="path/to/nanomotif/bin-motifs.tsv"
CONTIG_BIN="path/to/contig_bin.tsv"
OUT="path/to/output"
nanomotif detect_contamination --pileup $PILEUP --assembly $ASSEMBLY --bin_motifs $BIN_MOTIFS --contig_bins $CONTIG_BIN --out $OUT
```



```
usage: nanomotif detect_contamination [-h] --pileup PILEUP --assembly ASSEMBLY --bin_motifs BIN_MOTIFS --contig_bins CONTIG_BINS [-t THREADS]
                                      [--min_valid_read_coverage MIN_VALID_READ_COVERAGE] [--methylation_threshold METHYLATION_THRESHOLD] [--num_consensus NUM_CONSENSUS]
                                      [--force] [--write_bins] --out OUT [--contamination_file CONTAMINATION_FILE]

optional arguments:
  -h, --help            show this help message and exit
  -t THREADS, --threads THREADS
                        Number of threads to use for multiprocessing
  --min_valid_read_coverage MIN_VALID_READ_COVERAGE
                        Minimum read coverage for calculating methylation [used with methylation_util executable]
  --methylation_threshold METHYLATION_THRESHOLD
                        Filtering criteria for trusting contig methylation. It is the product of mean_read_coverage and N_motif_observation. Higher value means stricter
                        criteria. [default: 24]
  --num_consensus NUM_CONSENSUS
                        Number of models that has to agree for classifying as contaminant
  --force               Force override of motifs-scored-read-methylation.tsv. If not set existing file will be used.
  --write_bins          If specified, new bins will be written to a bins folder. Requires --assembly_file to be specified.
  --contamination_file CONTAMINATION_FILE
                        Path to an existing contamination file if bins should be outputtet as a post-processing step

Mandatory Arguments:
  --pileup PILEUP       Path to pileup.bed
  --assembly ASSEMBLY   Path to assembly file [fasta format required]
  --bin_motifs BIN_MOTIFS
                        Path to bin-motifs.tsv file
  --contig_bins CONTIG_BINS
                        Path to bins.tsv file for contig bins
  --out OUT             Path to output directory
```



### Include unbinned contigs

After motif identification, it is possible to assign unbinned contigs to bins using the `bin-motifs.tsv`, `assembly`, and `pileup`.
The `include_contigs` command assigns unbinned contigs in the assembly file to bins by training three supervised classifiers, random forest, linear discriminant analysis, and k-neighbors classifier. This will generate `include_contigs.tsv` specifying the contigs, assighned a new bin. See [here](https://nanomotif.readthedocs.io/en/latest/output.html#include-contigs-tsv) for detailed output description.

If decontamination should not be performed, the `include_contigs` can be run without the `--run_detect_contamination` flag or without the `--contamination_file` flag.

> Note: Assigning contigs based purely on methylation patterns can lead to errors as MAGs can share methylation patterns, which is especially problematic for unrecovered MAGs.

**QUICK START**
```shell
ASSEMBLY="path/to/assembly.fasta"
PILEUP="path/to/pileup.bed"
BIN_MOTIFS="path/to/nanomotif/bin-motifs.tsv"
CONTIG_BIN="path/to/contig_bin.tsv"
OUT="path/to/output"
nanomotif include_contigs --pileup $PILEUP --assembly $ASSEMBLY --bin_motifs $BIN_MOTIFS --contig_bins $CONTIG_BIN --run_detect_contamination --out $OUT
```



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
The installation requires conda to generate required environments. 
When the additional dependencies are installed you can run the workflow using `MTase-linker run`. See [here](https://nanomotif.readthedocs.io/en/latest/output.html#mtase-linker) for detailed output description.

**QUICK START**
```shell
ASSEMBLY="path/to/assembly.fasta"
CONTIG_BIN="path/to/contig_bin.tsv"
BIN_MOTIFS="path/to/nanomotif/bin-motifs.tsv"
ML_DEPENDICIES="path/to/ML_dependencies"
OUT="path/to/output"
nanomotif MTase-linker run --assembly $ASSEMBLY --contig_bin $CONTIG_BIN --bin_motifs $BIN_MOTIFS -d $ML_DEPENDICIES --out $OUT
```

```
usage: nanomotif MTase-linker run [-h] [-t THREADS] [--forceall FORCEALL] [--dryrun DRYRUN] --assembly ASSEMBLY --contig_bin CONTIG_BIN --bin_motifs BIN_MOTIFS [-d DEPENDENCY_DIR] [-o OUT] [--identity IDENTITY] [--qcovs QCOVS] [--minimum_motif_methylation MINIMUM_MOTIF_METHYLATION]

options:
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
  --minimum_motif_methylation MINIMUM_MOTIF_METHYLATION
                        Minimum fraction of motif occurrences in the bin that must be methylated for the motif to be considered. Default is 0.5
```
