# Usage

## Motif discovery

Motif discovery is meant for identifying motif at contig and bin level. It consist of three commands `find_motifs`, `score_motifs` & `bin_consensus`. We provide a wrapper command that executes these three commands togehter; `motif_discovery`.

We recommend always using `motif_discovery` unless there is a specific reason for using the seperate commands. 

QUICK START
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

COMING SOON


## MTase linking

COMING SOON