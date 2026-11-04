# Usage

---

## Motif discovery

The motif discovery process identifies motifs at bin levels. The command require the relation between bins and contigs to be specified. This can be done in three ways:
1) A TSV file specifying which bin contigs belong to (`-c` or `--contig_bin`).
2) A list of bin FASTA files with contig names as headers (`-f` or `--files`).
3) A directory containing bin FASTA files with contig names as headers (`-d` or `--directory`). The file extension of the bin FASTA files can be specified using the `--extension` flag.
The motif discovery process requires an assembly file and a modkit pileup file as input. The output will be written to the specified output folder, which will contain a `bin-motifs.tsv` file summarizing the identified motifs per bin. See [here](https://nanomotif.readthedocs.io/en/latest/output.html#motifs-tsv) for detailed output description.

The primary parameter to tune is the `--min_motif_score`, which determines the minimum score for a motif to be kept after identification. A lower value will result in more motifs being identified, but may also include more false positives. The default value is set to 1, which is a good starting point for most datasets. 
1) Minimum score of 0.5: Very sensitive, will identify many motifs, but also many false positives. In addition to small palindromic and non-palindromic motifs, most very degenerate or long bipartite motifs may be identified. In addition to motifs for systematic errors around other methylation types (e.g., C5mCGG for 6mA, giving rise to 6mA motif variant containing CCGG, e.g., 6mACCGG), motifs for systematic errors around non-methylated positions may be included, particularily 5mC motifs in high GC% organismns.
2) Minimum score of 1.0: Balanced sensitivity and specificity, suitable for most datasets. Will identify all small palindromic and non-palindromic motifs. Some very degenerate or long bipartite motifs may be missed. Motifs for systematic errors around other methylation type may be included (C5mCWGG for 6mA, giving rise to 6mA motif variant containing CCWGG, e.g. 6mANCCTGG).
3) Minimum score of 1.5: Very specific, will identify fewer motifs, but with high confidence. May miss some true motifs, particularly very degenerate or long bipartite motifs. Motifs for systematic errors around other methylation types or in high GC% organismns are unlikely to be included.


**QUICK START**
```shell
ASSEMBLY="path/to/assembly.fasta"
PILEUP="path/to/pileup.tsv"
BINS="path/to/bins"
OUT="path/to/output"
nanomotif motif_discovery $ASSEMBLY $PILEUP -d $BINS --out $OUT
```

```
usage: nanomotif motif_discovery (-c CONTIG_BIN | -f FILES [FILES ...] | -d DIRECTORY) [--extension EXTENSION] [--out OUT] [--methylation_threshold_low METHYLATION_THRESHOLD_LOW]
                                 [--methylation_threshold_high METHYLATION_THRESHOLD_HIGH] [--search_frame_size SEARCH_FRAME_SIZE] [--minimum_kl_divergence MINIMUM_KL_DIVERGENCE]
                                 [--min_motif_score MIN_MOTIF_SCORE] [--threshold_valid_coverage THRESHOLD_VALID_COVERAGE] [--min_motifs_bin MIN_MOTIFS_BIN] [-t THREADS] [-v]
                                 [--seed SEED] [-h]
                                 assembly pileup

positional arguments:
  assembly              path to the assembly file.
  pileup                path to the modkit pileup file.

contig bin arguments, use one of::
  -c CONTIG_BIN, --contig_bin CONTIG_BIN
                        TSV file specifying which bin contigs belong.
  -f FILES [FILES ...], --files FILES [FILES ...]
                        List of bin FASTA files with contig names as headers.
  -d DIRECTORY, --directory DIRECTORY
                        Directory containing bin FASTA files with contig names as headers.
  --extension EXTENSION
                        File extension of the bin FASTA files if using -d (DIRECTORY) argument. Default is '.fasta'.

Options:
  --out OUT             path to the output folder
  --methylation_threshold_low METHYLATION_THRESHOLD_LOW
                        A position is considered non-methylated if fraction of methylation is below this threshold. Default: 0.3
  --methylation_threshold_high METHYLATION_THRESHOLD_HIGH
                        A position is considered methylated if fraction of methylated reads is above this threshold. Default: 0.7
  --search_frame_size SEARCH_FRAME_SIZE
                        length of the sequnces sampled around confident methylation sites. Default: 40
  --minimum_kl_divergence MINIMUM_KL_DIVERGENCE
                        Minimum KL-divergence for a position to considered for expansion in motif search. Higher value means less exhaustive, but faster search. Default: 0.05
  --min_motif_score MIN_MOTIF_SCORE
                        Minimum score for a motif to be kept after identification. Default: 1
  --threshold_valid_coverage THRESHOLD_VALID_COVERAGE
                        Minimum valid base coverage (Nvalid_cov) for a position to be considered. Default: 5
  --min_motifs_bin MIN_MOTIFS_BIN
                        Minimum number of motif observations in a bin. Default: 20

general arguments:
  -t THREADS, --threads THREADS
                        Number of threads to use. Default is 1
  -v, --verbose         Increase output verbosity. (set logger to debug level)
  --seed SEED           Seed for random number generator. Default: 1
  -h, --help            show this help message and exit
```


---

## Bin improvement

The bin contamination (`detect_contamination`) and inclusion (`include_contigs`) commands have been removed and are no longer supported.

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
