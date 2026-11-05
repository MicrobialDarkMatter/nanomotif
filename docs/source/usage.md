# Usage

Nanomotif has three commands:

- `motif_discovery`: finds methylation motifs per bin (or genome).
- `MTase-linker`: links the motifs to methyltransferase genes ([below](#mtase-linker)).
- `check_installation`: runs a small test.

The bin contamination (`detect_contamination`) and inclusion (`include_contigs`) commands were removed
in version 1.2.0; for bin refinement from methylation, see
[EpiRefine](https://github.com/SorenHeidelbach/epirefine).

---

## Motif discovery

Motif discovery needs an assembly, a modkit pileup and the bin of each contig (see
[Required files](required_files.md)). The bins can be given in one of three ways:

1. a TSV file with the bin of each contig (`-c` / `--contig_bin`);
2. a list of bin FASTA files, one file per bin (`-f` / `--files`);
3. a directory of bin FASTA files (`-d` / `--directory`), with `--extension` for their file extension.

With FASTA files, the bin name is the file name without its extension.

**Quick start**
```shell
ASSEMBLY="path/to/assembly.fasta"
PILEUP="path/to/pileup.bed"
CONTIG_BIN="path/to/contig_bin.tsv"
OUT="path/to/output"
nanomotif motif_discovery $ASSEMBLY $PILEUP -c $CONTIG_BIN -t 8 --out $OUT
```

The motifs are written to `$OUT/bin-motifs.tsv`; see [Output](output.md). How the search works is
described in [How it works](how_it_works.md).

### Parameters

The defaults suit most datasets. The ones you may want to change:

- **`--min_motif_score`** (default 1.5): the score a motif needs to be reported. It is the main way to
  trade sensitivity for specificity:
  1. **0.5:** very sensitive. Finds small palindromic and non-palindromic motifs and most degenerate or
     long bipartite motifs, but also many false positives: motifs from systematic miscalls around other
     modification types (for example 6mA calls around C5mCGG, giving 6mA motifs containing CCGG), and
     around unmethylated positions, particularly 5mC motifs in high-GC genomes.
  2. **1.0:** balanced. Finds all small palindromic and non-palindromic motifs; may miss some very
     degenerate or long bipartite motifs, and may include motifs from miscalls around other modification
     types (for example 6mA motifs containing CCWGG).
  3. **1.5 (default):** specific. Fewer motifs, with high confidence. Very degenerate or long bipartite
     motifs may be missed; motifs from miscalls are unlikely.
- **`--methylation_threshold_high` / `--methylation_threshold_low`** (defaults 0.7 / 0.3): a position is
  methylated if at least 70% of its reads are called modified, unmethylated if at most 30% are.
  Positions in between are not used.
- **`--threshold_valid_coverage`** (default 5): positions need at least this many valid reads. Lower it
  for low-coverage data, at the cost of noisier calls. (Before version 1.2.0 this option was ignored and
  positions needed more than 5 reads; use 6 to reproduce that.)
- **`--min_motifs_bin`** (default 50): a motif is only reported if it has at least this many scored
  occurrences (`n_mod` + `n_nomod`) in the bin.
- **`--pileup_cache FILE`** (off by default): keep the filtered pileup in `FILE` and reuse it in later
  runs. Useful when you run motif discovery several times on the same data, for example with different
  search settings. The first run writes the file (one compressed file, about a tenth of the pileup's
  size, at little extra time); later runs with the same pileup, contig bins, `--threshold_valid_coverage` and
  `--methylation_threshold_high` read it instead of the pileup. If any of these differ, or the pileup
  file has changed, the cache is rebuilt and the log says why. Other options (for example
  `--min_motif_score` or `--search_frame_size`) can change between runs that share a cache.
- **`-t` / `--threads`**: bins are processed in parallel, so more threads help most with many bins. Memory
  use grows with the number of threads and the size of the largest bins.
- **`--seed`**: runs with the same seed and inputs give the same result.

```
usage: nanomotif motif_discovery (-c CONTIG_BIN | -f FILES [FILES ...] | -d DIRECTORY)
                                 [--extension EXTENSION] [--out OUT]
                                 [--methylation_threshold_low METHYLATION_THRESHOLD_LOW]
                                 [--methylation_threshold_high METHYLATION_THRESHOLD_HIGH]
                                 [--search_frame_size SEARCH_FRAME_SIZE]
                                 [--minimum_kl_divergence MINIMUM_KL_DIVERGENCE]
                                 [--min_motif_score MIN_MOTIF_SCORE]
                                 [--threshold_valid_coverage THRESHOLD_VALID_COVERAGE]
                                 [--pileup_cache FILE] [--min_motifs_bin MIN_MOTIFS_BIN] [-t THREADS]
                                 [-v] [--seed SEED] [-h]
                                 assembly pileup

positional arguments:
  assembly              Path to the assembly (FASTA).
  pileup                Path to the modkit pileup (plain or bgzipped).

contig bin arguments, use one of::
  -c CONTIG_BIN, --contig_bin CONTIG_BIN
                        Two-column TSV file (contig, bin) without header: the bin of each contig.
  -f FILES [FILES ...], --files FILES [FILES ...]
                        List of bin FASTA files with contig names as headers.
  -d DIRECTORY, --directory DIRECTORY
                        Directory containing bin FASTA files with contig names as headers.
  --extension EXTENSION
                        File extension of the bin FASTA files if using -d (DIRECTORY) argument. Default is
                        '.fasta'.

Options:
  --out OUT             Path to the output folder. Default: nanomotif
  --methylation_threshold_low METHYLATION_THRESHOLD_LOW
                        A position is considered unmethylated if its fraction of modified reads is at most
                        this. Default: 0.3
  --methylation_threshold_high METHYLATION_THRESHOLD_HIGH
                        A position is considered methylated if its fraction of modified reads is at least
                        this. Default: 0.7
  --search_frame_size SEARCH_FRAME_SIZE
                        Length of the sequences sampled around confidently methylated positions. Default: 40
  --minimum_kl_divergence MINIMUM_KL_DIVERGENCE
                        Minimum KL divergence for a position to be considered for expansion in the motif
                        search. A higher value means a less exhaustive, but faster search. Default: 0.05
  --min_motif_score MIN_MOTIF_SCORE
                        Minimum score for a motif to be kept after identification. Default: 1.5
  --threshold_valid_coverage THRESHOLD_VALID_COVERAGE
                        Minimum number of valid reads (Nvalid_cov) for a position to be used. Default: 5
  --pileup_cache FILE   Keep the filtered pileup in this file (one compressed Arrow file) and reuse it in
                        later runs with the same pileup, contig bins, --threshold_valid_coverage and
                        --methylation_threshold_high, skipping the pileup pass. Off by default (the
                        filtered pileup is temporary).
  --min_motifs_bin MIN_MOTIFS_BIN
                        Minimum number of scored motif occurrences (n_mod + n_nomod) in a bin for a motif to
                        be reported. Default: 50

general arguments:
  -t THREADS, --threads THREADS
                        Number of threads to use. Default is 1
  -v, --verbose         Increase output verbosity. (set logger to debug level)
  --seed SEED           Seed for random number generator. Default: 1
  -h, --help            show this help message and exit
```

---

## MTase-linker

This module links methylation motifs to the methyltransferase (MTase) genes that produce them and, when
present, to their entire restriction-modification (RM) system.

MTase-linker has additional dependencies that are not installed with Nanomotif; install them first with
`nanomotif MTase-linker install` (needs conda; see [Installation](installation.md#mtase-linker-installation)).
It then runs as a Snakemake workflow on the assembly, the contig-bin file and Nanomotif's
`bin-motifs.tsv`. See [Output](output.md#mtase-linker) for the output files.

**Quick start**
```shell
ASSEMBLY="path/to/assembly.fasta"
CONTIG_BIN="path/to/contig_bin.tsv"
BIN_MOTIFS="path/to/nanomotif/bin-motifs.tsv"
ML_DEPENDENCIES="path/to/ML_dependencies"
OUT="path/to/output"
nanomotif MTase-linker run -t 8 --assembly $ASSEMBLY --contig_bin $CONTIG_BIN --bin_motifs $BIN_MOTIFS \
    -d $ML_DEPENDENCIES --out $OUT
```

Leave out `--forceall` and `--dryrun` unless you want them: any value given to them, including `False`,
turns them on (see [Known issues](known_issues.md)).

```
usage: nanomotif MTase-linker run [-h] [-t THREADS] [--forceall FORCEALL] [--dryrun DRYRUN] --assembly
                                  ASSEMBLY --contig_bin CONTIG_BIN --bin_motifs BIN_MOTIFS
                                  [-d DEPENDENCY_DIR] [-o OUT] [--identity IDENTITY] [--qcovs QCOVS]
                                  [--minimum_motif_methylation MINIMUM_MOTIF_METHYLATION]

optional arguments:
  -h, --help            show this help message and exit
  -t THREADS, --threads THREADS
                        Number of threads to use. Default is 1
  --forceall FORCEALL   Flag for snakemake. Forcerun workflow regardless of already created output (default
                        = False)
  --dryrun DRYRUN       Flag for snakemake. Dry-run the workflow. Default is False
  --assembly ASSEMBLY   Path to assembly file.
  --contig_bin CONTIG_BIN
                        tsv file specifying which bin contigs belong.
  --bin_motifs BIN_MOTIFS
                        bin-motifs.tsv output from nanomotif.
  -d DEPENDENCY_DIR, --dependency_dir DEPENDENCY_DIR
                        Path to the ML_dependencies directory created during installation of the MTase-
                        linker module. Default is cwd/ML_dependencies
  -o OUT, --out OUT     Path to output directory. Default is cwd
  --identity IDENTITY   Identity threshold for motif prediction. Default is 80
  --qcovs QCOVS         Query coverage for motif prediction. Default is 80
  --minimum_motif_methylation MINIMUM_MOTIF_METHYLATION
                        Minimum fraction of motif occurrences in the bin that must be methylated for the
                        motif to be considered. Default is 0.5
```
