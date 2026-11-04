# How it works

Nanomotif finds the sequence motifs that a genome's methyltransferases (MTases) methylate. It works one
bin (or genome) and one modification type at a time: 6mA (`a`), 5mC (`m`) and 4mC (`21839`).

The idea: an MTase methylates its target motif at (nearly) every copy of it in the genome, and nowhere
else. Positions with a high fraction of methylated reads therefore share a sequence context, and that
context is the motif. Nanomotif searches for the shortest motif that explains a group of methylated
positions while matching few unmethylated ones.

[1. Reading the pileup](#1-reading-the-pileup) · [2. Methylated and unmethylated positions](#2-methylated-and-unmethylated-positions) ·
[3. How a motif is scored](#3-how-a-motif-is-scored) · [4. Searching for a motif](#4-searching-for-a-motif) ·
[5. Finding all motifs of a bin](#5-finding-all-motifs-of-a-bin) · [6. Clean-up](#6-clean-up) ·
[7. Output](#7-output) · [Performance](#performance) · [Limitations](#limitations)

## 1. Reading the pileup

The modkit pileup has one row per position, strand and modification type, with the number of reads at
the position (`Nvalid_cov`) and the percentage of them called modified. Nanomotif keeps the rows of
contigs that belong to a bin and filters them:

1. **Coverage:** positions need at least `--threshold_valid_coverage` valid reads (default 5).
2. **Enough methylation:** for each contig and modification type, the rows are kept only if more than 50
   positions have a methylated fraction of at least `--methylation_threshold_high`. Contigs without real
   methylation of a type then add no noise to the search for that type.
3. **Neighbouring calls:** of methylated positions within 8 bp of each other on the same strand, only the
   most methylated one is kept. Methylation calls often bleed into the bases next to a methylated base;
   this keeps one call per methylation event.

The pileup is read once, in a parallel streaming pass, and split into one file per bin and modification
type, which the worker processes read. Plain and bgzipped (`.bed.gz`) pileups are read the same way.

## 2. Methylated and unmethylated positions

Each remaining position is then called:

- **methylated** if its fraction of modified reads is at least `--methylation_threshold_high` (default
  0.7);
- **unmethylated** if it is at most `--methylation_threshold_low` (default 0.3);
- **undecided** in between. Undecided positions are not used to score motifs.

A motif's occurrences in the bin are counted at these positions, on both strands: `n_mod` methylated and
`n_nomod` unmethylated occurrences. A motif is a pattern of fixed bases (IUPAC codes allowed) around the
methylated base, for example `GATC` with the methylated A at position 1.

## 3. How a motif is scored

How methylated a motif is gets described by a Beta-Bernoulli model: the share of its occurrences that are
methylated has a Beta distribution, starting from Beta(5, 5) and updated with `n_mod` and `n_nomod`.

Candidate motifs are compared with their less specific parent (the motif with one fixed base less). The
score of a child motif measures whether its sites are predicted better by its own model than the sites
the parent has and the child does not, weighted by how much more methylated the child is than the
parent. In practice:

- a child that separates methylated from unmethylated occurrences scores high (for example `GATC` versus
  `ATC` in a genome with a GATC methylase);
- a child that only selects a random subset of its parent's sites scores near 0.

`--min_motif_score` (default 1.5) is the score a motif must reach to be reported.

## 4. Searching for a motif

The search starts from the methylated base alone (for example `A` for 6mA) and adds fixed bases one at a
time, best first:

1. **Windows.** Around every methylated position, Nanomotif takes the sequence window of
   `--search_frame_size` + 1 bases (41 by default), with the methylated base in the middle. Windows on the
   minus strand are reverse complemented. As a background, it samples windows centred on the same base
   type at random positions: 1% of each contig's length, at least 50 per contig.
2. **Where to extend.** For the windows that match the current motif, Nanomotif computes the base
   frequencies at every position (a position frequency matrix) and compares them with the background,
   using the Kullback-Leibler divergence. The position where they differ most is extended next, if its
   divergence is at least `--minimum_kl_divergence` (default 0.05).
3. **Which bases.** Each base that occurs in more than 15% of the matching windows at that position, and
   more often than half its background frequency, becomes a child motif.
4. **Order.** Children go into a priority queue, ranked by how much they reduce the unmethylated
   occurrences while keeping the methylated ones. Motifs with fewer than 10 occurrences, or more than 25
   positions, are not extended further.
5. **Stop.** The search stops when 30 extensions in a row have not found a better-scoring motif. The
   best-scoring motif found is the candidate.

## 5. Finding all motifs of a bin

A bin often has several MTases, so the search is repeated:

1. The candidate is **pruned**: each of its fixed bases is removed in turn, and a base that adds little to
   the score is left out. This removes bases that only came along by chance.
2. If the pruned motif scores at least `--min_motif_score`, it is kept. Either way, the windows it
   matches are removed, so the next search looks for what is still unexplained.
3. The search restarts on the remaining windows. It ends after 25 candidates that did not reach the
   score, or when less than 0.1% of the methylated windows remain.
4. High-scoring motifs that were explored during the searches but are not related to a kept motif are
   added as well.

## 6. Clean-up

The kept motifs are then cleaned up, in this order (the motifs after each step are saved in
`precleanup-motifs/`):

1. **Noise.** Motifs with isolated fixed bases (a fixed base with only unfixed positions within 3 bases
   around it) are removed; they are rarely real recognition sequences.
2. **Merging.** Similar motifs (at most two positions different) are merged into one motif with IUPAC
   codes, for example `GAGCTC` and `GGGCCC` into `GRGCYC`, when the merged motif describes their sites as
   well as the separate motifs do.
3. **Sub-motifs.** When one motif is a more specific form of another, the more specific one is kept only
   if it scores clearly better (score above 0.5 against the other).
4. **Complements.** A motif and its reverse complement (the same recognition site on the other strand,
   for example the two halves of a non-palindromic or bipartite motif) are reported on one line.

Finally, motifs with fewer than `--min_motifs_bin` scored occurrences (`n_mod` + `n_nomod`, default 50)
are dropped, and the rest are written to `bin-motifs.tsv`.

## 7. Output

`bin-motifs.tsv` has one line per motif and bin: the motif in IUPAC notation, the position and type of
the modified base, `n_mod` and `n_nomod`, the motif type, and the complement motif if found. See
[Output](output.md).

## Performance

Motif discovery is built to run on large metagenomes:

- **Fast motif counting.** For every bin and modification type, Nanomotif builds an index of the
  sequence around each methylated and unmethylated position (one bitset per window position and base).
  Counting a candidate motif's occurrences is then a few bit operations instead of a scan of the
  contigs.
- **Streaming, bounded memory.** The pileup is read once, in parallel, and written to compressed files per
  bin and modification type. No process holds the whole pileup, and each worker reads only its own bin.
- **Compiled kernels.** The pileup reader and the inner loops run in compiled code (`nanomotif-rs`, Rust),
  installed with Nanomotif as prebuilt wheels on Linux, macOS and Windows. Without it, the same steps run
  in Python, with the same results.

For example, on a 12-bin oral mock community (4 GB pileup, 8 threads) motif discovery takes about 20 s
and 2.3 GB of memory; a 490-bin fecal metagenome with a 330 GB pileup took 21 minutes.

## Limitations

- **Thresholds.** Positions are called methylated or unmethylated with fixed thresholds, so motifs that
  are only partially methylated (for example 50% of reads at most sites) can be missed or reported with
  low counts.
- **Coverage.** Contigs with fewer than about 50 confidently methylated positions of a type are not used
  for that type, so short and low-coverage contigs contribute little.
- **Caller artefacts.** Systematic miscalls of the basecaller can look like motifs, for example 6mA calls
  around 5mC motifs, or 5mC noise in high-GC genomes; see [Known issues](known_issues.md).
- **Mixed bins.** A contaminated bin mixes the methylation of several genomes. Motifs of a minority genome
  may be missed or reported with many unmethylated occurrences.
