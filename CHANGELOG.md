# Changelog

## 1.2.2 (2026-11-04)

### Fixed
- Works with polars 2 (installed on Python 3.10+): contig-bin files are read by column position
  (polars 2 names header-less columns from column_0), and `read_ipc` is called without `memory_map`.
  1.2.1 failed its tests on Python 3.10-3.12 because of this. Same motifs with polars 1 and 2 (checked
  on two gold-standard isolates); tests pass on Python 3.9 (polars 1) and 3.12 (polars 2).

## 1.2.1 (2026-11-04)

First release published to PyPI with the compiled kernels; motif results are identical to 1.2.0.

### Changed
- epimetheus-py is no longer a dependency. It ships wheels only for Python 3.9 on Linux x86_64, so on
  Python 3.10+ installing it needed a Rust compiler (and failed in CI). Bgzipped pileups are now read by
  nanomotif-rs (streamed, decompressed once and processed in parallel, whole contigs per thread), or
  without it through pysam's tabix reader.
- nanomotif-rs buffers rows compactly (contig index instead of a string per row).
- Site index: windows cut from a strided view (no int64 index matrix) and bitsets built in one Rust pass.

### Measured
- Gold-standard monocultures (11 isolates, 66 motifs, 8 threads): same motifs as 1.2.0; total time
  145 -> 122 s, largest peak memory 2,088 -> 1,571 MB.
- ZymoOral (8 threads): plain pileup 23.4 s / 3.86 GB -> 21.0 s / 2.27 GB; bgzipped pileup 35.1 s /
  2.68 GB (1.2.0, epimetheus) -> 31.2 s / 2.57 GB.

## 1.2.0 (2026-11-04)

### Removed
- `detect_contamination` and `include_contigs` (bin contamination and inclusion) are removed and no
  longer supported, with their tests, test data, docs and the scikit-learn and hdbscan dependencies.

### Fixed (motif discovery)
- `--threshold_valid_coverage` was ignored (coverage > 5 was always used); positions now need at least
  this many valid reads, as documented. With the default 5 this keeps positions with exactly 5 reads,
  which 1.1.2 dropped; use 6 for the 1.1.2 behaviour.
- One contig without methylated positions far enough from its ends discarded motif discovery for its
  whole bin and mod type; such contigs are now skipped.
- `--methylation_threshold_low/high` were ignored by motif merging and by two pileup filters (fixed
  0.3/0.7).
- The noisy-motif filter kept all motifs when all were noisy; the isolated-base check never looked at a
  motif's last position.
- IUPAC ambiguity codes in the assembly crashed discovery (KeyError).
- Results could differ between runs with the same `--seed` (contig order of background sampling).
- Methylated positions exactly `padding` bases from a contig start were dropped.
- The missed-candidate filter now drops candidates related to a best candidate in either direction.
- When filtering left no pileup data the run crashed instead of writing an empty bin-motifs.tsv.
- Progress bar finished early in the plain-pileup path; debug FASTAs in temp/ held only the last contig;
  `nm.datasets.geobacillus_plasmids_assembly()` failed.

### Faster, less memory
- Candidate motifs are counted from a per-bin index of the scored positions (bitsets per window column
  and base) instead of regex scans of every contig.
- Pileup loading and filtering vectorised; the adjacency filter compares only positions above the
  threshold.
- Workers no longer receive every bin's data with every task (on ZymoOral 3.3 GB pickled per task); they
  read their bin's pileup from lz4 Arrow IPC files and their contigs from the FASTA index.
- Search windows stored as 4-bit base masks; no array copy per search step; direct KL divergence.
- Rust extension `nanomotif-rs` (`rust/`): a streaming, parallel pass that filters and partitions a
  plain pileup in bounded memory, and a fused window filter + PSSM kernel. Installed with nanomotif as
  a prebuilt wheel on Linux (x86_64, aarch64), macOS and Windows (x86_64), so no compiler is needed;
  on other platforms the pure-Python code runs. Wheels are built by .github/workflows/release.yml,
  which replaces python-publish.yml and uploads nanomotif-rs before nanomotif.

ZymoOral (12 bins, 4 GB plain pileup, 8 threads): 18 min 45 s and 17.3 GB peak memory with 1.1.2 ->
24 s and 3.5 GB with nanomotif_rs installed (39 s and 12.7 GB without it). Motifs are those of the
bug-fixed search (42 vs 40 with 1.1.2: all 40 found, plus two 6mA motifs in one bin).
ZymoFecal (490 bins, 330 GB plain pileup, 12 threads): 21 min, 16 GB peak system memory.
