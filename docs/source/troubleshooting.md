# Troubleshooting

If you run into a problem that is not listed here, please open an issue on
[GitHub](https://github.com/MicrobialDarkMatter/nanomotif/issues).

## `bin-motifs.tsv` has only a header

No motifs passed the filters. Common reasons:

- **Low coverage.** Positions need at least `--threshold_valid_coverage` reads (default 5), and a contig
  needs more than 50 confidently methylated positions of a modification type to be used for it. Check the
  coverage, or lower `--threshold_valid_coverage`.
- **Contig names do not match.** The contig names in the contig-bin file, the assembly and the pileup must
  be the same (the first word of each FASTA header). Contigs of the contig-bin file that are not in the
  assembly are removed with a warning in the log.
- **No methylation calls.** The reads must have been basecalled with a modification model, and the
  pileup made from a mapping that keeps the `MM`/`ML` tags (see [Required files](required_files.md)).

## "pileup rows of contig ... are not contiguous"

The pileup's rows of each contig must be together in the file, as `modkit pileup` writes them. If the
file was concatenated or re-sorted, sort it by contig first:

```shell
sort -k1,1 -k2,2n pileup.bed > pileup.sorted.bed
```

## "Tabix index for ... not found"

A bgzipped pileup without `nanomotif-rs` installed is read through its tabix index. Create it with
`tabix -p bed pileup.bed.gz`, or install `nanomotif-rs` (`pip install nanomotif-rs`), which needs no index.

## `detect_contamination` / `include_contigs`: invalid choice

These commands were removed in version 1.2.0. For bin refinement from methylation, see
[EpiRefine](https://github.com/SorenHeidelbach/epirefine).

## Results differ from an earlier version

Version 1.2.0 fixed several bugs that changed results, most visibly in low-coverage bins (see the
[changelog](changelog.md)). In particular, `--threshold_valid_coverage` is now used: positions need at
least 5 reads by default, where earlier versions required more than 5. Use
`--threshold_valid_coverage 6` to compare with older results.

## "Pileup cache ... was made with a different ..."

`--pileup_cache` reuses its file only for the same pileup (path, size and modification time), contig bins,
`--threshold_valid_coverage`, `--methylation_threshold_high` and Nanomotif version. The message names
what differed; the cache is then rebuilt, which takes as long as a run without it.

## Running out of memory

Memory use grows with the number of threads (bins processed at the same time) and the size of the largest
bins. Use fewer threads (`-t`) if a run is killed. Without `nanomotif-rs`, a plain pileup is loaded into
memory as a whole; install `nanomotif-rs` (default with pip on supported platforms) to stream it instead.
