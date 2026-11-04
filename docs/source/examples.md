# Examples

## Finding motifs in a monoculture

For a monoculture, give the assembly itself as the only "bin" file, so all contigs are treated as one
genome:

```shell
nanomotif motif_discovery \
    assembly.fasta \
    pileup.bed \
    -f assembly.fasta \
    -t 8 \
    --out nanomotif_out
```

The bin is named after the file without its extension (here `assembly`), and `nanomotif_out/bin-motifs.tsv`
lists the genome's motifs.

## Finding motifs in a metagenome

With bins from a binner, give the bin of each contig as a TSV file (see
[Required files](required_files.md#contig-bin-relationship)):

```shell
nanomotif motif_discovery \
    assembly.fasta \
    pileup.bed \
    -c contig_bin.tsv \
    -t 16 \
    --out nanomotif_out
```

Or point to the folder of bin FASTA files:

```shell
nanomotif motif_discovery assembly.fasta pileup.bed -d bins/ --extension .fa -t 16 --out nanomotif_out
```

Bins are processed in parallel; with many bins, more threads give the largest speed-up.

## Using a bgzipped pileup

A bgzipped pileup is read the same way and saves disk space:

```shell
bgzip -@ 8 pileup.bed
nanomotif motif_discovery assembly.fasta pileup.bed.gz -c contig_bin.tsv -t 16 --out nanomotif_out
```

If `nanomotif-rs` is not installed (see [Installation](installation.md)), index the file first with
`tabix -p bed pileup.bed.gz`.

## Finding motifs in each contig

To find motifs separately for each contig, give every contig its own bin. For an assembly with contigs
`contig_1`, `contig_2` and `contig_3`:

```shell
contig_1   contig_1
contig_2   contig_2
contig_3   contig_3
```

A file like this can be made from the assembly with:

```shell
grep ">" assembly.fasta | sed 's/>//; s/ .*//' | awk '{print $1 "\t" $1}' > contig_contig.tsv
```

Then:

```shell
nanomotif motif_discovery assembly.fasta pileup.bed -c contig_contig.tsv -t 8 --out nanomotif_out
```

Short contigs often have too few methylated positions for motif discovery; their motifs are then
missing rather than wrong.

## Linking motifs to MTase genes

After motif discovery, MTase-linker links the motifs to methyltransferase genes (install its
dependencies first, see [Installation](installation.md#mtase-linker-installation)):

```shell
nanomotif MTase-linker run -t 8 \
    --assembly assembly.fasta \
    --contig_bin contig_bin.tsv \
    --bin_motifs nanomotif_out/bin-motifs.tsv \
    -d ML_dependencies \
    -o mtase_linker_out
```
