# Example

## Finding motifs in monocultures

When working with monoculture assemblies, you can directly use the contig FASTA file for motif discovery. Assuming your monoculture assembly is in `assembly.fasta`, you can run:

```shell
nanomotif motif_discovery \
    assembly.fasta \
    pileup.tsv \
    -f assembly.fasta \
    --out output_directory
```

Then Nanomotif will identify motifs in the full assembly, considering all contigs in the file as one bin. In the output file, the bin ID will be the same as the assembly file name.


## Finding motifs in metagenomes with binning

When working with metagenomic assemblies, you may want to identify motifs within specific bins. This requires a contig-bin relationship file to map contigs to their respective bins.

## Finding motifs in contigs

If you are interested in identifying motifs seperately for each contig, you can create a contig-bin relationship file where each contig is assigned to its own unique bin. For example, if your assembly contains contigs `contig_1`, `contig_2`, and `contig_3`, your contig-bin file would look like this:

```shell
contig_1   contig_1
contig_2   contig_2
contig_3   contig_3
```

Then provide this contig to contig relationship file to Nanomotif:

```shell
nanomotif motif_discovery \
    assembly.fasta \
    pileup.tsv \
    -c contig_bin.tsv \
    --out output_directory
```