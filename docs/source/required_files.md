# Generating Required files


## Required files

To identify methylated motifs, the following files are required: 
- Assembly
- methylation pileup
- Contig-bin relationship*

**Only for bin-consensus (and by extension, complete-workflow)*

#### Assembly
Assembly from any assembler can be used. Nanomotif was developed using assemblies made with [metaFlye](https://github.com/fenderglass/Flye)

#### Methylation pileup
Generated using ONT [modkit](https://github.com/nanoporetech/modkit/blob/master/book/src/advanced_usage.md#pileup). Assembly and basecalls with methylations are required.

*OBS: when demultiplexing, trimming of reads may result in errors downstream. We therefore recommend using untrimmed reads for mapping*
```shell
samtools fastq -T MM,ML {DORADO_BASECALLS.sam} | \
        minimap2 -ax map-ont -y {ASSEMBLY.fasta} - | \
        samtools view -bS | \
        samtools sort -o {ALIGNMENT.bam}
modkit pileup --only-tabs {ALIGNMENT.bam} {OUT.bed}
```

#### Contig-bin
File describing which contigs belongs to which bins. It should be headerless, tab-separated file with contig id in the first column and bin in the second column. If the bins are outputted in a folder with one fasta file pr. bin, the contig-bin file can be generated using the following snippet:
```
BINS=/path/to/bins    # Bin directory
EXT="fa"              # Filename extension
grep ">" ${BINS}/*.${EXT} | \
        sed "s/.*\///" | \
        sed "s/.${EXT}:>/\t/" | \
        awk -F'\t' '{print $2 "\t" $1}' > contig_bin.tsv
```
