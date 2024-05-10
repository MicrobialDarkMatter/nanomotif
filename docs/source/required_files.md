# Required File Preparation

To identify methylated motifs, the following files are required: 
- Assembly
- methylation pileup
- contig-bin relationship

## Assembly

Assembly file containing all the conting for evaluation.
Assembly from any assembler can be used. 
The assembly files should be in the standard fasta format. Nanomotif was developed using assemblies made with [metaFlye](https://github.com/fenderglass/Flye). 


## Methylation pileup
File containing the number of reads methylated at each position on a contig.
Gerated by mapping reads with methylation calls in the header to the assembly mentioned in [Assembly](####Assembly). Then, using ONT [modkit](https://github.com/nanoporetech/modkit/blob/master/book/src/advanced_usage.md#pileup) generate the methylation pileup. 


Code snippet for generating pileup file:
```shell
MODCALLS="path/to/reads/with/methylation/calls.bam"
ASSEMBLY="path/to/assembly.fa"
MAPPING="path/to/generated/mapping.bam"
PILEUP="path/to/generated/pileup.bed
samtools fastq -T MM,ML $MODCALLS | \
        minimap2 -ax map-ont -y $ASSEMBLY - | \
        samtools view -bS | \
        samtools sort -o $MAPPING
modkit pileup --only-tabs $MAPPING $PILEUP
```



Caling head on the generated pileup file should show a table similair to the one below:

|    |  |  |  |  |  |  |  | |  |  |  |  |  |  |  |  |  |
|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|
| contig_3 | 0     | 1   | m    | 133    | -      | 0      | 1    | 255,0,0   | 133  | 0.00  | 0    | 133    | 0       | 0          | 6    | 0          | 0     |
| contig_3 | 1     | 2   | a    | 174    | +      | 1      | 2    | 255,0,0   | 174  | 1.72  | 3    | 171    | 0       | 0          | 3    | 0          | 0     |
| contig_3 | 2     | 3   | a    | 172    | +      | 2      | 3    | 255,0,0   | 172  | 2.33  | 4    | 168    | 0       | 0          | 7    | 0          | 0     |
| contig_3 | 3     | 4   | a    | 178    | +      | 3      | 4    | 255,0,0   | 178  | 0.56  | 1    | 177    | 0       | 0          | 2    | 0          | 0     |
| contig_3 | 4     | 5   | a    | 177    | +      | 4      | 5    | 255,0,0   | 177  | 2.82  | 5    | 172    | 0       | 0          | 5    | 0          | 0     |
| contig_3 | 5     | 6   | a    | 179    | +      | 5      | 6    | 255,0,0   | 179  | 2.79  | 5    | 174    | 0       | 0          | 3    | 2          | 0     |
| contig_3 | 5     | 6   | m    | 1      | +      | 5      | 6    | 255,0,0   | 1    | 0.00  | 0    | 1      | 0       | 0          | 3    | 180        | 0     |
| contig_3 | 5     | 6   | a    | 1      | -      | 5      | 6    | 255,0,0   | 1    | 0.00  | 0    | 1      | 0       | 0          | 0    | 156        | 0     |
| contig_3 | 6     | 7   | m    | 183    | +      | 6      | 7    | 255,0,0   | 183  | 0.55  | 1    | 182    | 0       | 0          | 1    | 0          | 0     |
| contig_3 | 6     | 7   | a    | 4      | -      | 6      | 7    | 255,0,0   | 4    | 0.00  | 0    | 4      | 0       | 0          | 0    | 153        | 0     |



**Considerations**

- When demultiplexing, trimming of reads may result in errors downstream. We therefore recommend using untrimmed reads for mapping

- Running `modkit pileup` with default parameters results in modkit estimating the threshold for calling a methylation. This can result in very low methylation calling score threshold such as 0.6 or even lower. This is detrimental to Nanomotif motif identification, but can include a high degree of noise, and loss of some motifs. We generally recommend a `--filter-threshold` of 0.7. 

## Contig-bin
File describing which contigs belongs to which bins. It should be headerless, tab-separated file with contig id in the first column and bin in the second column. 


If the bins are outputted in a folder with one fasta file pr. bin, the contig-bin file can be generated using the following snippet:
```shell
BINS="/path/to/bins/fasta"    # Bins directory
BIN_EXT="fa"                  # Bins file extension
OUT="contig_bin.tsv"          # contig-bin output destination
grep ">" ${BINS}/*.${EXT} | \
        sed "s/.*\///" | \
        sed "s/.${EXT}:>/\t/" | \
        awk -F'\t' '{print $2 "\t" $1}' > $OUT
```
Caling head on the generated contig-bin file should show a table similair to the one below:

| | |
|-|-|
| contig_1 |        b1 |
| contig_2 |        b1 |
| contig_3 |        b1 |
| contig_4 |        b2 |
| contig_5 |        b2 |
| contig_6 |        b3 |
| contig_7 |        b3 |
| contig_8 |        b3 |
| contig_9 |        b3 |
| contig_10 |       b1 |