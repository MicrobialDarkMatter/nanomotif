# Required File Preparation

Before running nanomotif for motif detection and analysis, ensure that you have prepared the necessary input files. These include a genome assembly file, a methylation pileup file, and a contig-bin relationship file.

---

## Assembly

The assembly file should contain all contigs in FASTA format. Each header should have a unique contig identifier. The sequence should only include standard nucleotide or IUPAC characters (either upper or lower case). Nanomotif has been primarily developed and tested using assemblies generated by [Flye](https://github.com/fenderglass/Flye).

**Requirements:**
- Format: FASTA
- Contains all contigs for evaluation
- Contig ID in the FASTA header
- IUPAC-compliant characters only

---

## Methylation Pileup

The methylation pileup file indicates how many mapped reads at each position show evidence of methylation. To generate this file:

1. Map reads (with methylation calls) to the assembly.
2. Use [modkit pileup](https://github.com/nanoporetech/modkit/blob/master/book/src/advanced_usage.md#pileup) to create the pileup.

**Example commands:**
```shell
MODCALLS="path/to/reads/with/methylation/calls.bam"
ASSEMBLY="path/to/assembly.fa"
MAPPING="path/to/generated/mapping.bam"
PILEUP="path/to/generated/pileup.bed"

samtools fastq -T MM,ML $MODCALLS | \
    minimap2 -ax map-ont -y $ASSEMBLY - | \
    samtools view -bS | \
    samtools sort -o $MAPPING

modkit pileup --only-tabs $MAPPING $PILEUP
```


Expected format: The pileup file is a tab-delimited table where each row represents a position on a contig, including information about methylation status.

Running "head" on the pileup file should produce a table similar to the one below:

| contig_3 | 0  | 1 | m | 133 | - | 0 | 1 | 255,0,0 | 133 | 0.00 | 0 | 133 | 0 | 0 | 6 | 0 | 0 |
|----------|----|---|---|-----|---|---|---|---------|-----|------|---|-----|---|---|---|---|---|
| contig_3 | 1  | 2 | a | 174 | + | 1 | 2 | 255,0,0 | 174 | 1.72 | 3 | 171 | 0 | 0 | 3 | 0 | 0 |
| contig_3 | 2  | 3 | a | 172 | + | 2 | 3 | 255,0,0 | 172 | 2.33 | 4 | 168 | 0 | 0 | 7 | 0 | 0 |
| contig_3 | 3  | 4 | a | 178 | + | 3 | 4 | 255,0,0 | 178 | 0.56 | 1 | 177 | 0 | 0 | 2 | 0 | 0 |
| contig_3 | 4  | 5 | a | 177 | + | 4 | 5 | 255,0,0 | 177 | 2.82 | 5 | 172 | 0 | 0 | 5 | 0 | 0 |
| contig_3 | 5  | 6 | a | 179 | + | 5 | 6 | 255,0,0 | 179 | 2.79 | 5 | 174 | 0 | 0 | 3 | 2 | 0 |
| contig_3 | 5  | 6 | m | 1   | + | 5 | 6 | 255,0,0 | 1   | 0.00 | 0 | 1   | 0 | 0 | 3 | 180 | 0 |
| contig_3 | 5  | 6 | a | 1   | - | 5 | 6 | 255,0,0 | 1   | 0.00 | 0 | 1   | 0 | 0 | 0 | 156 | 0 |
| contig_3 | 6  | 7 | m | 183 | + | 6 | 7 | 255,0,0 | 183 | 0.55 | 1 | 182 | 0 | 0 | 1 | 0 | 0 |
| contig_3 | 6  | 7 | a | 4   | - | 6 | 7 | 255,0,0 | 4   | 0.00 | 0 | 4   | 0 | 0 | 0 | 153 | 0 |

Considerations:  
- Use untrimmed reads for mapping to avoid downstream errors.  
- Running modkit pileup with default parameters may set a low methylation threshold and introduce noise. A filter-threshold of 0.7 is recommended to reduce noise and improve motif detection quality.

---

## Contig-Bin Relationship

For analyses that require binning, you need contig-bin relationship. This maps each contig to its corresponding bin, which is essential for binning-based motif discovery.

This informaiton can be passed in one of three ways:
1. **Contig-Bin File**: A file that explicitly maps contigs to bins
2. **Bin FASTA Files**: A directory of bin FASTA files, where each file corresponds to a bin
3. **List of bin FASTAs**: A list of bin FASTA files, where each file corresponds to a bin


**Createing a Contig-Bin File**

This file links each contig to its corresponding bin. It is a tab-separated file with two columns and no header:
- Column 1: Contig ID  
- Column 2: Bin ID

If you have a folder of bin FASTA files (one file per bin), you can generate the contig-bin file by extracting contig IDs and their associated bin filenames, then formatting this information into a two-column TSV.
```shell
BINS="/path/to/bins/fasta"    # Bins directory
BIN_EXT="fa"                  # Bins file extension
OUT="contig_bin.tsv"          # contig-bin output destination

grep ">" ${BINS}/*.${BIN_EXT} | \
        sed "s/.*\///" | \
        sed "s/.${BIN_EXT}:>/\t/" | \
        awk -F'\t' '{print $2 "\t" $1}' > $OUT
```


Example output:

| contig_1  | bin1 |
|-----------|------|
| contig_2  | bin1 |
| contig_3  | bin1 |
| contig_4  | bin2 |
| contig_5  | bin2 |
| contig_6  | bin3 |
| contig_7  | bin3 |
| contig_8  | bin3 |
| contig_9  | bin3 |
| contig_10 | bin1 |

---
