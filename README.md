# Nanomotif

Nanomotif is a Python package that provides functionality for identifying methylation motifs on references using Nanopore sequencing.

## Installation

#### Local Environment

To install Nanomotif in a local Python environment:

```shell
python3 -m venv nanomotif
source nanomotif/bin/activate
pip install nanomotif
```

#### Conda Environment

If you prefer using Conda for managing your Python environments, you can create a new environment and install Nanomotif as follows:

```shell
conda create -n nanomotif python=3.9
conda activate nanomotif
python -m pip install nanomotif
```
#### Check installation
Once installed, the installation can be checked by running:
```shell
nanomotif check-installation
```
This runs a test run on a small dataset, ensuring everything works.

## Required files

To identify methylated motifs, the follwoing files are required: 
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

## Usage

```
usage: nanomotif [-h] [--version] {find-motifs,score-motifs,bin-consensus,complete-workflow,check-installation} ...

Motif identification and utilisation commands

positional arguments:
                        Command descriptions
    find-motifs         identifies motifs in contigs
    score-motifs        generate feature complete output (all identified motifs have an entry for all contigs)
    bin-consensus       generate consensus set of motif for each bin
    complete-workflow   run find-motifs, score-motifs and bin-consensus
```
#### Monoculture
If you are interested in finding methylated motifs in a monoculture sample, we recomment just running `find-motifs`. 
```
nanomotif find-motifs ASSEMBLY.fasta PILEUP.bed
```

#### Metagenomic sample
If you have metagenomic sample with multiple organims, we recommend running `complete-workflow`
```
nanomotif complete-workflow ASSEMBLY.fasta PILEUP.bed CONTIG_BIN.tsv
```

#### Help
Help page of `complete-workflow` contains description of all arguments
```
usage: nanomotif complete-workflow [-h] [--out OUT] [-t THREADS] [-v] [--seed SEED] [--threshold_methylation_general THRESHOLD_METHYLATION_GENERAL]
                                   [--search_frame_size SEARCH_FRAME_SIZE] [--threshold_methylation_confident THRESHOLD_METHYLATION_CONFIDENT]
                                   [--threshold_valid_coverage THRESHOLD_VALID_COVERAGE] [--minimum_kl_divergence MINIMUM_KL_DIVERGENCE]
                                   assembly pileup bins
positional arguments:
  assembly              path to the assembly file.
  pileup                path to the modkit pileup file.
  bins                  tsv file specifying which bin contigs belong.

optional arguments:
  -h, --help            show this help message and exit
  --out OUT             path to the output folder
  -t THREADS, --threads THREADS
                        number of threads to use. Default is 1
  -v, --verbose         increase output verbosity. (set logger to debug level)
  --seed SEED           seed for random number generator. Default: 1
  --threshold_methylation_general THRESHOLD_METHYLATION_GENERAL
                        minimum fraction of reads that must be methylated at a position for the position to be methylated. These position are used for counting number of
                        methylated position of a motif. Default: 0.6
  --search_frame_size SEARCH_FRAME_SIZE
                        length of the sequnces sampled around confident methylatyion sites. Default: 40
  --threshold_methylation_confident THRESHOLD_METHYLATION_CONFIDENT
                        minimum fraction of reads that must be methylated at a position for the position to be considered confiently methylated. These position are used to
                        search for candidate motifs. Default: 0.8
  --threshold_valid_coverage THRESHOLD_VALID_COVERAGE
                        minimum valid base coverage for a position to be considered. Default: 5
  --minimum_kl_divergence MINIMUM_KL_DIVERGENCE
                        minimum KL-divergence for a position to considered for expansion in motif search. Higher value means less exhaustive, but faster search. Default: 0.2
```
## Output description
### find-motifs and score-motifs

`find-motifs` outputs results to `motifs.tsv` and `score-motifs` outputs result to `motifs-scored.tsv`.

`find-motifs`and `score-motifs` follow the same output format, described in the table below:

| **Column**       | **Description**                                                                                       |
|------------------|-------------------------------------------------------------------------------------------------------|
| **contig**       | contig in which the motif was found                                                                |
| **motif**        | sequence of the detected motif in IUPAC format |
| **mod_position** | position within the motif where the methylation is located. 0-based index.                            |
| **mod_type**     | the type of modification [a (6mA) or m (5mC)]                                                                     |
| **n_mod**        | number of motif positions that are methylated in the contig                                                       |
| **n_nomod**         | number of motif positions that are not methylated in the contig                                                        |
| **motifs_type** | type of motif the sequence (palindrome, non-palindrome, bipartite or ambiguous) |
| **motif_complement**       | Sequence of the complement motif if present in IUPAC format.                                            |
| **mod_position_complement**           | Position within the complement motif where the methylation is located. 0-based index.                 |
| **n_mod_complement**       | Number of motif positions that are methylated in the contig.                               |
| **n_nomod_complement**     | Number of motif positions that are not methylated in the contig.                           |

Running `find-motifs` generates pre-cleanup folder, whihc contains motif that got removed in the postprocessing steps. The name of the file indicate which postprocessing steps have been run on the motifs.

### bin-consensus 
`bin-consensus` outputs results to `bin-motifs.tsv`
The format is almost identical the the output of find-motifs, except everything is aggregated to bin level and the contig column is replaced by a bin column
| **Column**       | **Description**                                                                                       |
|------------------|-------------------------------------------------------------------------------------------------------|
| **bin**          | bin to which the motif belong                                                                         |
| **motif**        | sequence of the detected motif in IUPAC format                                                        |
| **mod_position** | position within the motif where the methylation is located. 0-based index.                            |
| **mod_type**     | the type of modification [a (6mA) or m (5mC)]                                                                     |
| **n_mod**        | number of motif positions that are methylated in all contigs in the bin                                                       |
| **n_nomod**      | number of motif positions that are not methylated in all contigs in the bin                                                       |
| **motifs_type** | type of motif the sequence (palindrome, non-palindrome, bipartite or ambiguous) |
| **motif_complement**       | Sequence of the complement motif if present.                                            |
| **mod_position_complement**           | Position within the complement motif where the methylation is located. 0-based index.                 |
| **n_mod_complement**       | Number of motif positions that are methylated in all contigs in the bin                               |
| **n_nomod_complement**     | Number of motif positions that are not methylated in all contigs in the bin                           |

## License

Nanomotif is released under the [MIT License](https://github.com/your-username/nanomotif/blob/main/LICENSE). Feel free to use, modify, and distribute the package in accordance with the terms of the license.

## Acknowledgments

Nanomotif builds upon various open-source libraries and tools that are instrumental in its functionality. We would like to express our gratitude to the developers and contributors of these projects for their valuable work.


