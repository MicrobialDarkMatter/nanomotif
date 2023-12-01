# Nanomotif

Nanomotif is a Python package that provides functionality for identifying methylation motifs on references using Nanopore sequencing.

## Installation

#### Local Environment

To install Nanomotif in your local Python environment, follow these steps:

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

## Required files

To identify methylated motifs, two files are required: 
- [Modkit pileup](https://github.com/nanoporetech/modkit/blob/master/book/src/advanced_usage.md#pileup) output
- Assembly

If you want to identify bin consensus motifs, a file specifying contig bin relationship is required. It should be headerless and formatted as a TSV file with contig names in column 1 and respective bin in column 2.

*OBS: when demultiplexing, trimming of reads may result in errors downstream. We therefore recommend using untrimmed reads for mapping*
```shell
samtools fastq -T MM,ML {DORADO_BASECALLS.sam} > {READS.fastq}
minimap2 -ax map-ont -y {ASSEMBLY.fasta} {READS.fastq} |
        samtools view -bS |
        samtools sort > {ALIGNMENT.bam}
        samtools index {ALIGNMENT.bam}
modkit pileup {ALIGNMENT.bam} --only-tabs
```
## Example Usage

Nanomotif have several commands. The main command is complete-workflow, which runs motifs finding, motif scoring and bin consensus. 

```
usage: nanomotif [-h] [--version] {complete-workflow,find-motifs,score-motifs,bin-consensus} ...

Motif identification and utilisation commands

positional arguments:
  {complete-workflow,find-motifs,score-motifs,bin-consensus}
                        sub-command help
    complete-workflow   Run the complete workflow
    find-motifs         Find motifs in contigs of the assembly
    score-motifs        Find degree of methylation of all identified motifs for all contigs in the assembly
    bin-consensus       Find consensus motif in each bin

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
```

```
usage: nanomotif complete-workflow [-h] [-t THREADS] [-v] [--search_frame_size SEARCH_FRAME_SIZE]
                                   [--threshold_methylation_confident THRESHOLD_METHYLATION_CONFIDENT]
                                   [--threshold_methylation_general THRESHOLD_METHYLATION_GENERAL] [--threshold_valid_coverage THRESHOLD_VALID_COVERAGE]
                                   [--minimum_kl_divergence MINIMUM_KL_DIVERGENCE]
                                   assembly pileup bins output

positional arguments:
  assembly              Assembly file.
  pileup                Modkit pileup file.
  bins                  File specifying to which bin contigs belong. (tsv file with no header and columns: contig, bin)
  output                Output directory.

optional arguments:
  -h, --help            show this help message and exit
  -t THREADS, --threads THREADS
                        Number of threads to use. Default is 1.
  -v, --verbose         Increase output verbosity. (set logger to debug level)
  --search_frame_size SEARCH_FRAME_SIZE
                        Length of the sequnces sampled around methylatyion sites. Default: 40
  --threshold_methylation_confident THRESHOLD_METHYLATION_CONFIDENT
                        Minimum fraction of reads that must be methylated to be considered in motif search. Default: 0.8
  --threshold_methylation_general THRESHOLD_METHYLATION_GENERAL
                        Minimum fraction of reads that must be methylated for a position to be considered modified. Default: 0.5
  --threshold_valid_coverage THRESHOLD_VALID_COVERAGE
                        Minimum valid base coverage for a position to be considered. Default: 5
  --minimum_kl_divergence MINIMUM_KL_DIVERGENCE
                        Minimum KL-divergence for a position to considered for expansion in motif search. Default: 0.2
```

## Output description


### Find motifs output
The main output is `motifs.tsv`, which contain the identified motifs within a contigs. It contains the following:

| **Column**       | **Description**                                                                                       |
|------------------|-------------------------------------------------------------------------------------------------------|
| **contig**       | contig in which the motif was found                                                                |
| **mod_type**     | the type of modification [a (6mA) or m (5mC)]                                                                     |
| **motif**        | sequence of the detected motif in regex form (braces indicate multi base position and . indicate N positions) |
| **mod_position** | position within the motif where the methylation is located. 0-based index.                            |
| **alpha**        | number of motif positions that are methylated in the contig                                                       |
| **beta**         | number of motif positions that are not methylated in the contig                                                        |
| **motif_iupac**        | sequence of the detected motif in [IUPAC](https://www.bioinformatics.org/sms/iupac.html) form  |
| **motifs_type** | what kind of motif the sequence represents (palindrome, non-palindrome, bipartite or ambiguous)

In addition to `motifs.tsv`, a `precleanup-motifs` is outputted containing raw motifs. The formats of these files follow ``motifs.tsv

### Score motifs output

The output is exactly the same format as the output from find motifs. It only contain alpha and beta values for all identified motifs in all contigs for all contigs. 

### Bin consensus 

Identifies the consensus motifs with in bins and reports degree of methylation within the bin. The output is written to `bin-motfs.tsv`, which contain the following:

| **Column**       | **Description**                                                                                       |
|------------------|-------------------------------------------------------------------------------------------------------|
| **bin**       | the bin                                                                |
| **motif**        | sequence of the motif in [IUPAC](https://www.bioinformatics.org/sms/iupac.html) form  |
| **mod_position** | position within the motif where the methylation is located. 0-based index.                            |
| **mod_type**     | the type of modification [a (6mA) or m (5mC)]                                                                     |
| **contig_count**        | number contigs in the bin                                                       |
| **contigs_with_motif**         | number of contigs in the bin with mean methylation above 50%                                                        |
| **motifs_type** | what kind of motif the sequence represents (palindrome, non-palindrome, bipartite or ambiguous) |
| **methylated_count** | number of methylated positions in all contigs in the bin |
| **non_methylated_count** | number of non-methylated position in all contigs in the bin |
| **mean_methylation_bin** | mean methylation of all contigs in the bin |



## License

Nanomotif is released under the [MIT License](https://github.com/your-username/nanomotif/blob/main/LICENSE). Feel free to use, modify, and distribute the package in accordance with the terms of the license.

## Acknowledgments

Nanomotif builds upon various open-source libraries and tools that are instrumental in its functionality. We would like to express our gratitude to the developers and contributors of these projects for their valuable work.


