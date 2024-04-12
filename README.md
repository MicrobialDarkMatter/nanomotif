# Nanomotif

Nanomotif is a Python package designed to explore methylation in prokaryotic genomes using Nanopore sequencing. Nanomotif is a fast, scalable, and sensitive tool for identification and utilization of methylation motifs in monocultures and metagenomic samples.  Nanomotif offers de novo methylated motif identification, metagenomic bin contamination detection, bin association of unbinned contigs, and linking of MTase genes to methylation motifs.

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

## Usage
### Motif identification
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
                        length of the sequences sampled around confident methylation sites. Default: 40
  --threshold_methylation_confident THRESHOLD_METHYLATION_CONFIDENT
                        minimum fraction of reads that must be methylated at a position for the position to be considered confiently methylated. These position are used to
                        search for candidate motifs. Default: 0.8
  --threshold_valid_coverage THRESHOLD_VALID_COVERAGE
                        minimum valid base coverage for a position to be considered. Default: 5
  --minimum_kl_divergence MINIMUM_KL_DIVERGENCE
                        minimum KL-divergence for a position to considered for expansion in motif search. Higher value means less exhaustive, but faster search. Default: 0.2
```

#### Output description
##### find-motifs and score-motifs

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

Running `find-motifs` generates pre-cleanup folder, which contains motif that got removed in the postprocessing steps. The name of the file indicates which postprocessing steps have been run on the motifs.

##### bin-consensus 
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


### Bin contamination
After motif identification it is possible to identify contamination in bins using the bin-motifs, contig-bin and motif-scored files. 

```
usage: nanomotif detect_contamination [-h] --motifs_scored MOTIFS_SCORED --bin_motifs BIN_MOTIFS --contig_bins CONTIG_BINS [-t THREADS] [--mean_methylation_cutoff MEAN_METHYLATION_CUTOFF]
                                      [--n_motif_contig_cutoff N_MOTIF_CONTIG_CUTOFF] [--n_motif_bin_cutoff N_MOTIF_BIN_CUTOFF] [--ambiguous_motif_percentage_cutoff AMBIGUOUS_MOTIF_PERCENTAGE_CUTOFF]
                                      [--write_bins] [--assembly_file ASSEMBLY_FILE] --out OUT [--contamination_file CONTAMINATION_FILE]

optional arguments:
  -h, --help            show this help message and exit
  --motifs_scored MOTIFS_SCORED
                        Path to motifs-scored.tsv from nanomotif
  --bin_motifs BIN_MOTIFS
                        Path to bin-motifs.tsv file
  --contig_bins CONTIG_BINS
                        Path to bins.tsv file for contig bins
  -t THREADS, --threads THREADS
                        Number of threads to use for multiprocessing
  --mean_methylation_cutoff MEAN_METHYLATION_CUTOFF
                        Cutoff value for considering a motif as methylated
  --n_motif_contig_cutoff N_MOTIF_CONTIG_CUTOFF
                        Number of motifs that needs to be observed in a contig before it is considered valid for scoring
  --n_motif_bin_cutoff N_MOTIF_BIN_CUTOFF
                        Number of motifs that needs to be observed in a bin to be considered valid for scoring
  --ambiguous_motif_percentage_cutoff AMBIGUOUS_MOTIF_PERCENTAGE_CUTOFF
                        Percentage of ambiguous motifs defined as mean methylation between 0.05 and 0.40 in a bin. Motifs with an ambiguous methylation percentage of more than this value are removed from
                        scoring. Default is 0.40
  --write_bins          If specified, new bins will be written to a bins folder. Requires --assembly_file to be specified.
  --assembly_file ASSEMBLY_FILE
                        Path to assembly.fasta file
  --out OUT             Path to output directory
  --contamination_file CONTAMINATION_FILE
                        Path to an existing contamination file if bins should be outputtet as a post-processing step
```

If `detect_contamination` was run without the `--write_bins` flag, bins can be written as a post processing step if the `--contamination_file` is specified along with the `--write_bins` flag and the `--assembly_file` flag.

The output is a `bin_contamination.tsv` file with the following columns:

| **Column**       | **Description**                                                                                       |
|------------------|-------------------------------------------------------------------------------------------------------|
| **bin**          | bin to which the motif belong                                                                         |
| **bin_contig_compare** | The contig for which the methylation pattern is compared. The name is a concatenation of bin + "_" + contig. |
| **binary_methylation_missmatch_score** | Number of methylation missmatches between the contig and the bin_consensus pattern. |
| **non_na_comparisons** | Number of non-NA comparisons between the contig and the bin_consensus pattern. |
| **contig** | The contig for which the methylation pattern is compared. | 

### Include unbinned contigs
This module tries to assign contigs in the assembly file to bins by comparing the methylation pattern of the contig to the bin consensus. the contig must have a unique perfect match to the bin consensus to be assigned to a bin. Additionally, the `include_contigs` assigns all the contigs in the `bin_contamination.tsv` file to as unbinned. If decontamination should not be performed, the `include_contigs` can be run without the `--run_detect_contamination` flag or without the `--contamination_file` flag.

```
usage: nanomotif include_contigs [-h] --motifs_scored MOTIFS_SCORED --bin_motifs BIN_MOTIFS --contig_bins CONTIG_BINS [-t THREADS] [--mean_methylation_cutoff MEAN_METHYLATION_CUTOFF]
                                 [--n_motif_contig_cutoff N_MOTIF_CONTIG_CUTOFF] [--n_motif_bin_cutoff N_MOTIF_BIN_CUTOFF] [--ambiguous_motif_percentage_cutoff AMBIGUOUS_MOTIF_PERCENTAGE_CUTOFF]
                                 [--write_bins] [--assembly_file ASSEMBLY_FILE] --out OUT [--contamination_file CONTAMINATION_FILE | --run_detect_contamination]
                                 [--min_motif_comparisons MIN_MOTIF_COMPARISONS]

optional arguments:
  -h, --help            show this help message and exit
  --motifs_scored MOTIFS_SCORED
                        Path to motifs-scored.tsv from nanomotif
  --bin_motifs BIN_MOTIFS
                        Path to bin-motifs.tsv file
  --contig_bins CONTIG_BINS
                        Path to bins.tsv file for contig bins
  -t THREADS, --threads THREADS
                        Number of threads to use for multiprocessing
  --mean_methylation_cutoff MEAN_METHYLATION_CUTOFF
                        Cutoff value for considering a motif as methylated
  --n_motif_contig_cutoff N_MOTIF_CONTIG_CUTOFF
                        Number of motifs that needs to be observed in a contig before it is considered valid for scoring
  --n_motif_bin_cutoff N_MOTIF_BIN_CUTOFF
                        Number of motifs that needs to be observed in a bin to be considered valid for scoring
  --ambiguous_motif_percentage_cutoff AMBIGUOUS_MOTIF_PERCENTAGE_CUTOFF
                        Percentage of ambiguous motifs defined as mean methylation between 0.05 and 0.40 in a bin. Motifs with an ambiguous methylation percentage of more than this value are removed from
                        scoring. Default is 0.40
  --write_bins          If specified, new bins will be written to a bins folder. Requires --assembly_file to be specified.
  --assembly_file ASSEMBLY_FILE
                        Path to assembly.fasta file
  --out OUT             Path to output directory
  --contamination_file CONTAMINATION_FILE
                        Path to an existing contamination file to include in the analysis
  --run_detect_contamination
                        Indicate that the detect_contamination workflow should be run first
  --min_motif_comparisons MIN_MOTIF_COMPARISONS
                        Minimum number of non-NA motif comparisons required to include a contig in the analysis. Default is 5
```

The output is two files: `include_contigs.tsv` and `new_contig_bin.tsv`. The `include_contigs.tsv` file is the contigs that were assigned based on the methylation pattern and the `new_contig_bin.tsv` is the updated contig-bin file.

The `include_contigs.tsv` file has the following columns:

| **Column**       | **Description**                                                                                       |
|------------------|-------------------------------------------------------------------------------------------------------|
| **bin**          | bin to which the motif belong                                                                         |
| **bin_compare** | The contig for which the methylation pattern is compared. The name is a concatenation of bin + "_" + contig. |
| **binary_methylation_missmatch_score** | Number of methylation missmatches between the contig and the bin_consensus pattern. |
| **non_na_comparisons** | Number of non-NA comparisons between the contig and the bin_consensus pattern. |
| **contig** | The contig for which the methylation pattern is compared. | 

### Associating MTases and RM-systems to motifs
This module tries to link methylation motifs to their corresponding MTase and, when present, their entire RM system. 

The MTase-Linker module has additional dependencies that are not automatically installed with Nanomotif. Therefore, before using this module, you must manually install these dependencies using the `MTase-linker install`.

```
usage: nanomotif MTase-linker install [-h] [-d DEPENDENCY_DIR]

optional arguments:
  -h, --help            show this help message and exit
  -d DEPENDENCY_DIR, --dependency_dir DEPENDENCY_DIR
                        Path to the directory, where dependencies should be installed. A folder named
                        ML_dependencies will be generated. Default is cwd.
```

When the additional dependencies are installed you can run the workflow using `MTase-linker run`

To run the MTase-linker workflow, users must provide:
- A directory containing either the assembly file or bin files (files should have a .fa extension).
- A `contig_bin.tsv` file, which is necessary even when running the workflow on monoculturs. For monocultures, the bin column will contain the same entry for each contig.
- A `bin-motifs.tsv` file generated by the `nanomotif complete-workflow`
- Path to the `ML_dependencies` directory created during installation of the MTase-linker module.

```
usage: nanomotif MTase-linker run [-h] [-t THREADS] [--forceall FORCEALL] [--dryrun DRYRUN] --binsdir BINSDIR --contig_bin CONTIG_BIN --bin_motifs BIN_MOTIFS [-d DEPENDENCY_DIR] [-o OUTPUTDIR] [--identity IDENTITY] [--qcovs QCOVS]

optional arguments:
  -h, --help            show this help message and exit
  -t THREADS, --threads THREADS
                        Number of threads to use. Default is 1
  --forceall FORCEALL   Forcerun workflow regardless of already created output (default = False)
  --dryrun DRYRUN       Dry-run the workflow. Default is False
  --binsdir BINSDIR     Directory with bin files or assembly files. Needs to have the .fa extension.
  --contig_bin CONTIG_BIN
                        tsv file specifying which bin contigs belong.
  --bin_motifs BIN_MOTIFS
                        bin-motifs.tsv output from nanomotif.
  -d DEPENDENCY_DIR, --dependency_dir DEPENDENCY_DIR
                        Same as for installation. Path to directory of the ML_dependencies directory created during installation of the MTase-linker module. Default is cwd/ML_dependencies
  -o OUTPUTDIR, --outputdir OUTPUTDIR
                        Output directory. Default is cwd
  --identity IDENTITY   Identity threshold for motif prediction. Default is 80
  --qcovs QCOVS         Query coverage for motif prediction. Default is 80
```

**Example run**
```
nanomotif MTase-linker run -t 10 --binsdir bins_directory --contig_bin contig_bin.tsv --bin_motifs nanomotif/bin_motifs.tsv -d ML_dependencies -o mtase_linker
```
Running the nanomotif MTase-linker run command will generate output, containing two primary output files: mtase_assignment_table.tsv and nanomotif_assignment_table.tsv. The first file lists all predicted MTase genes in the genome along with their predicted methylation characteristics and whether the module was able to unambiguously assign any detected motifs to the MTase (`linked` = (True/False)).
The second file includes data from the bin-motifs.tsv of the nanomotif output with two additional columns `linked` and `candidate_genes`. The `linked` variable is a boolean indicator if the motif could be unambiguously linked to a MTase in the bin/genome (TRUE/FALSE). If True the gene_id of the MTase is provided in `candidate_gene`. If False, the `candidate_gene` variable lists feasible candidate facilitators of the modification based on motif type and modification type predictions.

MTase-linker also provides the raw outputs filtered outputs of dependency tools, which can be used for further analysis into cognate RM-system restriction enzymes and specificity subunits.

## License

Nanomotif is released under the [MIT License](https://github.com/your-username/nanomotif/blob/main/LICENSE). Feel free to use, modify, and distribute the package in accordance with the terms of the license.

## Acknowledgments

Nanomotif builds upon various open-source libraries and tools that are instrumental in its functionality. We would like to express our gratitude to the developers and contributors of these projects for their valuable work.


