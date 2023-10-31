# Nanomotif

Nanomotif is a Python package that provides functionality for identifying methylation motifs on references using Nanopore sequencing.

## Installation

### Local Environment

To install Nanomotif in your local Python environment, follow these steps:

```shell
python3 -m venv nanomitf
source nanomitf/bin/activate
pip install nanomotif
```

### Conda Environment

If you prefer using Conda for managing your Python environments, you can create a new environment and install Nanomotif as follows:

```shell
conda create -n nanomotif python=3.9
conda activate nanomotif
python -m pip install nanomotif
```

## Required files

The required files are modkit pileup output and assembly sequences. 


```shell
samtools fastq -T MM,ML {basecalls} > {fastq}

# Your prefered assembeler > assembly

minimap2 -ax map-ont -y {assembly} {fastq} |
        samtools view -bS |
        samtools sort > {alignment}
        samtools index {alignment}

modkit pileup {alignment} --only-tabs
```
## Example Usage

```shell
nanomotif [assembly] [modkit pileup] [output]
```


## Documentation [Not yet implemented]

For detailed documentation and examples of all available functionalities in Nanomotif, please refer to the [official documentation](https://nanomotif-docs/docs). It provides comprehensive information on the various classes, methods, and parameters, along with usage examples and explanations.


## Output description

Nanomotif output motif at different processing levels
- **motifs-raw.tsv** - contain all naively detected motifs
- **motifs-score-filtered.tsv**, minimum score filtration
- **motifs-score-sub-filtered**, above filtration and removal of motif that are a submotif of another motif within contigs
- **motifs-score-sub-noise-filtered.tsv**, above filtration and removal of motifs with too many isolated bases, e.g. ..G....C..A is not considered a motif
- **motifs.tsv**, above filtration and merging of motif based on edit distance. required edit distance is based on expected distanve vs. motif length (estimated from REBASE gold standar motifs)

Description of the columns in all of the above mentioned files
| **Column**       | **Description**                                                                                       |
|------------------|-------------------------------------------------------------------------------------------------------|
| **sequence**     | padded motif all of equal length                                                                      |
| **score**        | score used for selecting motif during search. -1 indicate a motif merging has occoured and scoring is not present |
| **contig**       | reference in which the motif was found                                                                |
| **mod_type**     | the type of modification [a or m]                                                                     |
| **motif**        | trimmed motif without padding. braces are used to indicate multi base match and . used to indicate any base (N) |
| **mod_position** | position within the motif where the modification is located, 0-based index.                            |
| **alpha**        | alpha parameter of the BetaBernoulli posterior                                                       |
| **beta**         | beta parameter of the BetaBernoulli posterior                                                        |


## Contributing

We welcome contributions to Nanomotif! If you encounter any issues, have suggestions for improvements, or would like to add new features, please open an issue or submit a pull request on the [Nanomotif GitHub repository](https://github.com/SorenHeidelbach/nanomotif). We appreciate your feedback and contributions to make Nanomotif even better.

## License

Nanomotif is released under the [MIT License](https://github.com/your-username/nanomotif/blob/main/LICENSE). Feel free to use, modify, and distribute the package in accordance with the terms of the license.

## Acknowledgments

Nanomotif builds upon various open-source libraries and tools that are instrumental in its functionality. We would like to express our gratitude to the developers and contributors of these projects for their valuable work.


