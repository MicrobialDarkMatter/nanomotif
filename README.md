# Nanomotif

Nanomotif is a Python package that provides functionality for identifying methylation motifs on references using Nanopore sequencing.

## Installation

### Local Environment

To install Nanomotif in your local Python environment, follow these steps:

```shell
python3 -m venv myenv
source myenv/bin/activate
git clone https://github.com/SorenHeidelbach/nanomotif.git
cd nanomotif
pip install -r requirements.txt
pip install .
```

### Conda Environment

If you prefer using Conda for managing your Python environments, you can create a new environment and install Nanomotif as follows:

```shell
conda create -n nanomotif-env python=3.9
conda activate nanomotif-env
git clone https://github.com/your-username/nanomotif.git
cd nanomotif
conda install --file requirements.txt
python -m pip install .
```

## Required files

The required files are modkit pileup output and assembly sequences. 


```shell
samtools fastq -T MM,ML {basecalls} > {fastq}

# Your prefered assembeler > assembly

minimap2 -ax map-ont -y {assembly} {fastq} |
        sed -e 's/;C+?,/;C+C?,/g' |
        samtools view -bS |
        samtools sort > {alignment}
        samtools index {alignment}

modkit pileup {alignment} --only-tabs
```
## Example Usage




## Documentation [Not yet implemented]

For detailed documentation and examples of all available functionalities in Nanomotif, please refer to the [official documentation](https://nanomotif-docs/docs). It provides comprehensive information on the various classes, methods, and parameters, along with usage examples and explanations.

## Contributing

We welcome contributions to Nanomotif! If you encounter any issues, have suggestions for improvements, or would like to add new features, please open an issue or submit a pull request on the [Nanomotif GitHub repository](https://github.com/SorenHeidelbach/nanomotif). We appreciate your feedback and contributions to make Nanomotif even better.

## License

Nanomotif is released under the [MIT License](https://github.com/your-username/nanomotif/blob/main/LICENSE). Feel free to use, modify, and distribute the package in accordance with the terms of the license.

## Acknowledgments

Nanomotif builds upon various open-source libraries and tools that are instrumental in its functionality. We would like to express our gratitude to the developers and contributors of these projects for their valuable work.


