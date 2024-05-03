# Installation

Nanomotif by default include the **motif discovery**, **bin contamination**, and **include unbinned** contigs modules. These modules can operate independently without the additional dependencies required for the **MTase-Linker module**. For details on the MTase-Linker setup, see [MTase-linker installation](#mtase-linker-installation).

## Conda Environment

Using Conda for managing your Python environments, you can create a new environment and install Nanomotif as follows:

```shell
conda create -n nanomotif  python=3.9
conda activate nanomotif
conda install -c bioconda nanomotif
```

## Local Environment

To install Nanomotif in a local Python environment:

```shell
python3 -m venv nanomotif
source nanomotif/bin/activate
pip install nanomotif
```

## Check installation
Once installed, the installation can be checked by running:
```shell
nanomotif check-installation
```
This runs a test run on a small dataset, ensuring everything works.


## MTase-linker installation
The MTase-Linker module has additional dependencies that are not automatically installed with Nanomotif. Therefore, before using this module, you must manually install these dependencies using the `MTase-linker install` command.
The `MTase-linker` module requires that conda is available on your system.

```shell
nanomotif MTase-linker install
```

This will create a folder named `ML_dependencies` in your current working directory, containing the required dependencies for the MTase-linker module. You can use the `--dependency_dir` flag to change the installation location of the `ML_dependencies` folder.

```
usage: nanomotif MTase-linker install [-h] [-d DEPENDENCY_DIR]

optional arguments:
  -h, --help            show this help message and exit
  -d DEPENDENCY_DIR, --dependency_dir DEPENDENCY_DIR
                        Path to the directory, where dependencies should be installed. A folder named
                        ML_dependencies will be generated. Default is cwd.
```