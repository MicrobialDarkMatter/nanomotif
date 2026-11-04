# Installation

Nanomotif runs on Python 3.9 or newer. The motif discovery module works on its own; the
[MTase-linker](#mtase-linker-installation) module needs extra dependencies, installed separately.

## With pip (recommended)

```shell
python3 -m venv nanomotif
source nanomotif/bin/activate
pip install nanomotif
```

This also installs `nanomotif-rs`, Nanomotif's compiled kernels (fast pileup reading and motif counting),
as a prebuilt wheel on:

- Linux x86_64 and aarch64
- macOS (Apple silicon and Intel)
- Windows x86_64

No compiler is needed. On other platforms, pip installs Nanomotif without `nanomotif-rs`, and the same
steps run in pure Python. Results are the same; it is slower and uses more memory for plain pileups.

To check whether the compiled kernels are in use:

```shell
python -c "import nanomotif_rs; print('nanomotif-rs installed')"
```

## With conda

```shell
conda create -n nanomotif python=3.12
conda activate nanomotif
conda install -c conda-forge -c bioconda nanomotif
```

The Bioconda package can lag behind PyPI, and it may not include `nanomotif-rs`. To add the compiled
kernels to a conda environment, run `pip install nanomotif-rs==<the installed nanomotif version>` in it.

## Check the installation

```shell
nanomotif check_installation
```

This runs motif discovery on a small test dataset that is shipped with Nanomotif.

---

## MTase-linker installation

The MTase-linker module has additional dependencies that are not installed with Nanomotif. Before using
it, install them with the `MTase-linker install` command. It needs conda on your system.

MTase-linker also needs snakemake. It is installed with Nanomotif on Python 3.11 and newer. On Python 3.9
and 3.10 it is not (there, one of its dependencies has to be compiled); install it yourself, for example
with `conda install -c conda-forge -c bioconda snakemake`, or use Python 3.11 or newer. Motif discovery
does not need snakemake.

```shell
nanomotif MTase-linker install
```

This creates a folder named `ML_dependencies` in the current directory, containing the dependencies of
MTase-linker. Use `--dependency_dir` to install it elsewhere. The installation takes a while: it creates
conda environments and runs the workflow on a small dataset to check it.

```
usage: nanomotif MTase-linker install [-h] [-d DEPENDENCY_DIR]

optional arguments:
  -h, --help            show this help message and exit
  -d DEPENDENCY_DIR, --dependency_dir DEPENDENCY_DIR
                        Path to installation location of dependencies. A folder named ML_dependencies
                        will be generated. Default is cwd.
```

---

## Installing from source

For development:

```shell
git clone https://github.com/MicrobialDarkMatter/nanomotif.git
cd nanomotif
pip install -e .
```

`pip install -e .` installs `nanomotif-rs` from PyPI on the supported platforms. To build the compiled
kernels from the repository instead (needs a Rust toolchain), run `pip install ./rust` first.
