from setuptools import setup, find_packages
from setuptools.command.install import install
import os
import shutil
import urllib.request
import sys
import platform

METHYLATION_UTILS_URL = {
    "Linux": "https://github.com/SebastianDall/methylation_utils/releases/download/v0.1.0/methylation_utils",
}


def download_methylation_utils(url, dest_path):
    try:
        print(f"Attempting to download binary from {url}...")
        urllib.request.urlretrieve(url, dest_path)
        print("Download completed successfully.")
        # Make the file executable for Unix-like systems
        if platform.system() != "Windows":
            os.chmod(dest_path, 0o755)
    except urllib.error.URLError as e:
        print(f"Failed to download binary from {url}. URL Error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred while downloading binary: {e}")
        sys.exit(1)


class InstallCommand(install):
    def run(self):
        # Determine the binary URL based on the platform
        system = platform.system()
        binary_url = METHYLATION_UTILS_URL.get(system)
        if not binary_url:
            sys.exit(f"Unsupported platform: {system}")
        
        # Download the binary to the correct location
        destination_path = os.path.join("nanomotif/bin", "methylation_utils")
        download_methylation_utils(binary_url, destination_path)

        # Continue with the rest of the installation
        super().run()

    

exec(open('nanomotif/_version.py').read())
setup(
    name='nanomotif',
    version=__version__,
    description='Identifying methlyation motifs in nanopore data',
    author='AAU_DarkScience',
    author_email='shei@bio.aau.com',
    license='MIT',
    packages=find_packages(),
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    include_package_data=True,
    package_data={'nanomotif': [
        'datasets/*',
        'mtase_linker/*.smk',
        '*.yaml',
        "mtase_linker/envs/*",
        "mtase_linker/src/*",
        "bin/methylation_utils"
    ]},
    zip_safe=False,
    install_requires=[
        "wheel",
        "requests",
        "numpy>=1.24.4",
        "pandas>=2.0.2",
        "polars>=0.19,<=0.20.23",
        "scipy>=1.10.1",
        "networkx>=3.1",
        "pyarrow>=15.0.2",
        "Bio>=1.6.2",
        "snakemake>=7.32.4",
        "progressbar2>=3.53.1"
    ],
    cmdclass = {
        "install": InstallCommand
    },
    entry_points={
            'console_scripts': [
                  'nanomotif = nanomotif.main:main'
            ]
    }
)
