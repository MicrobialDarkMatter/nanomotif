from setuptools import setup, find_packages
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
        "mtase_linker/src/*"
    ]},
    zip_safe=False,
    install_requires=[
        "wheel",
        "requests",
        "numpy>=1.24.4",
        "pandas>=2.0.2",
        "polars>=0.19",
        "scipy>=1.10.1",
        "networkx>=3.1",
        "pyarrow>=15.0.2",
        "Bio>=1.6.2",
        "snakemake>=7.32.4",
        "progressbar2>=3.53.1"
    ],
    entry_points={
            'console_scripts': [
                  'nanomotif = nanomotif.main:main'
            ]
    }
)
