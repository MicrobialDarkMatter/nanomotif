from setuptools import setup, find_packages

setup(
    name='nanomotif',
    version='0.1.0',
    description='Identifying methlyation motifs in nanopore data',
    author='AAU_DarkScience',
    author_email='shei@bio.aau.com',
    license='MIT',
    packages=find_packages(),
    zip_safe=False,
    install_requires=[
        "wheel",
        "requests",
        "numpy==1.24.4",
        "pandas==2.0.2",
        "polars==0.18.3",
        "seaborn==0.12.2",
        "scipy==1.10.1",
        "networkx==3.1"
    ],
    entry_points={
            'console_scripts': [
                  'nanomotif = nanomotif.main:main'
            ]
    }
)