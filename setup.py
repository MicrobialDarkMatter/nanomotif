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
    include_package_data=True,
    package_data={
        'nanomotif': ['datasets/*'],
    },
    zip_safe=False,
    install_requires=[
        "wheel",
        "requests",
        "numpy==1.24.4",
        "pandas==2.0.2",
        "polars>=0.19",
        "seaborn==0.12.2",
        "scipy==1.10.1",
        "networkx==3.1",
        "progressbar2==3.53.1"
    ],
    entry_points={
            'console_scripts': [
                  'nanomotif = nanomotif.main:main'
            ]
    }
)