# snakemake tutorial
Adapted from: https://snakemake.readthedocs.io/en/stable/tutorial/basics.html

Requirements:
- Linux or Mac
- Conda >= 4.9.0


After complete Conda installation you need to add the bioconda channel:

`conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge`

Create a conda environment for snakemake tutorial with the following command:

`conda env create -f envs/config.yaml -n snakemake-lab-training`


To activate the **snakemake-lab-training** environment, execute:
`conda activate snakemake-lab-training`

To install snakemake in a easy way, execute:
`mamba install -c conda-forge -c bioconda snakemake`

Everything is ready to start!!

## Basics: An example of variant calling workflow
