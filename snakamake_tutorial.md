# snakemake tutorial
Adapted from: https://snakemake.readthedocs.io/en/stable/tutorial/basics.html

## 1 Requirements:
- Snakemake is more suitable for Linux or Mac Systems but it can be applied for Windows (https://snakemake.readthedocs.io/en/stable/tutorial/setup.html)
- Conda >= 4.9.0 (https://docs.conda.io/en/latest/miniconda.html)


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

The genome of a living organism encodes its hereditary information. It serves as a blueprint for proteins, which form living cells, carry information and drive chemical reactions. Differences between species, populations or individuals can be reflected by differences in the genome. Certain variants can cause syndromes or predisposition for certain diseases, or cause cancerous growth in the case of tumour cells that have accumulated changes with respect to healthy cells. In order to recover the genome of the sample, one has to map these reads against a known reference genome (for example, the human one obtained during the famous human genome project). This task is called read mapping. Often, it is of interest where an individual genome is different from the species-wide consensus represented with the reference genome. Such differences are called variants. They are responsible for harmless individual differences (like eye color), but can also cause diseases like cancer.

### Step 1: Mapping reads
Create a new file called Snakefile with an editor of your choice. We propose to use the Atom editor, since it provides out-of-the-box syntax highlighting for Snakemake. In the Snakefile, define the following rule:
`rule bwa_map:
    input:
        "data/genome.fa",
        "data/samples/A.fastq"
    output:
        "mapped_reads/A.bam"
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"`

Show the execution plan:
`snakemake -np mapped_reads/A.bam`

Execute the workflow:
`snakemake --cores 1 mapped_reads/A.bam`

### Step 2: Generalizing the read mapping rule using WildCards

Snakemake allows generalizing rules by using named wildcards (https://en.wikipedia.org/wiki/Wildcard_character). Simply replace the A in the second input file and in the output file with the wildcard {sample}, leading to:

`rule bwa_map:
    input:
        "data/genome.fa",
        "data/samples/{sample}.fastq"
    output:
        "mapped_reads/{sample}.bam"
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"`
