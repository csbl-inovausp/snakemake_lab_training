# snakemake lab training tutorial

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

    rule bwa_map:
      input:
          "data/genome.fa",
          "data/samples/A.fastq"
      output:
          "mapped_reads/A.bam"
      shell:
          "bwa mem {input} | samtools view -Sb - > {output}"

Show the execution plan:

`snakemake -np mapped_reads/A.bam`

Show the execution plan visually of the resulting directed acyclic graph (DAG) of jobs:

`snakemake --dag mapped_reads/A.bam | dot -Tsvg | display`

it should looks like this:

![Example](./dags/dag_bwa.svg)

Execute the workflow:

`snakemake --cores 1 mapped_reads/A.bam`

It fails by first time because we dont have **bwa** and **samtools** installed in our environment. In high performace cluster systems (HPC), it can be preferable to use environment modules for deployment of optimized versions of certain standard tools. Snakemake allows to define environment modules per rule:

rule bwa_map:
  input:
      "data/genome.fa",
      "data/samples/A.fastq"
  output:
      "mapped_reads/A.bam"
  conda:
      "envs/bwa.yaml"    
  shell:
      "bwa mem {input} | samtools view -Sb - > {output}"

### Step 2: Generalizing the read mapping rule using WildCards

Snakemake allows generalizing rules by using named wildcards (https://en.wikipedia.org/wiki/Wildcard_character). Simply replace the A in the second input file and in the output file with the wildcard {sample}, leading to:


      rule bwa_map:
          input:
              "data/genome.fa",
              "data/samples/{sample}.fastq"
          output:
              "mapped_reads/{sample}.bam"
          conda:
              "envs/bwa.yaml"
          shell:
              "bwa mem {input} | samtools view -Sb - > {output}"



Show the execution plan visually of the resulting directed acyclic graph (DAG) of jobs:

`snakemake --dag mapped_reads/A.bam | dot -Tsvg | display`

![Example](./dags/dag_bwa_generalized.svg)

Question: What changed in the resulting DAG?

### Step 3: Sorting read alignments
For later steps, we need the read alignments in the BAM files to be sorted. This can be achieved with the samtools sort command. We add the following rule beneath the bwa_map rule:

    rule samtools_sort:
      input:
          "mapped_reads/{sample}.bam"
      output:
          "sorted_reads/{sample}.bam"
      conda:
          "envs/bwa.yaml"
      shell:
          "samtools sort -T sorted_reads/{wildcards.sample} "
          "-O bam {input} > {output}"          


Snakemake dependencies are resolved automatically by matching file names.

For example, by typying:

`snakemake -np sorted_reads/B.bam`

You will see how Snakemake wants to run first the rule **bwa_map** and then the rule **samtools_sort** to create the desired target file.

### Step 4: Indexing read alignments and visualizing the DAG of jobs

Next, we need to use samtools again to index the sorted read alignments so that we can quickly access reads by the genomic location they were mapped to. This can be done with the following rule:

    rule samtools_index:
      input:
          "sorted_reads/{sample}.bam"
      output:
          "sorted_reads/{sample}.bam.bai"
      conda:
          "envs/bwa.yaml"
      shell:
          "samtools index {input}"



Having three steps already, it is a good time to take a closer look at the resulting directed acyclic graph (DAG) of jobs. By executing:

`snakemake --dag sorted_reads/{A,B}.bam.bai | dot -Tsvg | display`

![Example](./dags/sorted_reads.svg)

### Step 5: Calling genomic variants

The next step in our workflow will aggregate the mapped reads from all samples and jointly call genomic variants on them (see Background). For the variant calling, we will combine the two utilities samtools and bcftools. Snakemake provides a helper function for collecting input files that helps us to describe the aggregation in this step.

we can add the following rule to our Snakefile:

    rule bcftools_call:
        input:
            fa="data/genome.fa",
            bam=expand("sorted_reads/{sample}.bam", sample=SAMPLES),
            bai=expand("sorted_reads/{sample}.bam.bai", sample=SAMPLES)
        output:
            "calls/all.vcf"
        conda:
            "envs/bwa.yaml"
        shell:
            "samtools mpileup -g -f {input.fa} {input.bam} | "
            "bcftools call -mv - > {output}"

Here, we use only the simple case of expand. We first let Snakemake know which samples we want to consider.Remember that Snakemake works backwards from requested output, and not from available input. Thus, it does not automatically infer all possible output from, for example, the fastq files in the data folder. Also remember that Snakefiles are in principle Python code enhanced by some declarative statements to define workflows. Hence, we can define the list of samples ad-hoc in plain Python at the top of the Snakefile:

`SAMPLES = ["A", "B"]`


### Step 6: Using custom scripts

Usually, a workflow not only consists of invoking various tools, but also contains custom code to for example calculate summary statistics or create plots. While Snakemake also allows you to directly write Python code inside a rule, it is usually reasonable to move such logic into separate scripts. For this purpose, Snakemake offers the script directive. Add the following rule to your Snakefile:

    rule plot_quals:
        input:
            "calls/all.vcf"
        output:
            "plots/quals.svg"
        conda:
            "envs/matplots.yaml"
        script:
            "scripts/plot-quals.py"



### Step 7: Adding a target rule

So far, we always executed the workflow by specifying a target file at the command line. Apart from filenames, Snakemake also accepts rule names as targets if the requested rule does not have wildcards. Hence, it is possible to write target rules collecting particular subsets of the desired results or all results. Moreover, if no target is given at the command line, Snakemake will define the first rule of the Snakefile as the target. Hence, it is best practice to have a rule all at the top of the workflow which has all typically desired target files as input files.

Here,  we add a rule at the top of the workflow:

    rule all:
        input:
            "plots/quals.svg"




From now, we can see the execution plan:

`snakemake -n`

Or visually by typying:

`snakemake --dag | dot -Tsvg | display `

![Example](./dags/final_dag.svg)

Finally, we can execute everything in a simple and concise command:

`snakemake --use-conda`


### BONUS: Generating a Snakemake report

From Snakemake 5.1 on, it is possible to automatically generate detailed self-contained HTML reports that encompass runtime statistics, provenance information, workflow topology and results. Let's update our previous rule by adding a simple report with the generated figure:

    rule plot_quals:
        input:
            "calls/all.vcf"
        output:
            report("plots/quals.svg", caption="report/quals.rst", category="Variant Calling")
        conda:
            "envs/matplots.yaml"
        script:
            "scripts/plot-quals.py"

Finally, to create a snakemake report simply run:

`snakemake --report report.html`


## Summary

In total, the resulting workflow should look like this:

    SAMPLES = ["A", "B"]


    rule all:
        input:
            "plots/quals.svg"


    rule bwa_map:
        input:
            "data/genome.fa",
            "data/samples/{sample}.fastq"
        output:
            "mapped_reads/{sample}.bam"
        conda:
            "envs/bwa.yaml"
        shell:
            "bwa mem {input} | samtools view -Sb - > {output}"


    rule samtools_sort:
        input:
            "mapped_reads/{sample}.bam"
        output:
            "sorted_reads/{sample}.bam"
        conda:
            "envs/bwa.yaml"
        shell:
            "samtools sort -T sorted_reads/{wildcards.sample} "
            "-O bam {input} > {output}"


    rule samtools_index:
        input:
            "sorted_reads/{sample}.bam"
        output:
            "sorted_reads/{sample}.bam.bai"
        conda:
            "envs/bwa.yaml"
        shell:
            "samtools index {input}"


    rule bcftools_call:
        input:
            fa="data/genome.fa",
            bam=expand("sorted_reads/{sample}.bam", sample=SAMPLES),
            bai=expand("sorted_reads/{sample}.bam.bai", sample=SAMPLES)
        output:
            "calls/all.vcf"
        conda:
            "envs/bwa.yaml"
        shell:
            "samtools mpileup -g -f {input.fa} {input.bam} | "
            "bcftools call -mv - > {output}"


    rule plot_quals:
        input:
            "calls/all.vcf"
        output:
            report("plots/quals.svg", caption="report/quals.rst", category="Variant Calling")
        conda:
            "envs/matplots.yaml"
        script:
            "scripts/plot-quals.py"
