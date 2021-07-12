# SNAKEMAKE
Make conda env and install as follows

Install mamba in base envs
    conda install -n base -c conda-forge mamba
    conda activate base

Create sankemake env
    mamba create -c conda-forge -c bioconda -n snakemake snakemake

After installing mamba, create snakmake env
    conda activate snakemake

Copy METEORE data using command below
    git clone https://github.com/comprna/METEORE.git
    cd METEORE/

Run the pipeline with your own dataset by replacing example folder in the data directory with your folder containing the fast5 files. You will use the fast5 folder name to specify your target output file in the Snakemake pipeline
Simply replace example in the output file with your fast5 folder name in the command line below.
place the reference genome file in .fasta format in a folder named data, and re-define the reference genome file within the Snakefile (Nanopolish, Deepsignal1, Tombo, Guppy) by replacing ecoli_k12_mg1655.fasta with your specified reference genome

## Nanopolish snakmake pipeline

To remove an env do as follows:
    conda env remove --name env_name

Create your conda envs
    mamba create -c conda-forge -c bioconda -n meteore snakemake
    conda activate meteore
    mamba install -c bioconda nanopolish samtools r-data.table r-dplyr r-plyr

As memory requirements are unknown, run in himem partition but don't ask for too much memory 
    screen -a
    srun --partition himem --time 0-12:00:00 --mem-per-cpu 40G --cpus-per-task 24 --pty bash

Basic running of programme
    snakemake -s Nanopolish nanopolish_results/example_nanopolish-freq-perCG.tsv --cores all
