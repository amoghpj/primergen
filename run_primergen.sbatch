#!/bin/bash
#SBATCH -n 1
#SBATCH -t 0-1:00
#SBATCH -p short
#SBATCH --mem=1G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=amogh_jalihal@hms.harvard.edu
#SBATCH -o %A.out
#SBATCH -e %A.err
snakemake --snakefile primergen.snakefile -p\
 --rerun-incomplete\
 --cores 4\
 --use-conda\
 --configfile config.yaml\
 -j 100 --executor slurm --default-resource slurm_account=springer_ms182 slurm_partition=short\
 --keep-incomplete --latency-wait 60
