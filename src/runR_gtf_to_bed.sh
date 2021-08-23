#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=R_bam_to_bigwig
#SBATCH -c 1
#SBATCH --mem=80G
#SBATCH --time=1-18:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu


echo "---------------  Running R script gtf_to_bed.R ---------------------"

module load R
R --slave -f R/gtf_to_bed.R

