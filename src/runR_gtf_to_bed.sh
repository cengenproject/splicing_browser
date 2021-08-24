#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=R_gtf_to_bed
#SBATCH -c 1
#SBATCH --mem=80G
#SBATCH --time=1-18:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu


echo "---------------  Running R script gtf_to_bed.R ---------------------"

module load R
R --slave -f R/gtf_to_bed.R


echo "---------------  Converting to bigBed ---------------------"
out_dir="outputs_visualization/references"


bedSort $out_dir/WS277_canonical_geneset.bed12 $out_dir/WS277_canonical_geneset.sorted.bed

bedToBigBed $out_dir/WS277_canonical_geneset.sorted.bed \
            $out_dir/chrom.WS277.sizes \
            $out_dir/WS277_canonical_geneset.bb
