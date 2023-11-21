#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=gen_browser_data
#SBATCH -c 12
#SBATCH --mem=210G
#SBATCH --time=23:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu

set -e

# this script takes the output from a "bsn" alignment pipeline, and calls R scripts to
# generate the data to visualize in the genome browser.

# parameters
WS="289"
out_version="231116"
sj_dir="/gpfs/gibbs/pi/hammarlund/CeNGEN/bulk/bulk_alignments/bsn12_junctions"
bams_dir="/gpfs/gibbs/pi/hammarlund/CeNGEN/bulk/bulk_alignments/bsn12_bams"

outliers_to_ignore="data/outliers_to_ignore.txt"


# Initializations

echo "---------------  Starting $(date) ---------------------"

module load R

out_dir=data/outs/${out_version}_browser

mkdir -p $out_dir





## Junction processing ----

mkdir -p $out_dir/sj $out_dir/sj/single_sample $out_dir/sj/single_neuron $out_dir/sj/global

Rscript R/sj_to_bed.R $WS $out_version $outliers_to_ignore $sj_dir

# convert to BigBed
chr_sizes="/home/aw853/project/references/WS"$WS"/chrom.sizes"
for file in $out_dir/sj/*/*.bed
do
  bedToBigBed $file $chr_sizes ${file%.bed}.bb
done




## Coverage processing

mkdir -p $out_dir/coverage $out_dir/coverage/raw_RLEs $out_dir/coverage/single_sample
mkdir -p $out_dir/coverage/single_neuron $out_dir/coverage/global

Rscript R/bam_to_bigwig.R $WS $out_version $bams_dir $outliers_to_ignore



# Prepare for export
cd data/outs/
tar czf ${out_version}_browser.tar.gz ${out_version}_browser/*/*/*.bb ${out_version}_browser/*/*/*.bw

echo "Send to vps with:"
echo "scp $(pwd)/${out_version}_browser.tar.gz cengen-vps:/var/www/public_data/splicing"

