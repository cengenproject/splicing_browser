#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=gen_browser_data
#SBATCH -c 1
#SBATCH --mem=70G
#SBATCH --time=18:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu

set -e

# this script takes the output from a "bsn" alignment pipeline, and calls R scripts to
# generate the data to visualize in the genome browser.

# parameters
WS="281"
out_version="220114"
bams_combined="/home/aw853/scratch60/2021-11-08_alignments"


# Initializations

echo "---------------  Starting $(date) ---------------------"

module load R

out_dir=data/outs/${out_version}_browser

mkdir -p $out_dir




# Junction processing

mkdir $out_dir/sj $out_dir/sj/single_sample $out_dir/sj/single_neuron $out_dir/sj/global

Rscript R/sj_to_bed.R $WS $out_version

# convert to BigBed
chr_sizes="/home/aw853/project/references/WS281/chrom.sizes"
for file in $out_dir/sj/*/*.bed
do
  bedToBigBed $file $chr_sizes ${file%.bed}.bb
done


# Coverage processing
# mkdir $out_dir/coverage $out_dir/coverage/raw_RLEs $out_dir/coverage/single_sample
# mkdir $out_dir/coverage/single_neuron $out_dir/coverage/global



# Prepare for export
cd data/outs/
tar czf ${out_version}_browser.tar.gz ${out_version}_browser/*/*/*.bb #${out_version}_browser/*/*/*.bw

echo "Send to vps with:"
echo "scp $(pwd)/${out_version}_browser.tar.gz cengen-vps:/var/www/public_data/splicing"
