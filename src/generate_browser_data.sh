#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=gen_browser_data
#SBATCH -c 7
#SBATCH --mem=110G
#SBATCH --time=18:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu

set -e

# this script takes the output from a "bsn" alignment pipeline, and calls R scripts to
# generate the data to visualize in the genome browser.

# parameters
WS="281"
out_version="220322"
bams_combined="/home/aw853/scratch60/2022-03-18_alignments"
bams_orig="/SAY/standard/mh588-CC1100-MEDGEN/bulk_alignments/bsn9_bams/"

outliers_to_ignore="data/outliers_to_ignore.txt"


# Initializations

echo "---------------  Starting $(date) ---------------------"

module load R

out_dir=data/outs/${out_version}_browser

mkdir -p $out_dir





## Junction processing ----

mkdir $out_dir/sj $out_dir/sj/single_sample $out_dir/sj/single_neuron $out_dir/sj/global

chr_sizes="/home/aw853/project/references/WS281/chrom.sizes"


Rscript R/sj_to_bed.R -w $WS -c $chr_sizes -o $out_version -i $outliers_to_ignore

# convert to BigBed
for file in $out_dir/sj/*/*.bed
do
  bedToBigBed $file $chr_sizes ${file%.bed}.bb
done




## Coverage processing

# merge the bams

if [ ! -d $bams_combined ]
then
  mkdir $bams_combined
  srun src/merge_bams.sh $bams_orig $bams_combined
elif [ -z "$(ls -A $bams_combined)" ]
then
  srun src/merge_bams.sh $bams_orig $bams_combined
else
   echo "bams already combined."
fi



mkdir $out_dir/coverage $out_dir/coverage/raw_RLEs $out_dir/coverage/single_sample
mkdir $out_dir/coverage/single_neuron $out_dir/coverage/global

Rscript R/bam_to_bigwig.R $WS $out_version $bams_combined $outliers_to_ignore



# Prepare for export
cd data/outs/
tar czf ${out_version}_browser.tar.gz ${out_version}_browser/*/*/*.bb ${out_version}_browser/*/*/*.bw

echo "Send to vps with:"
echo "scp $(pwd)/${out_version}_browser.tar.gz cengen-vps:/var/www/public_data/splicing"

