#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=bam_to_jx
#SBATCH -c 4
#SBATCH --mem=80G
#SBATCH --time=1-18:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu



#-------------------------------------------------------------------------------------------------------------
# Create workspace paths


# use the merged alignments in scratch60 (technical replicates are merged)
alig_dir="/home/aw853/scratch60/2021-08-18_alignments"

out_dir="outs/2021-09-04_jx"




#-------------------------------------------------------------------------------------------------------------

set -e

echo "Start megadepth quantif of junctions on $(date)"


## Check inputs

if [ ! -d $alig_dir ]
then
  echo "Error: bam directory does not exist: $alig_dir"
  exit 1
fi


if [ ! -d $out_dir ]
then
  echo "Error: destination directory does not exist. If this is a new pipeline version, you need to create it yourself: $out_dir"
  exit 1
fi




echo "--------------------------------------------------------------"
## Read sample list and remove trailing extension

echo "Reading samples from bsn5"
mapfile -t samplePath < <(ls $alig_dir/*.bam)


if [ ${#samplePath[@]} -lt 1 ]
  then
  echo "Error: failed to find samples."
  exit 1
fi
  





start_dir=$(pwd)
cd $out_dir


for((i=0; i<$nb_samples; i++))
do
  sample=$(basename -s .bam ${samplePath[$i]})
  
  
  megadepth ${samplePath[$i]} --all-junctions --prefix $sample --threads $SLURM_CPUS_PER_TASK
  
  echo
  echo
  
done
echo

cd $start_dir


echo "Done on $(date)"


