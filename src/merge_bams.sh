#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=merge_index_bam
#SBATCH -c 10
#SBATCH --mem=25G
#SBATCH --time=5-00:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu


echo "Do not use this script, prefer the updated script in stringtie_quantif/src/prep_alignments.sh which also removes multimappers."
exit 1


# Transfer the bam files merging the technical replicates. Index the results (bai)

module load SAMtools


alig_dir_orig=$1
alig_dir=$2


mapfile -t sampleList < <(ls $alig_dir_orig/*.bam \
                          | xargs basename -a -s .bam \
                          | cut -f1 -d't'\
                          | uniq)

echo ${#sampleList[@]}" samples"
  
  
echo "Merging."
for sample in ${sampleList[@]}
do
  samtools merge -@ $SLURM_CPUS_PER_TASK $alig_dir/$sample".bam" $(echo $alig_dir_orig/$sample"*.bam")
done


echo "BAM files merged. Indexing."
for sample in ${sampleList[@]}
do
  if [ -f "$alig_dir/$sample.bam.bai" ]; then
    echo "$sample already indexed."
  else 
    echo "Indexing $sample"
    samtools index -@ $SLURM_CPUS_PER_TASK $alig_dir/$sample".bam" $alig_dir/$sample".bam.bai"
  fi
done


echo
echo "All done."
