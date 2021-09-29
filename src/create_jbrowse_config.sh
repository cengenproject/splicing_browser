#!/bin/bash



##### Initializations

set -e

ref_dir='/var/www/public_data/references'
base_dir='/var/www/public_data/splicing/210824_coverage_bw'
jbrowse_dir='/var/www/html/jbrowse2'
config_file=$jbrowse_dir'/config.json'

fasta_name='c_elegans.PRJNA13758.WS277.genomic.fa'
gff_name='c_elegans.WS277.canonical_geneset.sorted.gff3.gz'



##### Add main tracks
cur_dir=$(pwd)
cd $jbrowse_dir
jbrowse add-assembly $ref_dir/$fasta_name \
    --load symlink
cd $cur_dir



jbrowse add-track $ref_dir/$gff_name \
    --out $jbrowse_dir \
    --load symlink \
    --category=Annotations \
    --config='{"displayMode": "compact"}'


##### Add each bw track

for bw in $base_dir/global/*.bw
do
  jbrowse add-track $bw \
    --load symlink \
    --category=Global \
    --out $jbrowse_dir
done


for bw in $base_dir/means/*.bw
do
  jbrowse add-track $bw \
    --load symlink \
    --category=Neurons \
    --out $jbrowse_dir
done


for bw in $base_dir/single_sample/*.bw
do
  jbrowse add-track $bw \
    --load symlink \
    --category='Individual Samples' \
    --out $jbrowse_dir
done

##### Create text index
jbrowse text-index --out $jbrowse_dir


##### Add logo
cat $config_file | sed 's/"configuration": {},/"configuration": {\n    "logoPath": {\n      "uri": "misc\/logo_resized.svg"\n    }\n  },/g' > ${config_file}


