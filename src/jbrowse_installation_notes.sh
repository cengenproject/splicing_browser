# Run on server to install and configure JBrowse2
# note this is a collection of commands, not a script.

# Will not work on Node v10 (Debian repos)
node -v
#> v14.18.0


# Genome tools: install version >= https://github.com/genometools/genometools/commit/64ccb7bf61bf5d3bc686c5818899e81f6699a388
# Because of issue https://github.com/genometools/genometools/issues/990

# Load fasta
sudo wget ftp://wormbase.org/pub/wormbase/releases/WS277/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS277.genomic.fa.gz
gunzip c_elegans.PRJNA13758.WS277.genomic.fa.gz
samtools faidx c_elegans.PRJNA13758.WS277.genomic.fa
jbrowse add-assembly /var/www/public_data/references/c_elegans.PRJNA13758.WS277.genomic.fa --out /var/www/html/jbrowse2 --load symlink


# Upload GFF file generated with R/gtf_to_gff3
mv ~/c_elegans.WS277.canonical_geneset.gff3 /var/www/public_data/references

/opt/gt-unstable/bin/gt gff3 -sortlines -tidy -retainids c_elegans.WS277.canonical_geneset.gff3 > c_elegans.WS277.canonical_geneset.sorted.gff3
bgzip c_elegans.WS277.canonical_geneset.sorted.gff3
tabix c_elegans.WS277.canonical_geneset.sorted.gff3.gz

jbrowse add-track /var/www/public_data/references/c_elegans.WS277.canonical_geneset.sorted.gff3.gz --load symlink

# Note: use create_jbrowse_config2.sh to automate this
jbrowse add-track /var/www/public_data/splicing/210824_coverage_bw/global/mean.bw --load symlink --category=global
jbrowse add-track /var/www/public_data/splicing/210824_coverage_bw/global/maximum.bw --load symlink --category=global
jbrowse add-track /var/www/public_data/splicing/210824_coverage_bw/global/minimum.bw --load symlink --category=global
jbrowse text-index
