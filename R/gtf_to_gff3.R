# GTF from Wormbase into a GFF3 compatible with JBrowse2
# Note the rtracklayer direct approach readGFF->export does not change the fields
# and wrongly uses GTF names instead of GFF names (e.g. transcript rather than mRNA, gene_id instead of ID etc).
# Loading into txdb first better, but gene names and IDs not what desired.


library(wbData)
library(tidyverse)

WS <- 281

gids <- wb_load_gene_ids(WS)



gtf_txdb <- GenomicFeatures::makeTxDbFromGFF(wb_get_gtf_path(WS))

gff_contents <- rtracklayer::asGFF(gtf_txdb)

# Remove the "exon" entries (only use CDS)
gff_contents <- gff_contents[gff_contents$type != "exon",]


# Correct it to have the fields we want
gff_contents$ID <- sapply(strsplit(gff_contents$ID, ":"), \(x) x[2])
gff_contents$Parent <- S4Vectors::endoapply(gff_contents$Parent, \(.x) gsub("TxID:|GeneID:", "", .x))
gff_contents$Name[gff_contents$type == "gene"] <- i2s(gff_contents$Name[gff_contents$type == "gene"],
                                                      gids,
                                                      warn_missing = TRUE)




rtracklayer::export(gff_contents, paste0("data/intermediates/c_elegans.WS",WS,".canonical_geneset.gff3"), format = "gff3")

