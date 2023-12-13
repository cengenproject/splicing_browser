# GTF from Wormbase into a GFF3 compatible with JBrowse2
# Note the rtracklayer direct approach readGFF->export does not change the fields
# and wrongly uses GTF names instead of GFF names (e.g. transcript rather than mRNA, gene_id instead of ID etc).
# Loading into txdb first better, but gene names and IDs not what desired.


library(wbData)
library(tidyverse)

WS <- 289

gids <- wb_load_gene_ids(WS)




gff_contents <- rtracklayer::asGFF(wb_load_TxDb(WS))



# Correct it to have the fields we want
gff_contents$ID <- sapply(strsplit(gff_contents$ID, ":"), \(x) x[2])
gff_contents$Parent <- S4Vectors::endoapply(gff_contents$Parent, \(.x) gsub("TxID:|GeneID:", "", .x))
gff_contents$Name[gff_contents$type == "gene"] <- i2s(gff_contents$Name[gff_contents$type == "gene"],
                                                      gids,
                                                      warn_missing = TRUE)


# sort
gff_contents_sorted <- GenomicRanges::sort(gff_contents, ignore.strand = TRUE)


rtracklayer::export(gff_contents_sorted, paste0("data/intermediates/231212_WS",WS,".canonical_geneset.gff3"), format = "gff3")





