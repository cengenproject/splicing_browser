library(wbData)
library(GenomicFeatures)

## Import ----
wb_gtf <- rtracklayer::import(wb_get_gtf_path(277))


## // Not working: save with rtracklayer function ----

# wb_gtf$score <- 1000
# 
# wb_fasta <- Biostrings::readDNAStringSet(wb_get_genome_path(277))
# wb_all_seqlengths <- Biostrings::width(wb_fasta) |>
#   setNames(names(wb_fasta))
# seqlengths(wb_gtf) <- wb_all_seqlengths[names(seqlengths(wb_gtf))]
# 
# export(wb_gtf, "data/WS277_canonical_geneset.bed")



## Split into separate GRangesLists ----
grList_by_tx <- split(wb_gtf, wb_gtf$transcript_id)

all_exons <- endoapply(grList_by_tx,
                       \(gr) gr[gr$type == "exon"])
all_transcripts <- endoapply(grList_by_tx,
                             \(gr) gr[gr$type == "transcript"])
# all_UTR <- endoapply(grList_by_tx,
#                              \(gr) gr[gr$type %in% c("three_prime_utr","five_prime_utr")])
all_start_codons <- endoapply(grList_by_tx,
                              \(gr) gr[gr$type == "start_codon"])
all_stop_codons <- endoapply(grList_by_tx,
                              \(gr) gr[gr$type == "stop_codon"])

# save(all_exons,all_transcripts,all_stop_codons,
#      file = "data/tmp_grangelists_for_gtf2bigbed.rda")

## Compute BED fields ----

chrom <- vapply(seqnames(all_transcripts), \(.x) as.character(.x@values), character(1L))
chromStart <- as.integer(start(all_transcripts) -1)
chromEnd <- as.integer(end(all_transcripts))
name <- names(all_transcripts)
score <- rep(1000, length(all_transcripts))
strand <- vapply(strand(all_transcripts), \(.x) as.character(.x@values), character(1L))

thickStart <- as.integer(start(all_start_codons))
thickEnd <- as.integer(end(all_stop_codons))
itemRgb <- rep(".", length(all_transcripts))

# block_counts <- vapply(all_exons, length, integer(1L))
blockCount <- integer(length(all_transcripts))
blockStarts <- character(length(all_transcripts))
blockSizes <- character(length(all_transcripts))
for(cur_transcr in seq_along(all_exons)){
  blockCount[[cur_transcr]] <- length(all_exons[[cur_transcr]])
  blockStarts[[cur_transcr]] <- paste0(start(ranges(all_exons[[cur_transcr]])) - chromStart[[cur_transcr]]
                                       , collapse=",")
  blockSizes[[cur_transcr]] <- paste0(width(ranges(all_exons[[cur_transcr]])), collapse=",")
}


## Check results ----
lapply(list(chrom, chromStart, chromEnd, name,score,strand,thickStart,
            thickEnd,itemRgb,blockCount,blockSizes,blockStarts),
       length) |>
  purrr::reduce(identical)



## Write to disk ----

paste(chrom, chromStart, chromEnd, name,score,strand,thickStart,
      thickEnd,itemRgb,blockCount,blockSizes,blockStarts,
      sep = "\t") |>
  readr::write_lines("data/WS277_canonical_geneset.bed12")

## For bigBed ----

wb_fasta <- Biostrings::readDNAStringSet(wb_get_genome_path(277))
wb_all_seqlengths <- data.frame(names(wb_fasta),
                                Biostrings::width(wb_fasta))
rm(wb_fasta)

readr::write_tsv(wb_all_seqlengths,
                 file = "data/chrom.WS277.sizes",
                 col_names = FALSE)
