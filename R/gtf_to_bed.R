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

if(file.exists( "outputs_visualization/references/tmp_grangelists_for_gtf2bigbed.rda")){
  cat("Loading cahed file...\n")
  load("outputs_visualization/references/tmp_grangelists_for_gtf2bigbed.rda")
} else{
  cat("split\n")
  grList_by_tx <- split(wb_gtf, wb_gtf$transcript_id)
  
  cat("filter GRangeLists\n")
  all_exons <- endoapply(grList_by_tx,
                         function(gr) gr[gr$type == "exon"])
  all_transcripts <- endoapply(grList_by_tx,
                               function(gr) gr[gr$type == "transcript"])
  # all_UTR <- endoapply(grList_by_tx,
  #                              \(gr) gr[gr$type %in% c("three_prime_utr","five_prime_utr")])
  all_start_codons <- endoapply(grList_by_tx,
                                function(gr) gr[gr$type == "start_codon"])
  all_stop_codons <- endoapply(grList_by_tx,
                               function(gr) gr[gr$type == "stop_codon"])
  
  save(all_exons,all_transcripts,all_stop_codons,all_start_codons,
       file = "outputs_visualization/references/tmp_grangelists_for_gtf2bigbed.rda")
  
}




## Compute BED fields ----
cat("Compute BED4\n")
chrom <- vapply(seqnames(all_transcripts),
                function(.x) as.character(.x@values), character(1L))
chromStart <- as.integer(start(all_transcripts) -1)
chromEnd <- as.integer(end(all_transcripts))
name <- names(all_transcripts)

cat("Compute BED9\n")
score <- rep(1000, length(all_transcripts))
strand <- vapply(strand(all_transcripts), 
                 function(.x) as.character(.x@values), character(1L))

# start codons: problem when splicing within codon
all_stcods <- start(all_start_codons)
lengths_stcods <- vapply(all_stcods, length, integer(1L))

spliced_stcods <- which(lengths_stcods == 2)
for(stcod in spliced_stcods){
  if(strand(all_start_codons[[stcod]])@values == "+"){
    all_stcods[[stcod]] <- min(all_stcods[[stcod]])
  } else{
    all_stcods[[stcod]] <- max(all_stcods[[stcod]])
  }
}

thickStart <- as.integer(all_stcods)

# stop codons, same problem
all_stcods <- end(all_stop_codons)
lengths_stcods <- vapply(all_stcods, length, integer(1L))

spliced_stcods <- which(lengths_stcods == 2)
for(stcod in spliced_stcods){
  if(strand(all_stop_codons[[stcod]])@values == "-"){
    all_stcods[[stcod]] <- max(all_stcods[[stcod]])
  } else{
    all_stcods[[stcod]] <- min(all_stcods[[stcod]])
  }
}

thickEnd <- as.integer(all_stcods)

itemRgb <- rep(".", length(all_transcripts))

cat("Compute BED12\n")
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
# lapply(list(chrom, chromStart, chromEnd, name,score,strand,thickStart,
#             thickEnd,itemRgb,blockCount,blockSizes,blockStarts),
#        length) |>
#   purrr::reduce(identical)



## Write to disk ----
cat("Write to disk\n")
readr::write_lines(paste(chrom, chromStart, chromEnd, name,score,strand,thickStart,
                         thickEnd,itemRgb,blockCount,blockSizes,blockStarts,
                         sep = "\t"),
                   "data/WS277_canonical_geneset.bed12")


## For bigBed ----
cat("Write chrom sizes\n")
wb_fasta <- Biostrings::readDNAStringSet(wb_get_genome_path(277))
wb_all_seqlengths <- data.frame(names(wb_fasta),
                                Biostrings::width(wb_fasta))
rm(wb_fasta)

readr::write_tsv(wb_all_seqlengths,
                 file = "outputs_visualization/references/chrom.WS277.sizes",
                 col_names = FALSE)
