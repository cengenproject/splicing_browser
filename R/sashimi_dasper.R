# To execute on cluster using wrapper.
# Produce sashimi plots of each gene in each sample


# Init ----
library(tidyverse)
library(dasper)
library(wbData)


ref <- wb_load_TxDb(277)
gids <- wb_load_gene_ids(277)


junction_dir <- "data/outs/2021-09-04_jx/"
coverage_dir <- "data/outs/210824_coverage_bw/single_sample/"


gene_tx_wormbase <- function(gene_tx_id){
  if (is.null(gene_tx_id) | length(gene_tx_id) != 1) {
    stop("gene_tx_id must be set and be of length 1")
  } else if (stringr::str_detect(gene_tx_id, "ENSG") | stringr::str_detect(gene_tx_id, "WBGene")) {
    gene_tx_type <- "gene_id"
  } else if (stringr::str_detect(gene_tx_id, "ENST")) {
    gene_tx_type <- "tx_name"
  } else {
    stop("gene_tx_id does not include an ENST or ENSG prefix")
  }
  
  # create named list of gene/tx id
  # for filtering txdb
  gene_tx_list <- list(gene_tx_id)
  names(gene_tx_list) <- gene_tx_type
  
  return(gene_tx_list)
}

junctions <- readRDS("data/intermediates/210906_junctions.rds")




# Plot ----


plot_sashimi(
  junctions,
  coverage_paths_case = "data/outs/210824_coverage_bw/single_sample/VDr81.bw",
  ref,
  case_id = list(samp_id = "VDr81"),
  gene_tx_id = s2i("vab-8", gids),
  count_label = TRUE,
  assay_name = "raw",
  gene_tx_func = gene_tx_wormbase,
  include_control = FALSE,
  binwidth = 2
  # ,
  # region = GRanges(seqnames = "X",
  #                  ranges = IRanges(start=3453942,
  #                                   end = 3462544   ),
  #                  strand = "-"),
  # region_strict = TRUE
)

plot_sashimi(
  junctions,
  coverage_paths_case = "data/outs/210824_coverage_bw/single_sample/AVGr43.bw",
  ref,
  case_id = list(samp_id = "AVGr43"),
  gene_tx_id = s2i("madd-4", gids),
  count_label = TRUE,
  assay_name = "norm",
  gene_tx_func = gene_tx_wormbase,
  include_control = FALSE,
  binwidth = 2
  # ,
  # region = GRanges(seqnames = "III",
  #                  ranges = IRanges(start=4214491,
  #                                   end = 4230252),
  #                  strand = "-"),
  # region_strict = TRUE
)



plot_sashimi(
  junctions,
  coverage_paths_case = "data/outs/210824_coverage_bw/single_sample/AIYr66.bw",
  ref,
  case_id = list(samp_id = "AIYr66"),
  gene_tx_id = s2i("daf-12", gids),
  count_label = TRUE,
  assay_name = "norm",
  gene_tx_func = gene_tx_wormbase,
  include_control = FALSE,
  binwidth = 2
  ,
  region = GRanges(seqnames = "III",
                   ranges = IRanges(start=4214491,
                                    end = 4230252),
                   strand = "-"),
  region_strict = TRUE
)

