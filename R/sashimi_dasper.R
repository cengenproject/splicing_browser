# To execute on cluster using wrapper.
# Produce sashimi plots of each gene in each sample


# Init ----
library(tidyverse)
library(dasper)
library(wbData)

options(wb_dir_cache = "/gpfs/ycga/project/ysm/hammarlund/aw853/references/WS277/")

ref <- wb_load_TxDb(277)
gids <- wb_load_gene_ids(277)

junction_dir <- "outs/2021-09-04_jx/"
coverage_dir <- "outs/210824_coverage_bw/single_sample/"

plot_out_dir <- "outs/210905_sashimi_plots"

genes_of_interest <- s2i(c("nrx-1", "unc-2", "unc-64", "unc-44", "slo-1", "slo-2",
                           "inx-1", "inx-18", "gar-1", "gar-2", "unc-32", "unc-13",
                           "unc-43", "tax-6", "dyn-1", "snt-1", "cla-1", "ric-4"),
                         gids)


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


samples_table <- tibble(samp_id = list.files(path = junction_dir,
                                               pattern = ".*\\.all_jxs\\.tsv\\.sjout") %>%
                          str_match("([A-Zef0-9]{2,4}r[0-9]{1,3})\\.all_jxs\\.tsv\\.sjout") %>%
                          .[,2])

stopifnot(all.equal(samples_table$samp_id,
                    list.files(coverage_dir) %>%
                      str_match("([A-Zef0-9]{2,4}r[0-9]{1,3})\\.bw") %>%
                      .[,2]))


samples_table <- samples_table %>%
  slice_head(n=4) %>%
  mutate(junction_path = paste0(junction_dir, samp_id,".all_jxs.tsv.sjout"),
         coverage_path = paste0(coverage_dir, samp_id,".bw"),
         neuron = str_match(samp_id, "^([A-Zef0-9]{2,4})r[0-9]{1,3}")[,2])

junctions <- junction_load(junction_paths = samples_table$junction_path,
                           metadata = samples_table) %>%
  junction_annot(ref) %>%
  junction_filter() %>%
  junction_norm()

all_plots <- expand_grid(samp_id = samples_table$samp_id,
            gene_id = genes_of_interest) %>%
  left_join(samples_table, by = "samp_id") %>%
  mutate(sashimi_plot = 
           map2(coverage_path, gene_id,
                ~plot_sashimi(junctions,
                              coverage_paths_case = .x,
                              ref,
                              gene_tx_id = .y,
                              count_label = TRUE,
                              assay_name = "raw",
                              gene_tx_func = gene_tx_wormbase,
                              include_control = FALSE
                )))
         
pwalk(list(all_plots$samp_id, all_plots$gene_id, all_plots$sashimi_plot),
        function(samp_id,gene_id,.p) ggsave(paste0(samp_id, ".png"),
                                         plot = .p, path = plot_out_dir))

bb <- plot_sashimi(
  jctions,
  coverage_paths_case = "data/DAr89.bw",
  ref,
  # case_id = list(samp_id = "samp_1"),
  gene_tx_id = c("WBGene00006798"),
  count_label = TRUE,
  assay_name = "raw",
  gene_tx_func = gene_tx_wormbase,
  include_control = FALSE
)
ggsave("test_out.png", plot = bb)
