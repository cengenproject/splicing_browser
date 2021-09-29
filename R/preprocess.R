# Note: these are parts of scripts that were used on the cluster or on local computer

library(wbData)
library(tidyverse)
library(dasper)

# 
# options(wb_dir_cache = "/gpfs/ycga/project/ysm/hammarlund/aw853/references/WS277/")

ref <- wb_load_TxDb(277)
gids <- wb_load_gene_ids(277)

# junction_dir <- "outs/2021-09-04_jx/"
# coverage_dir <- "outs/210824_coverage_bw/single_sample/"
# 
# plot_out_dir <- "outs/210905_sashimi_plots"




plot_out_dir <- "data/outs/210905_sashimi_plots"

genes_of_interest <- s2i(c("nrx-1", "unc-2", "unc-64", "unc-44", "slo-1", "slo-2",
                           "inx-1", "inx-18", "gar-1", "gar-2", "unc-32", "unc-13",
                           "unc-43", "tax-6", "dyn-1", "snt-1", "cla-1", "ric-4"),
                         gids)




# Load ----

samples_table <- tibble(samp_id = list.files(path = junction_dir,
                                             pattern = ".*\\.all_jxs\\.tsv\\.sjout") %>%
                          str_match("([A-Zef0-9]{2,4}r[0-9]{1,3})\\.all_jxs\\.tsv\\.sjout") %>%
                          .[,2])

stopifnot(all.equal(samples_table$samp_id,
                    list.files(coverage_dir) %>%
                      str_match("([A-Zef0-9]{2,4}r[0-9]{1,3})\\.bw") %>%
                      .[,2]))


samples_table <- samples_table %>%
  mutate(junction_path = paste0(junction_dir, samp_id,".all_jxs.tsv.sjout"),
         coverage_path = paste0(coverage_dir, samp_id,".bw"),
         neuron = str_match(samp_id, "^([A-Zef0-9]{2,4})r[0-9]{1,3}")[,2])


junctions <- junction_load(junction_paths = samples_table$junction_path,
                           metadata = samples_table) %>%
  junction_annot(ref) %>%
  junction_filter(count_thresh = c(raw = 6),
                  n_samp = c(raw = 3)) %>%
  junction_norm()

saveRDS(junctions, "data/intermediates/210906_junctions.rds")




junctions <- readRDS("data/intermediates/210906_junctions.rds")


# Generate all plots ----
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





# Global analysis of splicing ----

exons <- wb_load_exon_coords(277)

ex_per_gene <- exons |>
  group_by(gene_id) |>
  summarize(nb_exons = n(),
            biotype = first(gene_biotype))

ggplot(ex_per_gene) +
  geom_histogram(aes(x = nb_exons), bins = 100, color = "white") +
  theme_classic() +
  scale_x_continuous(trans = "log1p") +
  scale_y_continuous(trans = "log1p")


plot(eulerr::euler(table(ex_per_gene$nb_exons > 1, ex_per_gene$biotype == "protein_coding")),
     labels = c("Multi-exon","Protein-coding"))


table(ex_per_gene$biotype[ex_per_gene$nb_exons > 1])
mean(ex_per_gene$nb_exons)
ex_per_gene[which.max(ex_per_gene$nb_exons),]
