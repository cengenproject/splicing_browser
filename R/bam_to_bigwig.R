# Visualization of raw alignments: take in bams, make coverage (Rle) objects from them,
# save the Rle as well as useful vizualisations: a mean per neuron type, and a min and max across
# neuron type.


# Initializations ----
library(GenomicAlignments)
library(tidyverse)
# library(furrr)
# library(tictoc)




bam_dir <- "/home/aw853/scratch60/2021-08-18_alignments"

output_dir <- "outputs_visualization"




all_covs <- tibble(path = list.files(bam_dir, pattern = "\\.bam$", full.names = TRUE),
                   sample = stringr::str_split_fixed(basename(path), "\\.", 2)[,1],
                   neuron = stringr::str_split_fixed(sample, "r", 2)[,1],
                   replicate = stringr::str_split_fixed(sample, "r", 2)[,2])

# Already saved
plan(multicore, workers = 6)

cat("Read bams and save as RLE.\n")
tic <- proc.time()[["elapsed"]]
future_walk2(all_covs$path,
             all_covs$sample,
             ~ {
                  cur_bam <- readGAlignmentPairs(.x)
                  cur_rle <- coverage(cur_bam)
                  cur_rle <- 1e6 * cur_rle / length(cur_bam) # as CPM
                  saveRDS(cur_rle,
                          paste0(output_dir,"/raw_RLEs/",.y,".rds"))
                })
cat("    toc: ", proc.time()[["elapsed"]] - tic,"\n\n"); rm(tic)

# cat("Read RLEs.\n")
# all_covs <- all_covs %>%
#   mutate(coverage = map(sample, ~ readRDS(paste0(output_dir,"/raw_RLEs/",.x,".rds"))))


cat("Write the single samples\n")
tic <- proc.time()[["elapsed"]]
walk2(all_covs$sample, all_covs$coverage,
      ~ rtracklayer::export.bw(object = .y,
                               con = file.path(output_dir, "single_neur", paste0(.x, ".bw"))))
cat("    toc: ", proc.time()[["elapsed"]] - tic,"\n\n"); rm(tic)


cat("Compute and write the averages per neuron class\n")
pmean <- function(list_of_covs){
  reduce(list_of_covs, `+`) / length(list_of_covs)
}

tic <- proc.time()[["elapsed"]]
reduced_cov <- all_covs %>%
  group_by(neuron) %>%
  summarize(mean_coverage = list(pmean(coverage)))

walk2(reduced_cov$neuron, reduced_cov$mean_coverage,
      ~ rtracklayer::export.bw(object = .y,
                               con = file.path(output_dir, "means", paste0("mean_",.x, ".bw"))))
cat("    toc: ", proc.time()[["elapsed"]] - tic,"\n\n"); rm(tic)


cat("Write the global min and max\n")
tic <- proc.time()[["elapsed"]]
min_coverage <- do.call(pmin, reduced_cov$mean_coverage)
max_coverage <- do.call(pmax, reduced_cov$mean_coverage)

rtracklayer::export.bw(min_coverage, file.path(output_dir, "global", "minimum.bw"))
rtracklayer::export.bw(max_coverage, file.path(output_dir, "global", "maximum.bw"))
cat("    toc: ", proc.time()[["elapsed"]] - tic,"\n\n"); rm(tic)


cat("Write the global mean\n")
tic <- proc.time()[["elapsed"]]
partial_sum <- reduced_cov$mean_coverage[[1]]
for(cur_cov in 2:nrow(reduced_cov)){
  cat("--  ", cur_cov)
  partial_sum <- partial_sum + reduced_cov$mean_coverage[[cur_cov]]
  cat(" - toc: ", proc.time()[["elapsed"]] - tic,"\n")
}

mean_coverage <- partial_sum / nrow(reduced_cov)

# mean_coverage <- do.call(pmean, reduced_cov$mean_coverage)
rtracklayer::export.bw(mean_coverage, file.path(output_dir, "global", "mean.bw"))
cat("    toc: ", proc.time()[["elapsed"]] - tic,"\n\n"); rm(tic)

cat("\n\nAll done.\n")

