# Visualization of raw alignments: take in bams, make coverage (Rle) objects from them,
# save the Rle as well as useful vizualisations: a mean per neuron type, and a min and max across
# neuron type.


# Initializations ----
library(GenomicAlignments)
library(tidyverse)
library(furrr)
library(tictoc)




bam_dir <- "/home/aw853/scratch60/2021-08-18_alignments"

output_dir <- "outputs_visualization"




all_covs <- tibble(path = list.files(bam_dir, pattern = "\\.bam$", full.names = TRUE),
                   sample = stringr::str_split_fixed(basename(path), "\\.", 2)[,1],
                   neuron = stringr::str_split_fixed(sample, "r", 2)[,1],
                   replicate = stringr::str_split_fixed(sample, "r", 2)[,2])

## Already saved
# plan(multicore, workers = 6)
# 
# cat("Read bams and save as RLE.\n")
# tic()
# future_walk2(all_covs$path,
#              all_covs$sample,
#              ~ saveRDS(coverage(readGAlignmentPairs(.x)),
#                        paste0(output_dir,"/raw_RLEs/",.y,".rds")))
# toc()

cat("Read RLEs.\n")
all_covs <- all_covs %>%
  mutate(coverage = map(sample, ~ readRDS(paste0(output_dir,"/raw_RLEs/",.x,".rds"))))


# cat("Write the single neurons\n")
# walk2(all_covs$sample, all_covs$coverage,
#       ~ rtracklayer::export.bw(object = .y,
#                                con = file.path(output_dir, "single_neur", paste0(.x, ".bw"))))


cat("Write the averages per neuron class\n")
pmean <- function(list_of_covs){
  reduce(list_of_covs, `+`) / length(list_of_covs)
}

reduced_cov <- all_covs %>%
  group_by(neuron) %>%
  summarize(mean_coverage = list(pmean(coverage)))

walk2(reduced_cov$neuron, reduced_cov$mean_coverage,
      ~ rtracklayer::export.bw(object = .y,
                               con = file.path(output_dir, "means", paste0("mean_",.x, ".bw"))))


cat("Write the global min, max, and average\n")
reduced_cov %>%
  pull(mean_coverage) %>%
  pmin() %>%
  rtracklayer::export.bw(file.path(output_dir, "global", "minimum.bw"))

reduced_cov %>%
  pull(mean_coverage) %>%
  pmax() %>%
  rtracklayer::export.bw(file.path(output_dir, "global", "maximum.bw"))

reduced_cov %>%
  pull(mean_coverage) %>%
  pmean() %>%
  rtracklayer::export.bw(file.path(output_dir, "global", "mean.bw"))


cat("\n\nAll done.\n")

