# Visualization of raw alignments: take in bams, make coverage (Rle) objects from them,
# save the Rle as well as useful vizualisations: a mean per neuron type, and a min and max across
# neuron type.


# Initializations ----
library(GenomicAlignments)
library(tidyverse)


bam_dir <- "/home/aw853/scratch60/2021-08-18_alignments"
output_dir <- "outputs_visualization"

cat("Read in all files\n")
all_covs <- tibble(path = list.files(bam_dir, pattern = "\\.bam$"),
                   sample = stringr::str_split_fixed(path, "\\.", 2)[,1],
                   neuron = stringr::str_split_fixed(sample, "r", 2)[,1],
                   replicate = stringr::str_split_fixed(sample, "r", 2)[,2]) %>%
  mutate(coverage = map(path, ~ coverage(readGAlignmentPairs(.x))))

cat("Write the single neurons\n")
all_covs %>%
  walk2(sample, coverage,
        ~ rtracklayer::export.bw(object = .y,
                                 con = file.path(output_dir, "single_neur", paste0(.x, ".bw"))))

cat("Write the averages per neuron class\n")
pmean <- function(list_of_covs){
  reduce(list_of_covs, `+`) / length(list_of_covs)
}

reduced_cov <- all_covs %>%
  group_by(neuron) %>%
  summarize(mean_coverage = pmean(coverage))

reduced_cov %>%
  walk2(neuron, mean_coverage,
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

