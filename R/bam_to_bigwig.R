# Visualization of raw alignments: take in bams, make coverage (Rle) objects from them,
# save the Rle as well as useful vizualisations: a mean per neuron type, and a min and max across
# neuron type.


# Initializations ----
suppressPackageStartupMessages({
  library(GenomicAlignments)
  library(tidyverse)
  library(furrr)
  library(Rcpp)
})

source("R/Rle_utils.R",
       echo = FALSE)


cat("Starting\n\n")


bam_dir <- "/home/aw853/scratch60/2021-11-08_alignments"

output_dir <- "outs/211109_coverage_bw/"




all_covs <- tibble(path = list.files(bam_dir, pattern = "\\.bam$", full.names = TRUE),
                   sample = stringr::str_split_fixed(basename(path), "\\.", 2)[,1],
                   neuron = stringr::str_split_fixed(sample, "r", 2)[,1],
                   replicate = stringr::str_split_fixed(sample, "r", 2)[,2])




# Convert bams to RLEs on disk ----

plan(multicore, workers = 6)

cat("Read bams and save as RLE.\n")
tic <- proc.time()[["elapsed"]]
future_walk2(all_covs$path,
             all_covs$sample,
             ~ {
                  if(file.exists(paste0(output_dir,"raw_RLEs/",.y,".rds"))){
                    cat("      already exists: ",.y,"\n")
                  } else{
                    cat("      treating: ", .y,"\n")
                    cur_bam <- readGAlignmentPairs(.x)
                    cur_rle <- coverage(cur_bam)
                    chrom_names <- names(cur_rle)
                    cur_rle <- 1e6 * cur_rle / length(cur_bam) # as CPM
                    names(cur_rle) <- chrom_names
                    saveRDS(cur_rle,
                            paste0(output_dir,"raw_RLEs/",.y,".rds"))
                  }
                })
cat("----toc: ", proc.time()[["elapsed"]] - tic,"\n\n"); rm(tic)

plan(sequential)

# Load RLEs ----
cat("Read RLEs.\n")
all_covs <- all_covs %>%
  slice_head(n=4)%>%
  mutate(coverage = map(sample, ~ readRDS(paste0(output_dir,"raw_RLEs/",.x,".rds"))))


# Write single samples BW ----
cat("Write the single samples\n")
tic <- proc.time()[["elapsed"]]
walk2(all_covs$sample, all_covs$coverage,
      ~ rtracklayer::export.bw(object = .y,
                               con = file.path(output_dir,
                                               "single_sample",
                                               paste0(.x, ".bw"))))
cat("----toc: ", proc.time()[["elapsed"]] - tic,"\n\n"); rm(tic)



# Per neuron average ----


cat("Compute and write the averages per neuron class\n")
pmean <- function(list_of_covs){
  purrr::reduce(list_of_covs, `+`) / length(list_of_covs)
}

tic <- proc.time()[["elapsed"]]
reduced_cov <- all_covs %>%
  group_by(neuron) %>%
  summarize(mean_coverage = list(pmean(coverage)))

# save to disk
walk2(reduced_cov$neuron, reduced_cov$mean_coverage,
      ~ rtracklayer::export.bw(object = .y,
                               con = file.path(output_dir, "means", paste0("mean_",.x, ".bw"))))
cat("----toc: ", proc.time()[["elapsed"]] - tic,"\n\n"); rm(tic)



# Global metrics (all samples) ----

cat("Write the global metrics\n")
tic <- proc.time()[["elapsed"]]

global_mean <- pmean(reduced_cov$mean_coverage)
rtracklayer::export.bw(global_mean, file.path(output_dir, "global", "mean.bw"))
cat("----mean done: ", proc.time()[["elapsed"]] - tic,"\n\n")

global_median <- do.call(pmedian, reduced_cov$mean_coverage)
rtracklayer::export.bw(global_median, file.path(output_dir, "global", "median.bw"))
cat("----median done: ", proc.time()[["elapsed"]] - tic,"\n\n")

global_lower <- do.call(plower, reduced_cov$mean_coverage)
rtracklayer::export.bw(global_lower, file.path(output_dir, "global", "lower.bw"))
cat("----lower done: ", proc.time()[["elapsed"]] - tic,"\n\n")

global_higher <- do.call(phigher, reduced_cov$mean_coverage)
rtracklayer::export.bw(global_mean, file.path(output_dir, "global", "higher.bw"))
cat("----higher done: ", proc.time()[["elapsed"]] - tic,"\n\n")



cat("\n\nAll done.\n")

