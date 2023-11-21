# Visualization of raw alignments: take in bams, make coverage (Rle) objects from them,
# save the Rle as well as useful visualizations: a mean per neuron type, and a min and max across
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


cat("Starting: ", date(), "\n")

# check arguments ----
args <- commandArgs(TRUE)

WS <- args[[1]]
output_dir <- args[[2]]
in_dir <- args[[3]]

if(! WS %in% 230:300){
  stop("WS not recognized: ", WS)
}

stopifnot(dir.exists(in_dir))

outliers_to_ignore_file <- args[[4]]
stopifnot(file.exists(outliers_to_ignore_file))
outliers_to_ignore <- read_lines(outliers_to_ignore_file)

cat("Arguments, WS", WS, ", output dir ", output_dir, ", input dir: ", in_dir,
    ", ignore file: ",outliers_to_ignore, "\n")




all_covs <- tibble(path = list.files(in_dir, pattern = "\\.bam$", full.names = TRUE),
                   sample = stringr::str_split_fixed(basename(path), "\\.", 2)[,1],
                   neuron = stringr::str_split_fixed(sample, "r", 2)[,1],
                   replicate = stringr::str_split_fixed(sample, "r", 2)[,2]) %>%
  filter(! sample %in% outliers_to_ignore)




# Convert bams to RLEs on disk ----

plan(multicore, workers = 7)

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
                    qs::qsave(cur_rle,
                            paste0(output_dir,"raw_RLEs/",.y,".qs"))
                  }
                })
cat("----toc: ", proc.time()[["elapsed"]] - tic,"\n\n"); rm(tic)

plan(sequential)

# Load RLEs ----
cat("Read RLEs.\n")
all_covs <- all_covs %>%
  mutate(coverage = map(sample, ~ qs::qread(paste0(output_dir,"raw_RLEs/",.x,".qs"))))


# Write single samples BW ----
cat("Write the single samples\n")
tic <- proc.time()[["elapsed"]]
walk2(all_covs$sample, all_covs$coverage,
      ~ rtracklayer::export.bw(object = .y,
                               con = file.path(output_dir,
                                               "single_sample",
                                               paste0(.x, "_cov.bw"))))
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
                               con = file.path(output_dir, "single_neuron", paste0(.x, "_cov.bw"))))
cat("----toc: ", proc.time()[["elapsed"]] - tic,"\n\n"); rm(tic)



# Global metrics (all samples) ----

cat("Write the global metrics\n")
tic <- proc.time()[["elapsed"]]

global_mean <- pmean(reduced_cov$mean_coverage)
rtracklayer::export.bw(global_mean, file.path(output_dir, "global", "mean_cov.bw"))

min_cov <- do.call(pmin, reduced_cov$mean_coverage)
rtracklayer::export.bw(min_cov, file.path(output_dir, "global", "min_cov.bw"))

max_cov <- do.call(pmax, reduced_cov$mean_coverage)
rtracklayer::export.bw(max_cov, file.path(output_dir, "global", "max_cov.bw"))

cat("----global done: ", proc.time()[["elapsed"]] - tic,"\n\n")



cat("\n\nAll done.\n\nAny warnings:\n")

warnings()

cat("\n ---- \n")
