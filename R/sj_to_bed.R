# From the SJ.out.tab file output by STAR, to a bigbed file to display in browser.
# Meant to run on cluster.

# Initializations ----


cat("Starting: ", date(),"\n")


#~ packages ----
suppressPackageStartupMessages({
  library(tidyverse)
})


#~ check arguments ----
args <- commandArgs(TRUE)

WS <- args[[1]]
out_version <- args[[2]]

if(! WS %in% 230:300){
  stop("WS not recognized: ", WS)
}

outliers_to_ignore_file <- args[[3]]
stopifnot(file.exists(outliers_to_ignore_file))
outliers_to_ignore <- read_lines(outliers_to_ignore_file)


cat("Arguments, WS", WS, ", version ", out_version,", ignore file: ",
    outliers_to_ignore, "\n")


# Parameters
max_sj_length <- 25000

min_reads_min <- 1
min_reads_max <- 12
min_reads_sum <- 20




# Chromosome sizes ----

ref_cache <- paste0("/home/aw853/project/references/WS", WS)
path_chr_sizes <- file.path(ref_cache, "chrom.sizes")


cat("Init chr sizes.\n")
if(file.exists(path_chr_sizes)){
  chr_sizes <- read.delim(path_chr_sizes,
                          header = FALSE,
                          col.names=c("name","size"))
} else{
  pp <- wbData::wb_get_genome_path(WS,
                                   dir_cache = ref_cache)
  
  genome <- Biostrings::readDNAStringSet(pp)
  
  chr_sizes <- data.frame(name = names(genome),
                          size = Biostrings::width(genome))
  
  write.table(chr_sizes, path_chr_sizes,
              row.names = FALSE, col.names = FALSE, quote = FALSE,sep = "\t")
  rm(genome)
}

seqlengths <- setNames(chr_sizes$size, chr_sizes$name)






# Read SJ files ----


# ~ directories ----

sj_dir <- "/SAY/standard/mh588-CC1100-MEDGEN/bulk_alignments/bsn9_junctions/"

output_dir <- paste0("data/outs/",out_version,"_browser/sj/")


#~ functions ----

read_sj_file <- function(cur_path){
  read_table(cur_path,
             col_names=c("chr", "start","end","strand","motif", "annotated",
                         "count_unique","count_multimap", "max_overhang"),
             col_types = cols(
               chr = col_factor(levels = c("I","II","III","IV","MtDNA", "V", "X")),
               start = col_integer(),
               end = col_integer(),
               strand = col_factor(levels = as.character(0:2)),
               motif = col_factor(levels = as.character(0:6)),
               annotated = col_factor(levels = as.character(0:1)),
               count_unique = col_integer(),
               count_multimap = col_integer(),
               max_overhang = col_integer()
             )) %>%
    mutate(strand = recode_factor(strand,
                                  `0`="*",
                                  `1`="+",
                                  `2`="-")) %>%
    select(-max_overhang, -count_multimap)
}



# fn must be a function operating on the rows of a matrix (e.g. rowSums, rowMaxs, etc...)
combine_sj <- function(sj_file, fn){
  reduce(sj_file, full_join, by = c("chr", "start", "end", "strand", "motif", "annotated")) %>%
    mutate(across(starts_with("count"), replace_na, 0)) %>%
    mutate(count_unique = fn(as.matrix(select(., starts_with("count_unique"))))) %>%
    select(-starts_with("count_unique.")) %>%
    arrange(chr, start, end)
}


write_bed <- function(sj_file, out_path, min_reads = 2){
  
  gr <- sj_file %>%
    filter(count_unique > min_reads) %>%
    arrange(chr, start, end) %>%
    mutate(name = format(count_unique, big.mark = ","),
           score = round(1000*count_unique/max(count_unique))) %>%
    GenomicRanges::GRanges(seqlengths = seqlengths)
  
  rtracklayer::export(gr, out_path)
}

filter_sj <- function(sj_tibble, max_sj_length){
  sj_tibble |>
    mutate(width = end - start + 1) |>
    filter(width < max_sj_length) |>
    select(-width)
}


#~ read individual samples ----
cat("Read individual files.\n")

all_files <- tibble(path = list.files(sj_dir, full.names = TRUE),
                    replicate = stringr::str_split_fixed(basename(path), "\\.", 2)[,1],
                    sample = stringr::str_split_fixed(replicate, "t", 2)[,1]) %>%
  filter(! sample %in% outliers_to_ignore) %>%
  mutate(sj_file = map(path, read_sj_file)) %>%
  group_by(sample) %>%
  summarize(sj_file_combined = list(combine_sj(sj_file, rowSums))) %>%
  mutate(sj_file_combined = map(sj_file_combined, filter_sj, max_sj_length)) %>%
  mutate(out_path = paste0(output_dir, "single_sample/", sample, "_sj.bed"))



#~ save individual samples ----

walk2(all_files$sj_file_combined, all_files$out_path, write_bed)


#~ neuron summed ----
cat("Sum neurons.\n")
all_neurons <- all_files %>%
  select(-out_path) %>%
  mutate(neuron = stringr::str_split_fixed(sample, "r", 2)[,1]) %>%
  group_by(neuron) %>%
  summarize(sj_file_combined = list(combine_sj(sj_file_combined, rowSums))) %>%
  mutate(out_path = paste0(output_dir, "single_neuron/", neuron, "_sj.bed"))

walk2(all_neurons$sj_file_combined, all_neurons$out_path, write_bed)


#~ Global ----
cat("Global descriptions.\n")
# sum of all samples
all_samples <- all_neurons %>%
  ungroup() %>%
  summarize(sj_file_combined = list(combine_sj(sj_file_combined, rowSums)))

write_bed(all_samples$sj_file_combined[[1]],
          paste0(output_dir, "global/sum_sj.bed"),
          min_reads = min_reads_sum)

all_samples_max <- all_neurons %>%
  ungroup() %>%
  summarize(sj_file_combined = list(combine_sj(sj_file_combined, matrixStats::rowMaxs)))

write_bed(all_samples_max$sj_file_combined[[1]],
          paste0(output_dir, "global/max_sj.bed"),
          min_reads = min_reads_max)


all_samples_min <- all_neurons %>%
  ungroup() %>%
  summarize(sj_file_combined = list(combine_sj(sj_file_combined, matrixStats::rowMins)))

write_bed(all_samples_min$sj_file_combined[[1]],
          paste0(output_dir, "global/min_sj.bed"),
          min_reads = min_reads_min)


cat("done.\n")







