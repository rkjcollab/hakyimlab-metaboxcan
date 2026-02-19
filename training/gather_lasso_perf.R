#! /usr/bin/env Rscript
# the script gathers the r2 from run_gw_lasso.R output into a single file

# load modules
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(argparse)))

# SETUP ARGUMENTS
parser <- ArgumentParser()
parser$add_argument('--input_dir')
parser$add_argument('--output_prefix')


args <- parser$parse_args()

# set input dir and prefix
indir <- args$input_dir
outfile_pred <- paste0(args$output_prefix,"_r2_performance.txt")

# create an empty data frame
results <- data.frame(R2 = numeric(),
                      Pearson = numeric(),
                      Spearman = numeric(),
                      phenotype = character())

files <- list.files(path = indir, pattern = "\\.performance.tsv$")

for(file in files){
  infile <- (paste0(indir,"/",file))
  batch_res <- read.table(infile,header = T, sep = "\t")
  results <- rbind(results, batch_res)
}

# Defined functions
write_output <- function(results, fp) {
  write_tsv(results, fp)
}

# write output
write_output(results, outfile_pred)
