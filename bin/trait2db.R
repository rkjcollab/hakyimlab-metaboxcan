#! /usr/bin/env Rscript

# Load packages
suppressMessages(library(argparse))
suppressMessages(library(RSQLite))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

# SETUP ARGUMENTS
parser <- ArgumentParser()
parser$add_argument('--input_file')
parser$add_argument('--output_db')

args <- parser$parse_args()

# set input and output
infile = args$input_file
out_db = args$output_db

# Load the trait
trait_ss <- fread(infile, data.table = FALSE)

if (("rsid" %in% names(trait_ss)) & (! "variant_id" %in% names(trait_ss))) {
   trait_ss <- trait_ss %>% dplyr::rename(variant_id = rsid)
}
# create db
driver <- dbDriver('SQLite')
outconn <- dbConnect(drv = driver, out_db)
dbWriteTable(outconn, 'summary_stats', trait_ss, overwrite = TRUE)

# Create indes on chromosome and position
dbExecute(outconn, "CREATE INDEX chr ON summary_stats (chromosome)")
dbExecute(outconn, "CREATE INDEX pos ON summary_stats (position)")
dbExecute(outconn, "CREATE INDEX chr_pos ON summary_stats (chromosome,position)")
dbDisconnect(outconn)

