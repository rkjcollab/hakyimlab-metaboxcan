#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(RSQLite)
library(yaml)
library(argparse)

# SETUP ARGUMENTS
parser <- ArgumentParser()
parser$add_argument('--input_dir')
parser$add_argument('--bim_file')
parser$add_argument('--r2_file')
parser$add_argument('--out_db')

args <- parser$parse_args()

### FUNCTIONS
load_preds <- function(fp){
  predictions <- fread(fp,sep = "\t", header = T)
  class(predictions) <- "data.frame"
  return(predictions)
}

# Create an empty data frame
weights <- data.frame(snpid = character(),
                      weight = numeric(),
                      REF = character(),
                      ALT = character(),
                      CHR = numeric(),
                      phenotype = character())

indir <- args$input_dir
files <- list.files(path = args$input_dir, pattern = "\\.weights.tsv.gz$")

for (file in files){
  #print(i)
  infile <- paste0(args$input_dir,"/",file)
  batch_res <- fread(infile,header = T, sep = "\t")
  weights <- rbind(weights, batch_res)
}

weights <- weights %>%
  dplyr::rename(rsid=snpid, ref_allele=REF, eff_allele=ALT, gene=phenotype, chr=CHR)

# bim file
bim_file <- args$bim_file
bim_info <- fread(bim_file, header = F, sep = "\t", data.table=FALSE)
bim_info <- bim_info %>% dplyr::rename(pos=V4, rsid=V2)

# Add chromosome position and generate varID
weights <- weights %>%
    inner_join(bim_info[,c(2,4)], by = "rsid") %>%
    mutate(varID = paste0("chr",chr,"_",pos,"_",ref_allele,"_",eff_allele,"_b38")) %>%
    select(rsid,gene,weight,ref_allele,eff_allele,varID)

# Count snps in each model
snp_counts <-plyr::count(weights, "gene")

# extras table
# prediction performance
resid_pred <- args$r2_file

preds_df <- load_preds(resid_pred)

#merge the datasets
extra_df <- preds_df %>% 
   inner_join(snp_counts, by = c("phenotype" = "gene")) %>%
   #filter(Spearman > 0.1 & Variance >  0.01) %>% # filter metabolites
   rename(gene = phenotype, pred.perf.R2=R2,n.snps.in.model=freq) %>%
   mutate(genename = gene) %>%
   mutate(pred.perf.qval = NA) %>% mutate(pred.perf.pval = NA) %>%
   select(gene,genename,pred.perf.R2,pred.perf.pval,pred.perf.qval,n.snps.in.model)

# Filter weights df
# filter weights to include heritable metabs only 
weights <- weights %>% 
  filter(gene %in% extra_df$gene)

# sanity check
head(weights)
length(unique(weights$gene))
head(extra_df)
length(unique(extra_df$gene))

## Create database
out_db <- args$out_db
# create a db connection
driver <- dbDriver('SQLite')
# Create tables
conn <- dbConnect(drv = driver, out_db)
dbWriteTable(conn, 'extra', extra_df, overwrite = TRUE)
dbExecute(conn, "CREATE INDEX metab_model_summary ON extra (gene)")

dbWriteTable(conn, 'weights', weights, overwrite = TRUE)
dbExecute(conn, "CREATE INDEX weights_snpid ON weights (varID)")
dbExecute(conn, "CREATE INDEX weights_metabolite ON weights (gene)")
dbExecute(conn, "CREATE INDEX weights_snpid_metab ON weights (varID, gene)")
dbDisconnect(conn)


