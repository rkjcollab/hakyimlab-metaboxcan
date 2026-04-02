#! /usr/bin/env bash

set -euo pipefail

usage() {
  echo "Usage: $0 -p <plink_prefix> -o <out_dir>"
  exit 1
}

while getopts ":p:o:" opt; do
  case "${opt}" in
    p) plink_prefix=${OPTARG};;
    o) out_dir=${OPTARG};;
    :) echo "Option -$OPTARG requires an argument." >&2; usage;;
    \?) echo "Invalid option: -$OPTARG" >&2; usage;;
  esac
done


#TODO: before this:
  # do we need to clean genetic data (SNP & sample missingness?, HWE?)
  # if add this step, could also add filter to HapMap3 SNPs only


#TODO: for now, adding variable paths but keeping hardcoded file names to
# enable ease of questions to Festus if needed.
#TODO: maybe build plink file name into output folder, but keep files general?
#TODO: add thread_num as variable

#TODO: revisit how want to handle environment
module load gcc/6.2.0
module load bcftools
module load plink/1.90
module load plink/2.0
module load gcta
module load bgen/1.1.3

# Function to check file exists, use for commands that could fail silently
check_file_exists() {
   if [ ! -f "$1" ]; then
      echo "Error: Required file '$1' not found." >&2
      exit 1
   fi
}

# Make output directories if needed
mkdir -p "${out_dir}/data/train"
mkdir -p "${out_dir}/data/test"

## get individual list
#TODO: need to either check input is .fam or ?
check_file_exists "${plink_prefix}.fam"
# cut -f1,2 -d " " "${plink_prefix}.fam" > "${out_dir}/data/individuals.txt"
awk '{print $1, $2}' "${plink_prefix}.fam" > "${out_dir}/data/individuals.txt"

## random shuffle with a seed to select individuals randomly
# get ids to split individuals into training and testing
shuf --random-source=<(yes 149) "${out_dir}/data/individuals.txt" | \
  head -n 400 > "${out_dir}/data/test_indiv.txt"

grep -v -x -f "${out_dir}/data/test_indiv.txt" \
  "${out_dir}/data/individuals.txt" > "${out_dir}/data/train_indiv.txt"

#TODO: here the file name includes hamap3, but code doens't actually subset
# to these variants. README says to do so, so I think I will add it?
## split the data into training and testing
# train
echo -e "\nSplitting the data\n"
plink \
  --bfile "${plink_prefix}" \
  --keep "${out_dir}/data/train_indiv.txt" \
  --make-bed \
  --out "${out_dir}/data/train/train"

# test
plink \
  --bfile "$plink_prefix" \
  --keep "${out_dir}/data/test_indiv.txt" \
  --make-bed \
  --out "${out_dir}/data/test/test"


## generate other formats for the train set
## make a grm
echo -e "\nCalculating a grm matrix\n"

gcta \
  --bfile "${out_dir}/data/train/train" \
  --make-grm \
  --thread-num 8 \
  --out "${out_dir}/data/train/train_grm"

## Make a grm gz (for ridge)
echo -e "\nCalculating a grm matrix (gz)\n"

gcta \
  --bfile "${out_dir}/data/train/train" \
  --make-grm-gz \
  --thread-num 8 \
  --out "${out_dir}/data/train/train_grm"


## Calculate pca
echo -e "\nCalculating pcas\n"

gcta \
  --grm "${out_dir}/data/train/train_grm" \
  --pca 20 \
  --thread-num 8 \
  --out "${out_dir}/data/train/train_pca"

echo -e "\nGenerate a pgen format\n"

## make pgen format (uses plink 2.0) for lasso
plink2 \
  --bfile "${out_dir}/data/train/train" \
  --make-pgen vzs \
  --out "${out_dir}/data/train/train_pgen"

# formating the psam to work 
mv "${out_dir}/data/train/train_pgen.psam" "${out_dir}/data/train/train_pgen.psam.bak"

awk -F"\t" '{print $2"\t"$2"\t"$3}' "${out_dir}/data/train/train_pgen.psam.bak" > \
  "${out_dir}/data/train/train_pgen.psam"

## Dont use underscore in names, has a problem with gw_lasso
#TODO: resolve this in PLINK command instead?
sed -i  -e 's/^IID/#FID/' -e 's/_/-/g' "${out_dir}/data/train/train_pgen.psam"


## make bgen for gw ridge prediction

plink2 \
  --bfile "${out_dir}/data/test/test" \
  --export bgen-1.3 \
  --out "${out_dir}/data/test/test_bgen"

bgenix -g "${out_dir}/data/test/test_bgen.bgen -index

## make vcf for prediction to test the models
plink2 \
  --bfile "${out_dir}/data/test/test" \
  --recode vcf-iid \
  --out "${out_dir}/data/test/test_vcf"

echo -e "\nDONE!!"
