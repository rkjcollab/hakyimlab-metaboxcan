#! /usr/bin/env bash

set -euo pipefail

usage() {
  echo "Usage: $0 -p <plink_prefix> -o <out_dir> -t <thread_num> -f <test_prop>"
  exit 1
}

while getopts ":p:o:t:f:" opt; do
  case "${opt}" in
    p) plink_prefix=${OPTARG};;
    o) out_dir=${OPTARG};;
    t) thread_num=${OPTARG};;
    f) test_prop=${OPTARG};;
    :) echo "Option -$OPTARG requires an argument." >&2; usage;;
    \?) echo "Invalid option: -$OPTARG" >&2; usage;;
  esac
done

# Add path to PLINK2 file prefix made in step 0 (will be
# <out_dir>/gt_clean/gt_clean_hapmap3) after -p, output directory path after
# -o, thread number want to use after -t, and proportion of data (0-1) to hold
# out for testing after -f.

# TO NOTE: for now, added variable paths but keeping hardcoded file names to
# enable ease of questions to Festus if needed.

# Make output dirs if needed
mkdir -p "${out_dir}/train"
mkdir -p "${out_dir}/test"

# Make temp version of PLINK fileset to add constant FID if needed
plink2 --pfile "${plink_prefix}" \
  --make-bed --out "${out_dir}/tmp_plink_fid"

## get individual list
awk '{print $1, $2}' "${out_dir}/tmp_plink_fid.fam" > "${out_dir}/individuals.txt"

## random shuffle with a seed to select individuals randomly
# get ids to split individuals into training and testing
n_indiv=$(wc -l < "${out_dir}/individuals.txt")
n_test=$(awk -v n="$n_indiv" -v p="$test_prop" 'BEGIN {printf "%.0f", n*p}')
n_train=$(( n_indiv - n_test ))
echo "\nTotal individuals: ${n_indiv}."
echo "Test individuals: ${n_test}, based on test_prop: ${test_prop}."
echo "Train individuals: ${n_train}.\n"
shuf --random-source=<(yes 149) "${out_dir}/individuals.txt" | \
  head -n "${n_test}" > "${out_dir}/test_indiv.txt"

grep -v -x -f "${out_dir}/test_indiv.txt" \
  "${out_dir}/individuals.txt" > "${out_dir}/train_indiv.txt"

# train
echo -e "\nSplitting the data.\n"
plink \
  --bfile "${out_dir}/tmp_plink_fid" \
  --keep "${out_dir}/train_indiv.txt" \
  --keep-allele-order \
  --make-bed --out "${out_dir}/train/train"

# test
plink \
  --bfile "${out_dir}/tmp_plink_fid" \
  --keep "${out_dir}/test_indiv.txt" \
  --keep-allele-order \
  --make-bed --out "${out_dir}/test/test"

## generate other formats for the train set
## make a grm
echo -e "\nCalculating a grm matrix.\n"

gcta64 \
  --bfile "${out_dir}/train/train" \
  --make-grm \
  --thread-num "${thread_num}" \
  --out "${out_dir}/train/train_grm"

## Make a grm gz (for ridge)
echo -e "\nCalculating a grm matrix (gz)\n"

gcta64 \
  --bfile "${out_dir}/train/train" \
  --make-grm-gz \
  --thread-num "${thread_num}" \
  --out "${out_dir}/train/train_grm"

## Calculate pca
echo -e "\nCalculating PCA.\n"

gcta64 \
  --grm "${out_dir}/train/train_grm" \
  --pca 20 \
  --thread-num "${thread_num}" \
  --out "${out_dir}/train/train_pca"

echo -e "\nGenerate pgen format.\n"

## make pgen format (uses plink 2.0) for lasso
plink2 \
  --bfile "${out_dir}/train/train" \
  --make-pgen vzs \
  --out "${out_dir}/train/train_pgen"

# formating the psam to work 
mv "${out_dir}/train/train_pgen.psam" "${out_dir}/train/train_pgen.psam.bak"

awk -F"\t" '{print $2"\t"$2"\t"$3}' "${out_dir}/train/train_pgen.psam.bak" > \
  "${out_dir}/train/train_pgen.psam"

## Dont use underscore in names, has a problem with gw_lasso
sed -i.bak  -e 's/^IID/#FID/' -e 's/_/-/g' "${out_dir}/train/train_pgen.psam"

echo -e "\nGenerate bgen format.\n"

## make bgen for gw ridge prediction
plink2 \
  --bfile "${out_dir}/test/test" \
  --export bgen-1.3 \
  --out "${out_dir}/test/test_bgen"

bgenix -g "${out_dir}/test/test_bgen.bgen" -index -clobber

echo -e "\nGenerate VCF.\n"

## make vcf for prediction to test the models
plink2 \
  --bfile "${out_dir}/test/test" \
  --recode vcf-iid \
  --out "${out_dir}/test/test_vcf"

# Cleanup
rm ${out_dir}/tmp_*
