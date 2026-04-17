#! /usr/bin/env bash

set -euo pipefail

usage() {
  echo "Usage: $0 -p <plink_prefix> -o <out_dir> -t <thread_num>"
  exit 1
}

while getopts ":p:o:t:" opt; do
  case "${opt}" in
    p) plink_prefix=${OPTARG};;
    o) out_dir=${OPTARG};;
    t) thread_num=${OPTARG};;
    :) echo "Option -$OPTARG requires an argument." >&2; usage;;
    \?) echo "Invalid option: -$OPTARG" >&2; usage;;
  esac
done

# TO NOTE: for now, added variable paths but keeping hardcoded file names to
# enable ease of questions to Festus if needed.

# Make output directories if needed
mkdir -p "${out_dir}/train"
mkdir -p "${out_dir}/test"

# Make temp version of PLINK fileset to add constant FID if needed
plink2 --pfile "${plink_prefix}" \
  --make-bed --out "${out_dir}/tmp_plink_fid"

## get individual list
awk '{print $1, $2}' "${out_dir}/tmp_plink_fid.fam" > "${out_dir}/individuals.txt"

## random shuffle with a seed to select individuals randomly
# get ids to split individuals into training and testing
shuf --random-source=<(yes 149) "${out_dir}/individuals.txt" | \
  head -n 400 > "${out_dir}/test_indiv.txt"

grep -v -x -f "${out_dir}/test_indiv.txt" \
  "${out_dir}/individuals.txt" > "${out_dir}/train_indiv.txt"

# train
echo -e "\nSplitting the data\n"
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
echo -e "\nCalculating a grm matrix\n"

gcta \
  --bfile "${out_dir}/train/train" \
  --make-grm \
  --thread-num "${thread_num}" \
  --out "${out_dir}/train/train_grm"

## Make a grm gz (for ridge)
echo -e "\nCalculating a grm matrix (gz)\n"

gcta \
  --bfile "${out_dir}/train/train" \
  --make-grm-gz \
  --thread-num "${thread_num}" \
  --out "${out_dir}/train/train_grm"

## Calculate pca
echo -e "\nCalculating pcas\n"

gcta \
  --grm "${out_dir}/train/train_grm" \
  --pca 20 \
  --thread-num "${thread_num}" \
  --out "${out_dir}/train/train_pca"

echo -e "\nGenerate a pgen format\n"

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

## make bgen for gw ridge prediction
plink2 \
  --bfile "${out_dir}/test/test" \
  --export bgen-1.3 \
  --out "${out_dir}/test/test_bgen"

bgenix -g "${out_dir}/test/test_bgen.bgen" -index

## make vcf for prediction to test the models
plink2 \
  --bfile "${out_dir}/test/test" \
  --recode vcf-iid \
  --out "${out_dir}/test/test_vcf"

echo -e "\nDONE!!"

# Cleanup
rm ${out_dir}/tmp_*
