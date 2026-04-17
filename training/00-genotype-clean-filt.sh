#! /usr/bin/env bash

set -euo pipefail

usage() {
  echo "Usage: $0 -r <repo_dir> -p <plink_prefix> -o <out_dir>"
  exit 1
}

while getopts ":r:p:o:" opt; do
  case "${opt}" in
    r) repo_dir=${OPTARG};;
    p) plink_prefix=${OPTARG};;
    o) out_dir=${OPTARG};;
    :) echo "Option -$OPTARG requires an argument." >&2; usage;;
    \?) echo "Invalid option: -$OPTARG" >&2; usage;;
  esac
done

#TODO: make list more generic name instead of using hapmap3?

# Based on README, this script adds basic genetic data cleaning (SNP & sample
# missingness, HWE), and then filters to HapMap3 SNPs only using the list of
# SNPs at repo https://github.com/hakyimlab/Yanyus-misc-tools/tree/master/hapmap3_snps
# with MAF = 0.01 and build b37 (hg19).

# TO NOTE: this script requires PLINK2 formatted file inputs, assumes that
# input build is hg38, and assumes that given HapMap3 SNP list is hg19.
input_build_num="38"
list_build_num="19"

# Get HapMap SNP list if don't have already
if [ ! -f "${out_dir}/hapmap3_snps_maf0.01_hg19.tsv" ]; then
  wget -O "${out_dir}/hapmap3_snps_maf0.01_hg19.tsv.gz" \
    https://uchicago.box.com/shared/static/junrcgxwpuyck03r6gq88j9b18g18vf5
  gunzip "${out_dir}/hapmap3_snps_maf0.01_hg19.tsv.gz"
fi

# Remove SNPs with missingness > 5% & samples with missingness > 5%
plink2 --pfile "${plink_prefix}" \
  --geno 0.05 \
  --mind 0.05 \
  --make-pgen \
  --out "${out_dir}/gt_filt_miss"

# Apply HWE filter 1e-6 & 1e-20 in MHC region
# Get chr6 MHC region from Paul Norman's coordinates & apply HWE
start=$(head -n 1 "${repo_dir}/refs/mhc_extended_hg${input_build_num}.bed" | \
    awk -F':' '{gsub(/-.*/, "", $2); print $2}')
stop=$(tail -n 1 "${repo_dir}/refs/mhc_extended_hg${input_build_num}.bed" | \
    awk -F':' '{gsub(/-.*/, "", $2); print $2}')

plink2 --pfile "${out_dir}/gt_filt_miss" \
    --chr 6 --from-bp $start --to-bp $stop \
    --make-pgen --out "${out_dir}/tmp_mhc"

plink2 --pfile "${out_dir}/tmp_mhc" \
  --hwe 1e-20 --nonfounders \
  --write-snplist --out "${out_dir}/tmp_mhc_hwe_list"

plink2 --pfile "${out_dir}/tmp_mhc" \
  --extract "${out_dir}/tmp_mhc_hwe_list.snplist" \
  --make-bed --out "${out_dir}/tmp_mhc_hwe"

# Get all SNPs not in MHC region & apply HWE
awk '/^#/ {next} {print $3}' "${out_dir}/tmp_mhc.pvar" > \
    "${out_dir}/tmp_chr6_mhc_var_id_list.txt"

plink2 --pfile "${out_dir}/gt_filt_miss" \
    --exclude "${out_dir}/tmp_chr6_mhc_var_id_list.txt" \
    --make-pgen --out "${out_dir}/tmp_non_mhc"

plink2 --pfile "${out_dir}/tmp_non_mhc" \
  --hwe 1e-6 --nonfounders \
  --write-snplist --out "${out_dir}/tmp_non_mhc_hwe_list"

plink2 --pfile "${out_dir}/tmp_non_mhc" \
  --extract "${out_dir}/tmp_non_mhc_hwe_list.snplist" \
  --make-bed --out "${out_dir}/tmp_non_mhc_hwe"

# Merge back together
plink --bfile "${out_dir}/tmp_mhc_hwe" \
    --bmerge "${out_dir}/tmp_non_mhc_hwe" \
    --keep-allele-order --allow-no-sex \
    --make-bed --out "${out_dir}/gt_filt_miss_hwe"

# Make bed file from HapMap3 SNP list for liftover
awk -F'\t' '
BEGIN { OFS="\t" }

NR==1 {
  for (i=1; i<=NF; i++) {
    col[$i] = i
  }
  next
}

{
  print $col["chromosome"], \
        $col["position"], \
        $col["position"], \
        $col["id"]
}
' "${out_dir}/hapmap3_snps_maf0.01_hg19.tsv" > "${out_dir}/tmp_hapmap3_${list_build_num}.bed"

# Liftover HapMap3 SNPs from hg19 to hg38
CrossMap bed "${repo_dir}/refs/hg${list_build_num}ToHg${input_build_num}.over.chain" \
  "${out_dir}/tmp_hapmap3_${list_build_num}.bed"  \
  "${out_dir}/hapmap3_${input_build_num}.bed"

# Filter data to HapMap3 SNPs only, force write out of FID (constant 0 if not present)
plink2 --bfile "${out_dir}/gt_filt_miss_hwe" \
  --extract bed1 "${out_dir}/hapmap3_${input_build_num}.bed" \
  --make-pgen --out "${out_dir}/gt_clean_hapmap3"

# Cleanup
rm ${out_dir}/tmp_*
