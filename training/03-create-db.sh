#! /usr/bin/env bash

set -euo pipefail

usage() {
  echo "Usage: $0 -l <lasso_wt_dir> -o <out_prefix>"
  exit 1
}

while getopts ":r:l:p:o" opt; do
  case "${opt}" in
    r) repo_dir=${OPTARG};;
    l) lasso_wt_dir=${OPTARG};;
    p) plink_prefix=${OPTARG};;
    o) out_prefix=${OPTARG};;
    :) echo "Option -$OPTARG requires an argument." >&2; usage;;
    \?) echo "Invalid option: -$OPTARG" >&2; usage;;
  esac
done

# Add path to top level of local repository after -r, lasso weights output
# directory after -l, training PLINK file subset to HapMap3 SNPs after -p,
# and an output file path prefix after -o.

# Gather the cv
Rscript ${repo_dir}/gather_lasso_perf.R \
  --input_dir "${lasso_wt_dir}" \
  --output_prefix "${out_prefix}_cv"

# Create a databsae
Rscript ${repo_dir}/create_lasso_db.R \
  --input_dir "${lasso_wt_dir}" \
  --bim_file "${plink_prefix}.bim" \
  --r2_file "${out_prefix}_cv_r2_performance.txt" \
  --out_db "${out_prefix}.db"
