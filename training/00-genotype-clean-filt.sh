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

# Based on README, this script adds basic genetic data cleaning (SNP & sample
# missingness, HWE), and then filters to HapMap3 SNPs only.
#TODO: do we want any other QC here? MAF?