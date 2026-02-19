#! /usr/bin/env bash

module load gcc/6.2.0
module load bcftools
module load plink/1.90
module load gcta

## get individual list
cut -f1,2 -d " " genotype/metsim_info03_n6136_rsid.fam > data/individuals.txt

## random shuffle with a seed to select individuals randomly
# get ids to split individuals into training and testing
shuf --random-source=<(yes 149) data/individuals.txt | head -n 400 > data/test_indiv.txt

grep -v -x -f data/test_indiv.txt data/individuals.txt > data/train_indiv.txt

## split the data into training and testing
# test
echo -e "\nSplitting the data\n"
plink \
  --bfile genotype/metsim_info03_n6136_rsid \
  --keep data/train_indiv.txt \
  --make-bed \
  --out data/train/metsim_train-hapmap3

# train
plink \
  --bfile genotype/metsim_info03_n6136_rsid \
  --keep data/test_indiv.txt \
  --make-bed \
  --out data/test/metsim_test-hapmap3


## generate other formats for the train set
## make a grm
echo -e "\nCalculating a grm matrix\n"

gcta \
--bfile data/train/metsim_train-hapmap3 \
--make-grm \
--thread-num 8 \
--out data/train/metsim_train-hapmap3_grm

## Make a grm gz (for ridge)
echo -e "\nCalculating a grm matrix (gz)\n"

gcta \
--bfile data/train/metsim_train-hapmap3 \
--make-grm-gz \
--thread-num 8 \
--out data/train/metsim_train-hapmap3_grm


## Calculate pca
echo -e "\nCalculating pcas\n"

gcta \
--grm data/train/metsim_train-hapmap3_grm \
--pca 20 \
--thread-num 8 \
--out data/train/metsim_train-hapmap3_pca 

echo -e "\nGenerate a pgen format\n"

## make pgen format (uses plink 2.0) for lasso
module load plink/2.0

plink2 \
--bfile data/train/metsim_train-hapmap3 \
--make-pgen vzs \
--out data/train/metsim_train-hapmap3_pgen

# formating the psam to work 
mv data/train/metsim_train-hapmap3_pgen.psam data/train/metsim_train-hapmap3_pgen.psam.bak

awk -F"\t" '{print $2"\t"$2"\t"$3}' data/train/metsim_train-hapmap3_pgen.psam.bak > data/train/metsim_train-hapmap3_pgen.psam

## Dont use underscore in names, has a problem with gw_lasso
sed -i  -e 's/^IID/#FID/' -e 's/_/-/g' data/train/metsim_train-hapmap3_pgen.psam


## make bgen for gw ridge prediction
module load bgen/1.1.3

plink2 \
--bfile data/test/metsim_test-hapmap3 \
--export bgen-1.3 \
--out data/test/metsim_test-hapmap3_bgen

bgenix -g data/test/metsim_test-hapmap3_bgen.bgen -index

## make vcf for prediction to test the models
plink2 \
--bfile data/test/metsim_test-hapmap3 \
--recode vcf-iid \
--out data/test/metsim_test-hapmap3

echo -e "\nDONE!!"
