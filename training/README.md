## Model training
Here we describe our we train our metabolite prediction model

## Genotype Quality Control
To train the metabolite models we need a high quality data to start with, first we filter our genotype to remove SNPs and individuals that have a high missingness, we also filter out SNPs that depart from HWE. Finally we filter to retain [HapMap3 SNPs only](https://github.com/hakyimlab/Yanyus-misc-tools/tree/master/hapmap3_snps), this is because metabolites don't have cis regions and the models are trained genomewide. Selecting the HapMap3 set of SNPs ensures we mantain a manageable computation cost.

Save the pre-processed genotype in plink format `.bed` and split it into training and held out testing if you need to. Using the genotype calculate thr GRM matrix and PCAs using `gcta` Also generate a `.pgen` format to use for training the model. A sample workflow is in [`01-genotype-prep.sh`](./01-genotype-prep.sh)

## Metabolite Quality Control
First it is recommended to perform an exploratory data analysis to understand the distribution, rate of missingness across indivdiduals and metabolites. You can drop or impute the missing metabolites,for imputation we use the SoftImpute package in R. 

Inverse normalize the metabolite levels if they are not inverse normalized and save the file for the next step of the model training. Your final matrix show have individuals in the rows and metabolites in the column, the first column is the individual ids (IID) then you can save your matrix like this

```r
fwrite(train, file = "data/metabolites/metsim-invnorm-lasso.txt", sep = "\t")

# extract the phenolist to use in the downstram process
pheno_list = data.frame(pheno=names(train)[-c(1)])
fwrite(pheno_list, file = "data/metabolites/all_phenos.txt", 
        sep = "\t", col.names = F)
# also a map to use if estimating h2
metab_map = data.frame(Index=c(1:(ncol(train)-1)),
                       Metabolite=names(train)[-c(1)])
fwrite(metab_map, file= "data/metabolites/metsim-invnorm_map.txt", sep ="\t")
```
The phenotype list is use in the gw_lasso model training step.

## Training the model
We use the genomewide lasso to train the metabolite prediction model. More details about  how we implement the method are available [here](https://github.com/liangyy/ukb_idp_genetic_arch/tree/master/methods/gw_lasso). 

The genomewide lasso has issues with `_` in names make sure to replace it with `-`, once all the data is prepared you can run your gw_lasso as shown in [02-train-lasso.pbs](./02-train-lasso.pbs) job script. Before running gw_lasso ensure you have the snpney_config configured with the correct paths to `plink` and `zstdcat` tools.

If you have many metabolites (phenotypes) you can split them into smaller chunks to train the model faster and then gather the weights from all the chunks. In my analysis I split the phenotype into chunks of 100 and train the models then gather the results using [03-create-db.pbs](./03-create-db.pbs)

