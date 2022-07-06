## script for differential methylation analysis
## m-values used here are after TCA-deconvolution and batch correction
## linear models include the following with covariates: age at death and gender

## load libraries
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(minfi)
library(RColorBrewer)
library(qvalue)
library(vctrs)
library(dplyr)
library(DMRcate)
library(limma)
library(wateRmelon)
library(stringr)
library(RnBeads)
library(tidyr)
library(viridis)
library(ggpubr)
library(biomaRt)
library(ggplot2)

## set working directory
setwd("O:/02182022_backup_Lena/90plus")
baseDir <- getwd()

## get annotation data for EPIC
annEpic <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

## set up ensembl mart
ensembl <- useMart("ensembl")
ensemblHuman <- useDataset("hsapiens_gene_ensembl",mart = ensembl)

## specify biotypes, regions
biotypes <- c("promoters", "genes", "cpgislands", "cpgs")
all_regions <- c("MFG", "CG", "DEN", "CA1", "ERC", "LOC", "SN", "CBL")

## load functions for dm analysis
source("O:/02182022_backup_Lena/EPIC_R_Scripts/2022_final_scripts/functions/dm_continuous.R")
source("O:/02182022_backup_Lena/EPIC_R_Scripts/2022_final_scripts/functions/dm_categorical.R")

## load function to add cpgs (illumina cpg-ids) to the corresponding ENSEMBL-IDs, applied to protein coding results only
source("O:/02182022_backup_Lena/EPIC_R_Scripts/2022_final_scripts/functions/add_cpgs_to_protein_coding_results.R")

## run differential methylation analysis
dm_categorical(phenotype = "conf_diag", header = "episcore_categorical_clinical")
dm_continuous(phenotype = "niaaaascore", header = "episcore_continuous_niaaaascore")
dm_continuous(phenotype = "niaaabscore", header = "episcore_continuous_niaaabscore")
dm_continuous(phenotype = "niaaacscore", header = "episcore_continuous_niaaacscore")
dm_continuous(phenotype = "int_adseverityscore", header = "episcore_continuous_int_adseverityscore")
## analyses that were not included in the publication:
# dm_continuous(phenotype = "tdp43", header = "episcore_continuous_tdp43")
# dm_continuous(phenotype = "mvl_score", header = "episcore_continuous_mvl_score")
# dm_continuous(phenotype = "braak_for_lb_orig", header = "episcore_continuous_braak_for_lb_orig")

## add cpgs to protein coding results
add_cpgs_to_protein_coding_results(category = "continuous")
add_cpgs_to_protein_coding_results(category = "categorical")
