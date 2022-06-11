## code to get estimation of cell type proportions using EPISCORE package.

## downloaded package from figshare https://figshare.com/articles/software/EpiSCORE_R_package/14401340 and saved in working directory
# install.packages("devtools")
library(devtools)
# install.packages("cli")
# install.packages("processx")
# install.packages("ellipsis")
# devtools::install_github('immunogenomics/presto')
# devtools::install_local('EpiSCORE_1.0.0.tar.gz')
# install.packages("AnnotationDbi")
# install.packages("org.Hs.eg.db")

library(cli)
library(usethis)
library(processx)
library(ellipsis)
library(presto)
library(EpiSCORE)
library(missMethyl)
library(ggpubr)
library(ggplot2)

setwd("O:/02182022_backup_Lena/90plus/")
baseDir <- paste0("O:/02182022_backup_Lena/90plus/")

## run EpiSCORE cell type proportion prediction

## load methylation data as previously merged in script (2)_data_preparation_for_episcore
beta_file <- paste0("./preprocessed_date/merged_beta/merged_scaled_betas_separate_brainregions_cpgs_bulk.rds")
beta <- readRDS(beta_file)

## load reference methylation matrix for all brain cell types
ref_file <- paste0("./episcore_reference/reference_dfs/Brain_Mref.csv")
ref_df <- read.csv(ref_file)
rownames(ref_df) <- ref_df$ï..Genes
ref_df <- as.matrix(ref_df)

## map cpgs to entrez ids
entrez_ids <- constAvBetaTSS(beta, type = "850k")

## estimate proportion of celltypes in our methylation dataset
est_proportions <- wRPC(entrez_ids,ref = ref_df, useW = TRUE, wth = 0.4, maxit = 300)

## read in phenotype data
targets_file <- "./targets_incl_comorb.rds"
targets <- readRDS(targets_file)

proportions <- as.data.frame(est_proportions$estF)
rownames(targets) <- targets$sample_name

targets <- targets[order(rownames(targets)),]
proportions <- proportions[order(rownames(proportions)),]

## add proportions to phenotype data and save to file
targets <- cbind(targets, proportions)
saveRDS(targets, "./targets_incl_episcore.rds")



