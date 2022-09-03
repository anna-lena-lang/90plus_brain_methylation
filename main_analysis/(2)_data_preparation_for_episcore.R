## code to merge methylation beta-values from separate brain regions 
## uses the results from script (1)_data_preprocessing, 
## ..merges the betas from different brain regions and
## ..saves the combined dataset for usage in cell type estimation: 
## ..(3)_cell_type_proportion_estimation_episcore

library(dplyr)
setwd("~/90plus/")
baseDir <- getwd()

## function to merge multiple matrices like left joint
multimerge <- function (mylist) {
  unames <- unique(unlist(lapply(mylist, rownames)))
  n <- length(unames)
  out <- lapply(mylist, function(df) {
    tmp <-
      matrix(
        nr  =  n,
        nc  =  ncol(df),
        dimnames  =  list(unames, colnames(df))
      )
    tmp[rownames(df),] <- as.matrix(df)
    rm(df)
    gc()
    
    return(tmp)
  })
  
  stopifnot(all(sapply(out, function(x)
    identical(rownames(x), unames))))
  
  bigout <- do.call(cbind, out)
  # colnames(bigout) <-  paste(rep(names(mylist), sapply(mylist, ncol)), unlist(sapply(mylist, colnames)), sep  =  "")
  return(bigout)
}

## Load dfs with beta values for each brain region. These are the results from script (1)_data_preprocessing
DG <- readRDS(paste0(baseDir, "/data/beta_90plus_ewasfilt_bmiq_orig_850k_DG.rds"))
EC <- readRDS(paste0(baseDir,"/data/beta_90plus_ewasfilt_bmiq_orig_850k_EC.rds"))
CA1 <- readRDS(paste0(baseDir,"/data/beta_90plus_ewasfilt_bmiq_orig_850k_CA1.rds"))
MFG <- readRDS(paste0(baseDir,"/data/beta_90plus_ewasfilt_bmiq_orig_850k_MFG.rds"))
CG <- readRDS(paste0(baseDir,"/data/beta_90plus_ewasfilt_bmiq_orig_850k_CG.rds"))
CBM <- readRDS(paste0(baseDir,"/data/beta_90plus_ewasfilt_bmiq_orig_850k_CBM.rds"))
LC <- readRDS(paste0(baseDir,"/data/beta_90plus_ewasfilt_bmiq_orig_850k_LC.rds"))
SN <- readRDS(paste0(baseDir,"/data/beta_90plus_ewasfilt_bmiq_orig_850k_SN.rds"))

## merge into one df
beta <- multimerge(list(DG, LC, CBM, CG, MFG, SN, CA1, EC))

## removing all rows with any NA, keeping only probes that all brain region datasets share
beta <- na.omit(beta)
saveRDS(beta, file = paste0(getwd(), "/data/merged_beta/merged_scaled_betas_separate_brainregions_cpgs_bulk.rds"))


