## script for cel type deconvolution using Tensor Composition Analysis


library(TCA)
library(magrittr)
library(readxl)
library(stringr)
library(tidyverse)
library(RhpcBLASctl)
blas_set_num_threads(16) #limit number of threads for openblas to 16


setwd("~/ad_methylome/deconvolution/scripts/")

### things to input at runtime
### output_folder = where to output the tca results
### decomposition file = location of file containing the cellcount2 results
output_folder <- "../output/90plus_850k_bmiq_normalized_episcore/"
dir.create(output_folder, showWarnings=FALSE)

decomposition_file <- "../data/episcore/targets_incl_episcore.rds"



## function: load_data
## loads the methylation data and the phenotype data
load_data <- function(region){
    # get methylation and targets file
    methylation_file <- paste0("~/ad_methylome/methylome_project/Methylome_Cohorts_Data/DATASETS/90PLUS/90PLUS_850k/90PLUS_850k_normalized/90PLUS_ewastools_filt_normalized_850k/90PLUS_ewastools_filt_normalized_850k_separated_by_region/beta_90plus_ewasfilt_bmiq_orig_850k_", region, ".rds")
    targets_file <- paste0("~/ad_methylome/methylome_project/Methylome_Cohorts_Data/DATASETS/90PLUS/90PLUS_850k/90PLUS_850k_normalized/90PLUS_ewastools_filt_normalized_850k/90PLUS_ewastools_filt_normalized_850k_combined_all_regions/targets_incl_comorb.rds")


    # Load data
    # read the 90plus methylation data 
    all_meths <- readRDS(methylation_file)
    targets <- readRDS(targets_file)


    # grab only the phenotype data for this brain region
    targets_br <- targets %>%
        filter(brain_region == region)
    samples <- targets_br$sample_name

    # subset methylation
    meths_matched <- all_meths[,samples]


    # Prepare cell type composition estimates
    # read in episcore cell type estimates
    meth_age <- readRDS(decomposition_file) %>% 
        filter(brain_region == region)

    # filter meth age
    meth_age_matched <- meth_age %>%
        filter(sample_name %in% samples)
    head(meth_age_matched)


    dim(meth_age_matched)
    dim(targets_br)
    dim(meths_matched)
    
    # make sure dimensions match
    if (length(unique(c(nrow(meth_age_matched), nrow(targets_br), ncol(meths_matched)))) > 1){
        print(paste("ERROR: dimensions of data do not match!"))
        exit()
    }


    # combine the oligo and opc proportions as sum since they have low proportions
    head(meth_age_matched) 
    cell_props <- meth_age_matched %>%
        select(sample_name, neuron=Neuron, oligo=Oligo, astro=Astro, opc=OPC, endo=Endo, microglia=Microglia) %>%
        mutate(oligo_opc = oligo + opc) %>% 
        select(-oligo, -opc) %>% 
        mutate(total = neuron + oligo_opc + astro + endo + microglia) %>%
        rownames_to_column('rn') %>%
        select(-rn) %>%
        column_to_rownames('sample_name') 
    head(cell_props)

    # adjust proportions to sum to 1
    adjusted_props <- cell_props[1:5] / cell_props$total
    head(adjusted_props)

    # all sums should now be 1
    sums <- apply(adjusted_props, 1, sum)
    table(sums)

    # remove any cell types that are zero across all samples
    nonzero_celltypes <- adjusted_props[,apply(adjusted_props, 2, sum) != 0]
    head(nonzero_celltypes)

    # save the cell types that have at least one non-zero proportion
    kept_celltypes <- colnames(nonzero_celltypes)
    savefile <- paste0(output_folder, region, "celltypes_present.rds")
    saveRDS(kept_celltypes, savefile)

    # convert to matrix
    celltypes_matrix <- as.matrix(nonzero_celltypes)
    head(celltypes_matrix)
    
    betas_ordered <- meths_matched[,rownames(celltypes_matrix)] #reorder samples to match annotation data
    head(betas_ordered)

    dim(celltypes_matrix)
    dim(meths_matched)

    # make sure samples match in both matrices
    if (!identical(rownames(celltypes_matrix), colnames(betas_ordered))) stop("Samples do not match")

    data <- list(celltypes_matrix=celltypes_matrix,
                betas_ordered=betas_ordered)
}


## function: runTca
## run TCA to deconvolve cell types
runTca <- function(betas_ordered, celltypes_matrix, region){

    # filename for output, skip if output already exists
    savefile <- paste(output_folder, "tca_tensor_", region, ".rds", sep = "")
    if (file.exists(savefile)) return()

    # Figure out weights for making tensor
    print(paste("Starting tca for", region, Sys.time()))
    tca_out <- tca(betas_ordered, celltypes_matrix, parallel = TRUE, num_cores = 16)
    saveRDS(tca_out, paste(output_folder, "tca_", region, ".rds", sep = ""))

    # project from original mixed beta matrix to celltype-specific tensor
    print(paste("Building tensor for", region, Sys.time()))
    blas_set_num_threads(16) #limit number of threads for openblas to 16
    tca_tensor <- tensor(betas_ordered, tca_out, parallel = FALSE)
    saveRDS(tca_tensor, paste(output_folder, "tca_tensor_", region, ".rds", sep = ""))

    1
}


## function --
## main function to run through entire workflow
main <- function(){
    
    # brain regions
    all_regions <- c("CA1", "CBL", "CG", "DEN", "ERC", "LOC", "MFG", "SN")
    #region <- all_regions[3]


    # process each brain region
    process_region <- function(region){
        print(paste("Running region", region, Sys.time()))

        # load the methylation data
        data <- load_data(region)

        betas_ordered <- data$betas_ordered
        celltypes_matrix <- data$celltypes_matrix

        head(betas_ordered)
        head(celltypes_matrix)

        try(
            runTca(betas_ordered, celltypes_matrix, region)
        )
        gc()

        1
    }

    # process each of the brain regions
    lapply(all_regions, process_region)
}

