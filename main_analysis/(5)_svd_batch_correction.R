# this script uses SVA to correct for batch in the 90+ data. 
# Run this on mvals from the normalized and QC'ed beta values
# Need to pick one phenotype at a time to preserve signal

library(isva) # required for smartSVA
library(wateRmelon)
library(RnBeads)
library(RnBeads.hg19)
library(limma)
library(tidyverse)

setwd("~/90plus/batch_correction/svd/scripts/")


## parameters to control the regions, cell types, biotypes, and traits to process
brain_regions <- c("CA1", "CBL", "CG", "DEN", "ERC", "LOC", "MFG", "SN")
cell_types <- c("neuron", "oligo_opc", "astro", "endo", "microglia", "bulk")
biotypes <- c("genes", "promoters", "cpgislands", "cpgs")
protected_traits <- c('niaaaascore', 'niaaabscore', 'niaaacscore')

## parameters:
## correct_pval: TRUE to use p-value correction in correlation analysis
## reference_type: set to episcore to use episcore cell type deconvolution data
correct_pval <- TRUE
reference_type <- "episcore"




# ------------------------------------------------------------------------
## functions defined here --


# the run_all function
# runs the entire pipeline
run_all <- function(parameters){
    # check if final output exists
    # and skip run if it does
    corrected_dir <- "uncorrected_pvals/"
    if (correct_pval == TRUE){
        corrected_dir <- "corrected_pvals/"
    }

    savedir <- paste0(parameters$datadir, "cleaned_mvals/", parameters$trait, "/", parameters$cell_type, "/", corrected_dir)
    savefile <- paste(parameters$cell_type, parameters$biotype, parameters$brain_region, parameters$trait, "cleaned_mvals.rds", sep="_")
    savepath <- paste0(savedir, savefile)
    savepath

    if (file.exists(savepath)){
        print("Final output exists. Skipping region.")
        return(NA)
    }


    # save the general biotype mvals
    # only need to do this once 
    ret <- save_biotype_mvals(parameters)
    if (is.na(ret)) return(NA)


    # load the formatted and scaled mvals + the targets file
    res <- get_mvals(parameters)
    if (all(is.na(res))) return(NA)

    mvals <- res$mvals
    targets <- res$targets

    sum(is.na(mvals))

    print("Running SVD")
    run_svd(parameters, mvals)

    print("Correlate SVs with trait")
    correlate_svs(parameters)

    print("Remove SVs")
    limma_remove_svs(parameters, mvals, targets)

    1
} 



# gets the targets file and filter down to relevant samples
# for the run parameters
get_targets_file <- function(parameters){
    region <- parameters$brain_region
    cell_type <- parameters$cell_type
    biotype <- parameters$biotype
    cohort_type <- parameters$cohort_type
    analysis_type <- parameters$analysis_type

    print_update("Retrieving targets for", parameters)

    # read in the targets data
    targets_file <- paste0("~/90plus/data/targets_incl_comorb.rds") 
    targets <- readRDS(targets_file)
    dim(targets)
    head(targets)


    # add numerical values for niaa scores
    targets$fac_adseverityscore <- factor(targets$adseverityscore,
                                          levels=c('low', 'intermediate', 'high')
    )
    targets$int_adseverityscore <- as.integer(targets$fac_adseverityscore)
    table(targets$adseverityscore, targets$int_adseverityscore)

    # filter down brain region
    targets_br <- targets %>%
        filter(brain_region == region)
    dim(targets_br)
    head(targets_br)

    targets_br
}

# save the uncentered biotype mvals to save time
save_biotype_mvals <- function(parameters){
    names(parameters)
    region <- parameters$brain_region
    cell_type <- parameters$cell_type
    biotype <- parameters$biotype
    brain_region <- parameters$brain_region


    # define the output directory for m-values
    mvals_savedir <- paste0("../output/")
    dir.create(mvals_savedir, showWarnings=FALSE)
    mvals_savedir <- paste0("../output/", reference_type, "/")
    dir.create(mvals_savedir, showWarnings=FALSE)
    mvals_savedir <- paste0("../output/", reference_type, "/mvals/")
    dir.create(mvals_savedir, showWarnings=FALSE)
    mvals_savedir

    savefile <- paste(cell_type, biotype, region, "mvals.rds", sep="_")
    savepath <- paste0(mvals_savedir, savefile)
    savepath

    # load and skip if the mvalues already exist
    if (file.exists(savepath)){
        res <- readRDS(savepath)
        return(1)
    }
    

    print(paste("Getting scaled mvals for", cell_type, region, biotype))


    print("Loading data")
    # get the correct mvals depending on cell type
    if (cell_type == "bulk"){
        ## load in the data
        betas_file <- paste0("~/90plus/data/beta_90plus_ewasfilt_bmiq_orig_850k_", region, ".rds")
    
        # read in the files
        betas <- readRDS(betas_file)

    } else{
        # read in the kept celltypes
        celltypes_file <- paste0("~/90plus/deconvolution/output/90plus_850k_bmiq_normalized_episcore/", region, "celltypes_present.rds")
        kept_celltypes <- readRDS(celltypes_file)
        kept_celltypes

        # some cell types are not estimated in TCA
        if (!cell_type %in% kept_celltypes) return(NA)


        ## load in the data
        betas_file <- paste0("~/90plus/deconvolution/output/90plus_850k_bmiq_normalized_episcore/tca_tensor_", region, ".rds")
        betas_tensor <- readRDS(betas_file)

        length(betas_tensor)

        
        # grab the correct cell type
        i <- which(kept_celltypes == cell_type)
        i

        betas <- betas_tensor[[i]] # load neuron_pos data
    }


    # set up the data for the functions
    betas_matrix <- as.matrix(betas)


    # deal with betas for cell type data
    # remove <0, >1, 0, 1 values
    if (cell_type != "bulk"){
        print("scaling betas after TCA")
        betas_matrix[betas_matrix < 0] <- 0
        betas_matrix[betas_matrix > 1] <- 1
        betas_matrix <- (betas_matrix * (ncol(betas_matrix) - 1) + 0.5) / ncol(betas_matrix)

    }

    sum(is.na(betas_matrix))

    # get the mvals for the appropriate region type
    mvals <- NA
    if (biotype == 'cpgs'){
        print("convert to mvals")
        ## need to get mvals
        # use wateRmelon package to get mvals
        mvals <- Beta2M(betas_matrix)
        dim(mvals) # cpgs X samples 

    } else{
        # biotype is genes, promoters, or cpgislands
        # use RnBeads to convert to mvals within regions
        
        # sort the samples by their paper diagnosis
        pheno_temp <- data.frame(Sample = colnames(betas))
        
        # figure out which CpGs go with the biotype using RnBeads labels
        print(paste("Gathering RnBeadSet data for", region, cell_type, biotype, sep = " "))
        rnb_set <- RnBeadSet(pheno = pheno_temp, probes = rownames(betas_matrix), betas = betas_matrix, platform = "EPIC")

        # RnBeads averages methylatin within genomic regions by default
        print(paste("Computing mvals for", region, cell_type, biotype, sep = " "))
        rnb_biotype <- mval(rnb_set, type = biotype, row.names = TRUE) 
       

        if (biotype == "cpgislands"){
            ## get the data for the regions
            region_ranges <- rnb.annotation2data.frame(rnb.get.annotation(biotype))
            region_ranges['regInd'] <- 1:nrow(region_ranges)

            # keep only the ranges for cpgislands that we have in our dataset 
            keep_ranges <- region_ranges[rownames(rnb_biotype),]

            # replace rownames with the cpgisland ranges
            rownames(rnb_biotype) <- unite(keep_ranges[c('Chromosome', 'Start', 'End')], newcol, remove=FALSE)$newcol  

        }

        mvals <- rnb_biotype        
    }
    
    dim(mvals)
    sum(is.na(mvals)) # should be zero

    # save results for next time
    savepath
    saveRDS(mvals, savepath)

    1
}

# get the mvals for the region type
# format and scale them
# load the targets file as well
get_mvals <- function(parameters){
    names(parameters)
    region <- parameters$brain_region
    cell_type <- parameters$cell_type
    biotype <- parameters$biotype
    datadir <- parameters$datadir

    # read in the uncentered mvals
    mvals_savedir <- paste0("../output/", reference_type, "/mvals/")
    mvals_savedir

    mvals_file <- paste(cell_type, biotype, region, "mvals.rds", sep="_")
    mvals_path <- paste0(mvals_savedir, mvals_file)
    mvals_path


    print_update("Getting scaled mvals for", parameters)


    print("Loading data")
    # get the targets data
    targets <- get_targets_file(parameters)

    # get the correct mvals depending on cell type
    mvals <- readRDS(mvals_path)
    sum(is.na(mvals))

    dim(targets)
    dim(mvals)

    # get the targets that we care about
    head(mvals)
    head(targets)
    keep_samples <- targets$sample_name
    sub_mvals <- mvals[,keep_samples]
    length(keep_samples)
    dim(sub_mvals)


    dim(targets)
    dim(sub_mvals)
    sum(is.na(sub_mvals))

    if (nrow(targets) != ncol(sub_mvals)){
        print("ERROR: targets and mvals do not match!!")
        exit()
    }


    cbind(targets$sample_name, colnames(sub_mvals)) %>% head()


    head(sub_mvals)[1:10]

    # remove zero variance probes
    dim(sub_mvals)
    sum(apply(sub_mvals, 1, sd) == 0)

    var_mvals <- sub_mvals[apply(sub_mvals, 1, sd) != 0,]
    dim(sub_mvals)
    dim(var_mvals)
    head(var_mvals)


    # scale the data
    # scale centers and scales the columns
    # we want to scale across cpgs so we need to transpose our matrix
    mvals_centered_scaled <- t(scale(t(var_mvals), scale=TRUE))
    #mvals_centered_scaled <- apply(sub_mvals, 2, function(y) (y - mean(y) / sd(y) ^ as.logical(sd(y))))

    # save the mvals

    # number of subjects should match
    print(paste("targets:")) 
    print(dim(targets))
    print(paste("mvals:")) 
    print(dim(mvals_centered_scaled))

    sum(is.na(mvals_centered_scaled))

    res <- list(targets=targets,
                mvals=mvals_centered_scaled)

    return(res)
}

# ------------------------------------------------------------------------


# perform SVD to get SVs without protecting any traits
run_svd <- function(parameters, mvals){
    names(parameters)
    region <- parameters$brain_region
    cell_type <- parameters$cell_type
    trait <- parameters$trait
    biotype <- parameters$biotype
    datadir <- parameters$datadir


    # keep the sample names
    targets_ordered <- colnames(mvals)

    # estimate dimension of data using Random Matrix Theory
    # rows = features
    # cols = samples
    n.sv <- EstDimRMT(mvals, FALSE)$dim 
    #EstDimRMT(mvals, FALSE)

    # run SVD on the matrix
    svd.o <- svd(mvals)

    # grab the top SVs
    svs <- svd.o$v[,1:n.sv]
    
    # grab the proportion explained
    eig_sq <- svd.o$d ^ 2
    pct_var <- eig_sq / sum(eig_sq)
    subset_pct_var <- pct_var[1:n.sv]


    dim(svs)

    print("Saving files")
    sv_savefile <- paste(cell_type, biotype, region, "scaled_centered_bmiq_svd_svs.rds", sep="_")
    var_savefile <- paste(cell_type, biotype, region, "scaled_centered_bmiq_svd_pct_var.rds", sep="_")
    targets_savefile <- paste(cell_type, biotype, region, "scaled_centered_bmiq_samples_ordered.rds", sep="_")


    # store the SVs for this brain region
    sv_savedir <- paste0(datadir)
    dir.create(sv_savedir, showWarnings=FALSE)
    sv_savedir <- paste0(datadir, "bmiq_", cell_type, "/")
    dir.create(sv_savedir, showWarnings=FALSE)
    sv_savedir <- paste0(datadir, "bmiq_", cell_type, "/", biotype, "/")
    dir.create(sv_savedir, showWarnings=FALSE)
    
    sv_savepath <- paste0(sv_savedir, sv_savefile)
    var_savepath <- paste0(sv_savedir, var_savefile)
    targets_savepath <- paste0(sv_savedir, targets_savefile)

    saveRDS(svs, sv_savepath)
    saveRDS(subset_pct_var, var_savepath)
    saveRDS(targets_ordered, targets_savepath)
    1
}


# -------------------------------------------------------------------


# correlate SVs with trait
correlate_svs <- function(parameters){
    names(parameters)
    trait <- parameters$trait
    cell_type <- parameters$cell_type
    region <- parameters$brain_region
    biotype <- parameters$biotype
    datadir <- parameters$datadir

    targets <- get_targets_file(parameters)
    head(targets)
    dim(targets)


    # traits that we'll test against 
    names(targets)
    targets$sex <- targets$gender
    traits <- c('sample_name', 'niaaaascore', 'niaaabscore', 'niaaacscore', 'adseverityscore',
                'agedeath', 'int_adseverityscore', 'conf_diag', 'adss_group',
                'sex', 'batch', 'bisulfite_batch', 'paper_diag', 'path',
                'braak_for_lb_orig', 'tdp43', 'mvl_score'
                )

    # subset the targets
    target_traits <- targets[traits]
    head(target_traits)

    table(target_traits$batch, target_traits$bisulfite_batch)


    # read in the needed files
    sv_file <- paste(cell_type, biotype, region, "scaled_centered_bmiq_svd_svs.rds", sep="_")
    samples_file <- paste(cell_type, biotype, region, "scaled_centered_bmiq_samples_ordered.rds", sep="_")

    sv_dir <- paste0(datadir, "bmiq_", cell_type, "/", biotype, "/")
    
    sv_path <- paste0(sv_dir, sv_file)
    samples_path <- paste0(sv_dir, samples_file)

    svs <- readRDS(sv_path)
    ordered_samples <- readRDS(samples_path)


    # make sure samples are ordered appropriately    
    ordered_traits <- target_traits[match(ordered_samples, target_traits$sample_name),]
    head(ordered_traits)

    # set numeric and factor traits
    numerics <- c('agedeath', 'niaaaascore', 'niaaabscore', 'niaaacscore', 'int_adseverityscore')
    factors <- traits[-which(traits %in% numerics)]
    ordered_traits[numerics] <- sapply(ordered_traits[numerics], as.numeric) 
    ordered_traits[factors] <- sapply(ordered_traits[factors], as.factor) 
    


    # put svs in dataframe
    svs_df <- as.data.frame(svs)

    # get the number of svs
    n.sv <- ncol(svs)

    # in case there was only one sv
    if (is.null(n.sv)) n.sv <- 1
    n.sv

    n.traits <- length(traits) - 1

    # create matrix to hold pvalue from correlations
    pvals_mat <- matrix(nrow=n.sv, ncol=n.traits)
    cols <- rep(0,n.traits)


    # correlte the SVs with each of the traits
    sv_num <- 1
    trait_num <- 1
    for (sv_num in 1:n.sv){
        for (trait_num in 1:n.traits){
            trait <- traits[trait_num+1]

            if (class(ordered_traits[,trait]) != 'numeric'){
                pvals_mat[sv_num, trait_num] <- kruskal.test(svs_df[,sv_num] ~ ordered_traits[,trait])$p.value
                #kruskal.test(svs_df[,sv_num] ~ ordered_traits[,trait])$p.value
            } else{
                #pvals_mat[sv_num, trait_num] <- summary(lm(svs_df[,sv_num] ~ ordered_traits[,trait]))$coeff[2,4]
                #summary(lm(svs_df[,sv_num] ~ ordered_traits[,trait]))$coeff
                pvals_mat[sv_num, trait_num] <- cor.test(svs_df[,sv_num], ordered_traits[,trait], method='spearman',
                                                         exact=FALSE)$p.value
            }

            cols[trait_num] <- trait
        }
    }
    colnames(pvals_mat) <- cols
    pvals_mat

    # store a corrected and uncorrected version of pvals
    pvals_long <- pvals_mat %>%
        as.data.frame() %>%
        mutate(sv=paste("sv", row_number(), sep='_')) %>%
        gather(key='trait', value='pval', -sv) %>%
        group_by(trait) %>%
        mutate(bf_pval = p.adjust(pval, method='bonferroni'))
    head(pvals_long)


    # save the output
    sv_savefile <- paste(cell_type, biotype, region, "scaled_centered_bmiq_lm_pvals.rds", sep="_")

    # store the SVs for this brain region
    sv_savedir <- paste0(datadir)
    dir.create(sv_savedir, showWarnings=FALSE)
    sv_savedir <- paste0(datadir, "bmiq_", cell_type, "/")
    dir.create(sv_savedir, showWarnings=FALSE)
    sv_savedir <- paste0(datadir, "bmiq_", cell_type, "/", biotype, "/")
    dir.create(sv_savedir, showWarnings=FALSE)
    
    sv_savepath <- paste0(sv_savedir, sv_savefile)
    sv_savepath

    saveRDS(pvals_long, sv_savepath)
    1
}



# -----------------------------------------------------------------


# remove the SVs from the original data
# if they are correlated with trait of interest
# using limma's removeBatchEffect
limma_remove_svs <- function(parameters, mvals, targets){
    names(parameters)
    trait <- parameters$trait
    region <- parameters$brain_region
    cell_type <- parameters$cell_type
    biotype <- parameters$biotype
    datadir <- parameters$datadir

    # read in the necessary data
    # mvals, SVs, lm results

    # SVs
    sv_file <- paste(cell_type, biotype, region, "scaled_centered_bmiq_svd_svs.rds", sep="_")
    samples_file <- paste(cell_type, biotype, region, "scaled_centered_bmiq_samples_ordered.rds", sep="_")
    sv_dir <- paste0(datadir, "bmiq_", cell_type, "/", biotype, "/")
    sv_path <- paste0(sv_dir, sv_file)
    samples_path <- paste0(sv_dir, samples_file)

    svs <- readRDS(sv_path)
    ordered_samples <- readRDS(samples_path)

    # lm results
    lm_file <- paste(cell_type, biotype, region, "scaled_centered_bmiq_lm_pvals.rds", sep="_")
    lm_dir <- paste0(datadir, "bmiq_", cell_type, "/", biotype, "/")
    lm_path <- paste0(lm_dir, lm_file)

    lm_pvals <- readRDS(lm_path)


    # choose svs to remove
    head(svs)
    head(lm_pvals)



    trait
    to_keep <- c('bisulfite_batch', 'batch', as.character(trait))

    # select the correct pval column
    pval_col <- 'pval'
    if (correct_pval == TRUE){
        pval_col <- 'bf_pval'
    }

    # select the correct pval column
    # subset the lm pvals
    sub_lm_pvals <- lm_pvals[c('sv', 'trait', pval_col)] %>%
        filter(trait %in% to_keep) %>%
        spread(trait, value='bf_pval')
    head(sub_lm_pvals)


    # figure out which SVs to protect
    # protect SVs that are correlated with trait
    # and not correlated with either batch variable
    threshold <- 0.1
    protected_svs <- sub_lm_pvals %>%
        filter((batch >= threshold) & 
               (bisulfite_batch >= threshold) &
                (trait < threshold))
    protected_svs



    # we only have one sv
    if (is.null(dim(svs))){
        svs <- data.frame(sv_1 = svs)
    }

    # add column names to match the rm_svs labels
    colnames(svs) <- paste0("sv_", 1:ncol(svs))
    head(svs)

    svs_df <- as.data.frame(svs)
    head(svs_df)

    unprotected_svs <- NA
    if (nrow(protected_svs) == 0){
        unprotected_svs <- svs_df
    } else{
        unprotected_svs <- svs_df[, !(colnames(svs_df) %in% protected_svs$sv), drop=FALSE]
    }
    head(unprotected_svs)

    unprotcol = ncol(unprotected_svs)
    if (nrow(unprotected_svs) == 0) unprotcol = 0


    if ((nrow(protected_svs)+ unprotcol) != ncol(svs_df)) print("ERROR: PROTECTED SVS + UNPROTECTED SVS != SVS")
    unprotcol
    ncol(svs_df)
    
    protected_svs
    unprotected_svs


    
    # remove SVs from mvals
    
    # should match 
    dim(unprotected_svs)
    length(ordered_samples)

    # subset targets
    targets_trait <- targets[c('sample_name', trait)]
    targets_ordered <- targets_trait[match(ordered_samples, targets_trait$sample_name),]

    # order mvals according to targets
    mvals_ordered <- mvals[,targets_ordered$sample_name]
    head(mvals_ordered)

    dim(mvals_ordered)
    dim(targets_ordered)

    head(targets_ordered)
    head(unprotected_svs)

    # remove SVs
    cleaned_mvals <- removeBatchEffect(x=mvals_ordered, covariates=unprotected_svs) 
    head(cleaned_mvals)
    head(mvals_ordered)


    savedir <- paste0(datadir, "cleaned_mvals/")
    dir.create(savedir, showWarnings=FALSE)
    savedir <- paste0(datadir, "cleaned_mvals/", trait, "/")
    dir.create(savedir, showWarnings=FALSE)
    savedir <- paste0(datadir, "cleaned_mvals/", trait, "/", cell_type, "/")
    dir.create(savedir, showWarnings=FALSE)

    corrected_dir <- "uncorrected_pvals/"
    if (correct_pval == TRUE){
        corrected_dir <- "corrected_pvals/"
    }

    savedir <- paste0(datadir, "cleaned_mvals/", trait, "/", cell_type, "/", corrected_dir)
    dir.create(savedir, showWarnings=FALSE)

    savefile <- paste(cell_type, biotype, region, trait, "cleaned_mvals.rds", sep="_")
    savepath <- paste0(savedir, savefile)
    savepath

    saveRDS(cleaned_mvals, savepath)


    # store this for checking later
    savedir <- paste0(savedir, "protected_svs/")
    dir.create(savedir, showWarnings=FALSE)

    savefile <- paste(cell_type, biotype, region, trait, "protected_svs.txt", sep="_")
    savepath <- paste0(savedir, savefile)
    write.table(protected_svs, savepath, sep='\t', quote=FALSE, row.names=FALSE)

    1
}

# print update while running
print_update <- function(message, parameters){
    names(parameters)
    print(paste(message, 
                "analysis_type:", parameters$analysis_type,
                parameters$brain_region, parameters$biotype, parameters$cell_type, parameters$trait))
    1
}

# start_run




# ------------------------------------------------------------------------

# run the script from here -- 


## -- create all combinations of input parameters to run
runs_df <- expand.grid(
                       trait = protected_traits,
                       cell_type = cell_types,
                       biotype = biotypes,
                       brain_region = brain_regions
                    )
head(runs_df)
dim(runs_df)


# process each of the runs
res <- apply(runs_df, 1, start_run)


# start running the program 
row <- runs_df[1,] # for testing
start_run <- function(row){
    # get the parameters for this run
    trait = row[['trait']] %>% as.character()
    cell_type = row[['cell_type']] %>% as.character()
    biotype = row[['biotype']] %>% as.character()
    brain_region = row[['brain_region']] %>% as.character()
    analysis_type <- trait


    # define the data directory
    datadir <- paste0("../output/", reference_type, "/")
    dir.create(datadir, showWarnings=FALSE)
    datadir <- paste0("../output/", reference_type, "/", analysis_type, "/")
    dir.create(datadir, showWarnings=FALSE)
    datadir

    # create a parameters list to pass to functions
    parameters <- list(analysis_type = analysis_type,
                        brain_region = brain_region,
                        biotype = biotype,
                        cell_type = cell_type,
                        trait = trait,
                        datadir = datadir
                        )

    print(paste("*** Processing", analysis_type, cell_type, biotype, brain_region, trait))

    # run through the functions
    run_all(parameters)

    1
}


