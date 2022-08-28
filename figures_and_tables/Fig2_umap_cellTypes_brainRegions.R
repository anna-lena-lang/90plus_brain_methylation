## code to generate UMAP plot (Fig 1) with cell types and brain regions
## run pca on the cell type data + bulk data for visualization


library(isva)
library(wateRmelon)
library(uwot)
library(ggplot2)
library(tidyverse)


# set up save directory
setwd("~/90plus/deconvolution/scripts/")
savedir <- paste0("../output/02_pca/")
dir.create(savedir, showWarnings=FALSE)


## function --
## This function loads the bulk data for all brain regions
## region = name of the brain region to load
load_bulk <- function(region){
    print(paste("Processing", region))


    # get methylation and targets file -- update these to match file names of uploaded data **
    methylation_file <- paste0("~/90plus/data/beta_90plus_ewasfilt_bmiq_orig_850k_", region, ".rds")
    targets_file <- paste0("~/90plus/data/targets_incl_comorb.rds")


    # Load data
    # read the 90plus methylation data 
    all_meths <- readRDS(methylation_file)
    targets <- readRDS(targets_file)

    # use only the data for this brain region
    targets_br <- targets %>%
        filter(brain_region == region)
    samples <- targets_br$sample_name

    # subset methylation
    meths_matched <- all_meths[,samples]
    head(meths_matched)
    dim(meths_matched)

    meths_matched
}

## function --
## This function loads the cell type data for all brain regions
## Paramters;
## - region = name of the brain region to load
## - cell_type = name of the cell type to load
load_celltype <- function(region, cell_type){

    # do not include microglia unless they are from LOC, SN
    if (cell_type == 'microglia' & !(region %in% c('LOC', 'SN'))) return(NA)
    # do not include oligo_opc CBL
    if (cell_type == 'oligo_opc' & (region == 'CBL')) return (NA)
    
    print(paste("Loading cell type", cell_type, "in", region, Sys.time()))

    # read in the kept celltypes
    celltypes_file <- paste0("~/90plus/deconvolution/output/90plus_850k_bmiq_normalized_episcore/", region, "celltypes_present.rds")
    kept_celltypes <- readRDS(celltypes_file)
    kept_celltypes

     # some cell types are not estimated in TCA, skip these cases
    if (!cell_type %in% kept_celltypes) return(NA)


    ## load in the deconvovled cell type data from TCA output
    betas_file <- paste0("~/90plus/deconvolution/output/90plus_850k_bmiq_normalized_episcore/tca_tensor_", region, ".rds")
    betas_tensor <- readRDS(betas_file)


    # grab the correct cell type
    i <- which(kept_celltypes == cell_type)
    i

    print(paste("loading celltype", kept_celltypes[i]))

    betas <- betas_tensor[[i]] # select correct cell type betas

    # set up the data for the functions
    betas_matrix <- as.matrix(betas)
    head(betas_matrix)


    # deal with betas for cell type data
    # remove <0, >1, 0, 1 values
    print("scaling betas after TCA")
    betas_matrix[betas_matrix < 0] <- 0
    betas_matrix[betas_matrix > 1] <- 1
    betas_matrix <- (betas_matrix * (ncol(betas_matrix) - 1) + 0.5) / ncol(betas_matrix)

    sum(is.na(betas_matrix))

    celltype_meth <- data.frame(betas_matrix)
    head(celltype_meth)


    celltype_meth
}


## function --
## format the data for downstream PCA
## Parameters:
## - meth: list of methyltaion data frames
## -- each list item contains a matrix of methylation beta values for each cell type
## -- matrix rows = methylation probes
## -- matrix columns = samples
format_data <- function(meth){
    # remove NAs     
    meth <- meth[!is.na(meth)]

    # find the intersection of cpgs across all cell type matrices
    cpgs <- Reduce(intersect, (lapply(meth, function(x) rownames(x))))
    head(cpgs)
    length(cpgs) # 498,519 CpGs 

    # combine brain region data
    ordered_meth <- lapply(meth, function(x) x[cpgs,])
    formatted_meth <- do.call(cbind, ordered_meth)
    head(formatted_meth)
    dim(formatted_meth) # 495,519 probes x 321 samples

    # convert beta values to mvals
    mvals <- beta2m(formatted_meth)

    mvals
}

## function --
## run PCA on the formatted methylation data
## select the number of PCs using the Marcenko-Pastur distribution
## Arguments:
## - formatted meth: data frame containing formatted methylation data
## -- rows = methylation probes
## -- columns = samples
run_pca <- function(formatted_meth){

    # estimate dims and get sig pcs for this gene
    head(formatted_meth)
    dim(formatted_meth)

    # scale the mvals
    # transcpose to scale the rows (probes)
    scaled_meth <- t(scale(t(formatted_meth), center=TRUE, scale=TRUE))
    dim(scaled_meth)

    # estimate dimensions using Marcenko-Pastur distribution
    # catch errors in case something goes wrong
    # rows = features
    # columns = samples
    est_dim <- tryCatch(
        {
            EstDimRMT((scaled_meth), plot=FALSE)$dim
        },
        error=function(e){
            return(e)
        }
    )
    est_dim

    # check if there was an error
    if (!is.numeric(est_dim)){
        print("Estimated dim not found")
        error <- as.character(est_dim)
        est_dim <- NA
        return(NA)
    }


    # run pca
    # rows = samples
    # cols = features
    dim(scaled_meth)
    pca_res <- prcomp(t(scaled_meth))

    # get the PCs
    pcs <- pca_res$x
    head(pcs)
    
    # select the sigificant PCs based on the MP distribution
    # only keep these for further analysis
    sig_pcs <- pcs[,1:est_dim,drop=FALSE] %>%
        as.data.frame()
    head(sig_pcs)
    dim(sig_pcs)

    sig_pcs
}

## function --
## combine all of the formatted cell type data into one data frame
## Paramters:
## - all_celltypes: list of cell type methylation matrices with
## -- rows = methylation probes
## -- columns = samples
combine_celltypes <- function(all_celltypes){
    # get the intersection of cpgs
    cpgs <- Reduce(intersect, (lapply(all_celltypes, function(x) rownames(x))))
    head(cpgs)
    length(cpgs)

    # select the common CpGs and combine each of the data frames
    ordered_meth <- lapply(all_celltypes, function(x) x[cpgs,])
    formatted_meth <- do.call(cbind, ordered_meth)
    head(formatted_meth)
    dim(formatted_meth)

    # save this data
    savefile <- paste0(savedir, "combined_celltypes_mvals.rds")
    savefile
    saveRDS(formatted_meth, savefile)

    formatted_meth
}


## function --
## Compute axes for UMAP visualization 
## Parameters:
## - pcs: data frame containing the selected PCs that will be input to UMAP function
## -- rows = samples
## -- cols = PCs
## - cell_type: string indicating what cell type is being analysed, used for output filename
plot_umap <- function(pcs, cell_type){
    # plot umap on the pcs

    # get the brain region data
    targets_file <- paste0("~/90plus/data/targets_incl_comorb.rds")
    targets <- readRDS(targets_file)
    head(targets)
    dim(targets)


    # run UMAP
    # rows = observations (samples)
    # cols = features (cpgs)
    dim(pcs)
    head(pcs)
    umap_data <- umap(X = pcs, init = "pca",  a = 1.8956, b = 0.8006)

    colnames(umap_data) <- c('UMAP_1', 'UMAP_2')
    head(umap_data)
        
    
    # format the umap output
    head(umap_data)
    formatted_umap <- umap_data %>%
        as.data.frame() %>%
        rownames_to_column('sample_celltype') %>%
        separate(sample_celltype, c('sample', 'celltype'), sep='-')
    head(formatted_umap)

    # attach the phenotype data
    head(targets)
    sub_targets <- targets[targets$sample_name %in% unique(formatted_umap$sample),] %>%
        select(sample=sample_name, brain_region)
    head(sub_targets)

    formatted_umap <- formatted_umap %>%
        left_join(sub_targets)
    head(formatted_umap)

    savefile <- paste0(savedir, cell_type, "_umap_output.rds")
    savefile
    saveRDS(formatted_umap, savefile)

    # plot
    savefile <- paste0(savedir, cell_type, "_umap_plot.pdf")
    savefile

  
    # change brain region naming and remove bulk data
    formatted_dims <-  formatted_umap %>%
        mutate(brain_region = as.character(brain_region),
                brain_region = ifelse(brain_region == 'DEN',
                                     'DG', brain_region),
                brain_region = ifelse(brain_region == 'ERC',
                                     'EC', brain_region),
                brain_region = ifelse(brain_region == 'LOC',
                                     'LC', brain_region),
                brain_region = ifelse(brain_region == 'CBL',
                                     'CBM', brain_region),
        ) %>%
        mutate(celltype = ifelse(celltype == 'astro',
                         'astrocytes', celltype),
                celltype = ifelse(celltype == 'endo',
                                 'endothelial cells', celltype),
                celltype = str_to_title(celltype),
                celltype = ifelse(celltype == 'Oligo_opc',
                                 'Olig/OPC', celltype)
            ) %>%
        rename(`UMAP 1` = UMAP_1,
                `UMAP 2` = UMAP_2,
                `Cell type` = celltype) %>%
        filter(`Cell type` != 'Bulk')


    # set factor levels
    formatted_dims$brain_region = factor(formatted_dims$brain_region,
                               levels=c('MFG', 'CG', 'CA1',
                                        'DG', 'EC', 'LC', 
                                        'SN', 'CBM'))
    formatted_dims$`Cell type` = factor(formatted_dims$`Cell type`,
                              levels=c(#'Bulk', 
                                       'Neuron',
                                       'Astrocytes', 
                                       'Endothelial Cells',
                                       'Microglia', 
                                       'Olig/OPC'))
    head(formatted_dims)


    # set colors for plotting each cell type
    cols = c(#"#000000", 
    "#648FFF", "#FE6100", "#785EF0", "#DC267F", "#FFB000")
    colorlist = setNames(cols,  
               c(#"Bulk", 
                 "Neuron", 
                 "Astrocytes", 
                 "Endothelial Cells", 
                 "Microglia", 
                 "Olig/OPC"))
    colorlist

    # set shapes for brain regions
    shapelist = setNames(c(15, 16, 17, 18, 
                 8, 6, 0, 1), 
               c('MFG', 'CG', 'CA1',
                 'DG', 'EC', 'LC', 
                 'SN', 'CBM'))


    # remove bulk cell type and create ggplot
    p <- formatted_dims %>%
        filter(`Cell type` != 'Bulk') %>%
        ggplot() +
        geom_point(aes(x=`UMAP 1`, y=`UMAP 2`, color=`Cell type`, shape=brain_region),
           size=3) +
        # facet_wrap(. ~ brain_region, nrow=3) +
        scale_color_manual(values=colorlist) +
        scale_shape_manual(values=shapelist) +
        theme_bw() +
        theme(strip.background =element_rect(fill="white"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
        labs(shape='Brain region')

        savefile <- paste0(savedir, "full_umap.pdf")
        savefile
        ggsave(savefile, width=8, height=6)

    1
}

## function --
## main function handles the workflow of the script
main <- function(){

    # cell types and brain regions to process
    all_regions <- c("CA1", "CBL", "CG", "DEN", "ERC", "LOC", "MFG", "SN")
    cell_types <- c("neuron", "oligo_opc", "astro", "endo", "microglia")


    # -- BULK DATA -- #
    # start with processing the bulk dataset
    print(paste("Process bulk data"))

    # load the bulk data
    bulk_meth <- lapply(all_regions, load_bulk)
    names(bulk_meth) <- all_regions
    
    # format the data
    formatted_bulk <- format_data(bulk_meth)



    ## -- CELL TYPE DATA -- ##
    print(paste("Processing cell type data"))
    savefile <- paste0(savedir, "combined_celltypes_mvals.rds")
    if (file.exists(savefile)){
        combined_celltypes <- readRDS(savefile)
    } else{
        #region <- "CBL"
        #cell_type <- cell_types[2]
        # load all cell types 
        process_celltypes <- function(cell_type){
            print(paste("Processing", cell_type))

            # stick to one cell type for now
            # load celltype across all brain regions
            print(paste("Loading cell type data", Sys.time()))
            celltype_meth <- lapply(all_regions, load_celltype, cell_type)

            if (all(is.na(celltype_meth))) return(NA)


            # format the data
            print(paste("Formatting data", Sys.time()))
            formatted_celltype <- format_data(celltype_meth)


            # add cell type to sample name
            colnames(formatted_celltype) <- paste(colnames(formatted_celltype), cell_type, sep='-')
            head(formatted_celltype)

            # return the formatted methylation data
            formatted_celltype
        }

        # load and format the methylation data for each cell type
        all_celltypes <- lapply(cell_types, process_celltypes)


        # combine the cell type data
        combined_celltypes <- combine_celltypes(all_celltypes)
    } 


    # combine the bulk and cell type data to run PCA and UMAP on combined datasets
    head(formatted_bulk)
    head(combined_celltypes)

    # add the cell type to bulk dataset
    colnames(formatted_bulk) <- paste0(colnames(formatted_bulk), "-bulk")

    # find the intersection of CpGs across bulk and cell type data
    cpgs <- intersect(rownames(formatted_bulk), rownames(combined_celltypes))
    head(cpgs)

    # match the ordering of rows
    ordered_bulk <- formatted_bulk[cpgs,]
    ordered_celltypes <- combined_celltypes[cpgs,]

    head(ordered_bulk)
    head(ordered_celltypes)
    
    # combine the matched datasets
    bulk_celltypes <- cbind(ordered_bulk, ordered_celltypes)
    head(bulk_celltypes)


    # -- run pca
    formatted_meth <- bulk_celltypes
    pcs <- run_pca(formatted_meth)
    
    # save this output
    cell_type <- "celltypes_with_bulk"
    savefile <- paste0(savedir, cell_type, "_sig_pcs.rds")
    savefile
    saveRDS(pcs, savefile)

    # -- create umap
    plot_umap(pcs, cell_type)

    1
}

main()
