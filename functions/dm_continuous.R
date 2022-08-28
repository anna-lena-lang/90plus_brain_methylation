## function to run differential methylation analysis for continuous variables

dm_continuous <- function( header, phenotype) {
  
  ## targets file location
  targets_file <- "targets_incl_comorb.rds"
  
  ## load targets file
  all_targets <- readRDS(targets_file)
  
  ## dm analysis for each cell type, brain region and biotype combination:
  for (region in all_regions) {
    ## load file containing which cell types are present for each brain region
    celltypespresent <- readRDS(paste0(baseDir,"/cleaned_mvals/cleaned_mvals_episcore/episcore_celltypes_present/",region, "celltypes_present.rds"))
    ## add "bulk" to the list of cell types present
    celltypespresent <- c("bulk", celltypespresent)
    
    ## filter targets for this brain region
    targets <- all_targets %>% dplyr::filter(brain_region == region)
    
    for (celltype in celltypespresent) {

      ## change adseverityscore from categorical to integer
      if(header == "episcore_continuous_int_adseverityscore"){
      targets$int_adseverityscore <- targets$adseverityscore
      targets$int_adseverityscore <- sub("low", 1, targets$int_adseverityscore)
      targets$int_adseverityscore <- sub("intermediate", 2, targets$int_adseverityscore)
      targets$int_adseverityscore <- sub("high", 3, targets$int_adseverityscore)
      targets$int_adseverityscore <- factor(targets$int_adseverityscore, levels = c("1", "2", "3")) 
      }
      
      for (biotype in biotypes) {
        
        ## build empty dataframe to save dm results in
        sig_dms_per_brainregion <- data.frame(matrix(nrow = 1))
        colnames(sig_dms_per_brainregion)[1] <- "brain_region"
        sig_dms_per_brainregion["brain_region"] <- region
        
        ## get mvals
        mvals_file <- paste0(getwd(),"/cleaned_mvals/cleaned_mvals_episcore/", phenotype, "/",celltype,"/",celltype,"/corrected_pvals/",   celltype,"_" ,biotype,"_", region, "_", phenotype,  "_cleaned_mvals.rds")
        mvals <- readRDS(mvals_file)
        
        ## subset targets to match mvals df
        targets <- targets[targets$sample_name %in% colnames(mvals),]
        
        ## bring in right format
        targets[sapply(targets, is.character)] <- lapply(targets[sapply(targets, is.character)], as.factor)
        
        ## add new column for easier handling
        targets$predictor <- targets[[phenotype]]
        
        ## remove NAs and re-level factors and numerics
        targets_sub <- targets[!(is.na(targets$predictor)), ]
        targets_sub$agedeath <- as.character(targets_sub$agedeath) %>% as.numeric()
        targets_sub$predictor <- as.character(targets_sub$predictor)  %>% as.numeric()
        targets_sub$gender <- as.character(targets_sub$gender) %>% as.factor()
        targets_sub$sample_name <- as.character(targets_sub$sample_name)
        
        ## sort data
        mvals_sub <- mvals[, colnames(mvals) %in% targets_sub$sample_name]
        targets_sub <- targets_sub[order(targets_sub$sample_name),]
        mvals_sub <- mvals_sub[,order(colnames(mvals_sub))]
        
        if(all(colnames(mvals_sub) == targets_sub$sample_name) == FALSE){
          stop("ERROR: colnames mvals don't match targets")
        }
        
        ## ---- dm Analysis with limma ----
        ## create model
        model_content <-  paste("~ predictor + agedeath + gender", collapse = " + ")
        
        ## design model 
        design <- model.matrix(as.formula(model_content), data = targets_sub)
        
        ## run linear model
        lmfit <- lmFit(mvals_sub, design)

        ## load number of proteced svs
        number_svs_protected_file <- paste0(getwd(),"/cleaned_mvals/cleaned_mvals_episcore/",phenotype, "/",celltype,"/",celltype,"/corrected_pvals/protected_svs/",   celltype,"_" ,biotype,"_", region,"_", phenotype, "_protected_svs.txt")
        number_svs_protected <- read.table(number_svs_protected_file)
        number_svs_protected <- nrow(number_svs_protected) -1
        
        ## load number of total svs
        number_svs_total_file <- paste0(getwd(),"/svds/episcore/",phenotype, "/bmiq_",celltype,"/",biotype,"/", celltype,"_" ,biotype,"_", region, "_scaled_centered_bmiq_svd_svs.rds")
        number_svs_total <- readRDS(number_svs_total_file)
        number_svs_total <- as.data.frame(number_svs_total)
        number_svs_total <- ncol(number_svs_total)
        
        ## get number of svs that were removed by substracting protected svs from total svs
        remove <- number_svs_total - number_svs_protected
        
        print(paste0("number svs total: ", number_svs_total, ", protected: ", as.numeric(number_svs_protected), ", removed: ", remove))
        
        ## remove number of unprotected svs from degrees of freedom in lmfit
        lmfit$df.residual <- lmfit$df.residual - remove
        lm <- eBayes(lmfit)

        saveRDS(lm, file = paste(baseDir, "/dm_results/lm/",
          paste(
            "lm",
            header,
            celltype,
            biotype,
            phenotype,
            region,
            sep = "_"
          ),
          ".rds",
          sep = ""
        ))
        
        ## extract results from lm-file
        results <- topTable(
          lm,
          coef = 2, 
          number = Inf
        )
        results$ensembl_gene_id <- rownames(results)
        
        ## save raw results
        raw_results <- results
        write.csv(raw_results, file=paste0(getwd(), "/dm_results/continuous/raw/raw_",phenotype, "_", celltype, "_", biotype, "_", region, ".csv"))
        
        ## for promoters save protein coding raw results only
        if(biotype == "promoters"){
          ## get protein coding raw results only
          IDs <- raw_results$ensembl_gene_id
          genes_with_description <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","gene_biotype","description"), filters = 'ensembl_gene_id',values= IDs, mart=ensemblHuman)
          
          ## merge the two dfs
          raw_results <- merge(x = raw_results,y = genes_with_description, by = "ensembl_gene_id",all.x=TRUE)
          protein_coding <- raw_results %>% dplyr::filter(gene_biotype =="protein_coding")
          write.csv(protein_coding, file=paste0(getwd(), "/dm_results/continuous/raw/protein_coding/raw_protein_coding",phenotype, "_", celltype, "_", biotype, "_", region, ".csv"))
        }
        ## filter for significant ENS-IDs
        results_sig <- results[results$adj.P.Val < 0.05, ]
        results_sig <- results_sig[abs(results_sig$logFC) > 0.0001 , ]
        
               ## save significant results
        if(nrow(results_sig) > 0){
          ##print no of results 
          print(
            paste(
              "**** sig.results for",
              celltype,
              biotype,
              phenotype,
              region,
              nrow(results_sig),
              sep = " "
            )
          )
          ## add in gene description and gene biotype to ensemble id 
          if(biotype == "cpgs"){
            results_sig$external_gene_name <- "cpg"
            
          } else{
            IDs <- results_sig$ensembl_gene_id
            genes_with_description <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","gene_biotype","description"), filters = 'ensembl_gene_id',values= IDs, mart=ensemblHuman)
            
            ## merge the two dfs
            results_sig <- merge(x = results_sig,y = genes_with_description, by = "ensembl_gene_id",all.x=TRUE)
            write.csv(results_sig, file=paste0(getwd(), "/dm_results/continuous/sig/sig_",phenotype, "_", celltype, "_", biotype, "_", region, ".csv"))
      
            ## amongst all sig get top 50
            results_sig <- results_sig[order(abs(results_sig$logFC), decreasing = TRUE),] 
            if(nrow(results_sig) > 50) {
              results_top50 <- results_sig[1:50,] ## pick top 50 based on logFC
            } else{
              results_top50 <- results_sig
            }
            if(nrow(results_top50) > 0){
              write.csv(results_top50, file=paste0(getwd(), "/dm_results/continuous/top50/top50_",phenotype, "_", celltype, "_", biotype, "_", region, ".csv"))
            }
            
            ## get all protein coding sig results
           protein_coding <- results_sig %>% dplyr::filter(gene_biotype =="protein_coding")
            if(nrow(protein_coding) > 0){
            write.csv(protein_coding, file=paste0(getwd(), "/dm_results/continuous/sig/protein_coding/protein_coding_sig_",phenotype, "_", celltype, "_", biotype, "_", region, ".csv"))
            }
           
           ## amongst protein coding sig get top 50
           if(nrow(protein_coding) > 50) {
             protein_coding <- protein_coding[order(abs(protein_coding$logFC), decreasing = TRUE),] 
             protein_coding_top50 <- protein_coding[1:50,] ## pick top 50 based on logFC
           } else{
             protein_coding <- protein_coding[order(abs(protein_coding$logFC), decreasing = TRUE),] 
             protein_coding_top50 <- protein_coding
           }
           if(nrow(protein_coding_top50) > 0){
             write.csv(protein_coding_top50, file=paste0(getwd(), "/dm_results/continuous/top50/protein_coding/protein_coding_top50_",phenotype, "_", celltype, "_", biotype, "_", region, ".csv"))
           }
           }
        }
        
        
        ## create df only continaing number of significant hits for plotting later on
        name <- paste0(phenotype)
        sig_dms_per_brainregion[name] <- nrow(results_sig)
        

        ## save df containing number of sig. dms for each celltype, region, biotype
        saveRDS(sig_dms_per_brainregion,
                file = paste(baseDir, "/dm_results/number_sig_dms/",
                  paste(
                    "sig_dms_per_region",
                    header,
                    celltype,
                    biotype,
                    region,
                    sep = "_"
                  ),
                  ".rds",
                  sep = ""
                ))
      }

    }
  }
}
