## function to run differential methylation analysis for categorical variables

dm_categorical <- function(header, phenotype) {
  
  ## targets file location
  targets_file <- "targets_incl_comorb.rds"
  
  ## load targets file 
  all_targets <- readRDS(targets_file)

  ## dm analysis for each celltype, brain region and biotype combination:
  for (region in all_regions) {
    ## open file that knows which celltypes were present for each brain region
    celltypespresent <- readRDS(paste0(baseDir,"/cleaned_mvals/cleaned_mvals_episcore/episcore_celltypes_present/",region, "celltypes_present.rds"))
    celltypespresent <- c("bulk", celltypespresent)
    
    for (celltype in celltypespresent) {
      
      ## filter targets for this brain region
      targets <- all_targets %>% filter(brain_region == region)
      
      for (biotype in biotypes) {
        
        ## build empty dataframe to save dm results in
        sig_dms_per_brainregion <- data.frame(matrix(nrow = 1))
        colnames(sig_dms_per_brainregion)[1] <- "brain_region"
        sig_dms_per_brainregion["brain_region"] <- region
        
        ## get mvals
        mvals_file <- paste0(getwd(),"/cleaned_mvals/cleaned_mvals_episcore/", phenotype, "/",celltype,"/",celltype,"/corrected_pvals/",   celltype,"_" ,biotype,"_", region, "_", phenotype,  "_cleaned_mvals.rds")
        mvals <- readRDS(mvals_file)
        
        ## subset targets to match mvals df
        targets <- all_targets[all_targets$sample_name %in% colnames(mvals),]
        
        ## bring in right format
        targets[sapply(targets, is.character)] <- lapply(targets[sapply(targets, is.character)], as.factor)
        
        ## ---- dm Analysis with limma ----
        ## add new column for easier handling
        targets$predictor <- targets[[phenotype]]
        
        ## remove NAs and re-level factors and numerics
        targets_sub <- targets[!(is.na(targets$predictor)), ]
        targets_sub$agedeath <- as.character(targets_sub$agedeath) %>% as.numeric()
        targets_sub$predictor <- as.factor(targets_sub$predictor)
        targets_sub$gender <- as.character(targets_sub$gender) %>% as.factor()
        
        ## sort data
        mvals_sub <- mvals[, colnames(mvals) %in% targets_sub$sample_name]
        targets_sub <- targets_sub[order(targets_sub$sample_name),]
        mvals_sub <- mvals_sub[,order(colnames(mvals_sub))]
        
        if(all(colnames(mvals_sub) == targets_sub$sample_name) == FALSE){
          stop("ERROR")
        }
        
        ## create model
        model_content <-  paste("~0+ predictor + agedeath + gender", collapse = " + ")
        
        ## design model 
        design <- model.matrix(as.formula(model_content), data = targets_sub)
        
        if(phenotype == "adss_group"){
          ## for paper-diag
          colnames(design)[1:2] <- sub("predictor", "", colnames(design)[1:2])
          contMatrix <- makeContrasts(low_intermediate - high,
                                      levels = design)
        } else if(phenotype == "conf_diag"){
          colnames(design)[1:3] <- sub("predictor", "", colnames(design)[1:3])
          contMatrix <- makeContrasts(normal - CIND,
                                      normal - demented,
                                      CIND - demented,
                                      levels = design)
        }else if(phenotype == "niaaaascore"){
          colnames(design)[1:3] <- sub("predictor", "niaaaa_", colnames(design)[1:3])
          contMatrix <- makeContrasts(niaaaa_1 - niaaaa_2,
                                      niaaaa_1 - niaaaa_3,
                                      niaaaa_2 - niaaaa_3,
                                      levels = design)
        }else if(phenotype == "niaaabscore"){
          colnames(design)[1:2] <- sub("predictor", "niaaab_", colnames(design)[1:2])
          contMatrix <- makeContrasts(niaaab_2 - niaaab_3,
                                      levels = design)
        }else if(phenotype == "niaaacscore"){
          colnames(design)[1:3] <- sub("predictor", "niaaac_", colnames(design)[1:3])
          contMatrix <- makeContrasts(niaaac_1 - niaaac_2,
                                      niaaac_1 - niaaac_3,
                                      niaaac_2 - niaaac_3,
                                      levels = design)
        }
        ## run linear model
        fit <- lmFit(mvals_sub, design)
        lmfit <- contrasts.fit(fit, contMatrix)
        
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
        remove <- number_svs_total -number_svs_protected
        
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
        for (cont in colnames(lm$contrasts)) {
          results <- topTable(
            lm,
            coef = cont,
            number = Inf
          )
          
          results$ensembl_gene_id <- rownames(results)
          
          ## save raw results
          raw_results <- results
          write.csv(raw_results, file=paste0(getwd(), "/dm_results/categorical/raw/raw_",phenotype, "_", celltype, "_", biotype, "_", region, "_", cont, ".csv"))
          
          ## filter for significant ENS-IDs
          results_sig <- results[results$adj.P.Val < 0.05, ]
          results_sig <- results_sig[abs(results_sig$logFC) > 0.0001 , ]
          
          ## add in gene description and gene biotype to ensemble id 
      
          
          ## save significant results
          if(nrow(results_sig) >0){
            print(
              paste(
                "**** sig.results for",
                celltype,
                biotype,
                phenotype,
                region,
                "contrast:",
                cont,
                ":nrow =" ,
                nrow(results_sig),
                sep = " "
              ))
            if(biotype == "cpgs"){
              results_sig$external_gene_name <- NA
            } else{
            ## add annotation to ensemble-ids
              IDs <- results_sig$ensembl_gene_id
              genes_with_description <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","gene_biotype","description"), filters = 'ensembl_gene_id',values= IDs, mart=ensemblHuman)
              
              ## merge the two dfs
              results_sig <- merge(x = results_sig,y = genes_with_description, by = "ensembl_gene_id",all.x=TRUE)
              write.csv(results_sig, file=paste0(getwd(), "/dm_results/categorical/sig/sig_",phenotype, "_", celltype, "_", biotype, "_", region, "_", cont, ".csv"))
              
              ## amongst all get top 50
              results_sig <- results_sig[order(abs(results_sig$logFC), decreasing = TRUE),] 
              if(nrow(results_sig) > 50) {
                results_top50 <- results_sig[1:50,] ## pick top 50 based on logFC
              } else{
                results_top50 <- results_sig
              }
              if(nrow(results_top50) > 0){
                write.csv(results_top50, file=paste0(getwd(), "/dm_results/categorical/top50/top50_",phenotype, "_", celltype, "_", biotype, "_", region, "_", cont, ".csv"))
              }
              
              ## get all protein coding results
              protein_coding <- results_sig %>% filter(gene_biotype =="protein_coding")
              if(nrow(protein_coding) >0){
                write.csv(protein_coding, file=paste0(getwd(), "/dm_results/categorical/sig/protein_coding/protein_coding_sig_results_",phenotype, "_", celltype, "_", biotype, "_", region, "_", cont,".csv"))
              }
              
              ## amongst protein coding get top 50
              if(nrow(protein_coding) > 50) {
                protein_coding <- protein_coding[order(abs(protein_coding$logFC), decreasing = TRUE),] 
                protein_coding_top50 <- protein_coding[1:50,] ## pick top 50 based on logFC
              } else{
                protein_coding <- protein_coding[order(abs(protein_coding$logFC), decreasing = TRUE),] 
                protein_coding_top50 <- protein_coding
              }
              if(nrow(protein_coding_top50) > 0){
                write.csv(protein_coding_top50, file=paste0(getwd(), "/dm_results/categorical/top50/protein_coding/protein_coding_top50_",phenotype, "_", celltype, "_", biotype, "_", region, "_", cont, ".csv"))
              }
            }
          }
          ## create df only continaing number of significant hits for plotting later on
          name <- paste0(cont)
          sig_dms_per_brainregion[name] <- nrow(results_sig)
                    } 
        }
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
