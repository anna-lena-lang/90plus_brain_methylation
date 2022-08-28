## function to add number and exact cpgs (illumina cpg ids) to each ENSEMBL-ID
## will also add the hgnc symbol (external_gene_name), ...
## ...biotype (e.g. protein coding) and description for the hgnc symbol
## this is for promoter associated (biotype = promoter) and protein coding ENS-IDs only
library(tidyverse)
library(dplyr)
library(biomaRt)


regions <- c("DEN", "ERC", "CBL", "SN", "LOC", "CA1", "CG", "MFG")
celltypes <- c("neuron", "oligo_opc", "endo", "astro", "microglia")
biotype <- c("promoters")

## function to add number and exact cpgs to each ENSEMBL-ID from the protein coding Ensembl IDS.
## will also add the hgnc symbol (external_gene_name), biotype (e.g. protein coding) and description for the hgnc symbol
add_cpgs_to_protein_coding_results<- function(category){
  ## load df that displays ENS_ID and corresponding cpg-probes
  ENS_to_cpg <- readRDS(paste0(getwd(),"/ENS_to_cpg_lists/ENS_to_cpg_", biotype, ".rds")) %>% tibble()
  ## set up mart to access biomaRt annotation
  directory <- paste0(getwd(),"/dm_results/", category, "/sig/protein_coding/") 
  files <- list.files(directory, pattern = paste0("protein_coding_.+", ".+promoter.+", ".+csv"))
  if(identical(files, character(0)) == FALSE){
    
    for(filename in files){
     ## open df containing results of differentially methylated ENS-IDs
        results <- read.csv(paste0(directory, filename))
        results$ENS <- results$ensembl_gene_id
        sub_ENS <- subset(ENS_to_cpg, ENS %in% results$ensembl_gene_id)
        
        ## combine all cpgs corresponding to one ENS ID in one cell, separated by semicolon
        output <- sub_ENS %>% 
          group_by(ENS)  %>%
          summarise_all(list(~toString(unique(.))))
        
        ## left join to add the info with cpgs to the results table
        combo <- merge(x=results,y=output,by="ENS",all.x=TRUE)
        
        ## order based on logFC
        combo <- combo[order(abs(combo$logFC), decreasing = TRUE),]
        
        if(nrow(combo) > 0){
        write.csv(combo, paste0(getwd(),"/dm_results/", category, "/sig/protein_coding/protein_coding_incl_cpgs/incl_cpgs_", filename))
        }
    }
    }
}

