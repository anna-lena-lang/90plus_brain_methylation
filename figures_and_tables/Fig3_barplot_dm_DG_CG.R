## Figure 3
## code to plot barplots displaying number of significant promoter-associated
## differentially methylated Ensembl-IDs for dentate gyrus (DG) and cingulate gyrus (CG)

library(tidyr)
library(viridis)
library(ggpubr)
library(dplyr)
library(RColorBrewer)
library(ggplot2)

## set working directory
setwd("~/90plus")
savedir <- getwd()

## define phenotypes and cell types 
phenotypes <- c("int_adseverityscore", "niaaaascore", "niaaabscore", "niaaacscore") 
celltypes <- c('bulk','neuron', 'astro', 'endo', 'oligo_opc', 'microglia')

## function -----
## will plot barplot for cingulate gyrus (CG) and dentae gyrus (DG) 
## displaying number of significant differentially methlyated sites by neuropath
## colorcode the bars by cell types
dm_barplot_by_score_DG_CG <- function(biotype = "promoters", phenotypes, celltypes){
  big_df <- NULL
  all_regions<- c("DG", "CG")
  for (phenotype in phenotypes){
    for (celltype in celltypes){
      for (region in all_regions){
        filename <- c(paste(savedir, paste("/dm_results/continuous/sig/protein_coding/protein_coding_sig", phenotype, celltype, biotype, region,sep="_"),".csv", sep=""))
        if (file.exists(filename)){
          number_of_sig_dms <- ''
          brainregion <- ''
          df <- data.frame(number_of_sig_dms, brainregion)
          df_full <- read.csv(filename)
          df$number_of_sig_dms <- nrow(df_full)
          df$brainregion <- region
          df$phenotypes <- phenotype
          df$celltype <- celltype
          df$biotype <- biotype
          df$data_to_keep <- paste0(biotype, "_",celltype)
          df$data_to_keep <- paste0(df$brainregion, "_", df$data_to_keep)
          big_df <- rbind(big_df, df)
        } else{
          number_of_sig_dms <- ''
          brainregion <- ''
          df <- data.frame(number_of_sig_dms, brainregion)
          df$number_of_sig_dms <- 0
          df$brainregion <- region
          df$phenotypes <- phenotype
          df$celltype <- celltype
          df$biotype <- biotype
          df$data_to_keep <- paste0(biotype, "_",celltype)
          df$data_to_keep <- paste0(df$brainregion, "_", df$data_to_keep)
          big_df <- rbind(big_df, df)
        }
      }
    }
  }
  ## define colors
  n_colors <- 4
  bio <- viridis_pal(option = "D")(n_colors)
  big_df$brainregion <- factor(big_df$brainregion, levels =  c("CG", "DG"))
  big_df$phenotypes <- as.factor(big_df$phenotypes)
  big_df$celltype <- sub("bulk", "Bulk", big_df$celltype)
  big_df$celltype <- sub("neuron", "Neurons", big_df$celltype)
  big_df$celltype <- sub("astro", "Astrocytes", big_df$celltype)
  big_df$celltype <- sub("endo", "Endothelial cells", big_df$celltype)
  big_df$celltype <- sub("microglia", "Microglia", big_df$celltype)
  big_df$celltype <- sub("oligo_opc", "Olig/OPCs", big_df$celltype)
  
  ## create column for samples that have more than 0 significant dms so that not all of the barplots get labels
  big_df$dms_label <- NA
  big_df$dms_label[big_df$number_of_sig_dms > 0] <- big_df$number_of_sig_dms[big_df$number_of_sig_dms > 0]
  
  colorblind_celltype <- c("#000000", "#648FFF", "#FE6100", "#785EF0", "#FFB000")
  ## remove data where celltype proportions were low
  datatokeep <- read.csv("~/90plus/cleaned_mvals/cleaned_mvals_episcore/episcore_data_to_keep/data_to_keep.csv")
  big_df <- big_df[big_df$data_to_keep %in% datatokeep$keep,]
  ## rename
  rename <- c(niaaaascore = "NIA-AA A score", niaaabscore = "NIA-AA B score",niaaacscore = "NIA-AA C score",int_adseverityscore= "AD severity score", CG = "CG",DG = "DG")
  ## refactor
  big_df$celltype <- factor(big_df$celltype, levels = c("Bulk", "Neurons", "Astrocytes", "Endothelial cells", "Olig/OPCs") )
  ## plot
  p <- ggplot(big_df, aes(x = phenotypes, y = number_of_sig_dms, fill = celltype)) + 
    geom_bar(aes(fill= celltype), position = "dodge", stat = "identity", alpha=1, width = 40)+
    scale_fill_manual(values=colorblind_celltype, name = "Cell type")+
    ylab("number of DM promoters")+
    xlab("")+
    facet_grid(brainregion ~  phenotypes, labeller = as_labeller(rename))+
    scale_color_manual(values=colorblind_celltype, guide = "none")+
    geom_text(aes(label = big_df$dms_label,  color = celltype, size =8), position=position_dodge(width=40), size=3,vjust=0)+
    scale_y_continuous(expand = c(0, 0), limits = c(0, (max(big_df$number_of_sig_dms)+ (max(big_df$number_of_sig_dms)/8))))+ 
    theme_minimal()+
    theme(
      strip.background = element_rect(colour="black",
                                      fill="white"), 
      text = element_text(size = 12), 
      axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = element_text(size=8),
      panel.border = element_rect(colour = "black", fill=NA, size=0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())
  
  
  
  ggsave(
    paste0(savedir,"/plots/final_figures/Fig3_dm_barplot_protein_coding_promoters_CG_DG.png"),
    p,
    width =8,
    height = 3,
    dpi = 380
  )
  
  pdf(file = paste0(savedir,"/plots/final_figures/Fig3_dm_barplot_protein_coding_promoters_CG_DG.pdf"),width=10, height=5)
  print(p)
  dev.off()
}


## plot -----
dm_barplot_by_score_DG_CG(biotype = "promoters", phenotypes, celltypes)
