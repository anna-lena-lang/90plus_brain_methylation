## code to plot Supp_Fig_3 figure with overview barplot 
## visualizing number of significant differentially methylated ENS-IDs
## for all brain regions and cell types
## protein coding promoter regions only 

library(tidyr)
library(viridis)
library(ggpubr)
library(dplyr)
library(RColorBrewer)
library(ggplot2)

## set working directory
setwd("~/90plus")
baseDir <- getwd()
all_regions <- c("MFG", "CG", "DG", "CA1", "EC", "LC", "SN", "CBM")
phenotypes <- c("int_adseverityscore", "niaaaascore", "niaaabscore", "niaaacscore") 
celltypes <- c('bulk','neuron', 'astro', 'endo', 'oligo_opc', 'microglia')

##  function
dm_barplots_overview <- function(biotype = "promoters"){
  big_df <- NULL
  for (phenotype in phenotypes){
    header <- paste0("episcore_continuous_", phenotype)
    for (celltype in celltypes){
      for (region in all_regions){
        filename <- c(paste(baseDir, paste("/dm_results/continuous/sig/protein_coding/protein_coding_sig", phenotype, celltype, biotype, region,sep="_"),".csv", sep=""))
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
  # define colors
  n_colors <- 4
  bio <- viridis_pal(option = "D")(n_colors)
  big_df$brainregion <- factor(big_df$brainregion, levels =  c("MFG", "CG", "CA1", "DG", "EC", "LC", "SN", "CBM"))
  big_df$phenotypes <- sub("niaaaascore", "NIA-AA A score", big_df$phenotypes)
  big_df$phenotypes <- sub("niaaabscore", "NIA-AA B score", big_df$phenotypes)
  big_df$phenotypes <- sub("niaaacscore", "NIA-AA C score", big_df$phenotypes)
  big_df$phenotypes <- sub("int_adseverityscore", "AD severity score", big_df$phenotypes)
  
  
  # create column for samples that have more than 0 significant dms so that not all of the barplots get labels
  big_df$dms_label <- NA
  big_df$dms_label[big_df$number_of_sig_dms > 0] <- big_df$number_of_sig_dms[big_df$number_of_sig_dms > 0]
  colorBlind4 <- c("#E69F00", "#56B4E9", "#009E73",  "#CC79A7")
  # remove data where celltype proportions were low
  datatokeep <- read.csv("~/90plus/cleaned_mvals/cleaned_mvals_episcore/episcore_data_to_keep/data_to_keep.csv")
  big_df <- big_df[big_df$data_to_keep %in% datatokeep$keep,]
  rename <- c(bulk = "Bulk", astro = "Astrocytes",neuron = "Neurons",endo="Endothelial", oligo_opc = "Olig/OPCs", microglia = "Microglia",
              MFG = "MFG", CG = "CG", CA1= "CA1", DG = "DG", EC = "EC", LC = "LC", SN= "SN", CBM = "CBM")
  # remove tdp43, lb, mvlscore
  #big_df <- big_df[!(big_df$phenotypes %in% c("mvl_score", "braak_for_lb_orig", "tdp43")),]
  big_df$celltype <- factor(big_df$celltype, levels = c("bulk", "neuron", "astro", "endo", "oligo_opc", "microglia") )
  p <- ggplot(big_df, aes(x = phenotypes, y = number_of_sig_dms, fill = phenotypes)) + 
    geom_bar(aes(fill= phenotypes), position ="dodge", stat = "identity", alpha=0.5, width = 0.7)+
    scale_fill_manual(values=colorBlind4, name = "Score")+
    ylab("number of DM promoters")+
    xlab("")+
    facet_grid(brainregion ~  celltype, labeller = as_labeller(rename),
               scales = 'free', space = 'free')+
    scale_color_manual(values=colorBlind4, guide = "none")+
    scale_y_continuous(expand = c(0, 0), limits = c(0, (max(big_df$number_of_sig_dms)+ (max(big_df$number_of_sig_dms)/8))))+ 
    geom_text(aes(label = big_df$dms_label,  color = phenotypes, size =8), position=position_dodge(width=0.7), size=2, vjust=0)+
    theme_minimal()+
    theme(
      strip.background = element_rect(colour="black",
                                      fill="white"), 
      text = element_text(size = 10), 
      axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = element_text(size=8),
      panel.border = element_rect(colour = "black", fill=NA, size=0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())
  
  ggsave(
    paste0(baseDir,"/plots/final_figures/Supp_Fig_3_dm_barplot_overview.png"),
    p,
    width = 174,
    height = 160,
    dpi = 1200, 
    units = 'mm'
  )
}

## plot
dm_barplots_overview(biotype  = "promoters")
