## plot barplot Online_Resource4 displaying number of signigicant
## differentially methylated protein coding promoter regions 
## for clinical diagnosis only

library(limma)
library(biomaRt)
library(ggpubr)
library(RnBeads)

setwd("O:/02182022_backup_Lena/90plus/")
baseDir <- paste0("O:/02182022_backup_Lena/90plus/")

## specify biotypes, regions
biotype <- "promoters"
all_regions <- c("MFG", "CG", "DEN", "CA1", "ERC", "LOC", "SN", "CBL")
celltypes <- c("bulk"    ,  "neuron"  ,  "astro"    , "endo"  ,    "oligo_opc" ,"microglia")

## function 
dm_barplots_by_celltype_categorical <- function(biotype = "promoter", comparisons){
  big_df <- NULL
  for (celltype in celltypes){
    for (region in all_regions){
      for (comparison in comparisons){
        directory <- paste0(getwd(),"/dm_results/categorical/sig/protein_coding/") 
        filename <- list.files(directory, pattern = paste0("protein_coding_sig_results_", phenotype, ".+", celltype,".+", biotype,".+", region, ".+", comparison, ".+csv"))
        if(identical(filename, character(0)) == FALSE){
          filename <- paste0(directory, filename)
          number_of_sig_dms <- ''
          brainregion <- ''
          df <- data.frame(number_of_sig_dms, brainregion)
          df_full <- read.csv(filename)
          df$number_of_sig_dms <- nrow(df_full)
          df$brainregion <- region
          df$phenotypes <- phenotype
          df$celltype <- celltype
          df$biotype <- biotype
          df$comparison <- comparison
          df$data_to_keep <- paste0(biotype, "_",celltype)
          df$data_to_keep <- paste0(df$brainregion, "_", df$data_to_keep)
          big_df <- rbind(big_df, df)
        }else{
          number_of_sig_dms <- ''
          brainregion <- ''
          df <- data.frame(number_of_sig_dms, brainregion)
          df$number_of_sig_dms <- 0
          df$brainregion <- region
          df$phenotypes <- phenotype
          df$celltype <- celltype
          df$biotype <- biotype
          df$comparison <- comparison
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
  big_df$brainregion <- sub("DEN", "DG", big_df$brainregion)
  big_df$brainregion <- sub("LOC", "LC", big_df$brainregion)
  big_df$brainregion <- sub("CBL", "CBM", big_df$brainregion)
  big_df$brainregion <- sub("ERC", "EC", big_df$brainregion)
  big_df$brainregion <- factor(big_df$brainregion, levels =  c("MFG", "CG", "CA1", "DG", "EC", "LC", "SN", "CBM"))
  rename <- c(bulk = "Bulk", astro = "Astrocytes",neuron = "Neurons",endo="Endothelial Cells", oligo_opc = "Olig/OPCs", microglia = "Microglia",
              MFG = "MFG", CG = "CG", CA1= "CA1", DG = "DG", EC = "EC", LC = "LC", SN= "SN", CBM = "CBM")
  big_df$celltype <- factor(big_df$celltype, levels = c("bulk", "neuron", "astro", "endo", "oligo_opc", "microglia") )
  # create column for samples that have more than 0 significant dms so that not all of the barplots get labels
  big_df$dms_label <- NA
  big_df$dms_label[big_df$number_of_sig_dms > 0] <- big_df$number_of_sig_dms[big_df$number_of_sig_dms > 0]
  # remove data where celltype proportions were low
  datatokeep <- read.csv("O:/02182022_backup_Lena/90plus/cleaned_mvals/cleaned_mvals_episcore/episcore_data_to_keep/data_to_keep.csv")
  big_df <- big_df[big_df$data_to_keep %in% datatokeep$keep,]
  colorBlind4   <- c("#E69F00", "#56B4E9", "#009E73",  "#CC79A7")
  
  p <- ggplot(big_df, aes(x = comparison, y = number_of_sig_dms, fill = comparison)) + 
    geom_bar(aes(fill= comparison), position ="dodge", stat = "identity", alpha=0.5, width = 0.7)+
    ylab("number of DM promoters")+
    facet_grid(brainregion ~  celltype, labeller = as_labeller(rename))+
    scale_fill_manual(values=colorBlind4, name = "Contrast")+
    scale_color_manual(values=colorBlind4, guide = "none")+
    xlab("")+
    
    # facet_wrap(~celltype, nrow=2)+
    geom_text(aes(label = big_df$dms_label,  color = comparison, size =8), position=position_dodge(width=0.7), size=2, vjust=0)+
    scale_y_continuous(expand = c(0, 0), limits = c(0, (max(big_df$number_of_sig_dms)+ (max(big_df$number_of_sig_dms)/8))))+ 
    # coord_cartesian(ylim=c(0, 1000))+
    # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x))+
    theme_minimal()+
    theme(
      strip.background = element_rect(colour="black",
                                      fill="white"), 
      text = element_text(size = 10), 
      axis.text.x = element_text(angle = 60, vjust = 1, 
                                 size = 8, hjust = 1),
      axis.ticks = element_blank(),
      axis.text.y = element_text(size=5),
      panel.border = element_rect(colour = "black", fill=NA, size=0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())
  
  ggsave(
    paste0(baseDir,"/plots/final_figures/OR_4_episcore_categorical_conf_diag_barplot.png"),
    p,
    width = 8,
    height = 9,
    dpi = 800
  )
  pdf(file=paste0(baseDir,"/plots/final_figures/OR_4_episcore_categorical_conf_diag_barplot.pdf"),
      width = 8,
      height = 9)
  print(p)
  dev.off()
}

## plot 
dm_barplots_by_celltype_categorical(biotype = "promoters", comparisons = c("CIND - demented", "normal - demented", "normal - CIND"))
