## code to plot scatterplot displayng PEN-2 methylation m-values f
## or all four neurpath scores in DG neurons

library(cowplot)
library(magrittr)
library(tidyverse)
library(Cairo)
library(biomaRt)
library(dplyr)
library(tidyr)
library(ggpubr)
library(grid)
library(RnBeads)


setwd("~/90plus")

## set up ensembl mart
ensembl <- useMart("ensembl")
ensemblHuman <- useDataset("hsapiens_gene_ensembl",mart = ensembl)

## load sample sheet
targets <- readRDS("~/90plus/data/targets_incl_episcore.rds")
targets$int_adseverityscore <- targets$adseverityscore
targets$int_adseverityscore <- sub("intermediate", "2",targets$int_adseverityscore)
targets$int_adseverityscore <- sub("low", "1",targets$int_adseverityscore)
targets$int_adseverityscore <- sub("high", "3",targets$int_adseverityscore)
targets<- targets %>% dplyr::select(c("sample_name", "sample_source", "block_id","niaaaascore", "niaaabscore", "niaaacscore", "int_adseverityscore"))

## define ensembl ID of PEN2
PEN2 <- "ENSG00000205155"

## define comparisons and cell type
comparisons <- c('niaaaascore', 'niaaabscore', 'niaaacscore','int_adseverityscore')
celltype = 'neuron'

## plot
myplots <- list()
for(comparison in comparisons){
  # read in raw results
  raw_results <- read.csv(paste0("~/90plus/dm_results/continuous/raw/raw_",comparison,"_neuron_promoters_DEN.csv"))
  IDs <- raw_results$ensembl_gene_id
  genes_with_description <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","gene_biotype","description"), filters = 'ensembl_gene_id',values= IDs, mart=ensemblHuman)
  # merge the two dfs
  raw_results <- merge(x = raw_results,y = genes_with_description, by = "ensembl_gene_id",all.x=TRUE)
  PEN2 <- raw_results %>% filter(raw_results$external_gene_name == 'PSENEN')
  PEN2 <- PEN2 %>% dplyr::select("ensembl_gene_id", "adj.P.Val")
  neuron_promoters <- readRDS(paste0("~/90plus/cleaned_mvals/cleaned_mvals_episcore/", comparison, "/neuron/neuron/corrected_pvals/neuron_promoters_DEN_", comparison, "_cleaned_mvals.rds"))
  neuron_promoters <- rnb.mval2beta(neuron_promoters)
  neuron_promoters <- neuron_promoters[PEN2$ensembl_gene_id,] %>% as.data.frame()
  colnames(neuron_promoters)[1] <- 'beta'
  neuron_promoters[["sample_name"]] <- rownames(neuron_promoters)
  neuron_promoters$ensembl_gene_id <- PEN2$ensembl_gene_id
  neuron_promoters <- neuron_promoters  %>% as.data.frame() %>% left_join(PEN2, by=  "ensembl_gene_id")
  df <- left_join(neuron_promoters, targets, by="sample_name")
  df[[comparison]] <- df[[comparison]] %>% as.numeric
  colnames(df) <- c("beta", "sample_name", "ensembl_gene_id", "adj.P.Val", "NIA-AA A score", "NIA-AA B score", "NIA-AA C score", "AD severity score")
  if (comparison == "niaaaascore"){
    comparison <- "NIA-AA A score"
  }else if(comparison == "niaaabscore"){
    comparison <- "NIA-AA B score"
  }else if(comparison == "niaaacscore"){
    comparison <- "NIA-AA C score"
  }else if(comparison == "int_adseverityscore"){
    comparison <- "AD severity score"
  }
  
  p <- ggplot(df, aes(x = df[[comparison]], y = beta)) +
    stat_smooth(method = "lm",colour="grey") +
    geom_point(aes(alpha=0.3)) +
    theme_bw() +
    scale_color_manual(values="grey")+
    geom_text(
      size    = 3,
      color = "black",
      data    = df,
      mapping = aes(
        x = -Inf,
        y = -Inf,
        label = paste("FDR p < ", round(adj.P.Val, digits=3)),
        hjust   = -0.1,
        vjust   = -0.8
      )
    ) +
    ylim(0,1)+
    scale_y_continuous(breaks=seq(0,1,0.5), limits=c(0,1))+
    ylab("Methylation beta value")+
    xlab(paste0(comparison))+
    theme(legend.position = "none")+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))
  
  myplots[[comparison]] <- ggplotGrob(p)
}
dev.off()
png(width=6,height=5, units ='in', res=800,paste0("./plots/final_figures/OR_6_PEN2_DG_neuron_all_scores.png"))
plot<- 
  ggarrange(plotlist = myplots,
            labels = c("a", "b", "c", "d"),
            ncol = 2, nrow =2)
## add title
plot <- annotate_figure(plot, top = text_grob("PEN-2 methylation in DG neurons by AD neuropathology", 
                                      color = "black", size = 14))
print(plot)
dev.off()
