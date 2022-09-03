## Plot Online Ressource 7, Figure of neurons in the dentate gyrus, gamma secretase

library(cowplot)
library(magrittr)
library(tidyverse)
library(Cairo)
library(biomaRt)
library(RnBeads)

setwd("~/90plus")
targets <- readRDS("./data/targets_incl_comorb.rds")
PEN2 <- "ENSG00000205155"

## read in raw results
raw_results <- read.csv("./dm_results/continuous/raw/raw_niaaaascore_neuron_promoters_DG.csv")
IDs <- raw_results$ensembl_gene_id
genes_with_description <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","gene_biotype","description"), filters = 'ensembl_gene_id',values= IDs, mart=ensemblHuman)

## merge the two dfs
raw_results <- merge(x = raw_results,y = genes_with_description, by = "ensembl_gene_id",all.x=TRUE)
APH1A <- raw_results %>% filter(raw_results$external_gene_name == 'APH1A')
NCSTN <- raw_results %>% filter(raw_results$external_gene_name == 'NCSTN')
PSEN1<- raw_results %>% filter(raw_results$external_gene_name == 'PSEN1')
PSEN2 <- raw_results %>% filter(raw_results$external_gene_name == 'PSEN2')
PEN2 <- raw_results %>% filter(raw_results$external_gene_name == 'PSENEN')
gamma_secretase <- rbind(APH1A, PSEN1, PSEN2, PEN2, NCSTN)

## APH1A
neuron_promoters <- read_rds("./cleaned_mvals/cleaned_mvals_episcore/niaaaascore/neuron/neuron/corrected_pvals/neuron_promoters_DG_niaaaascore_cleaned_mvals.rds")
neuron_promoters_APH1A <- neuron_promoters[APH1A$ensembl_gene_id,]
neuron_promoters_APH1A <- rnb.mval2beta(neuron_promoters_APH1A)
neuron_promoters_APH1A <- tibble(sample_name = names(neuron_promoters_APH1A), Methylation = neuron_promoters_APH1A, gene = "APH-1A", adj_p = APH1A$adj.P.Val)  %>% 
  left_join(targets, by= "sample_name")
neuron_promoters_APH1A$niaaaascore %<>% as.numeric
neuron_promoters_APH1A$ensembl_gene_id <- APH1A$ensembl_gene_id

## NCSTN
neuron_promoters_NCSTN <- neuron_promoters[NCSTN$ensembl_gene_id,]
neuron_promoters_NCSTN <- rnb.mval2beta(neuron_promoters_NCSTN)
neuron_promoters_NCSTN <- tibble(sample_name = names(neuron_promoters_NCSTN), Methylation = neuron_promoters_NCSTN, gene = "NCSTN **",adj_p = NCSTN$adj.P.Val)  %>% 
  left_join(targets, by= "sample_name")
neuron_promoters_NCSTN$niaaaascore %<>% as.numeric
neuron_promoters_NCSTN$ensembl_gene_id <- NCSTN$ensembl_gene_id


## PSEN1
neuron_promoters_PSEN1 <- neuron_promoters[PSEN1$ensembl_gene_id,]
neuron_promoters_PSEN1 <- rnb.mval2beta(neuron_promoters_PSEN1)
neuron_promoters_PSEN1 <- tibble(sample_name = names(neuron_promoters_PSEN1), Methylation = neuron_promoters_PSEN1, gene = "PSEN-1",adj_p = PSEN1$adj.P.Val)  %>% 
  left_join(targets, by= "sample_name")
neuron_promoters_PSEN1$niaaaascore %<>% as.numeric
neuron_promoters_PSEN1$ensembl_gene_id <- PSEN1$ensembl_gene_id

## PSEN2
neuron_promoters_PSEN2 <- neuron_promoters[PSEN2$ensembl_gene_id,]
neuron_promoters_PSEN2 <- rnb.mval2beta(neuron_promoters_PSEN2)
neuron_promoters_PSEN2 <- tibble(sample_name = names(neuron_promoters_PSEN2), Methylation = neuron_promoters_PSEN2, gene = "PSEN-2",adj_p = PSEN2$adj.P.Val)  %>% 
  left_join(targets, by= "sample_name")
neuron_promoters_PSEN2$niaaaascore %<>% as.numeric
neuron_promoters_PSEN2$ensembl_gene_id <- PSEN2$ensembl_gene_id


## PEN2
neuron_promoters_PEN2 <- neuron_promoters[PEN2$ensembl_gene_id,]
neuron_promoters_PEN2 <- rnb.mval2beta(neuron_promoters_PEN2)
neuron_promoters_PEN2 <- tibble(sample_name = names(neuron_promoters_PEN2), Methylation = neuron_promoters_PEN2, gene = "PEN-2 **",adj_p = PEN2$adj.P.Val)  %>% 
  left_join(targets, by= "sample_name")
neuron_promoters_PEN2$niaaaascore %<>% as.numeric
neuron_promoters_PEN2$ensembl_gene_id <- PEN2$ensembl_gene_id


## combine all dataframes
joined_df <- rbind(neuron_promoters_PEN2,
                   neuron_promoters_PSEN1, 
                   neuron_promoters_PSEN2,
                   neuron_promoters_APH1A,
                   neuron_promoters_NCSTN)
joined_df$gene <- factor(joined_df$gene, levels = c("APH-1A",
                                                            "NCSTN **",
                                                            "PEN-2 **",
                                                            "PSEN-1",
                                                            "PSEN-2"))
joined_df <- joined_df  %>% as.data.frame() %>% left_join(gamma_secretase, by=  "ensembl_gene_id")

## plot 
p <- ggplot(joined_df, aes(niaaaascore, Methylation)) +
  stat_smooth(method = "lm",colour="grey") +
  geom_point() + 
  annotate(geom="text", x=2, y=2.5, label=  paste0("adj.p. = ", round(joined_df$adj_p,3)),
           color="black", size=4)+
  facet_wrap(~joined_df$gene)+
  scale_x_continuous(breaks = c(1, 2, 3)) +
  theme_classic() +
  theme(strip.background = element_rect(colour="black",
                                        fill="white"), 
        text = element_text(size = 16))+
  theme(panel.border = element_rect(fill = NA),
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  geom_text(
    size    = 3,
    color = "black",
    data    = joined_df,
    mapping = aes(
      x = -Inf,
      y = -Inf,
      label = paste("FDR p < ", round(adj.P.Val, digits=3)),
      hjust   = -0.1,
      vjust   = -0.8
    )
  )+
  xlab("NIA-AA A score") +
  ylab("Methylation beta value")+
  scale_y_continuous(breaks=seq(0,1,0.5), limits=c(0,1))

ggsave(
  paste0(baseDir,"/plots/final_figures/Supp_Fig_7_DG_neuron_gamma_secretase.png"),
  p,
  width =174,
  height = 100,
  dpi = 1200,
  units= 'mm'
)


