## plot PSENEN for all brain regions, neurons only
library(cowplot)
library(magrittr)
library(tidyverse)
library(Cairo)
library(RnBeads)
setwd("~/90plus")
baseDir <- getwd()

targets <- readRDS("~/90plus/data/targets_incl_comorb.rds")
psenen <- "ENSG00000205155"
## set up ensembl mart
ensembl <- useMart("ensembl")
ensemblHuman <- useDataset("hsapiens_gene_ensembl",mart = ensembl)

## load data
neuron_promoters_DG <- read_rds("./cleaned_mvals/cleaned_mvals_episcore/niaaaascore/neuron/neuron/corrected_pvals/neuron_promoters_DG_niaaaascore_cleaned_mvals.rds")
neuron_promoters_CA1 <- read_rds("./cleaned_mvals/cleaned_mvals_episcore/niaaaascore/neuron/neuron/corrected_pvals/neuron_promoters_CA1_niaaaascore_cleaned_mvals.rds")
neuron_promoters_EC <- read_rds("./cleaned_mvals/cleaned_mvals_episcore/niaaaascore/neuron/neuron/corrected_pvals/neuron_promoters_EC_niaaaascore_cleaned_mvals.rds")
neuron_promoters_MFG <- read_rds("./cleaned_mvals/cleaned_mvals_episcore/niaaaascore/neuron/neuron/corrected_pvals/neuron_promoters_MFG_niaaaascore_cleaned_mvals.rds")
neuron_promoters_CBM <- read_rds("./cleaned_mvals/cleaned_mvals_episcore/niaaaascore/neuron/neuron/corrected_pvals/neuron_promoters_CBM_niaaaascore_cleaned_mvals.rds")
neuron_promoters_LC <- read_rds("./cleaned_mvals/cleaned_mvals_episcore/niaaaascore/neuron/neuron/corrected_pvals/neuron_promoters_LC_niaaaascore_cleaned_mvals.rds")
neuron_promoters_SN <- read_rds("./cleaned_mvals/cleaned_mvals_episcore/niaaaascore/neuron/neuron/corrected_pvals/neuron_promoters_SN_niaaaascore_cleaned_mvals.rds")
neuron_promoters_CG <- read_rds("./cleaned_mvals/cleaned_mvals_episcore/niaaaascore/neuron/neuron/corrected_pvals/neuron_promoters_CG_niaaaascore_cleaned_mvals.rds")

## DG
neuron_promoters_DG <- neuron_promoters_DG[psenen,]
neuron_promoters_DG <- rnb.mval2beta(neuron_promoters_DG)
neuron_promoters_DG <- tibble(sample_name = names(neuron_promoters_DG), Methylation = neuron_promoters_DG, celltype = "PEN-2 DG Neuron **",sig_color= "sig") %>%
  left_join(targets)
raw_results <- read.csv(paste0("~/90plus/dm_results/continuous/raw/raw_niaaaascore_neuron_promoters_DG.csv"))
IDs <- raw_results$ensembl_gene_id
genes_with_description <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","gene_biotype","description"), filters = 'ensembl_gene_id',values= IDs, mart=ensemblHuman)
raw_results <- merge(x = raw_results,y = genes_with_description, by = "ensembl_gene_id",all.x=TRUE)
raw_PEN2 <- raw_results %>% filter(raw_results$external_gene_name == 'PSENEN')
raw_PEN2 <- raw_PEN2 %>% dplyr::select("ensembl_gene_id", "adj.P.Val")
neuron_promoters_DG$ensembl_gene_id <- raw_PEN2$ensembl_gene_id
neuron_promoters_DG <- neuron_promoters_DG  %>% as.data.frame() %>% left_join(raw_PEN2, by=  "ensembl_gene_id")
neuron_promoters_DG$niaaaascore %<>% as.numeric

## CA1
neuron_promoters_CA1 <- neuron_promoters_CA1[psenen,]
neuron_promoters_CA1 <- rnb.mval2beta(neuron_promoters_CA1)
neuron_promoters_CA1 <- tibble(sample_name = names(neuron_promoters_CA1), Methylation = neuron_promoters_CA1, celltype = "PEN-2 CA1 Neuron",sig_color= "not_sig") %>%
  left_join(targets)
raw_results <- read.csv(paste0("O:/02182022_backup_Lena/90plus/dm_results/continuous/raw/raw_niaaaascore_neuron_promoters_CA1.csv"))
raw_results <- merge(x = raw_results,y = genes_with_description, by = "ensembl_gene_id",all.x=TRUE)
raw_PEN2 <- raw_results %>% filter(raw_results$external_gene_name == 'PSENEN')
raw_PEN2 <- raw_PEN2 %>% dplyr::select("ensembl_gene_id", "adj.P.Val")
neuron_promoters_CA1$ensembl_gene_id <- raw_PEN2$ensembl_gene_id
neuron_promoters_CA1 <- neuron_promoters_CA1  %>% as.data.frame() %>% left_join(raw_PEN2, by=  "ensembl_gene_id")
neuron_promoters_CA1$niaaaascore %<>% as.numeric

## EC
neuron_promoters_EC <- neuron_promoters_EC[psenen,]
neuron_promoters_EC <- rnb.mval2beta(neuron_promoters_EC)
neuron_promoters_EC <- tibble(sample_name = names(neuron_promoters_EC), Methylation = neuron_promoters_EC, celltype = "PEN-2 EC Neuron",sig_color= "not_sig") %>%
  left_join(targets)
neuron_promoters_EC$niaaaascore %<>% as.numeric
raw_results <- read.csv(paste0("~/90plus/dm_results/continuous/raw/raw_niaaaascore_neuron_promoters_EC.csv"))
raw_results <- merge(x = raw_results,y = genes_with_description, by = "ensembl_gene_id",all.x=TRUE)
raw_PEN2 <- raw_results %>% filter(raw_results$external_gene_name == 'PSENEN')
raw_PEN2 <- raw_PEN2 %>% dplyr::select("ensembl_gene_id", "adj.P.Val")
neuron_promoters_EC$ensembl_gene_id <- raw_PEN2$ensembl_gene_id
neuron_promoters_EC <- neuron_promoters_EC  %>% as.data.frame() %>% left_join(raw_PEN2, by=  "ensembl_gene_id")
neuron_promoters_EC$niaaaascore %<>% as.numeric

## MFG
neuron_promoters_MFG <- neuron_promoters_MFG[psenen,]
neuron_promoters_MFG <- rnb.mval2beta(neuron_promoters_MFG)
neuron_promoters_MFG <- tibble(sample_name = names(neuron_promoters_MFG), Methylation = neuron_promoters_MFG, celltype = "PEN-2 MFG Neuron",sig_color= "not_sig") %>%
  left_join(targets)
neuron_promoters_MFG$niaaaascore %<>% as.numeric
raw_results <- read.csv(paste0("~/90plus/dm_results/continuous/raw/raw_niaaaascore_neuron_promoters_MFG.csv"))
raw_results <- merge(x = raw_results,y = genes_with_description, by = "ensembl_gene_id",all.x=TRUE)
raw_PEN2 <- raw_results %>% filter(raw_results$external_gene_name == 'PSENEN')
raw_PEN2 <- raw_PEN2 %>% dplyr::select("ensembl_gene_id", "adj.P.Val")
neuron_promoters_MFG$ensembl_gene_id <- raw_PEN2$ensembl_gene_id
neuron_promoters_MFG <- neuron_promoters_MFG  %>% as.data.frame() %>% left_join(raw_PEN2, by=  "ensembl_gene_id")
neuron_promoters_MFG$niaaaascore %<>% as.numeric

## SN
neuron_promoters_SN <- neuron_promoters_SN[psenen,]
neuron_promoters_SN <- rnb.mval2beta(neuron_promoters_SN)
neuron_promoters_SN <- tibble(sample_name = names(neuron_promoters_SN), Methylation = neuron_promoters_SN, celltype = "PEN-2 SN Neuron",sig_color= "not_sig") %>%
  left_join(targets)
raw_results <- read.csv(paste0("~/90plus/dm_results/continuous/raw/raw_niaaaascore_neuron_promoters_SN.csv"))
raw_results <- merge(x = raw_results,y = genes_with_description, by = "ensembl_gene_id",all.x=TRUE)
raw_PEN2 <- raw_results %>% filter(raw_results$external_gene_name == 'PSENEN')
raw_PEN2 <- raw_PEN2 %>% dplyr::select("ensembl_gene_id", "adj.P.Val")
neuron_promoters_SN$ensembl_gene_id <- raw_PEN2$ensembl_gene_id
neuron_promoters_SN <- neuron_promoters_SN  %>% as.data.frame() %>% left_join(raw_PEN2, by=  "ensembl_gene_id")
neuron_promoters_SN$niaaaascore %<>% as.numeric

## LC
neuron_promoters_LC <- neuron_promoters_LC[psenen,]
neuron_promoters_LC <- rnb.mval2beta(neuron_promoters_LC)
neuron_promoters_LC <- tibble(sample_name = names(neuron_promoters_LC), Methylation = neuron_promoters_LC, celltype = "PEN-2 LC Neuron",sig_color= "not_sig") %>%
  left_join(targets)
raw_results <- read.csv(paste0("~/90plus/dm_results/continuous/raw/raw_niaaaascore_neuron_promoters_LC.csv"))
raw_results <- merge(x = raw_results,y = genes_with_description, by = "ensembl_gene_id",all.x=TRUE)
raw_PEN2 <- raw_results %>% filter(raw_results$external_gene_name == 'PSENEN')
raw_PEN2 <- raw_PEN2 %>% dplyr::select("ensembl_gene_id", "adj.P.Val")
neuron_promoters_LC$ensembl_gene_id <- raw_PEN2$ensembl_gene_id
neuron_promoters_LC <- neuron_promoters_LC  %>% as.data.frame() %>% left_join(raw_PEN2, by=  "ensembl_gene_id")
neuron_promoters_LC$niaaaascore %<>% as.numeric

## CBM
neuron_promoters_CBM <- neuron_promoters_CBM[psenen,]
neuron_promoters_CBM <- rnb.mval2beta(neuron_promoters_CBM)
neuron_promoters_CBM <- tibble(sample_name = names(neuron_promoters_CBM), Methylation = neuron_promoters_CBM, celltype = "PEN-2 CBM Neuron",sig_color= "not_sig") %>%
  left_join(targets)
raw_results <- read.csv(paste0("~/90plus/dm_results/continuous/raw/raw_niaaaascore_neuron_promoters_CBM.csv"))
raw_results <- merge(x = raw_results,y = genes_with_description, by = "ensembl_gene_id",all.x=TRUE)
raw_PEN2 <- raw_results %>% filter(raw_results$external_gene_name == 'PSENEN')
raw_PEN2 <- raw_PEN2 %>% dplyr::select("ensembl_gene_id", "adj.P.Val")
neuron_promoters_CBM$ensembl_gene_id <- raw_PEN2$ensembl_gene_id
neuron_promoters_CBM <- neuron_promoters_CBM  %>% as.data.frame() %>% left_join(raw_PEN2, by=  "ensembl_gene_id")
neuron_promoters_CBM$niaaaascore %<>% as.numeric

## CG
neuron_promoters_CG <- neuron_promoters_CG[psenen,]
neuron_promoters_CG <- rnb.mval2beta(neuron_promoters_CG)
neuron_promoters_CG <- tibble(sample_name = names(neuron_promoters_CG), Methylation = neuron_promoters_CG, celltype = "PEN-2 CG Neuron",sig_color= "not_sig") %>%
  left_join(targets)
raw_results <- read.csv(paste0("~/90plus/dm_results/continuous/raw/raw_niaaaascore_neuron_promoters_CG.csv"))
raw_results <- merge(x = raw_results,y = genes_with_description, by = "ensembl_gene_id",all.x=TRUE)
raw_PEN2 <- raw_results %>% filter(raw_results$external_gene_name == 'PSENEN')
raw_PEN2 <- raw_PEN2 %>% dplyr::select("ensembl_gene_id", "adj.P.Val")
neuron_promoters_CG$ensembl_gene_id <- raw_PEN2$ensembl_gene_id
neuron_promoters_CG <- neuron_promoters_CG  %>% as.data.frame() %>% left_join(raw_PEN2, by=  "ensembl_gene_id")
neuron_promoters_CG$niaaaascore %<>% as.numeric

## join
joined_df <- rbind(neuron_promoters_MFG,
                   neuron_promoters_CG, 
                   neuron_promoters_CBM,
                   neuron_promoters_LC,
                   neuron_promoters_SN,
                   neuron_promoters_EC,
                   neuron_promoters_CA1,
                   neuron_promoters_DG)
joined_df$celltype <- factor(joined_df$celltype, levels = c("PEN-2 MFG Neuron",
                                                            "PEN-2 CG Neuron",
                                                            "PEN-2 DG Neuron **",
                                                            "PEN-2 CA1 Neuron",
                                                            "PEN-2 EC Neuron",
                                                            "PEN-2 LC Neuron",
                                                            "PEN-2 SN Neuron",
                                                            "PEN-2 CBM Neuron"))
## plot
p_PEN <- ggplot(joined_df, aes(niaaaascore, Methylation)) +
  stat_smooth(method = "lm",colour="grey") +
  geom_point() + 
  facet_wrap(~joined_df$celltype, nrow=2)+
  scale_x_continuous(breaks = c(1, 2, 3)) +
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
  ) +
  theme_classic() +
  theme(strip.background = element_rect(colour="black",
                                        fill="white"), 
        text = element_text(size = 16))+
  theme(panel.border = element_rect(fill = NA),
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  xlab("NIA-AA A score") +
  ylab("Methylation beta value")+
  scale_y_continuous(breaks=seq(0,1,0.5), limits=c(0,1))

ggsave(
  paste0(baseDir,"/plots/final_figures/Supp_Fig_5a_PEN2_niaaaascore_neuron_all_brain_regions.png"),
  p_PEN,
  width = 174,
  height = 60,
  units = 'mm',
  dpi = 1200
)

dev.off()
