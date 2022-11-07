
## code to plot SupFig 10c: Scatterplot for CD2AP in dentate gyrus (DG) 
## for all cell types

library(cowplot)
library(magrittr)
library(tidyverse)
library(Cairo)
library(RnBeads)
library(biomaRt)

setwd("~/90plus")
baseDir <- getwd()

## set up ensembl mart
# ensembl <- useMart("ensembl")
# ensemblHuman <- useDataset("hsapiens_gene_ensembl",mart = ensembl)

## read in sample sheet
targets <- readRDS("targets_incl_comorb.rds")
## define ensembl ID for CD2AP
CD2AP <- "ENSG00000198087"

## load data for each cell type
bulk_promoters <- read_rds("./cleaned_mvals/cleaned_mvals_episcore/niaaaascore/bulk/bulk/corrected_pvals/bulk_promoters_DEN_niaaaascore_cleaned_mvals.rds")
neuron_promoters <- read_rds("./cleaned_mvals/cleaned_mvals_episcore/niaaaascore/neuron/neuron/corrected_pvals/neuron_promoters_DEN_niaaaascore_cleaned_mvals.rds")
astro_promoters <- read_rds("./cleaned_mvals/cleaned_mvals_episcore/niaaaascore/astro/astro/corrected_pvals/astro_promoters_DEN_niaaaascore_cleaned_mvals.rds")
oligo_opc_promoters <- read_rds("./cleaned_mvals/cleaned_mvals_episcore/niaaaascore/oligo_opc/oligo_opc/corrected_pvals/oligo_opc_promoters_DEN_niaaaascore_cleaned_mvals.rds")
endo_promoters <- read_rds("./cleaned_mvals/cleaned_mvals_episcore/niaaaascore/endo/endo/corrected_pvals/endo_promoters_DEN_niaaaascore_cleaned_mvals.rds")
microglia_promoters <- read_rds("./cleaned_mvals/cleaned_mvals_episcore/niaaaascore/microglia/microglia/corrected_pvals/microglia_promoters_DEN_niaaaascore_cleaned_mvals.rds")

## prepare data for each cell type 

## bulk
bulk_promoters_CD2AP <- bulk_promoters[CD2AP,]
bulk_promoters_CD2AP <- rnb.mval2beta(bulk_promoters_CD2AP)
bulk_promoters_CD2AP_df <- tibble(sample_name = names(bulk_promoters_CD2AP), Methylation = bulk_promoters_CD2AP, celltype = "CD2AP DG Bulk",sig_color= "not_sig") %>%
  left_join(targets)
raw_results <- read.csv(paste0("~/90plus/dm_results/continuous/raw/raw_niaaaascore_bulk_promoters_DEN.csv"))
IDs <- raw_results$ensembl_gene_id
genes_with_description <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","gene_biotype","description"), filters = 'ensembl_gene_id',values= IDs, mart=ensemblHuman)
raw_results <- merge(x = raw_results,y = genes_with_description, by = "ensembl_gene_id",all.x=TRUE)
raw_CD2AP <- raw_results %>% filter(raw_results$external_gene_name == 'CD2AP')
raw_CD2AP <- raw_CD2AP %>% dplyr::select("ensembl_gene_id", "adj.P.Val")
bulk_promoters_CD2AP_df$ensembl_gene_id <- raw_CD2AP$ensembl_gene_id
bulk_promoters_CD2AP_df <- bulk_promoters_CD2AP_df  %>% as.data.frame() %>% left_join(raw_CD2AP, by=  "ensembl_gene_id")
bulk_promoters_CD2AP_df$niaaaascore %<>% as.numeric

## neuron
neuron_promoters_CD2AP <- neuron_promoters[CD2AP,]
neuron_promoters_CD2AP <- rnb.mval2beta(neuron_promoters_CD2AP)
neuron_promoters_CD2AP_df <- tibble(sample_name = names(neuron_promoters_CD2AP), Methylation = neuron_promoters_CD2AP, celltype = "CD2AP DG Neuron **", sig_color= "sig") %>%
  left_join(targets)
raw_results <- read.csv(paste0("~/90plus/dm_results/continuous/raw/raw_niaaaascore_neuron_promoters_DEN.csv"))
raw_results <- merge(x = raw_results,y = genes_with_description, by = "ensembl_gene_id",all.x=TRUE)
raw_CD2AP <- raw_results %>% filter(raw_results$external_gene_name == 'CD2AP')
raw_CD2AP <- raw_CD2AP %>% dplyr::select("ensembl_gene_id", "adj.P.Val")
neuron_promoters_CD2AP_df$ensembl_gene_id <- raw_CD2AP$ensembl_gene_id
neuron_promoters_CD2AP_df <- neuron_promoters_CD2AP_df  %>% as.data.frame() %>% left_join(raw_CD2AP, by=  "ensembl_gene_id")
neuron_promoters_CD2AP_df$niaaaascore %<>% as.numeric

## astrocytes
astro_promoters_CD2AP <- astro_promoters[CD2AP,]
astro_promoters_CD2AP <- rnb.mval2beta(astro_promoters_CD2AP)
astro_promoters_CD2AP_df <- tibble(sample_name = names(astro_promoters_CD2AP), Methylation = astro_promoters_CD2AP, celltype = "CD2AP DG Astrocytes",sig_color= "not_sig") %>%
  left_join(targets)
raw_results <- read.csv(paste0("~/90plus/dm_results/continuous/raw/raw_niaaaascore_astro_promoters_DEN.csv"))
raw_results <- merge(x = raw_results,y = genes_with_description, by = "ensembl_gene_id",all.x=TRUE)
raw_CD2AP <- raw_results %>% filter(raw_results$external_gene_name == 'CD2AP')
raw_CD2AP <- raw_CD2AP %>% dplyr::select("ensembl_gene_id", "adj.P.Val")
astro_promoters_CD2AP_df$ensembl_gene_id <- raw_CD2AP$ensembl_gene_id
astro_promoters_CD2AP_df <- astro_promoters_CD2AP_df  %>% as.data.frame() %>% left_join(raw_CD2AP, by=  "ensembl_gene_id")
astro_promoters_CD2AP_df$niaaaascore %<>% as.numeric

## oligo opcs
oligo_opc_promoters_CD2AP <- oligo_opc_promoters[CD2AP,]
oligo_opc_promoters_CD2AP <- rnb.mval2beta(oligo_opc_promoters_CD2AP)
oligo_opc_promoters_CD2AP_df <- tibble(sample_name = names(oligo_opc_promoters_CD2AP), Methylation = oligo_opc_promoters_CD2AP, celltype = "CD2AP DG Olig/OPCs",sig_color= "not_sig") %>%
  left_join(targets)
raw_results <- read.csv(paste0("~/90plus/dm_results/continuous/raw/raw_niaaaascore_oligo_opc_promoters_DEN.csv"))
raw_results <- merge(x = raw_results,y = genes_with_description, by = "ensembl_gene_id",all.x=TRUE)
raw_CD2AP <- raw_results %>% filter(raw_results$external_gene_name == 'CD2AP')
raw_CD2AP <- raw_CD2AP %>% dplyr::select("ensembl_gene_id", "adj.P.Val")
oligo_opc_promoters_CD2AP_df$ensembl_gene_id <- raw_CD2AP$ensembl_gene_id
oligo_opc_promoters_CD2AP_df <- oligo_opc_promoters_CD2AP_df  %>% as.data.frame() %>% left_join(raw_CD2AP, by=  "ensembl_gene_id")
oligo_opc_promoters_CD2AP_df$niaaaascore %<>% as.numeric

## endothelial cells
endo_promoters_CD2AP <- endo_promoters[CD2AP,]
endo_promoters_CD2AP <- rnb.mval2beta(endo_promoters_CD2AP)
endo_promoters_CD2AP_df <- tibble(sample_name = names(endo_promoters_CD2AP), Methylation = endo_promoters_CD2AP, celltype = "CD2AP DG Endothelial cells", sig_color= "not_sig") %>%
  left_join(targets)
endo_promoters_CD2AP_df$niaaaascore %<>% as.numeric
raw_results <- read.csv(paste0("~/90plus/dm_results/continuous/raw/raw_niaaaascore_endo_promoters_DEN.csv"))
raw_results <- merge(x = raw_results,y = genes_with_description, by = "ensembl_gene_id",all.x=TRUE)
raw_CD2AP <- raw_results %>% filter(raw_results$external_gene_name == 'CD2AP')
raw_CD2AP <- raw_CD2AP %>% dplyr::select("ensembl_gene_id", "adj.P.Val")
endo_promoters_CD2AP_df$ensembl_gene_id <- raw_CD2AP$ensembl_gene_id
endo_promoters_CD2AP_df <- endo_promoters_CD2AP_df  %>% as.data.frame() %>% left_join(raw_CD2AP, by=  "ensembl_gene_id")
endo_promoters_CD2AP_df$niaaaascore %<>% as.numeric

## microglia
microglia_promoters_CD2AP <- microglia_promoters[CD2AP,]
microglia_promoters_CD2AP <- rnb.mval2beta(microglia_promoters_CD2AP)
microglia_promoters_CD2AP_df <- tibble(sample_name = names(microglia_promoters_CD2AP), Methylation = microglia_promoters_CD2AP, celltype = "CD2AP DG Microglia",sig_color= "not_sig") %>%
  left_join(targets)
microglia_promoters_CD2AP_df$niaaaascore %<>% as.numeric
raw_results <- read.csv(paste0("~/90plus/dm_results/continuous/raw/raw_niaaaascore_microglia_promoters_DEN.csv"))
raw_results <- merge(x = raw_results,y = genes_with_description, by = "ensembl_gene_id",all.x=TRUE)
raw_CD2AP <- raw_results %>% filter(raw_results$external_gene_name == 'CD2AP')
raw_CD2AP <- raw_CD2AP %>% dplyr::select("ensembl_gene_id", "adj.P.Val")
microglia_promoters_CD2AP_df$ensembl_gene_id <- raw_CD2AP$ensembl_gene_id
microglia_promoters_CD2AP_df <- microglia_promoters_CD2AP_df  %>% as.data.frame() %>% left_join(raw_CD2AP, by=  "ensembl_gene_id")
microglia_promoters_CD2AP_df$niaaaascore %<>% as.numeric

## combine all dataframes
joined_df <- rbind(bulk_promoters_CD2AP_df,
                   neuron_promoters_CD2AP_df, 
                   astro_promoters_CD2AP_df,
                   endo_promoters_CD2AP_df,
                   microglia_promoters_CD2AP_df,
                   oligo_opc_promoters_CD2AP_df)
joined_df$celltype <- factor(joined_df$celltype, levels = c("CD2AP DG Bulk",
                                                            "CD2AP DG Neuron **",
                                                            "CD2AP DG Astrocytes",
                                                            "CD2AP DG Endothelial cells",
                                                            "CD2AP DG Microglia",
                                                            "CD2AP DG Olig/OPCs"))

## make a joint plot with factewrap
p <- ggplot(joined_df, aes(niaaaascore, Methylation)) +
  stat_smooth(method = "lm",colour="grey") +
  geom_point(alpha = 0.5) + 
  geom_text(
    size    = 3,
    color = "black",
    data    = joined_df,
    mapping = aes(
      x = -Inf,
      y = -Inf,
      label = paste("FDR p = ", round(adj.P.Val, digits=3)),
      hjust   = -0.1,
      vjust   = -0.8
    )
  ) +
  facet_wrap(~joined_df$celltype)+
  scale_x_continuous(breaks = c(1, 2, 3)) +
  theme_classic() +
  theme(strip.background = element_rect(colour="black",
                                        fill="white"), 
        text = element_text(size = 12))+
  theme(panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) +
  xlab("NIA-AA A score") +
  ylab("Methylation beta value")+
  scale_y_continuous(breaks=seq(0,1, 0.5), limits=c(0,1))

ggsave(
  paste0(baseDir,"/plots/final_figures/OR_fig_CD2AP_DG_all_cell_types.png"),
  p,
  width = 174,
  height = 80,
  dpi = 1200, 
  units = 'mm'
)

