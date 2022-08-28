## script to plot EpiSCORE results for Supp Figure 2
library(ggpubr)
library(ggplot2)
library(tidyr)
library(grid)

## set wd
setwd("~/90plus/")
baseDir <- paste0("~/90plus/")

## load function to plot legend inside empty facet when using facet_wrap for plotting
source("./functions/shift_legend.R")

## read in df with episcore cell type proportions
targets <- readRDS("./targets_incl_episcore.rds")

## define colors
col_brain <-
  c("#DDCC77",
    "#999933",
    "#CC6677",
    "#882255",
    "#AA4499",
    "#88CCEE",
    "#332288",
    "#117733"
  )
## combine olig and opc
targets$Olig_OPC <- targets$Oligo + targets$OPC
## rename brain regions for plotting
targets <- targets %>%
  mutate(brain_region = case_when(grepl("^ERC", brain_region) ~ "Entorhinal cortex",
                                  grepl("^DEN", brain_region) ~ "Dentate gyrus",
                                  grepl("^LOC", brain_region) ~ "Locus coeruleus",
                                  grepl("^CBL", brain_region) ~ "Cerebellar cortex",
                                  grepl("^MFG", brain_region) ~ "Middle frontal gyrus",
                                  grepl("^SN", brain_region) ~ "Substantia nigra",
                                  grepl("^CG", brain_region) ~ "Cingulate gyrus",
                                  grepl("^CA1", brain_region) ~ "Hippocampus CA1"
  ))
targets$brain_region <- factor(targets$brain_region, levels = c("Middle frontal gyrus",
                                                                "Cingulate gyrus",
                                                                "Dentate gyrus",
                                                                "Hippocampus CA1",
                                                                "Entorhinal cortex",
                                                                "Locus coeruleus",
                                                                "Substantia nigra",
                                                                "Cerebellar cortex"))

## plot EpiSCORE results with Olig/OPCs combined
targets <- targets[,c("brain_region", "Neuron", "Astro", "Microglia", "Endo","Olig_OPC" )]
colnames(targets) <- c("brain_region", "Neurons", "Astrocytes", "Microglia", "Endothelial cells","Olig/OPCs" )
targets_long <- pivot_longer(targets, cols=c("Neurons", "Astrocytes", "Microglia", "Endothelial cells","Olig/OPCs"), names_to = "cell_type", values_to = "proportion")
targets_long$cell_type <- factor(targets_long$cell_type, levels = c("Neurons", "Astrocytes", "Microglia", "Endothelial cells","Olig/OPCs"))

p <- ggplot(targets_long, aes(x=brain_region, y=proportion, color=brain_region))+
  geom_boxplot(width=0.6, outlier.shape = NA, alpha=0.05)+
  geom_jitter(width=0.1, size=1)+
  facet_rep_wrap(~cell_type, nrow=2,repeat.tick.labels = 'left')+
  ylab("Cell type proportion") +
  scale_color_manual(values=col_brain)+
  theme_classic()+
  theme(
    text = element_text(size = 12), 
    axis.title.y=element_text(size = 12),
    axis.text.y=element_text(size = 12),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    strip.text = element_text(size = 12))+
  xlab("")+  
  guides(color=guide_legend(title="Brain region"))

png(paste0("./plots/final_figures/Supp_Fig_2_episcore_proportions.png"),width=10,height=8, units ='in', res=800) 
grid.draw(shift_legend(p))
dev.off()

