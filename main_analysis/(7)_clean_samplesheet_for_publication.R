# code to clean samplesheet
setwd("~/90plus")
deconvo <- readRDS("targets_incl_episcore.rds")

#  rename
deconvo <- deconvo %>%
  mutate(brain_region = case_when(grepl("^ERC", brain_region) ~ "EC",
                                  grepl("^DEN", brain_region) ~ "DG",
                                  grepl("^LOC", brain_region) ~ "LC",
                                  grepl("^CBL", brain_region) ~ "CBM",
                                  grepl("^MFG", brain_region) ~ "MFG",
                                  grepl("^SN", brain_region) ~ "SN",
                                  grepl("^CG", brain_region) ~ "CG",
                                  grepl("^CA1", brain_region) ~ "CA1"
  ))
deconvo$Olig_opc <- deconvo$Oligo+deconvo$OPC

deconvo <- deconvo %>% dplyr::select(c("sample_name"    ,                   
                     "sample_source","brain_region" ,                     
                     "bisulfite_batch" ,                  
                     "batch", "gender", "education","genotype",
                     "conf_diag", "mmse_total", "cvlt_delay", "flanitot", 
                     "niaaaascore", "niaaabscore", "niaaacscore","adseverityscore",
                     "braak_for_lb_orig", "tdp43", "mvl_score", "Neuron", "Olig_opc", "Astro", "Endo", "Microglia"))

deconvo <- deconvo %>%
  mutate(case_when(grepl("^Neuron") ~ "proportion_neurons",
                                  grepl("^Olig_opc") ~ "proportion_olig_opc",
                                  grepl("^Astro") ~ "proportion_astrocytes",
                                  grepl("^Endo") ~ "proportion_endothelial",
                                  grepl("^Microglia") ~ "proportion_microglia"
  ))

deconvo <- dplyr::rename(deconvo, "proportion_neurons" = "Neuron")
deconvo <- dplyr::rename(deconvo, "proportion_olig_opc" = "Olig_opc") 
deconvo <- dplyr::rename(deconvo, "proportion_astrocytes" = "Astro")
deconvo <- dplyr::rename(deconvo, "proportion_endothelial" = "Endo")
deconvo <- dplyr::rename(deconvo, "proportion_microglia" = "Microglia") 

write.csv(deconvo, file = "~/90plus/Sup_File_cleaned_samplesheet.csv")
