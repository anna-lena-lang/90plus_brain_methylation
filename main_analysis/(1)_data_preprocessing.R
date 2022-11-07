## this is code to 
## - load all data from .idat files
## - run quality control on samples + plot density and snp-heatmaps 
## - remove low quality probes and normalize data for each brain region separately

library(dplyr)
library(IlluminaHumanMethylationEPICanno.ilm10b3.hg19) # BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b3.hg19")
library(minfi)
library(limma)
library(missMethyl)
library(Gviz)
library(DMRcate)
library(stringr)
library(grid)
library(wateRmelon)
library(methylumi)
library(ChAMP)
library(MethylAid)
library(ewastools)# devtools::install_github("hhhh5/ewastools")
library(reshape)
library(data.table)
library(rlang) 
library(minfiData) # needed for maxprobes
library(maxprobes) # remotes::install_github("markgene/maxprobes")
library(RColorBrewer)
library(Rtsne)
library(randomcoloR) 
library(plotly)
library(gplots)
library(ggplot2)
library(viridis)
## set working directory
setwd("~/90plus/idats/")
baseDir <- getwd()

## -----functions----------------------------------------------------------------
## clean targets sheet
clean_samplesheet <- function(targets){
  ## remove several batches that failed due to scanner problems
  targets <<- targets %>% 
    filter(!batch %in% c("2018_12_05","2018_12_11","2019_01_17", "2019_02_11", "2019_03_01")) %>%
    ## remove cases participant_48, participant_49, participant_50,participant_51,participant_52,participant_53,participant_54,participant_55
    filter(!sample_source %in% c("participant_48", "participant_49", "participant_50","participant_51","participant_52",
                                 "participant_53","participant_54","participant_55")) %>%
    ## participant_48, participant_53 not all FFPE blocks were available for processing - full cases excluded
    ## participant_49_CBM sample with high delta ct in qPCR, full case excluded
    ## participant_50 samples with high delta ct in qPCR, full case excluded
    ## participant_51 case with potential Parkinsons Disease and infarction, full case excluded
    ## participant_52 case with Parkinsons Disease, full case excluded
    ## participant_54 case with Glioma, full case excluded
    ## participant_55 case with massive hemorrhage, full case excluded
    ## only keep samples from 8 brain regions, remove all samples  from other projects
    filter(brain_region %in% c("MFG", "CG", "CA1", "DG", "EC", "LC", "SN", "CBM"))
  
  ## change NAs to 0 in biological_replicate
  targets$biological_replicate[is.na(targets$biological_replicate)] <- 0 
  targets[targets == "#NA"] <- NA
  targets[targets == "#N/A"] <- NA
  
  ## change empty rows to NA
  targets[targets == ""] <-  NA
  
  ## add braak for lb_orig 
  targets$braak_for_lb_orig <- targets$braak_for_lb
  
  ## change braak for lb scoring groups
  targets$braak_for_lb <- sub("0", "1", targets$braak_for_lb)
  targets$braak_for_lb <- sub("2", "1", targets$braak_for_lb)
  targets$braak_for_lb <- sub("3", "2", targets$braak_for_lb)
  targets$braak_for_lb <- sub("4", "2", targets$braak_for_lb)
  targets$braak_for_lb <- sub("5", "3", targets$braak_for_lb)
  targets$braak_for_lb <- sub("6", "3", targets$braak_for_lb)
  
  ## combine 0 and 1 within nia-score groups and ad-severity-score groups
  targets$niaaaascore <- sub("0", "1", targets$niaaaascore)
  targets$niaaacscore <- sub("0", "1", targets$niaaacscore)
  targets$adseverityscore <- sub("0", "1", targets$adseverityscore)
  
  ## change adseverityscore from numeric to character
  targets$adseverityscore <- sub("1", "low", targets$adseverityscore)
  targets$adseverityscore <- sub("2", "intermediate", targets$adseverityscore)
  targets$adseverityscore <- sub("3", "high", targets$adseverityscore)
  
  ## change gender 0/1 to Male/Female
  targets$gender <- factor(targets$gender, levels=c(0,1),labels=c("Male", "Female"))
  
  ## change character to factor for brain region
  targets$brain_region <- as.factor(targets$brain_region)
  targets$brain_region <- factor(targets$brain_region, levels =  c("MFG", "CG", "CA1", "DG", "EC", "LC", "SN", "CBM"))
  ## change character to factor forclinical diganosis
  targets$conf_diag <- as.factor(targets$conf_diag)
  targets$conf_diag <- factor(targets$conf_diag, levels=c(0, 1, 2),labels=c("normal", "CIND", "demented"))
  
  ## add new scores to targets: lbd yes/no, mvl yes/no, tdp43 yes/no
  #lb
  targets$lb_cat <- targets$braak_for_lb_orig
  targets$lb_cat <- sub("0", "no", targets$lb_cat)
  targets$lb_cat <- sub("1", "yes", targets$lb_cat)
  targets$lb_cat <- sub("2", "yes", targets$lb_cat)
  targets$lb_cat <- sub("3", "yes", targets$lb_cat)
  targets$lb_cat <- sub("4", "yes", targets$lb_cat)
  targets$lb_cat <- sub("5", "yes", targets$lb_cat)
  targets$lb_cat <- sub("6", "yes", targets$lb_cat)
  
  ## tdp43
  targets$tdp43_cat <- targets$tdp43
  targets$tdp43_cat <- sub("0", "no", targets$tdp43_cat)
  targets$tdp43_cat <- sub("1", "yes", targets$tdp43_cat)
  targets$tdp43_cat <- sub("2", "yes", targets$tdp43_cat)
  targets$tdp43_cat <- sub("3", "yes", targets$tdp43_cat)
  
  ## mvl
  targets$mvl_cat <- targets$mvl_score
  targets$mvl_cat <- sub("0", "no", targets$mvl_cat)
  targets$mvl_cat <- sub("1", "yes", targets$mvl_cat)
  targets$mvl_cat <- sub("2", "yes", targets$mvl_cat)
  targets$mvl_cat <- sub("3", "yes", targets$mvl_cat)
  

  print("saved cleaned targets sheet to targets_cleaned.rds, also assigned it to *targets* df in global environment")
  saveRDS(targets, file = "targets_cleaned.rds")
}

## bmiq normalization
bmiq_norm <- function(targets_90plus_ewasfilt_incl_R_850k, beta){
  ## BMIQ normalization 
  print(paste0("starting bmiq normalization for brain region " , brainregion))
  beta_90plus_ewasfilt_bmiq_incl_R_850k <<- champ.norm(beta = beta, method = "BMIQ", plotBMIQ = FALSE, arraytype = "EPIC", cores = 20)
  saveRDS(beta_90plus_ewasfilt_bmiq_incl_R_850k, file = paste0(baseDir, "/data/beta_90plus_ewasfilt_bmiq_incl_R_850k_", brainregion, ".rds"))
  targets_90plus_ewasfilt_bmiq_incl_R_850k <<- targets_90plus_ewasfilt_incl_R_850k
  saveRDS(targets_90plus_ewasfilt_bmiq_incl_R_850k, file = paste0(baseDir, "/data/targets_90plus_ewasfilt_bmiq_incl_R_850k_", brainregion, ".rds"))
  print("bmiq done and results saved")
}

## remove biological and technical replicates 
remove_reps <- function(targets, beta){
  ## remove technical and biological replicates
  print(paste0("removing tR and bR for brain region ",  brainregion))
  print(paste0("number of bR: " , sum(targets$biological_replicate == 1)))
  print(paste0("number of tR: " , sum(targets$technical_replicate == 1)))
  targets_90plus_ewasfilt_bmiq_orig_850k <<- targets[targets$original == 1,]
  beta_90plus_ewasfilt_bmiq_orig_850k <<- beta[ ,colnames(beta)%in% targets_90plus_ewasfilt_bmiq_orig_850k$sample_name]
  saveRDS(beta_90plus_ewasfilt_bmiq_orig_850k, file =  paste0(baseDir, "/data/beta_90plus_ewasfilt_bmiq_orig_850k_", brainregion, ".rds"))
  saveRDS(targets_90plus_ewasfilt_bmiq_orig_850k, file = paste0(baseDir, "/data/targets_90plus_ewasfilt_bmiq_orig_850k_", brainregion, ".rds"))
  print(paste0("tR and bR removed, new datasets saved. finished analysis for brain region ",  brainregion, "total count of samples: ", nrow(targets_90plus_ewasfilt_bmiq_orig_850k)))
}

## density plot for each batch
densityplots_bybatch <- function (rgset) {
  require(data.table)
  pd <- pData(rgset) %>%  as.data.frame() %>%
    dplyr::select(sample_name, batch, brain_region)
  pdf(paste0(baseDir, "/plots/density_plot/densityplots_by_batch.pdf"))
  for (i in unique(pd$batch)){
    pd2 <- pd[pd$batch==paste0(i),]
    keep <- colnames(rgset) %in% pd2$sample_name
    rgset2 <- rgset[,keep]
    ad <- getBeta(rgset2) %>% melt()
    colnames(ad) <- c('ProbeID', 'sample_name', 'beta')
    ad$sample_name <- as.character(ad$sample_name)
    ad <- as.data.table(ad)
    plot_data <- pd2 %>% inner_join(ad, by='sample_name')
    print(
      ggplot(plot_data , aes(x=beta, colour=as.factor(sample_name))) + 
        geom_density()+
        geom_hline(yintercept=0, colour="white", size=1)+#overrides the vertical line with white line
        stat_density(aes(x=beta, colour=as.factor(sample_name)), geom="line", position="identity")+
        theme(legend.position="right") +
        labs(colour='Sample ID')+
        ggtitle("density plot batch : ",i)+
        xlab(label = 'beta value')
    )
  }
  dev.off()
}

## interactive density plot using plotly, colored by batch
densityplots_plotly <- function(rgset){
  pd <- pData(rgset) %>%  as.data.frame() %>%
    dplyr::select(sample_name, batch, brain_region)
  ad <- getBeta(rgset) %>% melt()
  ad <- na.omit(ad)#remove NAs 
  colnames(ad) <- c('ProbeID', 'sample_name', 'beta')
  ad$sample_name <- as.character(ad$sample_name)
  ad <- as.data.table(ad)
  plot_data <- pd %>% inner_join(ad, by='sample_name')
  set.seed(1234)
  col <- distinctColorPalette(k = 27, altCol = FALSE, runTsne = FALSE)
  # plot_data2 <- plot_data[order(plot_data$batch),]
  ## plotly
  set.seed(1234)
  p <- ggplot(plot_data , aes(x=beta, colour=batch)) + 
    geom_density(aes(group = sample_name, color = batch),position="identity")+
    geom_hline(yintercept=0, colour="white", size=1)+#overrides the vertical line with white line
    theme(legend.position = "right") +
    labs(fill = 'batch') +
    scale_color_manual(values=col)+
    ggtitle("density Plot") +
    xlab(label = 'beta value')
  fig <- ggplotly(p)
  htmlwidgets::saveWidget(as_widget(fig), paste(baseDir,"/plots/density_plot/plotly_density_plot.html", sep=""))
}

## snp heatmap
snp_heatmap_bysource <- function(rgset) {
  pd <- pData(rgset) %>%  as.data.frame() %>% dplyr::select(sample_source, sample_name, batch, brain_region)
  my_col <- "heat.colors"## define color palette for heatmap
  pdf(paste0(baseDir, "/plots/SNP_heatmap/SNP_heatmaps_allsamples.pdf"), onefile=T, paper='A4r', width=9,height=7)
  for (i in unique(pd$sample_source)){
    pd_sub <- pd[pd$sample_source == i,]
    rgset_sub <- rgset[,colnames(rgset) %in% pd_sub$sample_name]
    snps <- minfi::getSnpBeta(rgset_sub) 
    snps <- snps[,order(colnames(snps))]
    print(
      heatmap.2(snps,scale = "none", trace = "none", margins = c(8,8),
                labRow = rownames(snps), labCol = colnames(snps), col = my_col,
                cexCol = 0.5, cexRow = 0.5, main = paste0("SNP heatmap, case ID:",i))
    )
  }
  dev.off()
}

## QC for samples will filter samples that: don't show beta distribution, with bisulfite conversion efficiency <80%, fail probe control metrics using ewastools,  are outliers in snp-heatmap (compared within each individual)
qualitycontrol_samples <- function(targets , rgset){
  ## make density plots and manually select the ones that don't show beta distribution (peak at 0 and 1)
  # densityplots_bybatch(rgset)
  # densityplots_plotly(rgset)
  
  ## make list containing all samples removed due to low quality
  low_q_samples <- list()
  ## define and remove samples that did not show beta value distribution in density plots
  print(paste0("sample number : ", nrow(targets), " Now removing samples with bad density plots")) ## printing numbers is just for my protocol + downstream overview.
  bad_density <- c("participant_7_Lx","participant_10_N_tR_04_samerun_otherchip", "participant_45_J2", targets$sample_name[targets$batch == "Batch_9"])
  rgset <- rgset[, !(colnames(rgset) %in% bad_density)]
  t_old <- targets
  targets <- targets[targets$sample_name %in% colnames(rgset), ]
  low_q_samples[["bad_density"]] <- nrow(t_old)- nrow(targets)
  ## bisulfite conversion efficiency
  print(paste0("sample number : ", nrow(targets), "Now removing samples with bisulfite conversion efficiency < 80%"))
  bsc <- bscon(rgset) ## checks bisulfite conversion efficiency
  bsc <- as.data.frame(bsc)
  bad_bsc <- rownames(bsc)[bsc$bsc < 80] ## remove samples with bc efficiency <80% = "participant_28_F"
  rgset <- rgset[ ,!(colnames(rgset) %in% bad_bsc)]
  t_old <- targets
  targets <- targets[targets$sample_name %in% colnames(rgset), ]
  print(paste0("sample number : ", nrow(targets)))
  low_q_samples[["bad_bsc"]] <- nrow(t_old)- nrow(targets)
  
  ## snp heatmap: manually removing snp outliers as seen in heatmaps from snp_heatmap_bysourc
  snp_heatmap_bysource(rgset) ## prints all heatmaps by sample source (=individual donor/patient) in one pdf
  print("removing samples that don't match snp heatmap pattern of samples from same donor")
  bad_snp <- c("participant_45_J2_redoBatch4fail", "participant_46_J2", "participant_6_TERx", "participant_8_DG", "participant_10_CA1", "participant_56_CA1_tR01")
  rgset <- rgset[, !(colnames(rgset) %in% bad_snp)]
  t_old <- targets
  targets <- targets[targets$sample_name %in% colnames(rgset), ]
  print(paste0("sample number : ", nrow(targets)))
  low_q_samples[["bad_snp"]] <- nrow(t_old)- nrow(targets)
  
  ## use contol_metrics function from ewastools to remove failed samples
  print("use contol_metrics function from ewastools to remove failed samples")
  meth <- read_idats(targets$Basename, quiet = FALSE) ## read in idats for the control_metrics function
  ctrls <- control_metrics(meth)
  targets$failed <- sample_failure(ctrls) ## check for failed samples
  ewastools_metrics_fail <- targets %>% filter(targets$failed == TRUE)
  rgset_90plus_raw_qcsfilt_incl_R_850k <<- rgset[, !(colnames(rgset) %in% targets[targets$failed == "TRUE",])] ## remove failed samples
  targets_90plus_raw_qcsfilt_incl_R_850k <<- targets[targets$sample_name %in% colnames(rgset), ]
  low_q_samples[["ewastools_metrics_fail"]] <- nrow(targets)- nrow(targets_90plus_raw_qcsfilt_incl_R_850k)
  print(paste0("sample number targets_90plus_raw_qcsfilt_incl_R_850k: ", nrow(targets_90plus_raw_qcsfilt_incl_R_850k)))
  
  ## save files
  print("sample filtering done. saving new datasets: rgset_90plus_raw_qcsfilt_incl_R_850k.rds, targets_90plus_raw_qcsfilt_incl_R_850k.rds")
  saveRDS(rgset_90plus_raw_qcsfilt_incl_R_850k, file = paste0(baseDir, "/data/rgset_90plus_raw_qcsfilt_incl_R_850k.rds"))
  saveRDS(targets_90plus_raw_qcsfilt_incl_R_850k, file = paste0(baseDir, "/data/targets_90plus_raw_qcsfilt_incl_R_850k.rds"))
  saveRDS(low_q_samples, file = paste0(baseDir, "/data/low_q_samples.rds"))
  }

## qc for probes, normalize using noob and bmiq and finally remove biological/technical replicates 
qc_probes_norm_remove_reps <- function(targets, rgset) {
  print(paste0("starting quality control of probes for brain region " , brainregion))
  print(paste0("probe number: ", nrow(rgset)))
  
  ## subset dataset by brain region
  targets <- targets[targets$brain_region == brainregion, ]
  rgset <- rgset[ ,colnames(rgset) %in% targets$sample_name]
  
  ## get detection p -values 
  detP_ewas <- ewastools::detectionP.minfi(rgset) ## calculating detP value using ewastools
  detP_ewas <- detP_ewas[, order(colnames(detP_ewas))]
  print(paste0("saving dataset with detection pvals, detP_ewas_90plus_raw_", brainregion, ".rds"))
  saveRDS(detP_ewas, file=paste0(baseDir, "/data/detP_ewas_90plus_raw_", brainregion, ".rds"))
  
  ## get beadcount for filtering later on
  beadcount <- beadcount(rgset)
  
  ## use noob (minfi::preprocessNoob) for dye bias Correction
  print(paste0("running noob dye bias correction for brain region ", brainregion))
  mset_90plus_noob_qcsfilt_incl_R_850k <- preprocessNoob(rgset, dyeCorr = TRUE, dyeMethod = "single")
  RSet_90plus_noob_qcsfilt_incl_R_850k <- ratioConvert(mset_90plus_noob_qcsfilt_incl_R_850k, what = "both", keepCN = TRUE)
  GRset_90plus_noob_qcsfilt_incl_R_850k <- mapToGenome(RSet_90plus_noob_qcsfilt_incl_R_850k)
  print(paste0("probe number: ", nrow(GRset_90plus_noob_qcsfilt_incl_R_850k)))
  
  ## 1) data includes males and females --> remove probes on the sex chromosomes
  print("removing probes on sex-chromosomes")
  keep <- !(featureNames(GRset_90plus_noob_qcsfilt_incl_R_850k) %in% annEpic$Name[annEpic$chr %in% c("chrX","chrY")])
  GRset_90plus_noob_qcsfilt_incl_R_850k_filt <- GRset_90plus_noob_qcsfilt_incl_R_850k[keep,]
  print(paste0("probe number: ", nrow(GRset_90plus_noob_qcsfilt_incl_R_850k_filt)))
  
  ## 2) remove probes with SNPs at CpG site
  print("removing probes associated with SNPs of minor allele frequency (MAF) 0.05")
  GRset_90plus_noob_qcsfilt_incl_R_850k_filt2genome <- mapToGenome(GRset_90plus_noob_qcsfilt_incl_R_850k_filt)
  GRset_90plus_noob_qcsfilt_incl_R_850k_filt3 <- dropLociWithSnps(GRset_90plus_noob_qcsfilt_incl_R_850k_filt2genome)
  print(paste0("probe number: ", nrow(GRset_90plus_noob_qcsfilt_incl_R_850k_filt3)))
  ## set up list that will save number of probes being remoed in each qc step
  low_qc_probes <- list()
  low_qc_probes[["snps_MAF"]] <-  nrow(GRset_90plus_noob_qcsfilt_incl_R_850k_filt2genome) - nrow(GRset_90plus_noob_qcsfilt_incl_R_850k_filt3)
  ## 3) finally exclude cross reactive probes (based on maxprobes github)
  print("removing cross reactive probes based on maxprobes::xreactive_probes")
  xreact <- xreactive_probes(array_type = "EPIC")
  keep2 <- !(featureNames(GRset_90plus_noob_qcsfilt_incl_R_850k_filt3) %in% xreact)
  GRset_90plus_noob_qcsfilt_incl_R_850k_filt4 <- GRset_90plus_noob_qcsfilt_incl_R_850k_filt3[keep2,]
  beta_90plus_noobfilt_incl_R_850k <- getBeta(GRset_90plus_noob_qcsfilt_incl_R_850k_filt4)
  print(paste0("probe number: ", nrow(beta_90plus_noobfilt_incl_R_850k)))
  ## include in list saving number of probes being removed in each qc step
  low_qc_probes[["x_reactive"]] <-  nrow(GRset_90plus_noob_qcsfilt_incl_R_850k_filt3) - nrow(GRset_90plus_noob_qcsfilt_incl_R_850k_filt4)
  
  ## 4) filter samples by detection P value from ewastools, using champ.filter 
  ## subset all datasets to match
  
  targets <- targets[order(targets$sample_name),]
  detP_90plus_ewasfilt_850k <- detP_ewas[rownames(detP_ewas)%in% rownames(beta_90plus_noobfilt_incl_R_850k),colnames(detP_ewas) %in% colnames(beta_90plus_noobfilt_incl_R_850k)]
  detP_90plus_ewasfilt_850k <- detP_90plus_ewasfilt_850k[order(rownames(detP_90plus_ewasfilt_850k)),]
  beta_90plus_noobfilt_incl_R_850k <- beta_90plus_noobfilt_incl_R_850k[order(rownames(beta_90plus_noobfilt_incl_R_850k)),]
  if(all(rownames(detP_90plus_ewasfilt_850k) == rownames(beta_90plus_noobfilt_incl_R_850k)) == FALSE){
    warning('rownames of dataframes dont match')
  }
  
  ## subset beadcount
  beadcount_90plus_850k <- beadcount[rownames(beta_90plus_noobfilt_incl_R_850k),]
  print("removing: probes with detP > 0.01, probes with <3 beads in at least 5% of samples, samples where > 10% of probes have detection p value > 0.01,
        probes with NA using champ.filter")
  
  ## filter based on beadcount and detP from ewastools using champ.filter
  champ <- champ.filter(beta = beta_90plus_noobfilt_incl_R_850k, 
                        pd = targets,
                        detP = detP_90plus_ewasfilt_850k,
                        beadcount = beadcount_90plus_850k,
                        autoimpute = FALSE,
                        ProbeCutoff = 0,
                        SampleCutoff = 0.1,
                        detPcut = 0.01,
                        filterDetP = TRUE,
                        beadCutoff = 0.05,
                        filterBeads = TRUE,
                        filterNoCG = FALSE,
                        filterSNPs = FALSE,
                        population = NULL,
                        filterMultiHit = FALSE,
                        filterXY = FALSE,
                        fixOutlier = FALSE,
                        arraytype = "EPIC")
  beta <- champ$beta
  targets_90plus_ewasfilt_incl_R_850k <<- champ$pd
  ## include in list saving number of probes being removed in each qc step
  name <- paste0("champ_ewastools_etc_", brainregion)
  low_qc_probes[[name]] <-  nrow(beta_90plus_noobfilt_incl_R_850k) - nrow(beta)
  
  ## get all samples that were being removed by champ filter (detP>0.01 in >10% of probes)
  champ_detP_low_q <- targets %>% filter(!(targets$sample_name %in% targets_90plus_ewasfilt_incl_R_850k$sample_name))
  ## add to vector containing all low quality samples 
  low_q_samples_file <- paste0(baseDir, "/data/low_q_samples.rds")
  low_q_samples <- readRDS(low_q_samples_file)
  name <- paste0("champ_detP_low_q_", brainregion)
  low_q_samples[[name]] <- champ_detP_low_q$sample_name

  print(paste0("probe number: ", nrow(beta)))
  
  print("removing probes with low variability (SD < 0.01) using rnb.execute.variability.removal")

  pheno_temp <- data.frame(Sample = colnames(beta))
  rnb_set <- RnBeadSet(pheno = pheno_temp, probes = rownames(beta), betas = beta, platform = "EPIC")
  rnb_set_filtered <- rnb.execute.variability.removal(rnb_set, 0.01)$dataset
  ## get mvals from rnb_set
  mvals <- mval(rnb_set_filtered, type = "sites" , row.names = TRUE)
  
  ## get beta
  beta_90plus_ewasfilt_incl_R_850k <<- lumi::m2beta(mvals)
  print(paste0("probe number: ", nrow(beta_90plus_ewasfilt_incl_R_850k)))
  
  name <- paste0("low_variability_", brainregion)
  low_qc_probes[[name]] <- nrow(beta) - nrow(beta_90plus_ewasfilt_incl_R_850k)
  
  ## save results
  saveRDS(beta_90plus_ewasfilt_incl_R_850k, file = paste0(baseDir, "/data/beta_90plus_ewasfilt_incl_R_850k_", brainregion, ".rds"))
  saveRDS(targets_90plus_ewasfilt_incl_R_850k, file = paste0(baseDir, "/data/targets_90plus_ewasfilt_incl_R_850k_", brainregion, ".rds"))
  print(paste0("Finished probe filtering, datasets incl replicates saved for ", brainregion))
  ## normalize
  bmiq_norm(targets_90plus_ewasfilt_incl_R_850k, beta_90plus_ewasfilt_incl_R_850k)
  ## remove replicates
  remove_reps(targets_90plus_ewasfilt_bmiq_incl_R_850k, beta_90plus_ewasfilt_bmiq_incl_R_850k)
  name <- paste0("replicates_", brainregion)
  ## get all samples that were being removed bc they are replicates
  reps <- targets_90plus_ewasfilt_bmiq_incl_R_850k %>% filter(!(targets_90plus_ewasfilt_bmiq_incl_R_850k$sample_name %in% targets_90plus_ewasfilt_bmiq_orig_850k$sample_name))
  low_q_samples[[name]] <- reps$sample_name
  ## save lists with numbers of qc control measures for flowchart on quality control
  saveRDS(low_q_samples, file = paste0(baseDir, "/data/low_q_samples.rds"))
  saveRDS(low_qc_probes, file = paste0(baseDir, "/data/low_qc_probes.rds"))
  } 

## combine all QC steps in one function
preprocess_bybrainregion <- function(targets, rgset){
  qualitycontrol_samples(targets, rgset) ## first running quality control on samples (checking bisuflite conversion, snp heatmaps and ewastools quality control to remove any outlier samples, this is done on full dataset
  for (brainregion in brainregions){ ## now for each brain region separately check quality of probes in each sample and remove probes with low quality (hight detection p-value detP > 0.01, beadcount, remove probes on sex chromosomes and remove cross-reactive probes. Finally remove do dye-bias-correction using noob and normalization using bmiq and remove all technical and biological replicates)
    brainregion <<- brainregion
    qc_probes_norm_remove_reps(targets_90plus_raw_qcsfilt_incl_R_850k, rgset_90plus_raw_qcsfilt_incl_R_850k)
  }
}


## ----- run all preprocessing steps---------------------------------------------
## get annotation data for epic array - needed for mapping probes to genome
annEpic <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b3.hg19)

## read in sample sheet containing phenotype and array information
targets <- read.metharray.sheet(baseDir, pattern = "mastersheet_all_samples.csv") 

## clean samplesheet
clean_samplesheet(targets) ## function below, changes nia-scores/braak_for_lb, all NAs in same format, remove some batches, saves sample sheet to targets_clean.rds

## read in raw idats
rgset_90plus_raw <- read.metharray.exp(targets = targets, extended = TRUE)
# saveRDS(rgset_90plus_raw, file = "rgset_90plus_raw.rds")

## assign column names
colnames(rgset_90plus_raw) <- rgset_90plus_raw@colData$sample_name 

## order dfs
rgset_90plus_raw <- rgset_90plus_raw[ ,order(colnames(rgset_90plus_raw))]
targets_90plus_raw <- targets[order(targets$sample_name), ]
all(colnames(rgset_90plus_raw) == targets_90plus_raw$sample_name) ## make sure that align

## set working directory to main folder
setwd("~/90plus/")
baseDir <- getwd()

## assign new samplenames to rgset
sampleNames(rgset_90plus_raw) <- targets_90plus_raw$sample_name 
saveRDS(targets_90plus_raw, file = "targets_90plus_raw.rds")

## run full analysis pipeline for each brain region separately
## this includes quality control of samples, probes, normalization and removal of replicates
brainregions <- c("MFG", "DG", "CG", "SN", "LC", "CBM", "CA1" , "EC")
preprocess_bybrainregion(targets_90plus_raw, rgset_90plus_raw)
