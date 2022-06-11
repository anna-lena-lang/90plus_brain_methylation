# code to plot Table 1
# only used to later copy paste values in excel sheet

library(boot) 
library(dplyr)
library(table1)

# set working directory
setwd("O:/02182022_backup_Lena/90plus")
brainregions <- c("MFG", "CG", "CA1", "DEN", "ERC", "LOC", "SN", "CBL")

# load samplesheet
targets <- readRDS("targets_incl_comorb.rds")

# format labels
targets$conf_diag <- targets$conf_diag %>% as.character() %>% as.factor()
targets$conf_diag <- sub("normal", "Normal",targets$conf_diag)
targets$conf_diag <- sub("demented", "Demented",targets$conf_diag)
targets$conf_diag <- factor(targets$conf_diag, levels= c("Normal", "CIND", "Demented"))

targets$agedeath <- as.numeric(targets$agedeath)
targets$mmse_total <- targets$mmse_total %>% as.character() %>% as.numeric()
targets$tdp43 <-  targets$tdp43 %>% as.character() %>% as.factor()
targets$mvl_score <-  targets$mvl_score %>% as.character() %>% as.factor()
targets$braak_for_lb_orig <- as.factor(targets$braak_for_lb_orig)
targets$adseverityscore <- targets$adseverityscore %>% as.character() %>% as.factor()
targets$adseverityscore <- sub("low", "Low",targets$adseverityscore)
targets$adseverityscore <- sub("intermediate", "Intermediate",targets$adseverityscore)
targets$adseverityscore <- sub("high", "High",targets$adseverityscore)
targets$adseverityscore <- factor(targets$adseverityscore, levels= c("Low", "Intermediate", "High"))

targets$niaaaascore <- targets$niaaaascore %>% as.character() %>% as.factor()
targets$niaaacscore <- targets$niaaacscore %>% as.character() %>% as.factor()
targets$niaaabscore <- targets$niaaabscore %>% as.character() %>% as.factor()

# functions to handle continuous and categorical data
my.render.cont <- function(x) {
  with(stats.apply.rounding(stats.default(x), digits=2), c("",
                                                           "Mean (SD)"=sprintf("%s (&plusmn; %s)", MEAN, SD)))
}
my.render.cat <- function(x) {
  c("", sapply(stats.default(x), function(y) with(y,
                                                  sprintf("%d (%.2f %%)", FREQ, PCT))))
}

# subset to keep one entry for each individual case
targets_47 <- targets[!(duplicated(targets$sample_source)),]

# add labels
labels <- list(
  variables=list(gender = "Sex",
                 agedeath = "Age at death",
                 genotype = "APOE genotype",
                 # mmse_total = "MMSE total",
                 adseverityscore= "AD severity score",
                 niaaaascore = "NIA-AA-A-Score",
                 niaaabscore = "NIA-AA-B-Score",
                 niaaacscore = "NIA-AA-C-Score",
                 conf_diag = "Clinical diagnosis"))
                 
                 # braak_for_lb_orig = "Braak for Lewy Bodies",
                 # tdp43 = "TDP 43 burden",
                 # mvl_score = "MVL burden"))

strata <-c(list(Total=targets_47))
table1(strata, labels,
       render.continuous=my.render.cont, render.categorical=my.render.cat)

# table including more neuropath info
labels <- list(
  variables=list(gender = "Sex",
                 agedeath = "Age at death",
                 mmse_total = "MMSE total",
                 conf_diag = "Clinical diagnosis",
                 genotype = "APOE genotype",
                 adseverityscore= "AD severity score",
                 niaaaascore = "NIA-AA-A-Score",
                 niaaabscore = "NIA-AA-B-Score",
                 niaaacscore = "NIA-AA-C-Score",
                 braak_for_lb_orig = "Braak for Lewy Bodies",
                 tdp43 = "TDP 43 burden",
                 mvl_score = "MVL burden"))

strata <-c(list(Total=targets_47))
table1(strata, labels,
       render.continuous=my.render.cont, render.categorical=my.render.cat)

# table split by apoe genotype
targets_47$genotype <- targets_47$genotype %>% replace_na("missing")
targets_47$apoe_combo <- targets_47$genotype

targets_47$apoe_combo <- sub("23", "4_absent", targets_47$apoe_combo)
targets_47$apoe_combo <- sub("33", "4_absent", targets_47$apoe_combo)

targets_47$apoe_combo <- sub("34", "4_present", targets_47$apoe_combo)
targets_47$apoe_combo <- sub("24", "4_present", targets_47$apoe_combo)

labels <- list(
  variables=list(gender = "Sex",
                 agedeath = "Age at death",
                 mmse_total = "MMSE total",
                 conf_diag = "Clinical diagnosis",
                 # genotype = "APOE genotype",
                 adseverityscore= "AD severity score",
                 niaaaascore = "NIA-AA-A-Score",
                 niaaabscore = "NIA-AA-B-Score",
                 niaaacscore = "NIA-AA-C-Score",
                 braak_for_lb_orig = "Braak for Lewy Bodies",
                 tdp43 = "TDP 43 burden",
                 mvl_score = "MVL burden"))
strata <-c(list(Total=targets_47), split(targets_47, targets_47$apoe_combo))
table1(strata, labels,
       render.continuous=my.render.cont, render.categorical=my.render.cat)
