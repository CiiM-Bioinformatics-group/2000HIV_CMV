rm(list = ls())

library(tidyverse)
library(MASS)
library(lmtest)
library(sandwich)
library("ggrepel")

# inverse ranking normaliation/ transformation (or can use the RNOmni package) 
inormal <- function(x) qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))

# load data =======================================================================
load("/vol/projects/BIIM/cohorts_CMV/CMV_2000HIV_NhanNguyen/processedDat/cohortDat.RData")

# run the rlm() model ------------------------------------------------
mebos <- names(cohortDat$allSample$mebo)

## discovery dataset ------------------------------------------------
top5_geneticPCs_discovery <- read.table(
  "/vol/projects/BIIM/cohorts_CMV/CMV_2000HIV_NhanNguyen/info_scripts_fromOthers/Discovery_allEthnicity_2000HIV.eigenvec",
  header = TRUE) %>% 
  as_tibble() %>% dplyr::select(IID, PC1, PC2, PC3, PC4, PC5)

inputDat_discovery <- cohortDat$donor_info %>%
  inner_join(top5_geneticPCs_discovery, by = c("Record.Id" = "IID")) %>%
  inner_join(cohortDat$Discovery$mebo %>% rownames_to_column("Record.Id")) %>%
  mutate(Institute.Abbreviation = as.factor(Institute.Abbreviation),
         SEX_BIRTH = as.factor(SEX_BIRTH), 
         CMV_IgG_Serology = as.factor(CMV_IgG_Serology))

rlm_res_discovery <- list()
for (mebo in mebos) {
  dat_temp <- inputDat_discovery %>% 
    dplyr::rename("valName" = mebo) %>% 
    mutate(valName = inormal(valName)) # inverse ranking transformation to have normal distribution
  
  ttest <- rlm(valName ~ CMV_IgG_Serology + AGE + SEX_BIRTH + BMI_BASELINE + 
                 Institute.Abbreviation + season_sin + season_cos +
                 PC1 + PC2 + PC3 + PC4 +PC5, 
               data = dat_temp, maxit=200)
  cf <- try(coeftest(ttest, vcov=vcovHC(ttest, type="HC0")))
  rlm_res_discovery[[mebo]] <- cf["CMV_IgG_Serology1", c("Estimate", "Std. Error", "z value", "Pr(>|z|)")]
  
}

outcome_discovery <- rlm_res_discovery %>% as_tibble() %>%
  t() %>%  as.data.frame() %>%
  dplyr::rename("effectSize" = 1, "StdError" = 2, "Z_Value" = 3, "pval" = 4) %>%
  mutate(padj = p.adjust(pval, method = "fdr"))


## validation dataset ------------------------------------------------
top5_geneticPCs_validation <- read.table(
  "/vol/projects/BIIM/cohorts_CMV/CMV_2000HIV_NhanNguyen/info_scripts_fromOthers/Validation_allEthnicity_2000HIV.eigenvec",
  header = TRUE) %>% 
  as_tibble() %>% dplyr::select(IID, PC1, PC2, PC3, PC4, PC5)

inputDat_validation <- cohortDat$donor_info %>%
  inner_join(top5_geneticPCs_validation, by = c("Record.Id" = "IID")) %>%
  inner_join(cohortDat$Validation$mebo %>% rownames_to_column("Record.Id")) %>%
  mutate(Institute.Abbreviation = as.factor(Institute.Abbreviation),
         SEX_BIRTH = as.factor(SEX_BIRTH), 
         CMV_IgG_Serology = as.factor(CMV_IgG_Serology))

sigMebos_discovery <- rownames(outcome_discovery %>% filter(padj < 0.05))

rlm_res_validation <- list()
#for (mebo in mebos) {
for (mebo in sigMebos_discovery) {
  dat_temp <- inputDat_validation %>% 
    dplyr::rename("valName" = mebo) %>% 
    mutate(valName = inormal(valName)) # inverse ranking transformation to have normal distribution
  
  ttest <- rlm(valName ~ CMV_IgG_Serology + AGE + season_sin + season_cos # +
                 #SEX_BIRTH + BMI_BASELINE + PC1 + PC2 + PC3 + PC4 +PC5
               , 
               data = dat_temp, maxit=200)
  cf <- try(coeftest(ttest, vcov=vcovHC(ttest, type="HC0")))
  rlm_res_validation[[mebo]] <- cf["CMV_IgG_Serology1", c("Estimate", "Std. Error", "z value", "Pr(>|z|)")]
}

outcome_validation <- rlm_res_validation %>% as_tibble() %>%
  t() %>%  as.data.frame() %>%
  dplyr::rename("effectSize" = 1, "StdError" = 2, "Z_Value" = 3, "pval" = 4) %>%
  mutate(padj = p.adjust(pval, method = "fdr"))

# save data ------------------------------------------------
rlm_res <- list("discovery" = outcome_discovery, "validation" = outcome_validation)

save(rlm_res, file = "processedDat/rlmRes_mebo.RData")
