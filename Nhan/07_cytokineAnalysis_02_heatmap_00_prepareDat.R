# for the heatmap, we need the side effects for all cytokine production in the validation cohort
# we do not have that in our cytokine production analysis
# so I wrote this code file to calculate the side effects for all cytokine production in the validation cohort

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
cytokines <- names(cohortDat$allSample$cytokine)

## discovery dataset ------------------------------------------------
top5_geneticPCs_discovery <- read.table(
  "/vol/projects/BIIM/cohorts_CMV/CMV_2000HIV_NhanNguyen/info_scripts_fromOthers/Discovery_allEthnicity_2000HIV.eigenvec",
  header = TRUE) %>% 
  as_tibble() %>% dplyr::select(IID, PC1, PC2, PC3, PC4, PC5)

inputDat_discovery <- cohortDat$donor_info %>%
  inner_join(top5_geneticPCs_discovery, by = c("Record.Id" = "IID")) %>%
  inner_join(cohortDat$Discovery$cytokine %>% rownames_to_column("Record.Id")) %>%
  mutate(Institute.Abbreviation = as.factor(Institute.Abbreviation),
         SEX_BIRTH = as.factor(SEX_BIRTH), 
         CMV_IgG_Serology = as.factor(CMV_IgG_Serology))

rlm_res_discovery <- list()
#rlm_cytokine_discovery <- list()
for (cytokine in cytokines) {
  dat_temp <- inputDat_discovery %>% 
    rename("valName" = cytokine) %>% 
    mutate(valName = inormal(valName)) %>% # inverse ranking transformation to have normal distribution
    mutate(isNAs = ifelse(is.na(valName), "NA", "nonNA"))
  
  ttest <- rlm(valName ~ CMV_IgG_Serology + AGE + SEX_BIRTH + BMI_BASELINE + 
                 Institute.Abbreviation + season_sin + season_cos +
                 PC1 + PC2 + PC3 + PC4 +PC5, 
               data = dat_temp, maxit=200)
  cf <- try(coeftest(ttest, vcov=vcovHC(ttest, type="HC0")))
  
  rlm_res_discovery[[cytokine]] <- cf["CMV_IgG_Serology1", c("Estimate", "Std. Error", "z value", "Pr(>|z|)")]
  #rlm_cytokine_discovery[[cytokine]] <- dat_temp %>% count(isNAs, CMV_IgG_Serology)
}

outcome_discovery <- rlm_res_discovery %>% as_tibble() %>%
  t() %>%  as.data.frame() %>%
  rename("effectSize" = 1, "StdError" = 2, "Z_Value" = 3, "pval" = 4) %>%
  mutate(padj = p.adjust(pval, method = "fdr")) 

# cytokine_discovery <- rlm_cytokine_discovery %>% 
#   bind_rows(.id = "cytokine" )

## validation dataset ------------------------------------------------
top5_geneticPCs_validation <- read.table(
  "/vol/projects/BIIM/cohorts_CMV/CMV_2000HIV_NhanNguyen/info_scripts_fromOthers/Validation_allEthnicity_2000HIV.eigenvec",
  header = TRUE) %>% 
  as_tibble() %>% dplyr::select(IID, PC1, PC2, PC3, PC4, PC5)

inputDat_validation <- cohortDat$donor_info %>%
  inner_join(top5_geneticPCs_validation, by = c("Record.Id" = "IID")) %>%
  inner_join(cohortDat$Validation$cytokine %>% rownames_to_column("Record.Id")) %>%
  mutate(SEX_BIRTH = as.factor(SEX_BIRTH), 
         CMV_IgG_Serology = as.factor(CMV_IgG_Serology))

rm_cytos <- c("pbmc_24h_tnf_rpmi", "pbmc_24h_il6_rpmi", "pbmc_7d_il17_rpmi", 
              "pbmc_7d_il22_rpmi") # cytokine only have measured values in CMV+ in validation cohort
selected_cytos <- cytokines[-which(cytokines %in% rm_cytos)]

rlm_res_validation <- list()
#rlm_cytokine_validation <- list()
for (cytokine in selected_cytos) {
  dat_temp <- inputDat_validation %>% 
    rename("valName" = cytokine) %>% 
    mutate(valName = inormal(valName)) %>% # inverse ranking transformation to have normal distribution
    mutate(isNAs = ifelse(is.na(valName), "NA", "nonNA"))
  
  ttest <- rlm(valName ~ CMV_IgG_Serology + AGE + SEX_BIRTH #+ 
               #  BMI_BASELINE + season_sin + season_cos + PC1 + PC2 + PC3 + PC4 +PC5
               , 
               data = dat_temp, maxit=200)
  cf <- try(coeftest(ttest, vcov=vcovHC(ttest, type="HC0")))
  
  rlm_res_validation[[cytokine]] <- cf["CMV_IgG_Serology1", c("Estimate", "Std. Error", "z value", "Pr(>|z|)")]
 # rlm_cytokine_validation[[cytokine]] <- dat_temp %>% count(isNAs, CMV_IgG_Serology)
}

outcome_validation <- rlm_res_validation %>% as_tibble() %>%
  t() %>%  as.data.frame() %>%
  rename("effectSize" = 1, "StdError" = 2, "Z_Value" = 3, "pval" = 4) %>%
  mutate(padj = p.adjust(pval, method = "fdr"))

# cytokine_validation <- rlm_cytokine_validation %>% 
#   bind_rows(.id = "cytokine" )

# save data ------------------------------------------------
rlm_res <- list("discovery" = outcome_discovery, "validation" = outcome_validation)
#cytokine_NAs <- list("discovery" = cytokine_discovery, "validation" = cytokine_validation)

save(rlm_res, #cytokine_NAs, 
     file = "processedDat/rlmRes_cytokine_forHeatmap.RData")
