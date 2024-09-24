rm(list = ls())

library(tidyverse)
library(openxlsx)
# load data =======================================================================
# cohort information
cohortDat <- list()
cohortDat$donor_info <- read.table("/vol/projects/CIIM/2000HIV/Phenotype/Phenotype_2000HIV_all_01.tsv", 
                          header = TRUE) # recent update by Xun

## protein and metabolite data - all availale sample ----------------------------
load("processedDat/RNAseqDat.RData")
load("processedDat/proteinDat.RData")
load("processedDat/meboDat.RData")
load("processedDat/cytokineDat.RData")
load("processedDat/cellcountDat.RData")

cohortDat$allSample$rna <- RNAseqDat %>%
  slice(which(rownames(RNAseqDat) %in% cohortDat$donor_info$Record.Id))

cohortDat$allSample$protein <- proteinDat %>%
  slice(which(rownames(proteinDat) %in% cohortDat$donor_info$Record.Id))

cohortDat$allSample$mebo <- meboDat %>%
  slice(which(rownames(meboDat) %in% cohortDat$donor_info$Record.Id))

cohortDat$allSample$cytokine <- cytokineDat %>%
  slice(which(rownames(cytokineDat) %in% cohortDat$donor_info$Record.Id))

cohortDat$allSample$cellcount <- cellcountDat %>%
  slice(which(rownames(cellcountDat) %in% cohortDat$donor_info$Record.Id))

## protein and metabolite data - split cohort to discovery and validation ----------------------------
### discovery cohort ----------------------------
cohortDat$Discovery$rna <- RNAseqDat %>% 
  slice(which(rownames(RNAseqDat) %in% 
                (cohortDat$donor_info %>% filter(Cohort == "Discovery"))$Record.Id))

cohortDat$Discovery$protein <- proteinDat %>% 
  slice(which(rownames(proteinDat) %in% 
                (cohortDat$donor_info %>% filter(Cohort == "Discovery"))$Record.Id))

cohortDat$Discovery$mebo <- meboDat %>% 
  slice(which(rownames(meboDat) %in% 
                (cohortDat$donor_info %>% filter(Cohort == "Discovery"))$Record.Id))

cohortDat$Discovery$cytokine <- cytokineDat %>% 
  slice(which(rownames(cytokineDat) %in% 
                (cohortDat$donor_info %>% filter(Cohort == "Discovery"))$Record.Id))

cohortDat$Discovery$cellcount <- cellcountDat %>% 
  slice(which(rownames(cellcountDat) %in% 
                (cohortDat$donor_info %>% filter(Cohort == "Discovery"))$Record.Id))

### validation cohort ----------------------------
cohortDat$Validation$rna <- RNAseqDat %>% 
  slice(which(rownames(RNAseqDat) %in% 
                (cohortDat$donor_info %>% filter(Cohort == "Validation"))$Record.Id))

cohortDat$Validation$protein <- proteinDat %>% 
  slice(which(rownames(proteinDat) %in% 
                (cohortDat$donor_info %>% filter(Cohort == "Validation"))$Record.Id))

cohortDat$Validation$mebo <- meboDat %>%
  slice(which(rownames(meboDat) %in% 
                (cohortDat$donor_info %>% filter(Cohort == "Validation"))$Record.Id))

cohortDat$Validation$cytokine <- cytokineDat %>%
  slice(which(rownames(cytokineDat) %in% 
                (cohortDat$donor_info %>% filter(Cohort == "Validation"))$Record.Id))

cohortDat$Validation$cellcount <- cellcountDat %>% 
  slice(which(rownames(cellcountDat) %in% 
                (cohortDat$donor_info %>% filter(Cohort == "Validation"))$Record.Id))

# save data ------------------------------------------------
save(cohortDat, file = "processedDat/cohortDat.RData")
