rm(list = ls())

library(tidyverse)
library(openxlsx)
# load data =======================================================================
raw_annotation <- read.xlsx("/vol/projects/CIIM/2000HIV/Metabolites/DATA_NORM.xlsx", sheet = "ions") %>%
  select(1:6)

# Used the final updated file listed in the README file (Nienke also used this file)
# It is the final matrix of raw metabolites consisting of 1902 samples and 1720 metabolites:
raw_ionMatrix <- read.csv("/vol/projects/CIIM/2000HIV/Metabolites/2000HIV_GeneralMetabolomics_raw_metabolome_data_outliers_and_pSS_excluded_1902samples_touse.csv")


raw_ionMatrix_v2 <- raw_ionMatrix  %>% select(-matches("ds")) %>% 
  column_to_rownames("Row.names") %>% t() %>% as.data.frame() %>% 
  rownames_to_column("ionMz") %>%
  mutate(ionMz = sprintf("%.4f",as.numeric(gsub("X", "", ionMz)))) # add missing 0s to the end of some ionMz values to match with the ionMz annotation

identical(raw_ionMatrix_v2$ionMz, raw_annotation$ionMz) # TRUE, so the ionMz in ionMatrix is matched with the ionMz in the annotation matrix

# Annotation: filter HMDB endogenous metabolites ------------------------------------------------
# HMDB endogenous metabolites
hmdb_endogenous <- read_csv("/vol/projects/CIIM/cohorts_CMV/CMV_NhanNguyen/reference/20221011_HMDB_endogenousMetabolites")

annotation_temp <- raw_annotation %>%
  slice(which(ionTopFormula %in% hmdb_endogenous$CHEMICAL_FORMULA))
length(unique(annotation_temp$ionIdx)) # 851 metabolites ionIdx
length(unique(annotation_temp$ionTopFormula)) # 851 metabolites formula (1 ionIdx - 1 formulas)

annotation <- annotation_temp %>%
  select(ionMz, ionTopFormula) %>% distinct()

# Covert the raw metabolites data to table format ------------------------------------------------
ionMatrix <- raw_ionMatrix_v2 %>% 
  inner_join(annotation) %>% select(-ionMz) %>%
  column_to_rownames("ionTopFormula") %>%
  t() %>% as.data.frame() # end with 851 formulas

# save data ------------------------------------------------
meboDat <- ionMatrix
save(meboDat, file = "processedDat/meboDat.RData")
