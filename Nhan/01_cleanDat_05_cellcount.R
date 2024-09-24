rm(list = ls())

library(tidyverse)
library(openxlsx)

# load data (lasted version on NOvember 2023) =======================================================================
path <- "/vol/projects/CIIM/2000HIV/2000HIV_Flowcytometry/UpdatedFlowcytometry/"


fileNames <- c("2000HIV_FLOW_ABS_panel123merged_QCed_untransformed(raw)data_1423samples_356vars_Nov062023.xlsx",  # final absolute counts data
               '2000HIV_FLOW_MFI_panel123merged_QCed_untransformed(raw)data_1434samples_521vars_Nov062023.xlsx',  # final mfi (mean fluorescence intensity?) data 
               '2000HIV_FLOW_MFI_panel123merged_QCed_untransformed(raw)data_1434samples_521vars_Nov062023.xlsx') # final percentage data

names(fileNames) <- c("ABS", "MFI", "PER")

raw_dat <- list()
for (fileName in names(fileNames)) {
  raw_dat[[fileName]] <- read.xlsx(
    paste0(path, fileNames[fileName])) 
}

intersect(names(raw_dat$ABS), names(raw_dat$MFI)) # only sample name column are the same
intersect(names(raw_dat$ABS), names(raw_dat$PER)) # only sample name column are the same


# # load data (old version, and should not use any more) =======================================================================
# # Use files based on the informaiton in the readme file
# path <- "/vol/projects/BIIM/2000HIV/2000HIV_Flowcytometry/Filtered2_merged_final/"
# 
# fileNames <- c("2000HIV_FLOW_ABS_panel123merged_QCed_untransformed(raw)data_1423samples_415vars.xlsx",  # final absolute counts data
#                '2000HIV_FLOW_MFI_panel123merged_QCed_untransformed(raw)data_1434samples_516vars.xlsx',  # final mfi (mean fluorescence intensity?) data 
#                '2000HIV_FLOW_PER_panel123merged_QCed_untransformed(raw)data_1434samples_403vars.xlsx') # final percentage data
# names(fileNames) <- c("ABS", "MFI", "PER")
# 
# raw_dat <- list()
# for (fileName in names(fileNames)) {
#   raw_dat[[fileName]] <- read.xlsx(
#     paste0(path, fileNames[fileName])) 
# }
# 
# intersect(names(raw_dat$ABS), names(raw_dat$MFI))
# intersect(names(raw_dat$ABS), names(raw_dat$PER)) # 404 colnames overlaps
# 
# cellcountDat <- raw_dat$ABS %>% column_to_rownames("SampleID")

# save data ------------------------------------------------
cellcountDat <- raw_dat$ABS %>% column_to_rownames("SampleID")
save(cellcountDat, file = "processedDat/cellcountDat.RData")

# cellcountDat_PER <- raw_dat$PER %>% column_to_rownames("SampleID")
# save(cellcountDat_PER, file = "processedDat/cellcountDat_PER.RData")
