rm(list = ls())

library(tidyverse)
library(openxlsx)
# load data =======================================================================
rawDat <- read_table("/vol/projects/CIIM/2000HIV/Olinkdata/updated_OLINK_Explore3072_1910samples_2367proteins_after_bridging_normalization_QCed_21nov2022.tsv")
proteinDat <- rawDat %>% column_to_rownames("SampleID")

# save data ------------------------------------------------
save(proteinDat, file = "processedDat/proteinDat.RData")
