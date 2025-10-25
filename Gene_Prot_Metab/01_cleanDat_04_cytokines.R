rm(list = ls())

library(tidyverse)
# load data =======================================================================
# from Xun: it is correct files (after remove the datapoint with max. titer), not transformation have done in this file yet
# then Nienke re-filtered and update the data.
cytokineDat <- read.table("/vol/projects/BIIM/2000HIV/Cytokines/merged_cytokines_cleaned.tsv", sep = "\t") %>%
  t() %>% as.data.frame()

# save data ------------------------------------------------
save(cytokineDat, file = "processedDat/cytokineDat.RData")

