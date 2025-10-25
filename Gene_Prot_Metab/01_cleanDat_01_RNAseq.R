rm(list = ls())

library(tidyverse)
library(gplots)
# load data =======================================================================
# from README filee: it is the "Normalized expression matrix of counts with genes filtered by clinical group"
raw_dat <- readRDS("/vol/projects/BIIM/2000HIV/RNAseq/2000HIV_bulk_transcriptomics_normalized_counts.RDS") 
identical(rownames(raw_dat), raw_dat$GENEID)

# check gene IDs which match with identical Symbols ------------------------------
# all cases of the same symbol, different gene IDs
genes_sameSymbol <- raw_dat %>% 
  filter(SYMBOL %in% raw_dat$SYMBOL[(duplicated(raw_dat$SYMBOL))]) %>%
  mutate(valName = paste0(SYMBOL, "_", GENEID)) %>%
  remove_rownames() %>% column_to_rownames("valName") 
genes_sameSymbol %>% dplyr::count(SYMBOL, sort = TRUE)

# Y_RNA case with 35 different gene IDs
genes_sameSymbol_yRNA <- genes_sameSymbol %>%
  filter(SYMBOL == "Y_RNA") %>% 
  select(-GENEID, -SYMBOL, -GENETYPE, -DESCRIPTION, -CHR) %>% 
  t()

# cases including Metazoa_SRP, U1, uc_338
genes_sameSymbol_noYrna_3cases <- genes_sameSymbol %>% 
  filter(SYMBOL %in% c("Metazoa_SRP", "U1", "uc_338")) %>% 
  select(-GENEID, -SYMBOL, -GENETYPE, -DESCRIPTION, -CHR) %>% 
  t()

# other cases
genes_sameSymbol_otherCase <- genes_sameSymbol %>% 
  filter(SYMBOL %in% c("BMS1P4", "COG8", "LINC01422", "LINC01481", "MATR3", 
                       "POLR2J4", "RGS5", "SNORA12", "SNORA31", "U3")) %>% 
  select(-GENEID, -SYMBOL, -GENETYPE, -DESCRIPTION, -CHR) %>% 
  t()

# cor plot for genes with the same symbol, different gene IDs
inputDat <- genes_sameSymbol_yRNA
inputDat <- genes_sameSymbol_noYrna_3cases
inputDat <- genes_sameSymbol_otherCase

corMat <- round(cor(inputDat), 2)
heatmap.2(corMat, margins = c(15, 15), 
          trace = "none", density.info = "none",
          col = "bluered")
# Note: gene with the same symbols but different gene IDs do not positively correlate to each others

# prepare the gene expression matrix ------------------------------
RNAseqDat_temp <- raw_dat %>% 
  mutate(valName = ifelse(SYMBOL %in% raw_dat$SYMBOL[(duplicated(raw_dat$SYMBOL))], 
                          paste0(SYMBOL, "_", GENEID), SYMBOL)) %>%
  remove_rownames() %>% column_to_rownames("valName") %>%
  select(-GENEID, -SYMBOL, -GENETYPE, -DESCRIPTION, -CHR) %>% 
  t()

# match sample IDs with patientID ------------------------------
sampleIDs <- readRDS("/vol/projects/BIIM/2000HIV/RNAseq/2000HIV_bulk_transcriptomics_sample_table.RDS")

RNAseqDat <- RNAseqDat_temp %>% 
  as.data.frame %>% rownames_to_column("ID") %>%
  full_join(sampleIDs %>% select(ID, DONOR_ID)) %>% 
  select(-ID) %>% 
  column_to_rownames("DONOR_ID")

# save data ------------------------------------------------
save(RNAseqDat, file = "processedDat/RNAseqDat.RData")

