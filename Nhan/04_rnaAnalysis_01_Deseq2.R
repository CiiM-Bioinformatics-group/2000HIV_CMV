rm(list = ls())

library(tidyverse)
library(DESeq2)
library(org.Hs.eg.db)

library("BiocParallel")
register(MulticoreParam(10))

# load data =======================================================================
load("/vol/projects/BIIM/cohorts_CMV/CMV_2000HIV_NhanNguyen/processedDat/cohortDat.RData")

## RNAseq raw data ------------------------------
raw_dat <- readRDS("/vol/projects/BIIM/2000HIV/RNAseq/2000HIV_bulk_transcriptomics_raw_counts.RDS") 
sampleIDs <- readRDS("/vol/projects/BIIM/2000HIV/RNAseq/2000HIV_bulk_transcriptomics_sample_table.RDS")

## convert Ensembl IDs into gene symbols
annots <- select(org.Hs.eg.db, 
                 keys= substring(rownames(raw_dat), 1, 15),
                 columns="SYMBOL", keytype="ENSEMBL")

## gene expression matrix
RNAseqDat <- raw_dat %>% t() %>% 
  as.data.frame %>% rownames_to_column("ID") %>% # convert sample id to donor id
  full_join(sampleIDs %>% dplyr::select(ID, DONOR_ID)) %>% dplyr::select(-ID) %>% 
  column_to_rownames("DONOR_ID") %>% t() %>% as.data.frame %>%
  rownames_to_column("ENSEMBL") %>% # convert ensembl IDs into gene symbols
  mutate(ENSEMBL = substring(ENSEMBL, 1, 15)) %>%  full_join(annots) %>%
  drop_na(SYMBOL) %>% dplyr::select(-ENSEMBL) %>%
  group_by(SYMBOL) %>% summarise_each(funs(mean)) %>% # calculate the average read count for duplicated gene symbols
  column_to_rownames("SYMBOL")

which(duplicated(RNAseqDat$SYMBOL) == TRUE)
unique(duplicated(RNAseqDat$SYMBOL)) # no duplication
# prepare dataset -------------------------------------------------------------
cohorts <- c("Discovery", "Validation")

## top 5 genetic PCs ----------------------------------------------------------
geneticPCs <- list()

for (cohort in cohorts) {
  geneticPCs[[cohort]] <- read.table(
    paste0("/vol/projects/BIIM/cohorts_CMV/CMV_2000HIV_NhanNguyen/info_scripts_fromOthers/",
           cohort, "_allEthnicity_2000HIV.eigenvec"),
    header = TRUE) 
}

top5_geneticPCs <- geneticPCs %>% 
  lapply(function(x) x %>% as_tibble() %>% dplyr::select(IID, PC1, PC2, PC3, PC4, PC5))

## prepare input data files ----------------------------------------------------------
metadata <- list()
overlapSamples <- list()
countDat <- list()
sampleInfo <- list()

for (cohort in cohorts) {
  metadata[[cohort]] <- cohortDat$donor_info %>% 
    dplyr::select(Record.Id, CMV_IgG_Serology, AGE, SEX_BIRTH, BMI_BASELINE, 
           Institute.Abbreviation, season_sin, season_cos) %>% 
    inner_join(top5_geneticPCs[[cohort]], by = c("Record.Id" = "IID")) %>%
    mutate(Institute.Abbreviation = as.factor(Institute.Abbreviation),
           SEX_BIRTH = as.factor(SEX_BIRTH), 
           CMV_IgG_Serology = as.factor(CMV_IgG_Serology)) %>% 
    drop_na(CMV_IgG_Serology)
  
  overlapSamples[[cohort]] <- intersect(metadata[[cohort]]$Record.Id, colnames(RNAseqDat))
  
  ## prepare input data
  countDat[[cohort]] <- RNAseqDat[, overlapSamples[[cohort]]]
  sampleInfo[[cohort]] <- (metadata[[cohort]] %>% column_to_rownames("Record.Id"))[overlapSamples[[cohort]], ]
  
}

# run DEseq2 -------------------------------------------------------------------------
identical(rownames(sampleInfo$Discovery), colnames(countDat$Discovery)) # TRUE, correct sample order -> can run DESeq2 now
identical(rownames(sampleInfo$Validation), colnames(countDat$Validation)) # TRUE, correct sample order -> can run DESeq2 now

DESeq2_res <- list()
DESeq2_resLFC <- list()

## discovery cohort ------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = round(countDat$Discovery), colData = sampleInfo$Discovery,
  design = ~ AGE + SEX_BIRTH + BMI_BASELINE + Institute.Abbreviation + 
    season_sin + season_cos + PC1 + PC2 + PC3 + PC4 +PC5 + 
    CMV_IgG_Serology) # compare CMV status

dds_v2 <- DESeq(dds)

DESeq2_res$discovery <- results(dds_v2, contrast=c("CMV_IgG_Serology", "1","0"), alpha=0.05)
# Shrinkage of effect size (LFC estimates), using the apeglm method (Zhu, Ibrahim, and Love 2018), which improves on the previous estimator.
DESeq2_resLFC$discovery <- lfcShrink(dds_v2, coef="CMV_IgG_Serology_1_vs_0", type="apeglm") 

## validation cohort ------------------------------------------
sigGene_discovery <- rownames(DESeq2_res$discovery %>% 
  as.data.frame %>% filter(padj < 0.05))
countDat$Validation_sigGene <- countDat$Validation[sigGene_discovery,]

dds <- DESeqDataSetFromMatrix(
 # countData = round(countDat$Validation), colData = sampleInfo$Validation,
  countData = round(countDat$Validation_sigGene), colData = sampleInfo$Validation,
  design = ~ AGE + season_sin + season_cos + 
    #SEX_BIRTH + BMI_BASELINE + PC1 + PC2 + PC3 + PC4 +PC5 + 
    CMV_IgG_Serology) # compare CMV status

dds_v2 <- DESeq(dds)

DESeq2_res$validation <- results(dds_v2, contrast=c("CMV_IgG_Serology", "1","0"), alpha=0.05)
# Shrinkage of effect size (LFC estimates), using the apeglm method (Zhu, Ibrahim, and Love 2018), which improves on the previous estimator.
DESeq2_resLFC$validation <- lfcShrink(dds_v2, coef="CMV_IgG_Serology_1_vs_0", type="apeglm") 

# save data ------------------------------------------------
save(DESeq2_res, DESeq2_resLFC, 
     file = "processedDat/DEseq2Res_rna.RData")
