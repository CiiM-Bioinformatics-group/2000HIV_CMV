rm(list = ls())

library(tidyverse)
library(venn)

# load data =======================================================================
# methylation 
#DMSasso_res <- read.table("info_scripts_fromOthers/cmv_season_BMI_SP_cor_eQTM_FDR_rep005_rlm_baconAdjusted.tsv") # old outcome form methylation analysis

DMSasso_res <- read.table("/vol/projects/BIIM/2000HIV/AnalysisDV/Discovery_MethylationDataBased/4EWAS_AllEthnicity/CMV_rlm/ForCMVPaper/eQTM.tsv", 
                          header = TRUE) %>% filter(FDR < 0.05)

## DEGs -------------------------------------------------------
load("processedDat/DEseq2Res_rna.RData")

discovery_padj <- DESeq2_res$discovery %>% as.data.frame %>% filter(padj < 0.05)
validation_pvalue <- DESeq2_res$validation %>% as.data.frame %>% filter(pvalue < 0.05)

DEGs_temp <- DESeq2_res %>% 
  lapply(function(x) x %>% as.data.frame %>% rownames_to_column("gene")) %>% 
  bind_rows(.id = "cohort") %>% 
  dplyr::select(-baseMean, -lfcSE, -stat) %>% 
  pivot_wider(names_from = "cohort", values_from = c("log2FoldChange", "pvalue", "padj")) %>%
  mutate(sig_bothCohorts = ifelse(padj_discovery < 0.05 & pvalue_validation < 0.05, "yes", "no")) %>%
  mutate(validated_DEGs = 
           ifelse(sig_bothCohorts == "yes" & log2FoldChange_discovery * log2FoldChange_validation  > 0, 
                  "yes", "no")) %>%
  filter(padj_discovery < 0.05)

DEGs <- (DEGs_temp %>% filter(validated_DEGs == "yes"))$gene

## DEPs -------------------------------------------------------
load("processedDat/rlmRes_protein.RData")

discovery_padj <- rlm_res$discovery %>% filter(padj < 0.05)
validation_pvalue <- rlm_res$validation %>% filter(pval < 0.05)

DEPs_temp <- rlm_res %>% 
  lapply(function(x) x %>% as.data.frame %>% rownames_to_column("gene")) %>% 
  bind_rows(.id = "cohort") %>% 
  dplyr::select(-StdError, -Z_Value) %>% 
  pivot_wider(names_from = "cohort", values_from = c("effectSize", "pval", "padj")) %>%
  mutate(sig_bothCohorts = ifelse(padj_discovery < 0.05 & pval_validation < 0.05, "yes", "no")) %>%
  mutate(delabel = ifelse(sig_bothCohorts == "yes", gene, NA)) %>%
  mutate(validated_DEPs = 
           ifelse(sig_bothCohorts == "yes" & effectSize_discovery * effectSize_validation  > 0, 
                  "yes", "no")) %>%
  filter(padj_discovery < 0.05)

DEPs <- (DEPs_temp %>% filter(validated_DEPs == "yes"))$gene

# venn diagram (with validated DEGs and DEPs) ----------------------------------------------------------------
valNames <- list( "DMSasso_adj" = DMSasso_res$Gene,
                  "validated_DEGs" = DEGs,
                 "validated_DEPs" = DEPs)

png("output/09_01_venn_sigGeneMethy_Gene_Protein.png", width = 528, height = 528)
venn(valNames, ilcs = 1, sncs = 1)
dev.off()

# list of overlapped between genes related to sig. Methylation sites and DEGs
intersect(DMSasso_res$Gene, DEGs) # 508 genes

# list of overlapped between DEGs and DEPs
intersect(DEGs, DEPs)

intersect(DMSasso_res$Gene, intersect(DEGs, DEPs))
