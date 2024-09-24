rm(list = ls())

library(tidyverse)
library(ggpubr)
library(ggrepel)

# load data =======================================================================
load("processedDat/cohortDat.RData")
compare_CMV <- list( c("CMV-", "CMV+"))
palette_colors = c("#A86A9D", "#E9C61D")

# cytokine ----------------------------------

## functions  --------------------------------------------------------------------
get.log2 <- function(dat) {
  library(tidyverse)
  library(rstatix)
  outcome <- dat %>% mutate(across(where(is.numeric), ~log2(.x))) %>%
    mutate(across(where(is.numeric), ~ifelse(is.infinite(.x), 0, .x)))
  return(outcome)
}

## prepare data  --------------------------------------------------------------------
inputDat <- cohortDat$donor_info %>% 
  right_join(cohortDat$allSample$cytokine %>% rownames_to_column("Record.Id")) %>%
  mutate(CMV_IgG_Serology = ifelse(CMV_IgG_Serology == 1, "CMV+",
                                   ifelse(CMV_IgG_Serology == 0, "CMV-", CMV_IgG_Serology))) %>%
  drop_na(CMV_IgG_Serology)

inputDat_log10 <- cohortDat$donor_info %>% 
  right_join(cohortDat$allSample$cytokine %>% mutate_all(~log10(.x)) %>% rownames_to_column("Record.Id")) %>%
  mutate(CMV_IgG_Serology = ifelse(CMV_IgG_Serology == 1, "CMV+",
                                   ifelse(CMV_IgG_Serology == 0, "CMV-", CMV_IgG_Serology))) %>%
  drop_na(CMV_IgG_Serology)

## get p-value from linear model ---------------------------
load("processedDat/rlmRes_cytokine.RData")
stat_rlm <- rlm_res %>% 
  lapply(function(x) x %>% as.data.frame %>% 
           rownames_to_column("cyto") %>%
           mutate(cyto = gsub("pbmc_", "", cyto))) %>% 
  bind_rows(.id = "cohort") %>% 
  mutate(group1 = "CMV+", group2 = "CMV-") %>%
  mutate(Cohort = ifelse(cohort == "discovery", "Discovery", "Validation")) %>% 
  mutate(pval = mapply(format, pval)) %>%
  dplyr::select(cyto, Cohort, group1, group2, pval)

## boxplot with p-value linear model ---------------------------

# individual boxplot for each cytokine
cytokine <- "pbmc_24h_il1b_cmv"
cytokine <- "pbmc_24h_il1ra_cmv"
cytokine <- "pbmc_24h_il8_cmv"   
cytokine <- "pbmc_24h_mcp1_cmv"

cytokine <- "pbmc_7d_il22_mtb" # this for Hua

# boxplot
inputDat %>% 
  ggboxplot(x = "CMV_IgG_Serology", y = cytokine,
            color = "CMV_IgG_Serology", palette = palette_colors,
            add = "jitter", add.params = list(size = 2, alpha = 0.5)) + 
  facet_wrap(~Cohort, nrow = 1) +
  stat_pvalue_manual(stat_rlm %>% filter(cyto %in% gsub("pbmc_", "", cytokine)), 
                     label = "pval", size = 5, 
                     y.position = max(inputDat[[cytokine]], na.rm = TRUE))  + 
  theme_classic(base_size = 16) +
  theme(plot.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.position = "none") + 
  ylab(paste0("log10( ", gsub("pbmc_", "", cytokine), ")"))


boxplot_cytokine <- inputDat_log10 %>% 
  ggboxplot(x = "CMV_IgG_Serology", y = cytokine,
            color = "CMV_IgG_Serology", palette = palette_colors,
            add = "jitter", add.params = list(size = 2, alpha = 0.5)) + 
  facet_wrap(~Cohort, nrow = 1) +
  stat_pvalue_manual(stat_rlm %>% filter(cyto %in% gsub("pbmc_", "", cytokine)), 
                     label = "pval", size = 6, 
                     y.position = max(inputDat_log10[[cytokine]], na.rm = TRUE))  + 
  theme_classic(base_size = 20) +
  theme(plot.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.position = "none") + 
  ylab(paste0("log10( ", gsub("pbmc_", "", cytokine), ")"))

boxplot_cytokine

# save the plot 
png(paste0("output/10_boxplot_cytokine_" ,cytokine, ".png"))
boxplot_cytokine
dev.off()

rm(cytokine_NAs, inputDat, inputDat_log10, rlm_res, stat_rlm, cytokine, get.log2)

# rna  --------------------------------------------------------------------

## prepare data 
inputDat <- cohortDat$donor_info %>% 
  right_join(cohortDat$allSample$rna %>% rownames_to_column("Record.Id")) %>%
  mutate(CMV_IgG_Serology = ifelse(CMV_IgG_Serology == 1, "CMV+",
                                   ifelse(CMV_IgG_Serology == 0, "CMV-", CMV_IgG_Serology))) %>%
  drop_na(CMV_IgG_Serology)

## get p-value from linear model ---------------------------
load("processedDat/DEseq2Res_rna.RData")

stat_DESeq2 <- DESeq2_res %>% 
  lapply(function(x) x %>% as.data.frame %>% 
           rownames_to_column("gene") ) %>% 
  bind_rows(.id = "cohort") %>% 
  mutate(group1 = "CMV+", group2 = "CMV-") %>%
  mutate(Cohort = ifelse(cohort == "discovery", "Discovery", "Validation")) %>% 
  mutate(pvalue = mapply(format, pvalue),
         padj = mapply(format, padj)) %>%
  dplyr::select(gene, Cohort, group1, group2, pvalue, padj)

## boxplot with p-value linear model ---------------------------

# individual boxplot for each gene
rna <- "FCRL6"
rna <- "GZMH"
rna <- "KLRD1"
rna <- "ADGRG1"
rna <- "GNLY"
rna <- "KIR2DL3"
rna <- "ITGAL"
rna <- "GBP1"
rna <- "CRTAM"
rna <- "DNMT3A" # for Xun
rna <- "KIR2DS4" #for Hua

boxplot_RNAseq <- inputDat %>% 
  ggboxplot(x = "CMV_IgG_Serology", y = rna,
            color = "CMV_IgG_Serology", palette = palette_colors,
            add = "jitter", add.params = list(size = 2, alpha = 0.5)) + 
  facet_wrap(~Cohort, nrow = 1) +
  stat_pvalue_manual(stat_DESeq2%>% filter(gene %in% rna), 
                     label = "pvalue", 
                     #label = "padj", 
                     size = 6, 
                     y.position = max(inputDat[[rna]], na.rm = TRUE))  + 
  theme_classic(base_size = 20) +
  theme(plot.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.position = "none")

boxplot_RNAseq

# save the plot 
png(paste0("output/10_boxplot_RNAseq_" ,rna, ".png"))
boxplot_RNAseq
dev.off()

rm(DESeq2_res, DESeq2_resLFC, stat_DESeq2, rna)

# protein  ----------------------------------

## prepare data 
inputDat <- cohortDat$donor_info %>% 
  right_join(cohortDat$allSample$protein %>% rownames_to_column("Record.Id")) %>%
  mutate(CMV_IgG_Serology = ifelse(CMV_IgG_Serology == 1, "CMV+",
                                   ifelse(CMV_IgG_Serology == 0, "CMV-", CMV_IgG_Serology))) %>%
  drop_na(CMV_IgG_Serology)

## get p-value from linear model for validated DEPs ---------------------------
load("processedDat/rlmRes_protein.RData")

stat_rlm <- rlm_res %>% 
  lapply(function(x) x %>% as.data.frame %>% 
           rownames_to_column("protein")) %>% 
  bind_rows(.id = "cohort") %>% 
  mutate(group1 = "CMV+", group2 = "CMV-") %>%
  mutate(Cohort = ifelse(cohort == "discovery", "Discovery", "Validation")) %>% 
  mutate(pval = mapply(format, pval)) %>%
  dplyr::select(protein, Cohort, group1, group2, pval)

### boxplot with p-value linear model ---------------------------

# individual boxplot with t-test
protein_name <- "FCRL6"
protein_name <- "GZMH"
protein_name <- "KLRD1"
protein_name <- "ADGRG1"
protein_name <- "GNLY"
protein_name <- "KIR2DL3"
protein_name <- "ITGAL"
protein_name <- "GBP1"
protein_name <- "CRTAM"
protein_name <- "KIR2DS4" #for Hua -> no p-value because it is non-validated DEPs

boxplot_protein <- inputDat %>% 
  ggboxplot(x = "CMV_IgG_Serology", y = protein_name,
            color = "CMV_IgG_Serology", palette = palette_colors,
            add = "jitter", add.params = list(size = 2, alpha = 0.5)) + 
  facet_wrap(~Cohort, nrow = 1) +
  stat_pvalue_manual(stat_rlm %>% filter(protein %in% protein_name), 
                     label = "pval", size = 6, 
                     y.position = max(inputDat[[protein_name]], na.rm = TRUE))  + 
  theme_classic(base_size = 20) +
  theme(plot.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.position = "none")

boxplot_protein

# save the plot 
png(paste0("output/10_boxplot_protein_" ,protein_name, ".png"))
boxplot_protein
dev.off()

rm(rlm_res, stat_rlm, protein_name)

## get p-value from linear model for non-validated DEPs ---------------------------
load("processedDat/rlmRes_protein_forBoxplot.RData")

stat_rlm <- rlm_res %>% 
  lapply(function(x) x %>% as.data.frame %>% 
           rownames_to_column("protein")) %>% 
  bind_rows(.id = "cohort") %>% 
  mutate(group1 = "CMV+", group2 = "CMV-") %>%
  mutate(Cohort = ifelse(cohort == "discovery", "Discovery", "Validation")) %>% 
  mutate(pval = mapply(format, pval)) %>%
  dplyr::select(protein, Cohort, group1, group2, pval)

### boxplot with p-value linear model ---------------------------

# individual boxplot with t-test
protein_name <- "KIR2DS4" # for Hua

boxplot_protein <- inputDat %>% 
  ggboxplot(x = "CMV_IgG_Serology", y = protein_name,
            color = "CMV_IgG_Serology", palette = palette_colors,
            add = "jitter", add.params = list(size = 2, alpha = 0.5)) + 
  facet_wrap(~Cohort, nrow = 1) +
  stat_pvalue_manual(stat_rlm %>% filter(protein %in% protein_name), 
                     label = "pval", size = 6, 
                     y.position = max(inputDat[[protein_name]], na.rm = TRUE))  + 
  theme_classic(base_size = 20) +
  theme(plot.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.position = "none")

boxplot_protein

# save the plot 
png(paste0("output/10_boxplot_protein_" ,protein_name, ".png"))
boxplot_protein
dev.off()

rm(inputDat, rlm_res, stat_rlm, protein_name)

