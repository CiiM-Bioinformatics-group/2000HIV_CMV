rm(list = ls())

library(tidyverse)
library(ggpubr)
library(ggrepel)

# load methylation data =======================================================================
# methyDat_sig <- read_tsv("info_scripts_fromOthers/EWAS_outputsOfValidatedSigCpG.tsv")
# identical(methyDat_sig$CpG, methyDat_sig$CpGsite_dis) # all the CpG sites are the same
# identical(methyDat_sig$CpG, methyDat_sig$CpGsite_rep) # all the CpG sites are the same

methyDat_allDiscovery <- readRDS("/vol/projects/BIIM/2000HIV/AnalysisDV/Discovery_MethylationDataBased/4EWAS_AllEthnicity/CMV_rlm/ForCMVPaper/corEWAS_dis_CMV_season_BMI_SP_cor_rlm_bacon.rds")

methyDat <- methyDat_allDiscovery %>%
  rename_with(~ paste0(.x, "_dis")) %>% rownames_to_column("CpG") #%>% full_join(methyDat_sig)


# load RNA data =======================================================================
load("processedDat/DEseq2Res_rna.RData")

rnaDat <- DESeq2_res %>% 
  lapply(function(x) x %>% as.data.frame %>% rownames_to_column("gene")) %>% 
  bind_rows(.id = "cohort") %>% 
  dplyr::select(-baseMean, -lfcSE, -stat) %>% 
  pivot_wider(names_from = "cohort", values_from = c("log2FoldChange", "pvalue", "padj"))

# load protein data =======================================================================
load("processedDat/rlmRes_protein.RData")

proDat <- rlm_res %>% 
  lapply(function(x) x %>% as.data.frame %>% rownames_to_column("gene")) %>% 
  bind_rows(.id = "cohort") %>% 
  dplyr::select(-StdError, -Z_Value) %>% 
  pivot_wider(names_from = "cohort", values_from = c("effectSize", "pval", "padj")) 

# load metabolite data =======================================================================
load("processedDat/rlmRes_mebo.RData")

meboDat <- rlm_res %>% 
  lapply(function(x) x %>% as.data.frame %>% rownames_to_column("gene")) %>% 
  bind_rows(.id = "cohort") %>% 
  dplyr::select(-StdError, -Z_Value) %>% 
  pivot_wider(names_from = "cohort", values_from = c("effectSize", "pval", "padj"))

# load cell count data =======================================================================
load("processedDat/rlmRes_cellcount.RData")

cellcountDat <- rlm_res %>% 
  lapply(function(x) x %>% as.data.frame %>% rownames_to_column("gene")) %>% 
  bind_rows(.id = "cohort") %>% 
  dplyr::select(-StdError, -Z_Value) %>% 
  pivot_wider(names_from = "cohort", values_from = c("effectSize", "pval", "padj"))

# load cytokine data =======================================================================
load("processedDat/rlmRes_cytokine.RData")

cytoDat <- rlm_res %>% 
  lapply(function(x) x %>% as.data.frame %>% rownames_to_column("gene")) %>% 
  bind_rows(.id = "cohort") %>% 
  dplyr::select(-StdError, -Z_Value) %>% 
  pivot_wider(names_from = "cohort", values_from = c("effectSize", "pval", "padj"))


# plot --------------------------------------------------------------------
plotDat <- methyDat %>% # methylation level
  mutate(varName = CpG, 
         effectSize_discovery = Estimate_dis, padj_discovery = FDR_dis,
         type = "methylation") %>% 
  full_join( # RNA / gene expression level
    rnaDat %>%
      mutate(varName = gene, 
             effectSize_discovery = log2FoldChange_discovery, 
             type = "gene expression")
  ) %>%
  full_join(proDat %>% mutate(varName = gene, type = "protein")) %>% # protein level
  full_join(meboDat %>% mutate(varName = gene, type = "metabolite")) %>% # metabolite level
  full_join(cellcountDat %>% mutate(varName = gene, type = "celltype")) %>% # cell count
  full_join(cytoDat %>% mutate(varName = gene, type = "cytokine")) %>% # cytokine production
  dplyr::select(type, varName, effectSize_discovery, padj_discovery) %>% 
  mutate(sig_discovery = ifelse(padj_discovery < 0.05, "yes", "no"),
         label = ifelse(sig_discovery == "yes", "FDR<0.05 in discovery cohort", "measured molecules")) %>%
  mutate(type = factor(type, levels = c("methylation", "gene expression", "protein", "metabolite", "celltype", "cytokine")))


pick <- function(condition){
  function(d) d %>% filter_(condition)
}

omicsOverview <- plotDat %>% 
  ggplot(aes(x = effectSize_discovery, y = type, sig_discovery, fill = label)) + 
  geom_jitter(data = pick(~sig_discovery == "no"), size = 3.5, colour = "#d3d3d3", alpha = 0.5, height = 0.3) + 
  geom_jitter(data = pick(~sig_discovery == "yes"), size = 3.5, colour = "#bb6f7c", alpha = 0.5, height = 0.3) + 
  theme_classic(base_size = 24) +
  theme(plot.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.position = "top")

omicsOverview

# save the plot 
png("output/10_omicsOverview.png", width = 816)
omicsOverview
dev.off()

png("output/10_omicsOverview_shortWidth.png", width = 672)
omicsOverview
dev.off()
