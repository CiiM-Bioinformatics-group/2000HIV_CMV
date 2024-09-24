rm(list = ls())

library(tidyverse)
library(venn)
library(gridExtra)
library(ggpubr)
library(ggrepel)

# cell count level, compare DEs between discovery and validation cohorts ==============
load("processedDat/rlmRes_cellcount.RData")
#load("processedDat/rlmRes_cellcount_addedConfounders.RData")

plotDat <- rlm_res %>% 
  lapply(function(x) x %>% as.data.frame %>% rownames_to_column("cellCounts")) %>% 
  bind_rows(.id = "cohort") %>% 
  dplyr::select(-StdError, -Z_Value) %>% 
  pivot_wider(names_from = "cohort", values_from = c("effectSize", "pval", "padj")) %>%
  mutate(sig_bothCohorts = ifelse(padj_discovery < 0.05 & pval_validation < 0.05, "yes", "no")) %>%
  mutate(delabel = ifelse(sig_bothCohorts == "yes", cellCounts, NA)) %>%
  mutate(validated_cellCounts = 
           ifelse(sig_bothCohorts == "yes" & effectSize_discovery * effectSize_validation  > 0, 
                  "yes", "no")) %>%
  filter(padj_discovery < 0.05)

## label the top 10 cell types -----------------------------------------
labelCells <- plotDat %>% top_n(10, effectSize_validation)
#labelCells <- plotDat %>% top_n(10, effectSize_discovery)

plotDat_v2 <- plotDat %>%
  mutate(delabel_v2 = ifelse(cellCounts %in% labelCells$cellCounts, cellCounts, NA)) %>%
  mutate(delabel_v3 = substring(delabel_v2, 13, nchar(delabel_v2)), 
         validated_cellCounts = factor(validated_cellCounts, levels = c("no", "yes")))

## scatter plot -----------------------------------------
comparePlot_cellcount <- ggplot(data = plotDat_v2, 
       aes(x = effectSize_discovery, y = effectSize_validation, 
           color = validated_cellCounts, label = delabel_v3)) + 
  geom_point(size = 3, alpha = 0.5) + 
  scale_color_manual(values=c("no" = "#d3d3d3", "yes" = "#bb6f7c"),
                     labels = c("no", "yes")) +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0)+
  geom_text_repel(size = 6, 
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 40),
                  box.padding = unit(0.3, "lines"), show.legend = FALSE) +
  theme_classic(base_size = 20) +
  theme(plot.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.position = "top")

comparePlot_cellcount

# save the plot 
png("output/10_comparePlot_cellcount.png", width = 600, height = 528)
comparePlot_cellcount
dev.off()

rm(labelCells, plotDat, plotDat_v2, rlm_res)

# RNA level, compare DEGs between discovery and validation cohorts ==============
load("processedDat/DEseq2Res_rna.RData")

plotDat <- DESeq2_res %>% 
  lapply(function(x) x %>% as.data.frame %>% rownames_to_column("gene")) %>% 
  bind_rows(.id = "cohort") %>% 
  dplyr::select(-baseMean, -lfcSE, -stat) %>% 
  pivot_wider(names_from = "cohort", values_from = c("log2FoldChange", "pvalue", "padj")) %>%
  mutate(sig_bothCohorts = ifelse(padj_discovery < 0.05 & pvalue_validation < 0.05, "yes", "no")) %>%
  mutate(validated_DEGs = 
           ifelse(sig_bothCohorts == "yes" & log2FoldChange_discovery * log2FoldChange_validation  > 0, 
                  "yes", "no")) %>%
  filter(padj_discovery < 0.05)

## label DEGs which also eQTM from validated sig. CpG sites at methylation levels -----------------------------------------
# label_DEGs <- c("LTBP3", "FCRL6", "ZNF418", "PATL2", "ZNF256",  "ZNF814", 
#                 "GSTM2","TESPA1", "ADAMTS1", "SPG11", "MXRA7", "NLRC5", 
#                 "WNT10B", "ZFP28", "GFI1", "ZNF677" ) # based on old methylation analysis

label_DEGs <- c("KIR2DL3", "ITGAL", "FCRL6", "ADGRG1", "GZMH", "GBP1", "GNLY", "KLRD1" ) # results based on updated methylation analysis

plotDat_v2 <- plotDat %>%
  mutate(label_DEGs = ifelse(gene %in% label_DEGs, gene, NA))

## scatter plot -----------------------------------------
comparePlot_RNAseq <- ggplot(data = plotDat_v2, 
       aes(x = log2FoldChange_discovery, y = log2FoldChange_validation,
           color = validated_DEGs, label = label_DEGs)) + 
  geom_point(size = 4, alpha = 0.5) + 
  scale_color_manual(values=c("no" = "#d3d3d3", "yes" = "#bb6f7c"),
                     labels = c("no", "yes")) +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0)+
  geom_label_repel(size = 6,
                   max.overlaps = getOption("ggrepel.max.overlaps", default = 40),
                   box.padding = unit(0.6, "lines"), show.legend = FALSE) +
  theme_classic(base_size = 20) +
  theme(plot.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.position = "top")

comparePlot_RNAseq

# save the plot 
png("output/10_comparePlot_RNAseq.png", width = 528)
comparePlot_RNAseq
dev.off()

rm(DESeq2_res, DESeq2_resLFC, plotDat, plotDat_v2, label_DEGs)

# protein level, compare DEGs between discovery and validation cohorts ==============
load("processedDat/rlmRes_protein.RData")

plotDat <- rlm_res %>% 
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

## label DEPs which are also validated DEGs -----------------------------------------
# label_DEPs <- c("ADGRG1", "CRTAM", "FCRL6", "GNLY", "GZMH", "ITGAL", "KIR2DL3", "KLRD1") # based on old analysis

label_DEPs <- c("KIR2DL3", "ITGAL", "FCRL6", "ADGRG1", "GZMH", "GBP1", "GNLY", "KLRD1" ) # results based on updated analysis

plotDat_v2 <- plotDat %>%
  mutate(label_DEPs = ifelse(gene %in% label_DEPs, gene, NA))

## scatter plot -----------------------------------------
comparePlot_protein <- ggplot(data = plotDat_v2, 
       aes(x = effectSize_discovery, y = effectSize_validation, 
           color = validated_DEPs, label = label_DEPs)) + 
  geom_point(size = 3, alpha = 0.5) + 
  scale_color_manual(values=c("no" = "#d3d3d3", "yes" = "#bb6f7c"),
                     labels = c("no", "yes")) +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0)+
  geom_label_repel(size = 6,
                   max.overlaps = getOption("ggrepel.max.overlaps", default = 40),
                   box.padding = unit(0.6, "lines"), show.legend = FALSE) +
  theme_classic(base_size = 20) +
  theme(plot.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.position = "top")

comparePlot_protein

# save the plot 
png("output/10_comparePlot_protein.png", width = 528)
comparePlot_protein
dev.off()

rm(plotDat, plotDat_v2, rlm_res, label_DEPs)
