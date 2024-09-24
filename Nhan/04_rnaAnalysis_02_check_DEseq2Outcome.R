rm(list = ls())

library(tidyverse)
library(venn)
library(gridExtra)
library(ggpubr)
library(ggrepel)
# load data =======================================================================
load("processedDat/DEseq2Res_rna.RData")

# venn diagram ---------------------------------
discovery_padj <- DESeq2_res$discovery %>% as.data.frame %>% 
  filter(padj < 0.05)
validation_pvalue <- DESeq2_res$validation %>% as.data.frame %>% 
  filter(pvalue < 0.05)

valNames <- list("discovery_padj" = rownames(discovery_padj),
                 "validation_pvalue" = rownames(validation_pvalue))
venn(valNames)

selected_vals <- intersect(rownames(discovery_padj), rownames(validation_pvalue))

# gene name and effected size ---------------------------------------------
rna_effects <- discovery_padj %>% 
  dplyr::select(log2FoldChange) %>% rownames_to_column("gene") %>% 
  filter(gene %in% selected_vals)

write.table(rna_effects, file = "processedDat/sigRNAs_log2FoldChange.txt", row.names = FALSE, quote = FALSE)


# ploting to check data, these plots are not in the manuscript ----------------------------------

# ## compare plot ----------------------------------
# plotDat <- DESeq2_res %>% 
#   lapply(function(x) x %>% as.data.frame %>% rownames_to_column("gene")) %>% 
#   bind_rows(.id = "cohort") %>% 
#   dplyr::select(-baseMean, -lfcSE, -stat) %>% 
#   pivot_wider(names_from = "cohort", values_from = c("log2FoldChange", "pvalue", "padj")) %>%
#   mutate(sig_bothCohorts = ifelse(padj_discovery < 0.05 & pvalue_validation < 0.05, "yes", "no")) %>%
#   mutate(validated_DEGs = 
#            ifelse(sig_bothCohorts == "yes" & log2FoldChange_discovery * log2FoldChange_validation  > 0, 
#                   "yes", "no")) %>%
#   filter(padj_discovery < 0.05)
# 
# ggplot(data = plotDat, 
#        aes(x = log2FoldChange_discovery, y = log2FoldChange_validation, color = validated_DEGs)) + 
#   geom_point(alpha = 0.8) + 
#   scale_color_manual(values=c("grey","red")) +
#   geom_point(data = plotDat %>% filter(validated_DEGs == "yes"),
#              aes(x = log2FoldChange_discovery, y = log2FoldChange_validation), 
#              alpha = 0.8, colour = "red") + 
#   geom_vline(xintercept = 0) + geom_hline(yintercept = 0)+
#   theme_classic()
# 
# plotDat %>% group_by(validated_DEGs) %>% summarise(countNum = n()) %>%
#   mutate(groupPercent = (countNum/sum(countNum))*100)
# 
# 
# ## volcano plot ----------------------------------
# outcome <- DESeq2_res$discovery %>% as.data.frame
# 
# plotDat <- outcome %>%
#   rownames_to_column("valName") %>%
#   mutate(sig = ifelse(padj < 0.05, "padj<0.05", "Not Sig"),
#          delabel = ifelse(sig == "Not Sig", NA, valName))
# 
# plotDat %>%
#   ggplot(aes(log2FoldChange, -log10(padj), col = sig, label = delabel)) +
#   geom_point() +
#   scale_color_manual(values=c("black", "red")) +
#   ggtitle("Differential expression between CMV+ vs. CMV-") +
#   geom_hline(yintercept=-log10(0.05), col="red") +
#   #geom_text_repel(max.overlaps = getOption("ggrepel.max.overlaps", default = 30)) +
#   theme_bw()
# 
# ## boxplot ----------------------------------
# load("processedDat/cohortDat.RData")
# 
# inputDat <- cohortDat$donor_info %>% 
#   right_join(cohortDat$allSample$rna %>% rownames_to_column("Record.Id")) %>%
#   mutate(CMV_IgG_Serology = ifelse(CMV_IgG_Serology == 1, "CMV+",
#                                    ifelse(CMV_IgG_Serology == 0, "CMV-", CMV_IgG_Serology))) %>%
#   drop_na(CMV_IgG_Serology)
# 
# # set up comparison
# compare_CMV <- list( c("CMV-", "CMV+"))
# 
# # individual boxplot
# rna <- "FCRL6"
# rna <- "ITGAL"
# rna <- "DNAJA4"
# rna <- "DNMT3A"
# rna <- "CRABP1" # NOT FOUND?
# rna <- "CTNND2" # NOT FOUND?
# rna <- "GZMH"
# rna <- "CRTAM"
# rna <- "GNLY"
# rna <- "GBP1"
# rna <- "KLRD1"
# rna <- "ADGRG1"
# rna <- "KIR2DL3"
# rna <- "DNMT3A"
# 
# rna <- "SPG11"   
# rna <-"NLRC5"  
# rna <-"ADAMTS1" 
# rna <-"GFI1"    
# rna <-"TAP1"    
# rna <-"WNT10B"  
# rna <-"MXRA7"   
# rna <-"ZNF418"  
# rna <-"TAP2"    
# rna <-"GSTM2"   
# rna <-"PATL2"  
# 
# inputDat %>% 
#   ggboxplot(x = "CMV_IgG_Serology", y = rna,
#             paletter = "jco", add = "jitter") + facet_wrap(~Cohort, nrow = 1) +
#   stat_compare_means(comparisons = compare_CMV, method = "t.test")
# 
# # polish plot for publication, v1
# inputDat %>% 
#   ggboxplot(x = "CMV_IgG_Serology", y = rna,
#             color = "CMV_IgG_Serology", palette = c("#0047ab", "#00743f"),
#             add = "jitter", add.params = list(size = 1, alpha = 0.5)) + 
#   # geom_point(alpha = 0.5) + 
#   facet_wrap(~Cohort, nrow = 1) +
#   stat_compare_means(comparisons = compare_CMV, method = "t.test") + 
#   theme_classic(base_size = 12) +
#   theme(plot.title = element_text(face = "bold", size = 14),
#         axis.title = element_text(face = "bold"),
#         legend.position = "none")
# 
# # plot all rnas
# plotList <- list()
# for (rna in selected_vals) {
#   plotList[[rna]] <- inputDat %>% 
#     ggboxplot(x = "CMV_IgG_Serology", y = rna,
#               paletter = "jco", add = "jitter") + facet_wrap(~Cohort, nrow = 1) +
#     stat_compare_means(comparisons = compare_CMV, method = "t.test")
# }
