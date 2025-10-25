rm(list = ls())

library(tidyverse)
library(venn)
library(gridExtra)
library(ggpubr)
library(ggrepel)
# load data =======================================================================
load("processedDat/rlmRes_cellcount.RData")

# check sig. cellcounts ---------------------------------
discovery_padj <- rlm_res$discovery %>% filter(padj < 0.05)
validation_pval <- rlm_res$validation %>% filter(pval < 0.05)

valNames <- list("discovery_padj" = rownames(discovery_padj),
                 "validation_pval" = rownames(validation_pval))
venn(valNames)

selected_vals <- intersect(rownames(discovery_padj), rownames(validation_pval))

# ploting to check data, these plots are not in the manuscript ----------------------------------
# ## volcano plot ----------------------------------
# outcome <- rlm_res$discovery
# 
# plotDat <- outcome %>%
#   rownames_to_column("valName") %>%
#   mutate(sig = ifelse(padj < 0.05, "padj<0.05", "Not Sig"),
#          #delabel = ifelse(sig == "Not Sig", NA, valName)
#          delabel = ifelse(valName %in% selected_vals, valName, NA))
# 
# 
# plotDat %>%
#   ggplot(aes(effectSize, -log10(padj), col = sig, label = delabel)) +
#   geom_point() +
#   scale_color_manual(values=c("black", "red")) +
#   ggtitle("Differential association between CMV+ vs. CMV-") +
#   geom_hline(yintercept=-log10(0.05), col="red") +
#   geom_text_repel(max.overlaps = getOption("ggrepel.max.overlaps", default = 30)) + theme_bw()
# 
# 
# ## boxplot ----------------------------------
# load("processedDat/cohortDat.RData")
# 
# inputDat <- cohortDat$donor_info %>% 
#   right_join(cohortDat$allSample$cellcount %>% rownames_to_column("Record.Id")) %>%
#   mutate(CMV_IgG_Serology = ifelse(CMV_IgG_Serology == 1, "CMV+",
#                                    ifelse(CMV_IgG_Serology == 0, "CMV-", CMV_IgG_Serology))) %>%
#   drop_na(CMV_IgG_Serology)
# 
# # set up comparison
# compare_CMV <- list( c("CMV-", "CMV+"))
# 
# # individual boxplot
# cellcount <- "Panel1_1045_TCRyd"
# cellcount <- "Panel1_1050_NK-T.like" 
# cellcount <- "Panel2_2289_NKcells_HLA-DR+"
# 
# inputDat %>% 
#   ggboxplot(x = "CMV_IgG_Serology", y = cellcount,
#             paletter = "jco", add = "jitter") + facet_wrap(~Cohort, nrow = 1) +
#   stat_compare_means(comparisons = compare_CMV, method = "t.test")
# 
# # plot all cellcounts
# plotList <- list()
# for (cellcount in selected_vals) {
#   plotList[[cellcount]] <- inputDat %>% 
#     ggboxplot(x = "CMV_IgG_Serology", y = cellcount,
#               paletter = "jco", add = "jitter") + facet_wrap(~Cohort, nrow = 1) +
#     stat_compare_means(comparisons = compare_CMV, method = "t.test")
# }


