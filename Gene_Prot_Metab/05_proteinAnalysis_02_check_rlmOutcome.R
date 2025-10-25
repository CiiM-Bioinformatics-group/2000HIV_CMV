rm(list = ls())

library(tidyverse)
library(venn)
library(gridExtra)
library(ggpubr)
library(ggrepel)
# load data =======================================================================
load("processedDat/rlmRes_protein.RData")

# check sig. proteins ---------------------------------
discovery_padj <- rlm_res$discovery %>% filter(padj < 0.05)
validation_pvalue <- rlm_res$validation %>% filter(pval < 0.05)

valNames <- list("discovery_padj" = rownames(discovery_padj),
                 "validation_pvalue" = rownames(validation_pvalue))
venn(valNames)

selected_vals <- intersect(rownames(discovery_padj), rownames(validation_pvalue))

# gene name and effected size ---------------------------------------------
protein_effects <- discovery_padj %>% 
  select(effectSize) %>% rownames_to_column("gene") %>% 
  filter(gene %in% selected_vals)

write.table(protein_effects, file = "processedDat/sigProteins_effect.txt", row.names = FALSE, quote = FALSE)

# ploting to check data, these plots are not in the manuscript ----------------------------------
# ## compare plot ----------------------------------
# plotDat <- rlm_res %>% 
#   lapply(function(x) x %>% as.data.frame %>% rownames_to_column("gene")) %>% 
#   bind_rows(.id = "cohort") %>% 
#   dplyr::select(-StdError, -Z_Value) %>% 
#   pivot_wider(names_from = "cohort", values_from = c("effectSize", "pval", "padj")) %>%
#   mutate(sig_bothCohorts = ifelse(padj_discovery < 0.05 & pval_validation < 0.05, "yes", "no")) %>%
#   mutate(delabel = ifelse(sig_bothCohorts == "yes", gene, NA)) %>%
#   mutate(validated_DEPs = 
#            ifelse(sig_bothCohorts == "yes" & effectSize_discovery * effectSize_validation  > 0, 
#                   "yes", "no")) %>%
#   filter(padj_discovery < 0.05)
# 
# ggplot(data = plotDat, 
#        aes(x = effectSize_discovery, y = effectSize_validation, 
#            color = validated_DEPs, label = delabel)) + 
#   geom_point(alpha = 0.8) + 
#   scale_color_manual(values=c("grey","red")) +
#   geom_point(data = plotDat %>% filter(validated_DEPs == "yes"),
#              aes(x = effectSize_discovery, y = effectSize_validation), 
#              alpha = 0.8, colour = "red") + 
#   geom_vline(xintercept = 0) + geom_hline(yintercept = 0)+
#   geom_text_repel(max.overlaps = getOption("ggrepel.max.overlaps", default = 30)) +
#   theme_classic()
# 
# plotDat %>% group_by(validated_DEPs) %>% summarise(countNum = n()) %>%
#   mutate(groupPercent = (countNum/sum(countNum))*100)
# 
# ## volcano plot ----------------------------------
# outcome <- rlm_res$discovery
# 
# plotDat <- outcome %>%
#   rownames_to_column("valName") %>%
#   mutate(sig = ifelse(padj < 0.05, "padj<0.05", "Not Sig"),
#          #delabel = ifelse(sig == "Not Sig", NA, valName)
#          #delabel = ifelse(valName %in% selected_vals, valName, NA)
#          delabel = ifelse(valName == "FCRL6", valName, NA))
# 
# plotDat %>%
#   ggplot(aes(effectSize, -log10(padj), col = sig, label = delabel)) +
#   geom_point() +
#   scale_color_manual(values=c("black", "red")) +
#   ggtitle("Differential association between CMV+ vs. CMV-") +
#   geom_hline(yintercept=-log10(0.05), col="red") +
#   geom_text_repel(max.overlaps = getOption("ggrepel.max.overlaps", default = 30)) + theme_bw()
# 
# ## boxplot ----------------------------------
# load("processedDat/cohortDat.RData")
# 
# inputDat <- cohortDat$donor_info %>% 
#   right_join(cohortDat$allSample$protein %>% rownames_to_column("Record.Id")) %>%
#   mutate(CMV_IgG_Serology = ifelse(CMV_IgG_Serology == 1, "CMV+",
#                                    ifelse(CMV_IgG_Serology == 0, "CMV-", CMV_IgG_Serology))) %>%
#   drop_na(CMV_IgG_Serology)
# 
# # set up comparison
# compare_CMV <- list( c("CMV-", "CMV+"))
# 
# # individual boxplot
# protein <- "FCRL6"
# protein <- "ITGAL"
# protein <- "DNAJA4"
# protein <- "GZMH"
# protein <- "CRTAM"
# protein <- "GNLY"
# protein <- "GBP1"
# protein <- "KLRD1"
# protein <- "ADGRG1"
# protein <- "KIR2DL3"
# 
# 
# inputDat %>% 
#   ggboxplot(x = "CMV_IgG_Serology", y = protein,
#             paletter = "jco", add = "jitter") + facet_wrap(~Cohort, nrow = 1) +
#   stat_compare_means(comparisons = compare_CMV, method = "t.test")
# 
# # polish plot for publication, v1
# inputDat %>% 
#   ggboxplot(x = "CMV_IgG_Serology", y = protein,
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
# # plot all proteins
# plotList <- list()
# for (protein in selected_vals) {
#   plotList[[protein]] <- inputDat %>% 
#     ggboxplot(x = "CMV_IgG_Serology", y = protein,
#               paletter = "jco", add = "jitter") + facet_wrap(~Cohort, nrow = 1) +
#     stat_compare_means(comparisons = compare_CMV, method = "t.test")
# }
# 
# pdf(file = "processedDat/proPadj_boxplot.pdf", height = 40, width = 15, onefile = TRUE)
# #do.call("grid.arrange", c(plotList, nrow = 3, ncol = 3))
# cowplot::plot_grid(plotlist = plotList, ncol = 3)
# dev.off()
# 
# ## boxplot with padj value from rlm() model -----------------------
# # make the statistic table of rlm() model to add into boxplot
# y.position <- inputDat %>% select(selected_vals) %>% 
#   summarise_if(is.numeric, max, na.rm = TRUE) %>% 
#   t() %>% as.data.frame %>% rownames_to_column(".y.") %>% rename("y.position" = V1)
# 
# stat.test <- discovery_padj %>% 
#   mutate(Cohort = "Discovery") %>% rownames_to_column(".y.") %>%
#   full_join(validation_padj %>% 
#               mutate(Cohort = "Validation") %>% rownames_to_column(".y.")) %>%
#   mutate(group1 = "CMV+", group2 = "CMV-", padj = round(padj, digits = 20)) %>%
#   full_join(y.position) %>%
#   rename("p.adj" = padj)  %>% as.data.frame()
# 
# # boxplot with rlm() padj value
# inputDat %>% 
#   ggboxplot(x = "CMV_IgG_Serology", y = protein,
#             paletter = "jco", add = "jitter") + facet_wrap(~Cohort, nrow = 1) +
#   stat_pvalue_manual(stat.test %>% filter(.y. == protein), label = "p.adj")
# 
# # plot all proteins with rlm() padj value
# plotList <- list()
# for (protein in selected_vals) {
#   plotList[[protein]] <- inputDat %>% 
#     ggboxplot(x = "CMV_IgG_Serology", y = protein,
#               paletter = "jco", add = "jitter") + facet_wrap(~Cohort, nrow = 1) +
#     stat_pvalue_manual(stat.test %>% filter(.y. == protein), label = "p.adj")
# }
# 
# pdf(file = "processedDat/proPadj_boxplot_rlm.pdf", height = 30, width = 15, onefile = TRUE)
# #do.call("grid.arrange", c(plotList, nrow = 3, ncol = 3))
# cowplot::plot_grid(plotlist = plotList, ncol = 4)
# dev.off()
# 
