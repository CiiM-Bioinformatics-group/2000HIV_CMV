rm(list = ls())

library(tidyverse)
library(venn)
library(gridExtra)

library(ggpubr)

get.log2 <- function(dat) {
  library(tidyverse)
  library(rstatix)
  outcome <- dat %>% mutate(across(where(is.numeric), ~log2(.x))) %>%
    mutate(across(where(is.numeric), ~ifelse(is.infinite(.x), 0, .x)))
  return(outcome)
}

# load data =======================================================================
load("processedDat/rlmRes_mebo.RData")

# check sig. metabolites ---------------------------------
discovery_padj <- rlm_res$discovery %>% filter(padj < 0.05)
validation_pvalue <- rlm_res$validation %>% filter(pval < 0.05)

valNames <- list("discovery_padj" = rownames(discovery_padj),
                 "validation_pvalue" = rownames(validation_pvalue))
venn(valNames)

selected_vals <- intersect(rownames(discovery_padj), rownames(validation_pvalue))

# ploting to check data, these plots are not in the manuscript ----------------------------------
# ## volcano plot ----------------------------------
# outcome <- rlm_res$discovery
# 
# plotDat <- outcome %>%
#   rownames_to_column("valName") %>%
#   mutate(sig = ifelse(padj < 0.05, "padj<0.05", "Not Sig"),
#          delabel = ifelse(sig == "Not Sig", NA, valName))
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
#   right_join(cohortDat$allSample$mebo %>% rownames_to_column("Record.Id")) %>%
#   mutate(CMV_IgG_Serology = ifelse(CMV_IgG_Serology == 1, "CMV+",
#                                    ifelse(CMV_IgG_Serology == 0, "CMV-", CMV_IgG_Serology))) %>%
#   drop_na(CMV_IgG_Serology)
# 
# inputDat_log2 <- cohortDat$donor_info %>% 
#   right_join(cohortDat$allSample$mebo %>% mutate_all(~log2(.x)) %>% rownames_to_column("Record.Id")) %>%
#   mutate(CMV_IgG_Serology = ifelse(CMV_IgG_Serology == 1, "CMV+",
#                                    ifelse(CMV_IgG_Serology == 0, "CMV-", CMV_IgG_Serology))) %>%
#   drop_na(CMV_IgG_Serology)
# 
# # set up comparison
# compare_CMV <- list( c("CMV-", "CMV+"))
# 
# # individual boxplot
# mebo <- "C5H5NO"
# mebo <- "C9H9NO2"
# mebo <- "C7H8O5S"
# mebo <- "C8H10O5S"
# mebo <- "C12H14O4"
# 
# inputDat %>% 
#   ggboxplot(x = "CMV_IgG_Serology", y = mebo,
#             paletter = "jco", add = "jitter") + facet_wrap(~Cohort, nrow = 1) +
#   stat_compare_means(comparisons = compare_CMV, method = "t.test")
# 
# # plot all mebos
# plotList <- list()
# for (mebo in selected_vals) {
#   plotList[[mebo]] <- inputDat %>% 
#     ggboxplot(x = "CMV_IgG_Serology", y = mebo,
#               paletter = "jco", add = "jitter") + facet_wrap(~Cohort, nrow = 1) +
#     stat_compare_means(comparisons = compare_CMV, method = "t.test")
# }
# 
