rm(list = ls())

library(tidyverse)

# load data for sig. cytokine =======================================================================
load("processedDat/rlmRes_cytokine.RData") # to get the sig. cytokines

sigCytos_discovery <- rownames(rlm_res$discovery %>% filter(padj < 0.05))
sigCytos_validation <- rownames(rlm_res$validation %>% filter(pval < 0.05))

#rm(rlm_res, cytokine_NAs)
rm(rlm_res)

# load data for heatmap s=======================================================================
load("processedDat/rlmRes_cytokine_forHeatmap.RData") # to get the effect size for all cytokine production

# get t-statistic ---------------------------------
tstat_temp <- rlm_res %>% 
  lapply(function(x) x %>% rownames_to_column("cyto")) %>% 
  bind_rows(.id = "cohort")


cytos_bothCohort <- tstat_temp$cyto[duplicated(tstat_temp$cyto)] #cytokine have effect size in both cohorts 
selected_cytos <-  cytos_bothCohort[-grep("rpmi", cytos_bothCohort)] # remove the cytokine measure at baseline

tstat <- tstat_temp %>% filter(cyto %in% selected_cytos )

tstat_pVal <- tstat %>% 
  mutate(pSig = ifelse(cohort == "discovery" & cyto %in% sigCytos_discovery, "yes",
                       ifelse(cohort == "validation" & cyto %in% sigCytos_validation, "yes", NA)))

# heatmap  ----------------------------------------------------
tstat_pVal %>%
  ggplot(aes(x = cohort, y = cyto, fill = effectSize)) + 
  geom_tile() +
  geom_text(aes(label = ifelse(pSig == "yes", "*", NA))) +
  scale_fill_gradient2(low = "#3C5488FF", mid = "white", high = "#DC0000FF") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1))

# plot heatmap with order
plotDat_wide <- tstat_pVal %>% dplyr::select(cyto, cohort, effectSize) %>%
  pivot_wider(names_from = "cyto", values_from =  "effectSize") %>% 
  column_to_rownames("cohort") %>% t()
plotDat_dendrogram <- as.dendrogram(hclust(d = dist(plotDat_wide)))
plotDat_denOrder <- order.dendrogram(plotDat_dendrogram)

plotDat_order <-  tstat_pVal%>%
  mutate(cyto = factor(cyto, levels = cyto[plotDat_denOrder], ordered = TRUE))

plotDat_order %>%
  ggplot(aes(x = cohort, y = cyto, fill = effectSize)) + 
  geom_tile() +
  geom_text(aes(label = ifelse(pSig == "yes", "*", NA))) +
  scale_fill_gradient2(low = "#8491B4FF", mid = "white", high = "#DC0000FF") + 
  theme_bw() + 
  theme(text=element_text(size=16))

heatmap_plot <- plotDat_order %>%
  ggplot(aes(x = cohort, y = cyto, fill = effectSize)) + 
  geom_tile() +
  geom_text(aes(label = ifelse(pSig == "yes", "*", NA))) +
  scale_fill_gradientn(
    limits = c(-2, 2),
    colours = c("#3C5488FF", "white", "#DC0000FF")) + 
  theme_bw() + 
  theme(text=element_text(size=16))

heatmap_plot

# save the plot ------------------------------------------------------------------
png("output/07_01_cytokine_heatmap.png",height = 1008)
heatmap_plot
dev.off()
