rm(list = ls())

library(tidyverse)
library(MASS)
library(lmtest)
library(sandwich)
library("ggrepel")

# inverse ranking normaliation/ transformation (or can use the RNOmni package) 
inormal <- function(x) qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))

# load data =======================================================================
load("dat/cohortDat.RData")


## prepare down sampling dataset for CMV+ ---------------------------------------
cohorts <- c("Discovery", "Validation")

# downSampling <- list()
# for (cohort in cohorts) {
#   for (idx in c(1:100)) {
#     downSampling[[cohort]][[idx]] <- read.csv(
#       paste0("CMV_downsample_2025-03/", tolower(cohort), "/raw/seed_", idx, ".csv"),
#       header = TRUE)
#   }
# }
# 
# save(downSampling, file = "dat/downSampling.RData")
load("dat/downSampling.RData")

## prepare down sampling dataset for CMV- ---------------------------------------
samples_CMVneg <- list()

for (cohort in cohorts) {
  samples_CMVneg[[cohort]] <- (cohortDat$donor_info %>% 
             filter(CMV_IgG_Serology == 0, Cohort == cohort))$Record.Id
}

# run the rlm() model ------------------------------------------------
proteins <- names(cohortDat$allSample$protein)

## discovery dataset ------------------------------------------------
top5_geneticPCs_discovery <- read.table(
  "Discovery_allEthnicity_2000HIV.eigenvec", header = TRUE) %>% 
  as_tibble() %>% dplyr::select(IID, PC1, PC2, PC3, PC4, PC5)

inputDat_discovery <- cohortDat$donor_info %>%
  inner_join(top5_geneticPCs_discovery, by = c("Record.Id" = "IID")) %>%
  inner_join(cohortDat$Discovery$protein %>% rownames_to_column("Record.Id")) %>%
  mutate(Institute.Abbreviation = as.factor(Institute.Abbreviation),
         SEX_BIRTH = as.factor(SEX_BIRTH), 
         CMV_IgG_Serology = as.factor(CMV_IgG_Serology))

rlm_res_discovery <- list()
for (idx in c(1:100)) {
  rlm_res_temp <- list()
  for (protein in proteins) {
    dat_temp <- inputDat_discovery %>% 
      filter(Record.Id %in% c(colnames(downSampling$Discovery[[idx]]), 
                              samples_CMVneg$Discovery)) %>%
      dplyr::rename("valName" = protein) %>% 
      mutate(valName = inormal(valName)) # inverse ranking transformation to have normal distribution
    
    ttest <- rlm(valName ~ CMV_IgG_Serology + AGE + SEX_BIRTH + BMI_BASELINE + 
                   Institute.Abbreviation + season_sin + season_cos +
                   PC1 + PC2 + PC3 + PC4 +PC5, 
                 data = dat_temp, maxit=200)
    cf <- try(coeftest(ttest, vcov=vcovHC(ttest, type="HC0")))
    rlm_res_temp[[protein]] <- cf["CMV_IgG_Serology1", c("Estimate", "Std. Error", "z value", "Pr(>|z|)")]
  }
  rlm_res_discovery[[idx]] <-  rlm_res_temp %>% as_tibble() %>%
    t() %>%  as.data.frame() %>%
    dplyr::rename("effectSize" = 1, "StdError" = 2, "Z_Value" = 3, "pval" = 4) %>%
    mutate(padj = p.adjust(pval, method = "fdr"))
  
  print(paste0("Finished seed ", idx))
}

save(rlm_res_discovery, file = "output/rlm_proteinDiscovery.RData")

# compare the downsampling with original analysis outcomes in discovery cohort ----------------
## padj in reanalysis is calculated in protein-wide data (all proteins) ----------------
# rm(list = ls())
# 
# library(tidyverse)
# library(gridExtra)
# library(ggpubr)
# library(ggrepel)
# 
# load("rlmRes_protein.RData")
# cohort <- "discovery"
# 
# inputDat <- list()
# inputDat$original_discovery <- rlm_res[[cohort]]
# 
# load("output/rlm_proteinDiscovery.RData")
# 
# plotDat <- list()
# plots_list <- list()
# for (idx in c(1:100)) {
#   inputDat$downsampled_discovery <- rlm_res_discovery[[idx]]
#   
#   label_DEPs <- c("FCRL6") 
#   
#   # prepare the plot data
#   plotDat[[idx]] <- inputDat %>% 
#     lapply(function(x) x %>% as.data.frame %>% rownames_to_column("gene")) %>% 
#     bind_rows(.id = "cohort") %>% 
#     dplyr::select(-StdError, -Z_Value) %>% 
#     pivot_wider(names_from = "cohort", values_from = c("effectSize", "pval", "padj")) %>%
#     filter(padj_original_discovery < 0.05) %>%
#     mutate(sig_bothCohorts = ifelse(padj_original_discovery < 0.05 & padj_downsampled_discovery < 0.05, "yes", "no")) %>%
#     mutate(validated_DEPs = 
#              ifelse(sig_bothCohorts == "yes" & effectSize_original_discovery * effectSize_downsampled_discovery  > 0, 
#                     "yes", "no")) %>% 
#     mutate(label_DEPs = ifelse(gene %in% label_DEPs & validated_DEPs == "yes", gene, NA))
#   
#   validated_percentage <- (plotDat[[idx]] %>% 
#                              group_by(validated_DEPs)  %>% summarise(n = n()) %>% mutate(per= n/sum(n)) %>%
#                              dplyr::select(validated_DEPs, per) %>% 
#                              pivot_wider(names_from = "validated_DEPs", values_from = "per"))$yes %>% round(., digits = 4) * 100 
#   
#   # create plot 
#   comparePlot_protein <- ggplot(data = plotDat[[idx]], 
#                                aes(x = effectSize_original_discovery, 
#                                    y = effectSize_downsampled_discovery,
#                                    color = validated_DEPs, label = label_DEPs)) + 
#     geom_point(size = 4, alpha = 0.5) + 
#     scale_color_manual(name = "validated_DEPs",
#                        values=c("no" = "#d3d3d3", "yes" = "#bb6f7c"),
#                        labels = c("no", "yes")) +
#     geom_vline(xintercept = 0) + geom_hline(yintercept = 0)+
#     geom_label_repel(size = 4,
#                      max.overlaps = getOption("ggrepel.max.overlaps", default = 40),
#                      box.padding = unit(0.6, "lines"), show.legend = FALSE) +
#     ggtitle(label = paste0('Downsampled_', idx, ", ", validated_percentage, "%")) +
#     theme_classic(base_size = 12) +
#     guides(colour = guide_legend(position = "inside")) +
#     theme(legend.position.inside = c(0.8, 0.1),
#           axis.title.x = element_blank(), axis.title.y = element_blank())
#   
#   plots_list[[idx]] <- comparePlot_protein 
# }
# 
# save(plotDat, plots_list, file = "output/rlmProteinRes_plotDat.RData")
# 
# # check the validated DEPs percentage
# plotDat[[2]] %>% group_by(sig_bothCohorts)  %>% summarise(n = n()) %>% mutate(per= n/sum(n))
# 
# validated_DEPs <- plotDat %>% 
#   lapply(function(x) x %>% 
#            group_by(validated_DEPs)  %>% summarise(n = n()) %>% mutate(per= n/sum(n)) %>%
#            dplyr::select(validated_DEPs, per) %>% 
#            pivot_wider(names_from = "validated_DEPs", values_from = "per")) %>% 
#   bind_rows(.id = "downsampled_seed")
# 
# range(validated_DEPs$yes) # 0.1578947 0.7894737
# hist(validated_DEPs$yes)
# summary(validated_DEPs$yes) # Mean = 0.3487
# 
# # create plot group and save
# #cowplot::plot_grid(plots_list[c(1:50)], nrow = 5, ncol = 10)
# 
# plotSW_1 <- do.call(grid.arrange, c(plots_list[c(1:50)], ncol = 10))
# plotSW_2 <- do.call(grid.arrange, c(plots_list[c(51:100)], ncol = 10))
# 
# png("reAnalysis_outcome/protein/10_revision_proteinPart1_discoveryCohort.png", width = 2400, height = 960)
# plot(plotSW_1)
# dev.off()
# 
# png("reAnalysis_outcome/protein/10_revision_proteinPart2_discoveryCohort.png", width = 2400, height = 960)
# plot(plotSW_2)
# dev.off()

## padj in reanalysis is calculated in DEPs from original analysis (38 proteins) ----------------
rm(list = ls())

library(tidyverse)
library(gridExtra)
library(ggpubr)
library(ggrepel)

load("dat/rlmRes_protein.RData")
cohort <- "discovery"

inputDat <- list()
inputDat$original_discovery <- rlm_res[[cohort]]

load("output/rlm_proteinDiscovery.RData")

plotDat <- list()
plots_list <- list()
for (idx in c(1:100)) {
  inputDat$downsampled_discovery <- rlm_res_discovery[[idx]]
  
  label_DEPs <- c("FCRL6") 
  
  # prepare the plot data
  plotDat[[idx]] <- inputDat %>% 
    lapply(function(x) x %>% as.data.frame %>% rownames_to_column("gene")) %>% 
    bind_rows(.id = "cohort") %>% 
    dplyr::select(-StdError, -Z_Value) %>% 
    pivot_wider(names_from = "cohort", values_from = c("effectSize", "pval", "padj")) %>%
    filter(padj_original_discovery < 0.05) %>%
    mutate(padj_downsampled_discovery_v2 = p.adjust(pval_downsampled_discovery, method = "fdr")) %>% # calculate padj in downsampled data based on 38 DEPs in original analysis
    mutate(sig_bothCohorts = ifelse(padj_original_discovery < 0.05 & padj_downsampled_discovery_v2 < 0.05, "yes", "no")) %>%
    mutate(validated_DEPs = 
             ifelse(sig_bothCohorts == "yes" & effectSize_original_discovery * effectSize_downsampled_discovery  > 0, 
                    "yes", "no")) %>% 
    mutate(label_DEPs = ifelse(gene %in% label_DEPs & validated_DEPs == "yes", gene, NA)) %>%
    mutate(validated_DEPs_v2 = ifelse(is.na(label_DEPs), validated_DEPs, "yes_highlight"))
  
  validated_percentage <- (plotDat[[idx]] %>% 
                             group_by(validated_DEPs)  %>% summarise(n = n()) %>% mutate(per= n/sum(n)) %>%
                             dplyr::select(validated_DEPs, per) %>% 
                             pivot_wider(names_from = "validated_DEPs", values_from = "per"))$yes %>% round(., digits = 4) * 100 
  
  # create plot 
  comparePlot_protein <- ggplot(data = plotDat[[idx]]  %>% arrange(validated_DEPs_v2), 
                               aes(x = effectSize_original_discovery, 
                                   y = effectSize_downsampled_discovery,
                                   color = validated_DEPs_v2, label = label_DEPs)) + 
    geom_point(size = 4, alpha = 0.5) + 
    scale_color_manual(name = "validated_DEPs",
                       values=c("no" = "#d3d3d3", "yes" = "#bb6f7c", "yes_highlight" = "#660000"),
                       labels = c("no", "yes", "FCRL6")) +
    geom_vline(xintercept = 0) + geom_hline(yintercept = 0)+
    geom_label_repel(size = 4,
                     max.overlaps = getOption("ggrepel.max.overlaps", default = 40),
                     box.padding = unit(0.6, "lines"), show.legend = FALSE) +
    ggtitle(label = paste0('Downsampled_', idx, ", ", validated_percentage, "%")) +
    theme_classic(base_size = 12) +
    guides(colour = guide_legend(position = "inside")) +
    theme(legend.position.inside = c(0.8, 0.1),
          axis.title.x = element_blank(), axis.title.y = element_blank())
  
  plots_list[[idx]] <- comparePlot_protein 
}

save(plotDat, plots_list, file = "output/rlmProteinRes_plotDat_v2.RData")

# check the validated DEPs percentage
plotDat[[2]] %>% group_by(sig_bothCohorts)  %>% summarise(n = n()) %>% mutate(per= n/sum(n))

validated_DEPs <- plotDat %>% 
  lapply(function(x) x %>% 
           group_by(validated_DEPs)  %>% summarise(n = n()) %>% mutate(per= n/sum(n)) %>%
           dplyr::select(validated_DEPs, per) %>% 
           pivot_wider(names_from = "validated_DEPs", values_from = "per")) %>% 
  bind_rows(.id = "downsampled_seed")

range(validated_DEPs$yes) # 0.4473684 1.0000000
hist(validated_DEPs$yes) 
summary(validated_DEPs$yes) #  Mean = 0.8013
# create plot group and save
#cowplot::plot_grid(plots_list[c(1:50)], nrow = 5, ncol = 10)

# plotSW_1 <- do.call(grid.arrange, c(plots_list[c(1:50)], ncol = 10))
# plotSW_2 <- do.call(grid.arrange, c(plots_list[c(51:100)], ncol = 10))
# 
# png("reAnalysis_outcome/protein/10_revision_proteinPart1_discoveryCohort_v2.png", width = 2400, height = 960)
# plot(plotSW_1)
# dev.off()
# 
# png("reAnalysis_outcome/protein/10_revision_proteinPart2_discoveryCohort_v2.png", width = 2400, height = 960)
# plot(plotSW_2)
# dev.off()

library(patchwork)

plotSW_1 <- wrap_plots(plots_list[1:50], ncol = 10) +
  plot_layout(guides = "collect") &
  theme(legend.position = "none")

plotSW_2 <- wrap_plots(plots_list[51:100], ncol = 10) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom") &
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 5),
                               theme = theme(legend.text = element_text(size=18)),
                               byrow = TRUE))

png("reAnalysis_outcome/protein/10_revision_proteinPart1_discoveryCohort_v3.png", width = 2400, height = 960)
plot(plotSW_1)
dev.off()

png("reAnalysis_outcome/protein/10_revision_proteinPart2_discoveryCohort_v3.png", width = 2400, height = 980)
plot(plotSW_2)
dev.off()
