# I modify the 04_rnaAnalysis_01_Deseq2.R code 
rm(list = ls())

library(tidyverse)
library(DESeq2)
library(org.Hs.eg.db)

#library("BiocParallel")
#register(MulticoreParam(10))

# load data =======================================================================
load("dat/cohortDat.RData")

# # RNAseq raw data ------------------------------
# raw_dat <- readRDS("2000HIV_bulk_transcriptomics_raw_counts.RDS")
# sampleIDs <- readRDS("2000HIV_bulk_transcriptomics_sample_table.RDS")
# 
# ## convert Ensembl IDs into gene symbols
# annots <- select(org.Hs.eg.db,
#                  keys= substring(rownames(raw_dat), 1, 15),
#                  columns="SYMBOL", keytype="ENSEMBL")
# 
# ## gene expression matrix
# RNAseqDat <- raw_dat %>% t() %>%
#   as.data.frame %>% rownames_to_column("ID") %>% # convert sample id to donor id
#   full_join(sampleIDs %>% dplyr::select(ID, DONOR_ID)) %>% dplyr::select(-ID) %>%
#   column_to_rownames("DONOR_ID") %>% t() %>% as.data.frame %>%
#   rownames_to_column("ENSEMBL") %>% # convert ensembl IDs into gene symbols
#   mutate(ENSEMBL = substring(ENSEMBL, 1, 15)) %>%  full_join(annots) %>%
#   drop_na(SYMBOL) %>% dplyr::select(-ENSEMBL) %>%
#   group_by(SYMBOL) %>% summarise_each(funs(mean)) %>% # calculate the average read count for duplicated gene symbols
#   column_to_rownames("SYMBOL")
# 
# save(RNAseqDat, file = "RNAseqDat.RData") 
# Note: this time get 36494 genes, but last time (original analysis) get 36000 genes. 
# Some genes in the original analysis with 36000 genes, such as A2MP1, are not in the new gene conversion.

load("dat/RNAseqDat.RData")

which(duplicated(RNAseqDat$SYMBOL) == TRUE)
unique(duplicated(RNAseqDat$SYMBOL)) # no duplication

# prepare dataset -------------------------------------------------------------

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
  temp <- (cohortDat$donor_info %>% 
             filter(CMV_IgG_Serology == 0, Cohort == cohort))$Record.Id
  samples_CMVneg[[cohort]] <- intersect(temp, colnames(RNAseqDat))
}

## top 5 genetic PCs ----------------------------------------------------------
geneticPCs <- list()

for (cohort in cohorts) {
  geneticPCs[[cohort]] <- read.table(
    paste0(cohort, "_allEthnicity_2000HIV.eigenvec"),
    header = TRUE) 
}

top5_geneticPCs <- geneticPCs %>% 
  lapply(function(x) x %>% as_tibble() %>% dplyr::select(IID, PC1, PC2, PC3, PC4, PC5))

## prepare input data files ----------------------------------------------------------
countDat <- list()
sampleInfo <- list()

for (idx in c(1:100)) {
  metadata <- list()
  overlapSamples <- list()
  
  for (cohort in cohorts) {
    metadata[[cohort]] <- cohortDat$donor_info %>% 
      dplyr::select(Record.Id, CMV_IgG_Serology, AGE, SEX_BIRTH, BMI_BASELINE, 
                    Institute.Abbreviation, season_sin, season_cos) %>% 
      inner_join(top5_geneticPCs[[cohort]], by = c("Record.Id" = "IID")) %>%
      mutate(Institute.Abbreviation = as.factor(Institute.Abbreviation),
             SEX_BIRTH = as.factor(SEX_BIRTH), 
             CMV_IgG_Serology = as.factor(CMV_IgG_Serology)) %>% 
      drop_na(CMV_IgG_Serology)
    
    overlapSamples[[cohort]] <- intersect(metadata[[cohort]]$Record.Id, 
                                          c(colnames(downSampling[[cohort]][[idx]]), 
                                            samples_CMVneg[[cohort]]))
    
    ## prepare input data
    countDat[[cohort]][[idx]] <- RNAseqDat[, overlapSamples[[cohort]]]
    sampleInfo[[cohort]][[idx]] <- (metadata[[cohort]] %>% column_to_rownames("Record.Id"))[overlapSamples[[cohort]], ]
    
  }
}

# run DEseq2 -------------------------------------------------------------------------
identical(rownames(sampleInfo$Discovery[[1]]), colnames(countDat$Discovery[[1]])) # TRUE, correct sample order -> can run DESeq2 now
identical(rownames(sampleInfo$Validation[[1]]), colnames(countDat$Validation[[1]])) # TRUE, correct sample order -> can run DESeq2 now

rm(cohortDat, downSampling, geneticPCs, metadata, overlapSamples, RNAseqDat, samples_CMVneg, top5_geneticPCs)

for (idx in c(1:100)) {
  DESeq2_res <- list()
  DESeq2_resLFC <- list()
  
  ## discovery cohort ------------------------------------------
  dds <- DESeqDataSetFromMatrix(
    countData = round(countDat$Discovery[[idx]]), colData = sampleInfo$Discovery[[idx]],
    design = ~ AGE + SEX_BIRTH + BMI_BASELINE + Institute.Abbreviation + 
      season_sin + season_cos + PC1 + PC2 + PC3 + PC4 +PC5 + 
      CMV_IgG_Serology) # compare CMV status
  
  dds_v2 <- DESeq(dds)
  
  DESeq2_res$discovery <- results(dds_v2, contrast=c("CMV_IgG_Serology", "1","0"), alpha=0.05)
  # Shrinkage of effect size (LFC estimates), using the apeglm method (Zhu, Ibrahim, and Love 2018), which improves on the previous estimator.
  DESeq2_resLFC$discovery <- lfcShrink(dds_v2, coef="CMV_IgG_Serology_1_vs_0", type="apeglm") 
  
  ## validation cohort ------------------------------------------
  sigGene_discovery <- rownames(DESeq2_res$discovery %>% 
                                  as.data.frame %>% filter(padj < 0.05))
  countDat$Validation_sigGene[[idx]] <- countDat$Validation[[idx]][sigGene_discovery,]
  
  dds <- DESeqDataSetFromMatrix(
    # countData = round(countDat$Validation[[idx]]), colData = sampleInfo$Validation[[idx]],
    countData = round(countDat$Validation_sigGene[[idx]]), colData = sampleInfo$Validation[[idx]],
    design = ~ AGE + season_sin + season_cos + 
      #SEX_BIRTH + BMI_BASELINE + PC1 + PC2 + PC3 + PC4 +PC5 + 
      CMV_IgG_Serology) # compare CMV status
  
  dds_v2 <- DESeq(dds)
  
  DESeq2_res$validation <- results(dds_v2, contrast=c("CMV_IgG_Serology", "1","0"), alpha=0.05)
  # Shrinkage of effect size (LFC estimates), using the apeglm method (Zhu, Ibrahim, and Love 2018), which improves on the previous estimator.
  DESeq2_resLFC$validation <- lfcShrink(dds_v2, coef="CMV_IgG_Serology_1_vs_0", type="apeglm") 
  
  ## save data ------------------------------------------------
  save(DESeq2_res, DESeq2_resLFC, 
       file = paste0("output/DEseq2Res_seed_", idx, ".RData"))
  
}

# compare the downsampling with original analysis outcomes in discovery cohort ----------------
## padj in reanalysis is calculated in genome-wide data (all genes) ----------------
# rm(list = ls())
# 
# library(tidyverse)
# library(gridExtra)
# library(ggpubr)
# library(ggrepel)
# 
# load("dat/DEseq2Res_rna.RData")
# cohort <- "discovery"
# 
# inputDat <- list()
# inputDat$original_discovery <- DESeq2_res[[cohort]]
# rm(DESeq2_res, DESeq2_resLFC)
# 
# plotDat <- list()
# plots_list <- list()
# for (idx in c(1:100)) {
#   load(paste0("output/DEseq2Res_seed_", idx, ".RData"))
#   inputDat$downsampled_discovery <- DESeq2_res[[cohort]]
#   rm(DESeq2_res, DESeq2_resLFC)
#   
#   label_DEGs <- c("FCRL6") 
#   
#   # prepare the plot data
#   plotDat[[idx]] <- inputDat %>% 
#     lapply(function(x) x %>% as.data.frame %>% rownames_to_column("gene")) %>% 
#     bind_rows(.id = "cohort") %>% 
#     dplyr::select(-baseMean, -lfcSE, -stat) %>% 
#     pivot_wider(names_from = "cohort", values_from = c("log2FoldChange", "pvalue", "padj")) %>%
#     filter(padj_original_discovery < 0.05) %>%
#     mutate(sig_bothCohorts = ifelse(padj_original_discovery < 0.05 & padj_downsampled_discovery < 0.05, "yes", "no")) %>%
#     mutate(validated_DEGs = 
#              ifelse(sig_bothCohorts == "yes" & log2FoldChange_original_discovery * log2FoldChange_downsampled_discovery  > 0, 
#                     "yes", "no")) %>% 
#     mutate(validated_DEGs_v2 = ifelse(is.na(validated_DEGs), "no", validated_DEGs)) %>%
#     mutate(label_DEGs = ifelse(gene %in% label_DEGs & validated_DEGs_v2 == "yes", gene, NA))
#   
#   validated_percentage <- (plotDat[[idx]] %>% 
#                              group_by(validated_DEGs_v2)  %>% summarise(n = n()) %>% mutate(per= n/sum(n)) %>%
#                              dplyr::select(validated_DEGs_v2, per) %>% 
#                              pivot_wider(names_from = "validated_DEGs_v2", values_from = "per"))$yes %>% round(., digits = 4) * 100 
#   
#   # create plot 
#   comparePlot_RNAseq <- ggplot(data = plotDat[[idx]], 
#                                aes(x = log2FoldChange_original_discovery, 
#                                    y = log2FoldChange_downsampled_discovery,
#                                    color = validated_DEGs_v2, label = label_DEGs)) + 
#     geom_point(size = 4, alpha = 0.5) + 
#     scale_color_manual(name = "validated_DEGs",
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
#   plots_list[[idx]] <- comparePlot_RNAseq 
#   
#   rm(DESeq2_res, DESeq2_resLFC)
# }
# 
# save(plotDat, plots_list, file = "output/DEseq2Res_plotDat.RData")
# 
# # check the validated DEGs percentage
# plotDat[[2]] %>% group_by(sig_bothCohorts)  %>% summarise(n = n()) %>% mutate(per= n/sum(n)) # has NA genes
# plotDat[[2]] %>% group_by(validated_DEGs_v2)  %>% summarise(n = n()) %>% mutate(per= n/sum(n)) # has no NA genes, only yes and no
# 
# validated_DEGs <- plotDat %>% 
#   lapply(function(x) x %>% 
#            group_by(validated_DEGs_v2)  %>% summarise(n = n()) %>% mutate(per= n/sum(n)) %>%
#            dplyr::select(validated_DEGs_v2, per) %>% 
#            pivot_wider(names_from = "validated_DEGs_v2", values_from = "per")) %>% 
#   bind_rows(.id = "downsampled_seed")
# 
# range(validated_DEGs$yes) # 0.09708738 0.58183079
# hist(validated_DEGs$yes)
# summary(validated_DEGs$yes) # Mean = 0.30307 
# # create plot group and save
# #cowplot::plot_grid(plots_list[c(1:50)], nrow = 5, ncol = 10)
# 
# plotSW_1 <- do.call(grid.arrange, c(plots_list[c(1:50)], ncol = 10))
# plotSW_2 <- do.call(grid.arrange, c(plots_list[c(51:100)], ncol = 10))
# 
# png("reAnalysis_outcome/RNA/10_revision_RNApart1_discoveryCohort.png", width = 2400, height = 960)
# plot(plotSW_1)
# dev.off()
# 
# png("reAnalysis_outcome/RNA/10_revision_RNApart2_discoveryCohort.png", width = 2400, height = 960)
# plot(plotSW_2)
# dev.off()

## padj in reanalysis is calculated in DEGs from original analysis (1442 genes) ----------------
rm(list = ls())

library(tidyverse)
library(gridExtra)
library(ggpubr)
library(ggrepel)

load("dat/DEseq2Res_rna.RData")
cohort <- "discovery"

inputDat <- list()
inputDat$original_discovery <- DESeq2_res[[cohort]]
rm(DESeq2_res, DESeq2_resLFC)

plotDat <- list()
plots_list <- list()
for (idx in c(1:100)) {
  load(paste0("output/DEseq2Res_seed_", idx, ".RData"))
  inputDat$downsampled_discovery <- DESeq2_res[[cohort]]
  rm(DESeq2_res, DESeq2_resLFC)
  
  label_DEGs <- c("FCRL6") 

  # prepare the plot data
  plotDat[[idx]] <- inputDat %>% 
    lapply(function(x) x %>% as.data.frame %>% rownames_to_column("gene")) %>% 
    bind_rows(.id = "cohort") %>% 
    dplyr::select(-baseMean, -lfcSE, -stat) %>% 
    pivot_wider(names_from = "cohort", values_from = c("log2FoldChange", "pvalue", "padj")) %>%
    filter(padj_original_discovery < 0.05) %>%
    mutate(padj_downsampled_discovery_v2 = p.adjust(pvalue_downsampled_discovery, method = "fdr")) %>% # calculate padj in downsampled data based on 1442 DEGs in original analysis
    mutate(sig_bothCohorts = ifelse(padj_original_discovery < 0.05 & padj_downsampled_discovery_v2 < 0.05, "yes", "no")) %>%
    mutate(validated_DEGs = 
             ifelse(sig_bothCohorts == "yes" & log2FoldChange_original_discovery * log2FoldChange_downsampled_discovery  > 0, 
                    "yes", "no")) %>% 
    mutate(validated_DEGs_v2 = ifelse(is.na(validated_DEGs), "no", validated_DEGs)) %>%
    mutate(label_DEGs = ifelse(gene %in% label_DEGs & validated_DEGs_v2 == "yes", gene, NA)) %>%
    mutate(validated_DEGs_v3 = ifelse(is.na(label_DEGs), validated_DEGs_v2, "yes_highlight"))

  validated_percentage <- (plotDat[[idx]] %>% 
    group_by(validated_DEGs_v2)  %>% summarise(n = n()) %>% mutate(per= n/sum(n)) %>%
    dplyr::select(validated_DEGs_v2, per) %>% 
    pivot_wider(names_from = "validated_DEGs_v2", values_from = "per"))$yes %>% round(., digits = 4) * 100 
  
  # create plot 
  comparePlot_RNAseq <- ggplot(data = plotDat[[idx]] %>% arrange(validated_DEGs_v3), 
                               aes(x = log2FoldChange_original_discovery, 
                                   y = log2FoldChange_downsampled_discovery,
                                   color = validated_DEGs_v3, label = label_DEGs)) + 
    geom_point(size = 4, alpha = 0.5) + 
    scale_color_manual(name = "validated_DEGs",
                       values=c("yes_highlight" = "#660000", "no" = "#d3d3d3", "yes" = "#bb6f7c"),
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
  
  plots_list[[idx]] <- comparePlot_RNAseq 
  
  rm(DESeq2_res, DESeq2_resLFC)
}

save(plotDat, plots_list, file = "output/DEseq2Res_plotDat_v2.RData")

# check the validated DEGs percentage
plotDat[[2]] %>% group_by(sig_bothCohorts)  %>% summarise(n = n()) %>% mutate(per= n/sum(n)) # has NA genes
plotDat[[2]] %>% group_by(validated_DEGs_v2)  %>% summarise(n = n()) %>% mutate(per= n/sum(n)) # has no NA genes, only yes and no

validated_DEGs <- plotDat %>% 
  lapply(function(x) x %>% 
           group_by(validated_DEGs_v2)  %>% summarise(n = n()) %>% mutate(per= n/sum(n)) %>%
           dplyr::select(validated_DEGs_v2, per) %>% 
           pivot_wider(names_from = "validated_DEGs_v2", values_from = "per")) %>% 
  bind_rows(.id = "downsampled_seed")

range(validated_DEGs$yes) # 0.4889043 0.8654646
hist(validated_DEGs$yes)
summary(validated_DEGs$yes) # Mean = 0.7149
# create plot group and save
#cowplot::plot_grid(plots_list[c(1:50)], nrow = 5, ncol = 10)

# plotSW_1 <- do.call(grid.arrange, c(plots_list[c(1:50)], ncol = 10))
# plotSW_2 <- do.call(grid.arrange, c(plots_list[c(51:100)], ncol = 10))
# 
# png("reAnalysis_outcome/RNA/10_revision_RNApart1_discoveryCohort_v2.png", width = 2400, height = 960)
# plot(plotSW_1)
# dev.off()
# 
# png("reAnalysis_outcome/RNA/10_revision_RNApart2_discoveryCohort_v2.png", width = 2400, height = 960)
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

png("reAnalysis_outcome/RNA/10_revision_RNApart1_discoveryCohort_v3.png", width = 2400, height = 960)
plot(plotSW_1)
dev.off()

png("reAnalysis_outcome/RNA/10_revision_RNApart2_discoveryCohort_v3.png", width = 2400, height = 980)
plot(plotSW_2)
dev.off()

# compare the downsampling with original analysis outcomes in validation cohort ----------------
# load("dat/DEseq2Res_rna.RData")
# 
# cohort <- "validation"
# 
# inputDat <- list()
# inputDat$original_validation <- DESeq2_res[[cohort]]
# rm(DESeq2_res, DESeq2_resLFC)
# 
# plotDat <- list()
# plots_list <- list()
# for (idx in c(1:100)) {
#   load(paste0("output/DEseq2Res_seed_", idx, ".RData"))
#   inputDat$downsampled_validation <- DESeq2_res[[cohort]]
#   rm(DESeq2_res, DESeq2_resLFC)
#   
#   # prepare the plot data
#   plotDat[[idx]] <- inputDat %>% 
#     lapply(function(x) x %>% as.data.frame %>% rownames_to_column("gene")) %>% 
#     bind_rows(.id = "cohort") %>% 
#     dplyr::select(-baseMean, -lfcSE, -stat) %>% 
#     pivot_wider(names_from = "cohort", values_from = c("log2FoldChange", "pvalue", "padj")) %>%
#     mutate(sig_bothCohorts = ifelse(pvalue_original_validation < 0.05 & pvalue_downsampled_validation < 0.05, "yes", "no")) %>%
#     mutate(validated_DEGs = 
#              ifelse(sig_bothCohorts == "yes" & log2FoldChange_original_validation * log2FoldChange_downsampled_validation > 0, 
#                     "yes", "no")) %>%
#     filter(pvalue_original_validation < 0.05)
#   
#   # create plot 
#   comparePlot_RNAseq <- ggplot(data = plotDat[[idx]], 
#                                aes(x = log2FoldChange_original_discovery, 
#                                    y = log2FoldChange_downsampled_discovery,
#                                    color = validated_DEGs)) + 
#     geom_point(size = 4, alpha = 0.5) + 
#     scale_color_manual(values=c("no" = "#d3d3d3", "yes" = "#bb6f7c"),
#                        labels = c("no", "yes")) +
#     geom_vline(xintercept = 0) + geom_hline(yintercept = 0)+
#     theme_classic(base_size = 12) +
#     theme(legend.position = "top")
#   
#   plots_list[[idx]] <- comparePlot_RNAseq 
#   
#   rm(DESeq2_res, DESeq2_resLFC)
# }
# 
# # check the validated DEGs percentage
# plotDat[[1]] %>% group_by(sig_bothCohorts) %>% summarise(n = n())
# plotDat[[2]] %>% group_by(sig_bothCohorts) %>% summarise(n = n())
# # note: there are around 809 sig. gene in validation cohort (p-value < 0.05) at the original analysis 
# # but in the downsampled analysis, only around 380-450 genes are passed the sig. threshold at discovery cohort (p-adj < 0.05), 
# # so the number of sig. gene in validation (p-value < 0.05) at the downsampled analysis is even much lower.
# 
# validated_DEGs <- plotDat %>% 
#   lapply(function(x) x %>% 
#            group_by(validated_DEGs)  %>% summarise(n = n()) %>% mutate(per= n/sum(n)) %>%
#            dplyr::select(validated_DEGs, per) %>% 
#            pivot_wider(names_from = "validated_DEGs", values_from = "per")) %>% 
#   bind_rows(.id = "downsampled_seed")
# 
# range(validated_DEGs$yes)
# hist(validated_DEGs$yes)
