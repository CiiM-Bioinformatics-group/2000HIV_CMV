# check pathway
rm(list = ls())
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(enrichplot)
library(stringr)
library(DOSE)
library(aPEAR)

## load data -------------------------------------------------------
load("processedDat/DEseq2Res_rna.RData")

dat <- DESeq2_res %>% 
  lapply(function(x) x %>% as.data.frame %>% rownames_to_column("gene")) %>% 
  bind_rows(.id = "cohort") %>% 
  dplyr::select(-baseMean, -lfcSE, -stat) %>% 
  pivot_wider(names_from = "cohort", values_from = c("log2FoldChange", "pvalue", "padj")) %>%
  mutate(sig_bothCohorts = ifelse(padj_discovery < 0.05 & pvalue_validation < 0.05, "yes", "no")) %>%
  mutate(validated_DEGs = 
           ifelse(sig_bothCohorts == "yes" & log2FoldChange_discovery * log2FoldChange_validation  > 0, 
                  "yes", "no")) %>%
  filter(padj_discovery < 0.05)

DEGs <- (dat %>% filter(validated_DEGs == "yes"))$gene

DEGs_up <- dat %>% filter(validated_DEGs == "yes", log2FoldChange_discovery > 0) %>%
  dplyr::select(gene) %>% unlist()

DEGs_down <- dat %>% filter(validated_DEGs == "yes", log2FoldChange_discovery < 0) %>%
  dplyr::select(gene) %>% unlist()

## Adapt code from Xun's pathway analysis ---------------------------------
# only using GO analysis because KEGG have no result since no gene can be mapped ???
geneList <- data.frame(bitr(geneID = DEGs, 
                            fromType = 'SYMBOL', toType = c('ENTREZID', 'ENSEMBL'),
                            OrgDb = org.Hs.eg.db, drop = T))

pathwayAnalysis <- function(genenames, pval = 0.05, qval = 0.05) {
  geneList <- data.frame(bitr(geneID = genenames, 
                              fromType = 'SYMBOL', toType = c('ENTREZID', 'ENSEMBL'),
                              OrgDb = org.Hs.eg.db, drop = T))
  
  enr.GO.BP <- enrichGO(gene = geneList$ENTREZID,
                        OrgDb = org.Hs.eg.db,
                        keyType = "ENTREZID",
                        ont = "BP", # Biological process
                        pAdjustMethod = "fdr", #none or fdr
                        pvalueCutoff = pval,
                        qvalueCutoff = qval,
                        readable = TRUE) 
  
  enr.GO.MF <- enrichGO(gene = geneList$ENTREZID,
                        OrgDb = org.Hs.eg.db,
                        keyType = "ENTREZID",
                        ont = "MF", # molecular function
                        pAdjustMethod = "fdr", #none or fdr
                        pvalueCutoff = pval,
                        qvalueCutoff = qval,
                        readable = TRUE)
  
  res = list('GO.BP' = enr.GO.BP, 'GO.MF' = enr.GO.MF)
  return(res)
}

## get pathway results -----------------------------------------------------------
res <- list("up" = pathwayAnalysis(genenames = DEGs_up), 
            "down" = pathwayAnalysis(genenames = DEGs_down)) %>% 
  rlist::list.flatten()


#res_sigPathways <- res %>% lapply(function(x) x@result %>% filter(qvalue < 0.05))

res_sigPathways <- res %>% 
  lapply(function(x) x@result %>% filter(qvalue < 0.05)) %>% 
  bind_rows(.id = "compare") %>% 
  rename("pathwayGenes" = "geneID") %>% 
  mutate( GeneRatio_v2 = sapply(GeneRatio, function(x) eval(parse(text = x))))

res_sigPathways_plot <- res_sigPathways %>% 
  mutate(group = ifelse(compare %in% c("down.GO.BP", "down.GO.MF"), "downDEGs", "upDEGs"))

res_sigPathways_plot %>% as.data.frame() %>% dplyr::count(group)

sigPathwayNumbers <- res_sigPathways_plot %>% 
  ggplot(aes(x = group)) + geom_bar() + 
  labs(y = "Number of pathways (qvalue < 0.05)") + ylim(0, 160) +
  theme_classic(base_size = 20) +
  theme(plot.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.position = "none") + coord_flip()

# save the plot
png("output/04_03_rna_sigPathwayNumbers.png", width = 576, height = 288)
sigPathwayNumbers 
dev.off()

# save data 
save(res, res_sigPathways_plot, 
     file = "processedDat/DEseq2Res_rna_pathwayAnalysis.RData")

## plot pathway network-----------------------------------------------------------
pathwayNetwork <- enrichmentNetwork(res_sigPathways, 
                  colorBy = "qvalue", 
                  colorType = "pval", 
                  nodeSize = "Count", 
                  # nodeSize = "GeneRatio_v2",
                  verbose = TRUE, 
                  drawEllipses = TRUE, 
                  fontSize = 6, repelLabels = TRUE) #+ scale_color_gradientn(colours = c("#D47B5F", "#5FB7D4"))
pathwayNetwork

# save the plot
png("output/04_03_rna_pathwayNetwork.png", width = 816, height = 816)
pathwayNetwork
dev.off()

pathwayNetwork_upDEGs <- enrichmentNetwork(res_sigPathways %>% filter(compare %in% c("up.GO.BP", "up.GO.MF")), 
                  colorBy = "qvalue", 
                  colorType = "pval", 
                  nodeSize = "Count", 
                  # nodeSize = "GeneRatio_v2",
                  verbose = TRUE, 
                  drawEllipses = TRUE, 
                  fontSize = 6, repelLabels = TRUE) #+ scale_color_gradientn(colours = c("#D47B5F", "#5FB7D4"))

pathwayNetwork_upDEGs

# save the plot
png("output/04_03_rna_pathwayNetwork_upDEGs.png", width = 720, height = 720)
pathwayNetwork_upDEGs
dev.off()

# plot pathway tree --------------------------------------------------------------------
options(enrichplot.colours = c("white","blue"))

nodeid.tbl_tree <- utils::getFromNamespace("nodeid.tbl_tree", "tidytree")
rootnode.tbl_tree <- utils::getFromNamespace("rootnode.tbl_tree", "tidytree")
offspring.tbl_tree <- utils::getFromNamespace("offspring.tbl_tree", "tidytree")
offspring.tbl_tree_item <- utils::getFromNamespace(".offspring.tbl_tree_item", "tidytree")
child.tbl_tree <- utils::getFromNamespace("child.tbl_tree", "tidytree")
parent.tbl_tree <- utils::getFromNamespace("parent.tbl_tree", "tidytree")

plotDat <- pairwise_termsim(res$up.GO.BP) %>% 
  treeplot(nWords = 5,  nCluster = 5, showCategory = 35, 
           cex_category = 0.7, 
           hexpand = 0.15, 
           offset = rel(3), offset_tiplab = 0.8,
          xlim = c(0,20),
          
           #hilight = FALSE, 
           # offset = 9,
             # legend_n = 3,
           frontsize = 60,
          label_format = 20,
           group_color = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442")
  ) + theme(legend.position = "bottom") 

plotDat$layers[[7]]$aes_params$size <- 4
plotDat 

ggsave(paste0("fig/pathwayTree.png"), 
       plotDat, device = "png")

# save the plot
png("output/04_03_rna_pathwayTree_upDEGs_GOBP.png", width = 624, height = 720)
plotDat
dev.off()


plotDat_leftLegend <- pairwise_termsim(res$up.GO.BP) %>% 
  treeplot(nWords = 5,  nCluster = 5, showCategory = 35, 
           cex_category = 0.7, 
           hexpand = 0.15, 
           offset = rel(3), offset_tiplab = 0.8,
           xlim = c(0,20),
           
           #hilight = FALSE, 
           # offset = 9,
           # legend_n = 3,
           frontsize = 60,
           label_format = 20,
           group_color = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442")
  ) 

# save the plot
png("output/04_03_rna_pathwayTree_upDEGs_GOBP_leftLegend.png", width = 720, height = 720)
plotDat_leftLegend
dev.off()

