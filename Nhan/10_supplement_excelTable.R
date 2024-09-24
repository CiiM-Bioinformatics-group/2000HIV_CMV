rm(list = ls())

library(openxlsx)
library(tidyverse)

# cytokine production level, compare DEs between discovery and validation cohorts =====================
load("processedDat/rlmRes_cytokine.RData")

res_cyto <- rlm_res

#rm(cytokine_NAs, rlm_res)
rm(rlm_res)

# cell count level, compare DEs between discovery and validation cohorts ============================
load("processedDat/rlmRes_cellcount.RData")

res_cellsubset <- rlm_res
rm(rlm_res)

# methylation (from Xun analysis) ===================================
# cis_eQTM <- list()
# 
# cis_eQTM$discovery <- read.table("/vol/projects/BIIM/2000HIV/AnalysisDV/Discovery_MethylationDataBased/4EWAS_AllEthnicity/CMV_rlm/ForCMVPaper/cis-eQTM_CMV-validated_DMS.tsv", header = TRUE)
# cis_eQTM$validation <- read.table("info_scripts_fromOthers/cmv_season_BMI_SP_cor_eQTM_FDR_rep005_rlm_baconAdjusted.tsv")
# 
# methy_mediationTest <- list()
# methy_mediationTest$toCytokine <- read.table("/vol/projects/BIIM/2000HIV/AnalysisDV/Discovery_MethylationDataBased/4EWAS_AllEthnicity/CMV_rlm/ForCMVPaper/CMV2DMS2CytokineResponse.tsv", header = TRUE)
# methy_mediationTest$toGene <- read.table("/vol/projects/BIIM/2000HIV/AnalysisDV/Discovery_MethylationDataBased/4EWAS_AllEthnicity/CMV_rlm/ForCMVPaper/CMV_ValidateDMS_Gene_mediationTest.tsv", header = TRUE)


eQTM <- read.table("/vol/projects/BIIM/2000HIV/AnalysisDV/Discovery_MethylationDataBased/4EWAS_AllEthnicity/CMV_rlm/ForCMVPaper/eQTM.tsv", 
                   header = TRUE)

# RNA level, compare DEGs between discovery and validation cohorts ===================================
load("processedDat/DEseq2Res_rna.RData")

res_geneExpression <- DESeq2_res %>% 
  lapply(function(x) x %>% as.data.frame %>% rownames_to_column("gene"))
rm(DESeq2_res, DESeq2_resLFC)

## RNA level, pathway analysis  -------------------------------------
load("processedDat/DEseq2Res_rna_pathwayAnalysis.RData")

res_RNAseq_sigPathways <- res_sigPathways_plot
rm(res_sigPathways_plot)

# protein level, compare DEGs between discovery and validation cohorts ============================
load("processedDat/rlmRes_protein.RData")

res_proteinLevels <- rlm_res
rm(rlm_res)

# metabolite level, compare DEGs between discovery and validation cohorts ============================
load("processedDat/rlmRes_mebo.RData")

res_metaboliteLevels <- rlm_res
rm(rlm_res)

# save as excel files ====================================================================================
write.xlsx(list(cytokines_discovery = res_cyto$discovery,
                cytokines_validation = res_cyto$validation,
                cellsubsets_discovery = res_cellsubset$discovery,
                cellsubsets_validation = res_cellsubset$validation,
                # methylation_eQTM_discovery = cis_eQTM$discovery, # old results from Methylation analysis
                # methylation_eQTM_validation = cis_eQTM$validation, # old results 
                # mediationTest_cmv_toDMS_toGene = methy_mediationTest$toGene, # old results 
                # mediationTest_cmv_toDMS_toCytokine = methy_mediationTest$toCytokine, # old results 
                # methylation_eQTM = eQTM, # update results from Methylation analysis
                geneExpressions_discovery = res_geneExpression$discovery,
                geneExpressions_validation = res_geneExpression$validation,
                geneExpressions_sigPathways = res_RNAseq_sigPathways,
                proteinLevels_discovery = res_proteinLevels$discovery,
                proteinLevels_validation = res_proteinLevels$validation,
                metaboliteLevels_discovery = res_metaboliteLevels$discovery,
                metaboliteLevels_validation = res_metaboliteLevels$validation), 
           "output/suppTables.xlsx", rowNames = TRUE)
