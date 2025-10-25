rm(list = ls())

library(tidyverse)
library(caret)
library(reshape2)

get.pca_plot <- function(pca, metadat, groupType) {
  library(tidyverse)
  library(caret)
  
  ggplot(data.frame(pca$x), aes(PC1, PC2, color = metadat[, groupType])) +
    geom_point(alpha=.5) +
    theme_classic() +
    theme(aspect.ratio = 1, plot.title = element_text(hjust = .5)) +
    ggsci::scale_color_nejm() +
    labs(x = paste0("PC1 [",
                    round(100*summary(pca)$importance["Proportion of Variance", "PC1"],
                          digits = 2),
                    "%]"),
         y = paste0("PC2 [",
                    round(100*summary(pca)$importance["Proportion of Variance", "PC2"],
                          digits = 2),
                    "%]"),
         color = groupType) #+ stat_ellipse()
  
}

# load data =======================================================================
load("/vol/projects/BIIM/cohorts_CMV/CMV_2000HIV_NhanNguyen/processedDat/cohortDat.RData")

# some covariates suggested by previous analysis from the group (other teammates, and the covariet data come from the "Phenotype_2000HIV_all_01.tsv" file)
phenoDat_temp <- cohortDat$donor_info %>% 
  mutate(center.EMC = ifelse(Institute.Abbreviation == "EMC", 1, 0),
         center.ETZ = ifelse(Institute.Abbreviation == "ETZ", 1, 0),
         center.OLV = ifelse(Institute.Abbreviation == "OLV", 1, 0)) %>%
  dplyr::select(-ETHNICITY.Pacific_Islander, # only have 0 and NA
         -ETHNICITY, -cplv, -EC_combined, -Institute.Abbreviation) %>% # have >2 groups
  dplyr::select(Record.Id, AGE, SEX_BIRTH, BMI_BASELINE, matches("center|season|ETHNICITY"))

# some additional covariates suggested by collaborators (which comes later)
subjectInfo_additional <- read_csv2("/vol/projects/CIIM/2000HIV/Phenotype/221110_2000hiv_study_export_processed_2.0_SIMPLIFIED.csv")

subjectInfo_selected <- subjectInfo_additional %>% 
  dplyr::select(ID, TIMETOLAB, COVID_VACC, COVID19, PANDEMIC_BEFOREAFTER) %>% 
  rename("Record.Id" = "ID")

phenoDat_prep <- phenoDat_temp %>% left_join(subjectInfo_selected)

# prepare data ------------------------------------------------------------------------
# input data (except cytokine)
datType <- "rna"
datType <- "protein"
datType <- "mebo"
datType <- "cellcount"

inputDat <- cohortDat$allSample[[datType]]

# input data for cytokine
datType <- "cytokine"
inputDat  <- cohortDat$allSample$cytokine[, which(colMeans(!is.na(cohortDat$allSample$cytokine)) > 0.5)] # use 67 out of 85 cytokines that have <50% NAs across samples

# Covariates - PC calculation  ------------------------------------------------------------------------

## PCA calculation ------------------------------------------------------
metadat <- cohortDat$donor_info %>% filter(Record.Id %in% rownames(inputDat)) %>%
  mutate(SEX_BIRTH = as.character(SEX_BIRTH), CMV_IgG_Serology = as.character(CMV_IgG_Serology))

inputDat_pca <-  predict(preProcess(x=inputDat, method = "knnImpute"), # input the NA value (in case)
                         inputDat) %>% prcomp()

## PCA plot and Cumulative Proportion plot ------------------------------------------------------
# # PCA plot
# get.pca_plot(inputDat_pca, metadat, "Cohort")
# get.pca_plot(inputDat_pca, metadat, "CMV_IgG_Serology")
# 
# # PCs Cumulative Proportion plot
# pca_PCs <- summary(inputDat_pca)$importance %>% 
#   t() %>% as.data.frame() %>% 
#   rownames_to_column("PCs") %>% 
#   mutate(PC = gsub("PC", "", PCs,)%>% as.numeric())
# 
# pca_PCs %>%
#   ggplot(aes(x = PC, y = `Cumulative Proportion`)) + 
#   geom_line()+ theme_bw()

## calculate the accumulated variance
pca_var <- inputDat_pca$sdev^2
pca_var_per <- (pca_var / sum(pca_var) * 100)
pca_var_per[1:30]

acumullated_variance <- sum(pca_var_per[1:30])

## correlations PCs - variables (adapt from Manoj codes) ------------------------
n_princomps <- 30 # PCs 1-30
princomps <- inputDat_pca$x[, 1:n_princomps] 

phenoDat <- phenoDat_prep %>% 
  slice(which(Record.Id %in% rownames(princomps))) %>% 
  arrange(match(Record.Id, rownames(princomps))) %>%
  column_to_rownames("Record.Id")

identical(rownames(phenoDat), rownames(princomps))

# calculate the correlation PCs - variable 
corr.pvalues <- matrix(nrow = ncol(princomps), 
                       ncol = ncol(phenoDat), 
                       dimnames = list(colnames(princomps), colnames(phenoDat)))

for (var1 in colnames(phenoDat)) {
  for ( pc in colnames(princomps)) {
    print(pc)
    print(var1)
    #res <- cor.test(princomps[, pc], pheno2[, var1])
    #pval <- -1 * log10(res$p.value)
    #pval <- res$p.value
    #corr.pvalues[pc, var1] <- pval
    
    lm_res <- lm(princomps[, pc] ~ phenoDat[, var1] ) # lm() similar to cor.test, but can compare categories
    lm_outcome <- data.frame(summary(lm_res)$coefficients[2,]) # if only take row 2 --> lm() for 2 categories also not work
    pvali <- lm_outcome[4,]
    corr.pvalues[pc, var1] <- pvali
  }
}

## prepare the plot data -----------------------------------------------------------------
df.res <- melt(corr.pvalues) %>%
  rename("PC" = "Var1", "cov" = "Var2", "pval" = "value") %>%
  mutate(significance = cut(pval, breaks=c(-Inf, 1e-10,1e-5, 0.01, 0.05, Inf), 
                            label=c("P<1e-10","P<1e-5", "P<0.01", "P<0.05", "NS"))) %>%
  mutate(PC = factor(PC, levels= colnames(inputDat_pca$x[, 1:n_princomps])),
         cov = factor(cov, levels = colnames(phenoDat)))


## plot the figure ---------------------------------------------------------
plot_pcaCovariates <- ggplot(data = df.res, aes(x = PC, y = cov, fill = significance)) +
  geom_tile() + scale_fill_manual(values = c("#E31A1C","#FC4E2A","#FEB24C","#FFF3B2","white")) +
  ylab("Covariates") + xlab("Principal component")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0),
        legend.position = "bottom", 
        text = element_text(size = 18)) +
  ggtitle(paste0("Accumulated variance: ", round(acumullated_variance, 2), "%, ", datType))

plot_pcaCovariates

# pca_PCs %>% filter(PC <= 31) %>%
#   ggplot(aes(x = PC, y = `Proportion of Variance`)) + 
#   geom_bar(stat = "identity") + theme_bw()

# save the plot ------------------------------------------------------------------
png(paste0("output/03_02_pcaCovariates_" , datType, ".png"), width = 1008)
plot_pcaCovariates
dev.off()





