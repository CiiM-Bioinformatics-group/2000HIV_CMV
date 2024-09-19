#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(ggplot2)

library(future)
library(future.apply)

library(TwoSampleMR) # mr
library(ieugwasr)
library(plinkbinr)
#library(LDlinkR)

root <- "/vol/projects/CIIM/meta_cQTL/out"
cohort <- "2000HIV-EU-discovery"
cov <- "main"

pheno2 <- "cmv"   # outcome
pheno2 <- "cytokines"   # outcome
# /vol/project/BIIM/2000HIV/RareVariants/outputs/associations/gwas/summary_statistic/cmv_infection
file <- "/home/nvanunen/nvanunen/extra_projects/hua_gwas_mr/out_new/FCRL6_Inflammation.tsv"
file <- "/home/nvanunen/nvanunen/extra_projects/hua_gwas_mr/out_new/ENSG00000181036.tsv"

pheno1 <- "expression"
trait1 <- "ENSG00000243772"

args <- commandArgs(trailingOnly=TRUE)
pheno1 <- args[1]
trait1 <- args[2]
pheno2 <- args[3] # cmv or cytokines




outpath <- paste0("/home/nvanunen/nvanunen/extra_projects/hua_gwas_mr/out_new")
dir.create(outpath, showWarnings = FALSE, recursive = TRUE)
message(outpath)

# clumping is done ahead with genomewide_clumping.R
trait1 <- file %>% basename() %>% str_remove("\\.tsv")
message(trait1)

exposure_dat <- fread(file)

# Get unique chr.exposure values from exposure_dat
unique_chr_exposure <- unique(exposure_dat$chr.exposure)

# Function to process data for each unique chr.exposure


if (pheno2 == "cytokines") {
    process_data <- function(chr_exposure) {
        cmd <- paste0("awk -F '\t' 'NR == 1 || ", 
                    paste0("$1 == \"", exposure_dat$ID[exposure_dat$chr.exposure == chr_exposure], collapse = "\" || "),
                    "\"' ", root, "/", cohort, "/", pheno2, "/mapping/", cov, "/", chr_exposure, ".tsv")
        tmp <- fread(cmd = cmd)
        message("Done with ", chr_exposure)
        tmp
    }

    # then get all the same snps from the other data 
    message("Getting all snps from other data ...")
    plan(multisession)
    outcome_all_list <- future_lapply(unique_chr_exposure, process_data)
    outcome_all_raw <- do.call(rbind, outcome_all_list)
    nrow(outcome_all_raw) # 661

    # fix some things
    message("Fixing some things ...")
    outcome_all <- outcome_all_raw %>%
        separate(SNP, c("chr", "location", "effect_allele", "other_allele", "SNP"), sep = "[:;]") %>%
        filter(!is.na(SNP)) %>%  #remove any where snp = NA
        rename(pval = "p-value") %>%
        mutate(se = abs(beta) / `t-stat`, location = as.numeric(location)) %>%
        inner_join(exposure_dat %>% 
                    select(SNP, eaf.exposure), by = "SNP") %>%
        mutate(pheno = pheno2)
} else if (pheno2 == "cmv") {
    gwas <- fread("/vol/projects/BIIM/2000HIV/RareVariants/outputs/associations/cmv_infection/gwas/summary_statistic/eur_discovery/gwas_sumstat.cmv_infection.glm.logistic.hybrid.gz")
    head(gwas)
    gwas <- gwas %>% filter(ID %in% exposure_dat$SNP)
    outcome_all <- gwas %>%
        mutate(beta = log(OR)) %>%
        select(pval = "P", chr = "#CHROM", location = "POS", SNP = "ID", effect_allele = "REF", other_allele = "ALT", beta, se =  `LOG(OR)_SE`) %>%
        mutate(gene = pheno2) %>%
        inner_join(exposure_dat %>% 
                    select(SNP, eaf.exposure), by = "SNP") %>%
        mutate(pheno = "GWAS")
}

fwrite(outcome_all, paste0(outpath, "/", trait1, "_", pheno2, ".tsv"), sep = "\t")
