#!/usr/bin/env Rscript

# mr between FCRL6 and CMV GWAS

library(tidyverse)
library(data.table)

library(future)
library(future.apply)

library(ieugwasr)
library(plinkbinr)
#library(LDlinkR)

#options(future.globals.maxSize = 8000 * 1024^2)

root <- "/vol/projects/CIIM/meta_cQTL/out"
cohort <- "2000HIV-EU-discovery"
cov <- "main"
bfile <- paste0(root, "/", cohort, "/genotype/allchr")
pval <- "1e-5" 

pheno <- "proteins"
trait1 <- "FCRL6_Inflammation"

pheno <- "expression"
trait1 <- "ENSG00000181036"
trait1 <- "ENSG00000243772" # Chromosome 19: 54,738,513-54,753,052 

args <- commandArgs(trailingOnly=TRUE)
pheno <- args[1]
trait1 <- args[2]

message("pheno ", pheno, " trait1 ", trait1)

# ADAMTS1 ENSG00000154734 Chromosome 21: 26,835,755-26,845,409
# ARHGEF12 ENSG00000196914 Chromosome 11: 120,336,413-120,489,937
# FCRL6 ENSG00000181036 Chromosome 1: 159,800,511-159,816,257
# GFI1
# GSTM2 ENSG00000213366 Chromosome 1: 109,668,022-109,709,551
# ITGB1 ENSG00000150093 Chromosome 10: 32,887,273-33,005,792
# LTBP3 ENSG00000168056 Chromosome 11: 65,538,559-65,558,930
# MXRA7 ENSG00000182534 Chromosome 17: 76,672,551-76,711,004
# NLRC5 ENSG00000140853 Chromosome 16: 56,989,485-57,083,531 
# PATL2 ENSG00000229474 Chromosome 15: 44,665,732-44,711,323 
# SPG11 ENSG00000104133 Chromosome 15: 44,554,818-44,663,688
# TAP1 ENSG00000168394 Chromosome 6: 32,845,209-32,853,816 
# TAP2
# WNT10B
# ZFP28
# ZNF256
# ZNF418
# ZNF677
# ZNF814. 

# ADGRG1
# CRTAM
# FCRL6
# GNLY
# GZMH
# ITGAL
# KIR2DL3
# KLRD1

genes <- fread("/home/nvanunen/nvanunen/extra_projects/cmv_revamp/genes.txt")
gene <- genes %>% dplyr::filter(ensembl_gene_id == trait1)
print(gene)
#chr <- 19
#min <- 54738513 - 1e6
#max <- 54753052  + 1e6
chr <- gene$chromosome_name
min <- gene$start_position - 1e6
max <- gene$end_position + 1e6

outpath <- paste0("/home/nvanunen/nvanunen/extra_projects/cmv_revamp/out/clump_cis")
dir.create(outpath, recursive = TRUE, showWarnings = FALSE)


cmd <- paste0("awk -F '[:;\\t]' 'NR == 1 || $6 == \"", trait1, "\" && $2 >= ", min, " && $2 <= ", max, "' ", root, "/", cohort, "/", pheno, "/mapping/main/nominal/chr", chr, ".tsv")
message(cmd)
exposure_snps <- fread(cmd = cmd) %>%
    setNames(c("SNP", "gene", "beta", "t-stat", "p-value")) %>% 
    separate(SNP, c("chr", "location", "ref", "alt", "snp"), sep = "[:;]", fill = "right", remove = FALSE) %>%
    dplyr::rename("pval" = `p-value`) %>%
    filter(!is.na(snp)) %>%  #remove any where snp = NA
    mutate(se = abs(beta) / `t-stat`, location = as.numeric(location)) %>%
    mutate(pheno = pheno)
nrow(exposure_snps) # 1137 prot - 497 expr

tmp_path <- file.path(root, cohort, "imputation", "freq", paste0("chr", chr, ".afreq"))
all_eaf_dat <- fread(tmp_path, header = TRUE) %>%
    separate(ID, c("other", "snp"), sep = ";", fill = "right") %>% 
    mutate(MAF = ifelse(ALT_FREQS < 0.5, ALT_FREQS, 1 - ALT_FREQS)) %>%
    dplyr::rename(eaf = ALT_FREQS) %>%
    dplyr::select(snp, eaf)

exposure_snps <- exposure_snps %>% inner_join(all_eaf_dat, by = "snp")
nrow(exposure_snps) # 787 prot - 497 expr


exposure_snps_filter <- exposure_snps %>% dplyr::filter(pval < 1e-3)
dim(exposure_snps_filter) # 169 prot - 27 expr

# clump
exposure_snps_clump <- exposure_snps_filter %>% 
    dplyr::rename(rsid = "SNP") %>%
    ld_clump(clump_kb = 10000, clump_r2 = 0.001, clump_p = 1, 
        bfile = bfile, plink_bin = plinkbinr::get_plink_exe()) #genetics.binaRies::get_plink_binary())
nrow(exposure_snps_clump) # 

exposure_dat <- exposure_snps_clump %>%
    dplyr::rename(ID = rsid, SNP = snp, pval.exposure = "pval", beta.exposure = "beta",
    chr.exposure = "chr", pos.exposure = "location", effect_allele.exposure = "ref",
    other_allele.exposure = "alt", se.exposure = "se", eaf.exposure = "eaf", 
    exposure = "gene", id.exposure = "pheno") %>%
    arrange(ID)


#if (nrow(exposure_dat) >= 3) {
# write to file
fwrite(exposure_dat, file.path(outpath, paste0(trait1, ".tsv")), sep = "\t")
#}