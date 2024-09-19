#!/usr/bin/env Rscript

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

#pheno <- "proteins"
#trait1 <- "FCRL6_Inflammation"

pheno <- "expression"
trait1 <- "ENSG00000243772"

args <- commandArgs(trailingOnly=TRUE)
pheno <- args[1]
trait1 <- args[2]

message("pheno ", pheno, " trait1 ", trait1)

outpath <- paste0("/home/nvanunen/nvanunen/extra_projects/cmv_revamp/out/clump")
dir.create(outpath, recursive = TRUE, showWarnings = FALSE)


exposure_snps <- fread(cmd = paste0("cat ", root, "/", cohort, "/", pheno, "/mapping/main_", pval, ".tsv | grep ", trait1), header = TRUE) %>%
    setNames(c("SNP", "gene", "beta", "t-stat", "p-value")) %>% 
    separate(SNP, c("chr", "location", "ref", "alt", "snp"), sep = "[:;]", fill = "right", remove = FALSE) %>%
    dplyr::rename("pval" = `p-value`) %>%
    filter(!is.na(snp)) %>%  #remove any where snp = NA
    mutate(se = abs(beta) / `t-stat`, location = as.numeric(location)) %>%
    mutate(pheno = pheno)
nrow(exposure_snps) # 193

# gather the effect allele frequency for each snp ;; using future
plan(multisession)
eaf_data_list <- future_lapply(seq_len(22), function(chr) {
    tmp_path <- file.path(root, cohort, "imputation", "freq", paste0("chr", chr, ".afreq"))
    tmp_dat <- fread(tmp_path, header = TRUE) %>%
        separate(ID, c("other", "snp"), sep = ";", fill = "right") %>% 
        mutate(MAF = ifelse(ALT_FREQS < 0.5, ALT_FREQS, 1 - ALT_FREQS)) %>%
        dplyr::rename(eaf = ALT_FREQS) %>%
        dplyr::select(snp, eaf)
    message("Processed chr ", chr)
    tmp_dat
})
all_eaf_dat <- do.call(rbind, eaf_data_list)
nrow(all_eaf_dat) # 8680983


exposure_snps <- exposure_snps %>% inner_join(all_eaf_dat, by = "snp")
nrow(exposure_snps) # 193 

# clump
exposure_snps_clump <- exposure_snps %>%
    dplyr::rename(rsid = "SNP") %>%
    ld_clump(clump_kb = 10000, clump_r2 = 0.001, clump_p = 1, 
        bfile = bfile, plink_bin = plinkbinr::get_plink_exe()) #genetics.binaRies::get_plink_binary())
nrow(exposure_snps_clump) # 17

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