#!/usr/bin/env Rscript

# mr between FCRL6 and CMV GWAS

library(tidyverse)
library(data.table)

library(future)
library(future.apply)

library(ieugwasr)
library(plinkbinr)


n_samples <- 1064

# load in eqtl cis and eqtl trans snp and perform MR
path <- "/home/nvanunen/nvanunen/extra_projects/cmv_revamp/out"
outpath <-  paste0(path, "/tier1")
dir.create(outpath, recursive = TRUE, showWarnings = FALSE)

pheno <- "expression"
trait <- "ENSG00000005844"

args <- commandArgs(trailingOnly=TRUE)
pheno <- args[1]
trait <- args[2]

file1 <- paste0(path, "/clump_cis/", trait, ".tsv")
file2 <- paste0(path, "/clump/", trait, ".tsv")
qtl_file <- paste0("/vol/projects/CIIM/meta_cQTL/out/2000HIV-EU-discovery/", pheno, "/mapping/main_1e-5.tsv")
qtls <- fread(qtl_file)

# load in
eqtl_cis <- fread(file1) %>% dplyr::select(-id)
dim(eqtl_cis) # 1
eqtl_trans <- fread(file2) %>% dplyr::select(-id)
dim(eqtl_trans) # 31
eqtl <- rbind(eqtl_cis, eqtl_trans) %>% unique()
dim(eqtl) # 32

# for snp in eqtl
rets <- data.frame()
i <- 2
for (i in 1:nrow(eqtl)) {
    message(i)
    row <- eqtl[i,]
    # find if SNP has more significant associations
    ret <- qtls %>% dplyr::filter(SNP == row$ID)
    min <- ret %>% dplyr::filter(`p-value` == min(ret$`p-value`))
    min$tot <- nrow(ret)
    rets <- rbind(rets, min)
}
uwu <- rets %>% 
    separate(SNP, c("chr", "pos", "ref", "alt", "snp"), sep = "[:;]", remove = FALSE) %>%
    mutate(chr = as.numeric(gsub("chr", "", chr))) %>%
    arrange(chr, pos)
tier1 <- uwu %>% dplyr::filter(tot <= 4)

# save tier1
write.table(tier1, paste0(outpath, "/", trait, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
