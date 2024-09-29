#!/usr/bin/env Rscript
# File: cmv_gwas.r
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Jan 05, 2024
# Updated: Jan 08, 2024

options(stringsAsFactors = FALSE, datatable.verbose = FALSE, datatable.showProgress = FALSE)
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)

  library(AnnotationHub)
  library(ensembldb)

  library(DESeq2)
  library(DEGreport)
})

# Annotations
annhub <- AnnotationHub::AnnotationHub()
ann_db <- AnnotationHub::query(annhub, pattern = c("GRCh38", "EnsDb"))

ensembl_db <- ann_db[["AH104864"]]
ensembl_genes <- ensembldb::genes(ensembl_db)
GenomeInfoDb::seqlevelsStyle(ensembl_genes) <- "UCSC"

# Create DDS of DESeq2
dds_save_to <- "~/Documents/projects/wp_2000hiv/outputs/rna_seq/DESeq2/objects/2000HIV.bulk_rnaseq.dds_object.RDS"
if (file.exists(dds_save_to)) {
  cat("Loading from disk ...", dds_save_to, "\n")
  dds <- readRDS(dds_save_to)
} else {
  meta_data_general <- readRDS("/vol/projects/BIIM/2000HIV/RNAseq/2000HIV_bulk_transcriptomics_sample_table.RDS") %>%
    dplyr::mutate(FID = DONOR_ID, IID = DONOR_ID) %>%
    dplyr::select(-DONOR_ID)
  kept_samples <- rownames(meta_data_general)

  read_counts <- readRDS("/vol/projects/BIIM/2000HIV/RNAseq/2000HIV_bulk_transcriptomics_raw_counts.RDS")
  feature_tab <- rownames(read_counts) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(gene_id = stringr::str_extract(value, "ENSG[0-9]+")) %>%
    dplyr::mutate(kept = !stringr::str_detect(value, "_PAR_Y") & gene_id %in% ensembl_genes$gene_id)
  kept_features <- dplyr::pull(feature_tab, kept)
  read_counts <- read_counts[kept_features, kept_samples]
  rownames(read_counts) <- dplyr::filter(feature_tab, kept) %>% dplyr::pull(gene_id)

  dds <- DESeqDataSetFromMatrix(countData = read_counts, colData = meta_data_general, design = ~ 1)
  dds@rowRanges <- (ensembl_genes[rownames(dds), ] %>% split(.$gene_id))[rownames(dds)]
  dds$AGE <- as.numeric(dds$AGE)
  dds$CENTER <- as.factor(dds$CENTER)
  dds$SEX_BIRTH <- as.factor(dds$SEX_BIRTH)
  dds$ETHNICITY <- as.factor(dds$ETHNICITY)
  mean_cart_duration <- mean(dds$CART_DURATION, na.rm = TRUE)
  dds$CART_DURATION <- dds$CART_DURATION %>% purrr::map_dbl(~ifelse(is.na(.x), mean_cart_duration, .x))

  mcols(dds)$geneID <- rownames(dds)
  mcols(dds)$symbol <- ensembl_genes[rownames(dds), ]$symbol

  cat("Dumping into disk ...", dds_save_to, "\n")
  saveRDS(dds, dds_save_to)
}

# Obtain normalized gene expression matrix
expmat_save_to <- "~/Documents/projects/wp_2000hiv/outputs/rna_seq/DESeq2/2000HIV.bulk_rnaseq.normalized_expr.csv"
if (file.exists(expmat_save_to)) {
  cat("Loading from disk ...", expmat_save_to, "\n")
  norm_counts <- fread(expmat_save_to)
} else {
  dds <- estimateSizeFactors(dds)
  norm_counts <- counts(dds, normalized = TRUE)
  old_names <- colData(dds) %>% as.data.frame() %>% dplyr::mutate(sample_index = rownames(.)) %>% dplyr::pull(FID, sample_index)
  colnames(norm_counts) <- old_names[colnames(norm_counts)]
  norm_counts <- norm_counts %>%
    as.data.frame() %>%
    dplyr::mutate(ensembl_id = rownames(.)) %>%
    dplyr::mutate(gene_symbol = ensembl_genes[ensembl_id, ]$symbol) %>%
    dplyr::filter(gene_symbol != "") %>%
    dplyr::mutate(`#chrom` = seqnames(ensembl_genes[ensembl_id, ]) %>% as.character()) %>%
    dplyr::mutate(start = start(ensembl_genes[ensembl_id, ]) %>% as.integer()) %>%
    dplyr::mutate(end = end(ensembl_genes[ensembl_id, ]) %>% as.integer()) %>%
    dplyr::mutate(strand = strand(ensembl_genes[ensembl_id, ]) %>% as.character()) %>%
    dplyr::filter(stringr::str_starts(`#chrom`, "chr")) %>%
    dplyr::select(-dplyr::starts_with("FAM")) %>%
    dplyr::relocate(`#chrom`, start, end, gene_symbol, ensembl_id, strand) %>%
    dplyr::arrange(`#chrom`, start, end) %>%
    as.data.table()

  cat("Dumping into disk ...", expmat_save_to, "\n")
  fwrite(norm_counts, expmat_save_to)
}


# Covariates
covmat_save_to <- "~/Documents/projects/wp_2000hiv/outputs/rna_seq/DESeq2/2000HIV.bulk_rnaseq.covariates.tsv"
if (file.exists(covmat_save_to)) {
  covar_mat <- fread(covmat_save_to)
} else {
  selected_donors <- colnames(norm_counts) %>% purrr::keep(~!.x %in% c("#chrom", "start", "end", "ensembl_id", "gene_symbol", "strand"))

  pca_eigen_tab <- "~/Documents/projects/wp_2000hiv/outputs/variants/gsa/pca/2000HIV.all.pca.eigenvec" %>%
    data.table::fread() %>%
    (function(tab) { colnames(tab) <- c("FID", "IID", paste0("PC", 1:(ncol(tab) - 2))); tab }) %>%
    dplyr::select(-c(IID, PC11:PC20)) %>%
    dplyr::filter(FID %in% selected_donors)

  covar_mat <- colData(dds) %>%
    as.data.frame() %>%
    dplyr::filter(FID %in% selected_donors) %>% #, ETHNICITY %in% c("White", "Hispanic")) %>%
    dplyr::select(FID, AGE, CENTER, SEX_BIRTH, ETHNICITY, HIV_DURATION, CART_DURATION, CD4_NADIR) %>%
    dplyr::left_join(pca_eigen_tab, by = "FID") %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), as.character)) %>%
    tidyr::pivot_longer(cols = -FID, names_to = "id", values_to = "value") %>%
    tidyr::pivot_wider(names_from = FID, values_from = value) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~dplyr::if_else(is.na(.x), "NA", .x))) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~stringr::str_remove_all(.x, " "))) %>%
    as.data.table()

  covar_mat %>% fwrite(covmat_save_to, sep = "\t")
}
