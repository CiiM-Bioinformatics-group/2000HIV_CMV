#!/usr/bin/env Rscript
# File: cmv_gwas.r
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Aug 30, 2023
# Updated: May 29, 2024
options(stringsAsFactors = FALSE, datatable.showProgress = FALSE, datatable.verbose = FALSE)
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(topr)
  library(ggsci)
})

#' Calulate inflation factor
inflation_factor <- function(pval) { median(qchisq(1 - pval, 1), na.rm = TRUE) / qchisq(0.5, 1) }


proj_dir <- "~/Documents/projects/wp_2000hiv"

partition_vec <- c("eur_discovery", "eur_validation")# , "eur_all")
measurement_vec <- c("cmv_infection")
pheno_list <- list(c("cmv_infection")) %>% purrr::set_names(measurement_vec)

# Plot Manhattan and QQ plots
p_threshold <- signif(5e-8 / Reduce(sum, lapply(pheno_list, length)), 3)
sumstat_save_to <- file.path(proj_dir, "outputs/associations/cmv_infection/gwas/overview/cmv_infection/eur.summary_statistic_tab.significant_5e-8.csv") 
if (!file.exists(sumstat_save_to)) {
  all_sumstat <- lapply(measurement_vec, function(per_ms) {
    # per_ms <- "percentage"
    base_save_to <- file.path(proj_dir, "outputs/associations/cmv_infection/gwas/overview/cmv_infection", per_ms)
    if (!dir.exists(base_save_to)) dir.create(base_save_to, recursive = TRUE)

    cat("[I]: Processing measurement:", per_ms, "...\n")
    lapply(pheno_list[[per_ms]], function(per_pheno) {
      # per_pheno <- pheno_list[[per_ms]][1]
      cat("[I]: Processing phenotype: ", per_pheno, " of ", per_ms, " measurement ...\n", sep = "")

      per_file <- paste0("gwas_sumstat.", per_pheno, ".glm.logistic.hybrid.gz")
      sumstat_list <- lapply(partition_vec, function(x) {
        sumstat_path <- file.path(proj_dir, "outputs/associations/cmv_infection/gwas/summary_statistic/", per_pheno, x, per_file)
        sumstat_tab <- fread(sumstat_path) %>%
          dplyr::rename(c("CHROM" = "#CHROM")) %>%
          dplyr::filter(0.05 < A1_FREQ, A1_FREQ < 0.95, nchar(REF) == 1, nchar(ALT) == 1)
      }) %>%
        purrr::set_names(partition_vec)

      # Manhattan and QQ plots
      qq_save_to <- file.path(base_save_to, paste0("cmv_infection.", per_pheno, ".eur_discovery_validation.qq_plot.pdf"))
      mp_dis_save_to <- file.path(base_save_to, paste0("cmv_infection.", per_pheno, ".eur_discovery.manhattan_plot.pdf"))
      mp_disval_save_to <- file.path(base_save_to, paste0("cmv_infection.", per_pheno, ".eur_discovery_validation.manhattan_plot.pdf"))
      if (all(file.exists(qq_save_to, mp_dis_save_to, mp_disval_save_to))) {
        cat("[I]: Skipping Manhattan and QQ plots. Working on", per_pheno, "...\n")
      } else {
        # Y axis max
        y_max <- (lapply(sumstat_list, function(x) max(-log10(x$P))) %>% unlist() %>% max(., -log10(p_threshold), na.rm = TRUE)) + 2

        # Manhattan plot
        p <- sumstat_list[[1]] %>%
          dplyr::filter(P < 0.05) %>%
          manhattan(legend_labels = "Discovery", color = c("#A86A9D"), sign_thresh = p_threshold, annotate = 5e-08, ntop = 1,, annotate_with = "ID", region_size = 100000000, ymax = y_max)
        ggsave(mp_dis_save_to, p, width = 16, height = 6)

        # Manhattan plot, with validation
        p <- lapply(sumstat_list, function(tab) dplyr::filter(tab, P < 0.05)) %>%
          manhattan(
            legend_labels = c("Discovery", "Validation"), sign_thresh = p_threshold, annotate = 5e-08, ntop = 1,
            region_size = 100000000, highlight_genes_ypos = -0.5, ymax = y_max
          )
        ggsave(mp_disval_save_to, p, width = 16, height = 6)

        # QQ plot
        lambda <- inflation_factor(sumstat_list[["eur_discovery"]]$P)
        p <- lapply(sumstat_list, function(x) dplyr::select(x, CHROM, POS, P)) %>%
          qqtopr(legend_position = "bottom", legend_labels = c("Discovery", "Validation"), size = 0.75, legend_name = "Dataset") +
          annotate("text", x = 1, y = y_max, label = paste0("Lambda = ", round(lambda, 2)))
        ggsave(qq_save_to, p, width = 8, height = 4)
      }

      # Save genome-wide significant SNPs, which are validated in validation cohort.
      sigvar_save_to <- file.path(base_save_to, paste0("cmv_infection.", per_pheno, ".significant_variants.csv"))
      dis_sig_tab <- sumstat_list$eur_discovery %>% dplyr::filter(P < 5e-8)
      val_sig_tab <- sumstat_list$eur_validation %>% dplyr::filter(ID %in% dis_sig_tab$ID, P < 0.05)

      sumtab <- NULL
      if (nrow(dis_sig_tab) > 0 && nrow(val_sig_tab) > 0) {
        cat("[I]: Found", nrow(dis_sig_tab), "significant SNPs in discovery cohort and", nrow(val_sig_tab), "in validation cohort.\n")
        sumtab <- list(dplyr::mutate(dis_sig_tab, partition = "Discovery"), dplyr::mutate(val_sig_tab, partition = "Validation")) %>%
          Reduce(rbind, .) %>%
          dplyr::mutate(Measurement = per_ms)
      }
    }) %>%
    Reduce(rbind, .)
  }) %>%
  Reduce(rbind, .)

  fwrite(all_sumstat, sumstat_save_to)
} else {
  all_sumstat <- fread(sumstat_save_to)
}


# Target SNP
target_snp <- all_sumstat %>%
  dplyr::filter(ID == "rs7180928") %>%
  dplyr::select(chrom = CHROM, pos = POS, snp_id = ID, ref = REF, alt = ALT, a1 = A1, or = OR, log_or_se = `LOG(OR)_SE`, p = `P`, partition) %>%
  dplyr::mutate(lower = exp(log(or) - 1.96 * log_or_se), upper = exp(log(or) + 1.96 * log_or_se)) %>%
  dplyr::mutate(partition = factor(partition, levels = c("Discovery", "Validation")))

p_gwas <- ggplot(data = target_snp) +
  geom_pointrange(aes(x = log2(or), y = partition, xmin = log2(lower), xmax = log2(upper), color = partition)) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_y_discrete(limits=rev) +
  labs(x = "Log2(OR)", y = NULL, color = "Cohort") +
  theme_classic() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

save_to <- file.path(proj_dir, "outputs/associations/gwas/function_analysis/cmv_infection/cmv_cancidate_snp.eur.gwas.effect_size.pdf")
ggsave(save_to, plot = p_gwas, width = 3, height = 2)


# Check the locus
per_pheno <- "cmv_infection"
per_file <- paste0("gwas_sumstat.", per_pheno, ".glm.logistic.hybrid.gz")
sumstat <- file.path(proj_dir, "outputs/associations/cmv_infection/gwas/summary_statistic/cmv_infection/eur_discovery/", per_file) %>%
  fread() %>%
  dplyr::rename(c("CHROM" = "#CHROM")) %>%
  dplyr::filter(0.05 < A1_FREQ, A1_FREQ < 0.95, nchar(REF) == 1, nchar(ALT) == 1)

rp_save_to <- file.path(proj_dir, "outputs/associations/gwas/function_analysis/cmv_infection/cmv_cancidate_snp.eur.gwas.regionplot.pdf")
pdf(rp_save_to, width = 9, height = 6)
regionplot(sumstat, gene = "IREB2", annotate_with_vline = 5e-08, gene_padding = 2e5, annotate = 5e-09, unit_ratios = "1:8:2.5")
dev.off()

inflation_factor(sumstat$P)


# eQTL effect of the to SNP
selected_genes <- c("ADAMTS1", "ADGRG1", "ARHGEF12", "CRTAM", "FCRL6", "GFI1", "GNLY", "GSTM2", "GZMH", "ITGAL", "ITGB1", "KIR2DL3", "KLRD1", "LTBP3", "MXRA7", "NLRC5", "PATL2", "SPG11", "TAP1", "TAP2", "WNT10B", "ZFP28", "ZNF256", "ZNF418", "ZNF677", "ZNF814")

eqtl_sumstat <- file.path(proj_dir, "outputs/associations/cmv_infection/xqtl/summary_statistics/2000HIV.bulk_rnaseq.all.cmv_gwas_genes.chr15_75355207_80252213.cis_eqtl.nominal.txt") %>%
    fread(col.names = c("feature", "feature_chrom", "feature_start", "feature_end", "feature_strand", "n_SNPs", "dist", "rsid", "snp_chrom", "snp_pos", "snp_pos_xx", "p_value", "r_squared", "z_score", "std_err", "is_top"))
# for (tar_feature in c("WDR61", "CIB2", "DNAJA4", "ACSBG1", "CRABP1", "IREB2", "CHRNB4", "GOLGA6GP")) {
for (tar_feature in c("CHRNB4")) {
  per_eqtl_sumstat <- eqtl_sumstat %>%
    dplyr::filter(feature == tar_feature) %>%
    dplyr::select(CHROM = snp_chrom, POS = snp_pos, ID = rsid, P = p_value)

  rp_save_to <- file.path(proj_dir, "outputs/associations/cmv_infection/gwas/function_analysis/cmv_infection", paste0("cmv_cancidate_snp.eur.", tar_feature, "_eqtl_effect.regionplot.pdf"))
  pdf(rp_save_to, width = 9, height = 6)
  regionplot(list(per_eqtl_sumstat, sumstat), color = c("#E9C61DB1", "#A86A9D"), legend_labels = c(tar_feature, "CMV GWAS"), gene = tar_feature, annotate_with_vline = 5e-08, gene_padding = 5e5, annotate = 5e-09, unit_ratios = "1:8:2.5", protein_coding_only = TRUE)
  dev.off()
}

eqtl_sumstat %>% dplyr::filter(rsid == "rs7180928", p_value < 5e-2)


# Check the expression difference between genotypes
tar_snp <- "rs7180928"
genotype_input_cmd <- paste("zgrep -w -e CHROM -e", tar_snp, "/vol/projects/BIIM/2000HIV/RareVariants/outputs/variants/gsa/vcf/2000HIV.chr5_15.annotated.snps.maf05.vcf.gz")
genotype_tab <- fread(cmd = genotype_input_cmd, skip = "#CHROM") %>%
  dplyr::select(-dplyr::starts_with("FAM"), -c(QUAL, FILTER, FORMAT, INFO)) %>%
  dplyr::rename("CHROM" = "#CHROM") %>%
  dplyr::mutate(dplyr::across(-c(CHROM, ID, POS, REF, ALT), ~ stringr::str_extract(.x, "[01]\\|[01]"))) %>%
  tidyr::pivot_longer(-c(CHROM, ID, POS, REF, ALT), values_to = "Genotype01", names_to = "DonorID") %>%
  dplyr::mutate(GenotypeRA = dplyr::case_when(
    Genotype01 == "0|0" ~ paste0(REF, REF),
    Genotype01 %in% c("1|0", "0|1") ~ paste0(REF, ALT),
    Genotype01 == "1|1" ~ paste0(ALT, ALT)
  ))

expr_file <- "/vol/projects/BIIM/2000HIV/RareVariants/outputs/rna_seq/DESeq2/2000HIV.bulk_rnaseq.normalized_expr.csv"
expression_tab <- fread(expr_file)

for (tar_feature in c("CHRNB4", "GOLGA6GP", "CRABP1", "ADAMTS7", "SH2D7")) {
  plot_tab <- expression_tab %>%
    dplyr::filter(gene_symbol == tar_feature) %>%
    tidyr::pivot_longer(-c(`#chrom`, start, end, gene_symbol, ensembl_id, strand), names_to = "DonorID", values_to = "Expression") %>%
    dplyr::select(DonorID, Expression) %>%
    dplyr::inner_join(genotype_tab, by = "DonorID") %>%
    dplyr::select(DonorID, Genotype01, GenotypeRA, Expression)

  p <- ggplot() +
    geom_violin(data = plot_tab, aes(x = GenotypeRA, y = log10(Expression + 1), color = GenotypeRA)) +
    geom_point(data = plot_tab, aes(x = GenotypeRA, y = log10(Expression + 1), color = GenotypeRA), position = "jitter") +
    labs(x = "Genotype", y = "Log10(Expression + 1)") +
    theme_classic() +
    theme(legend.position = "none")
  boxp_save_to <- file.path(proj_dir, "outputs/associations/gwas/function_analysis/cmv_infection", paste0("cmv_cancidate_snp.eur.", tar_snp, ".", tar_feature, "_eqtl_effect.violin_plot.pdf"))
  ggsave(boxp_save_to, plot = p, width = 4, height = 3)
}


# Obtain QTLs
# QTL summary statistic path: /vol/projects/CIIM/meta_cQTL/out/2000HIV-EU-discovery
top_snp_id <- all_sumstat %>% dplyr::slice_min(P) %>% dplyr::pull(ID)
top_snp_chrom <- all_sumstat %>% dplyr::filter(ID == top_snp_id) %>% dplyr::pull(CHROM) %>% unique()
qtl_sumstat_path <- tibble::tribble(
  ~QTL_type, ~Trait_type, ~Summary_statistic_path,
  "protein_qtl", "Protein", "/vol/projects/BIIM/meta_cQTL/out/2000HIV-EU-discovery/proteins/mapping/main_nominal.tsv",
  "cytokine_qtl", "Cytokine", "/vol/projects/BIIM/meta_cQTL/out/2000HIV-EU-discovery/cytokines/mapping/main_nominal.tsv",
  "metabolite_qtl", "Metabolite", "/vol/projects/BIIM/meta_cQTL/out/2000HIV-EU-discovery/metabolites/mapping/main_nominal.tsv",
  # "expression_qtl", "Expression", file.path("/vol/projects/BIIM/meta_cQTL/out/2000HIV-EU-discovery/expression/mapping/main", paste0("chr", top_snp_chrom, ".tsv"))
)

base_save_to <- file.path(proj_dir, "outputs/associations/gwas/function_analysis/cmv_infection")
.tmp <- apply(qtl_sumstat_path, 1, function(vec) {
  # vec <- c(QTL_type = "protein_qtl", Trait_type = "Protein", Summary_statistic_path = "/vol/projects/BIIM/meta_cQTL/out/2000HIV-EU-discovery/proteins/mapping/main_nominal.tsv")
  # vec <- c(QTL_type = "cytokine_qtl", Trait_type = "Cytokine", Summary_statistic_path = "/vol/projects/BIIM/meta_cQTL/out/2000HIV-EU-discovery/cytokines/mapping/main_nominal.tsv")
  # vec <- c(QTL_type = "metabolite_qtl", Trait_type = "Metabolite", Summary_statistic_path = "/vol/projects/BIIM/meta_cQTL/out/2000HIV-EU-discovery/metabolites/mapping/main_nominal.tsv")
  # vec <- c(QTL_type = "expression_qtl", Trait_type = "Expression", Summary_statistic_path = file.path("/vol/projects/BIIM/meta_cQTL/out/2000HIV-EU-discovery/expression/mapping/main", paste0("chr", top_snp_chrom, ".tsv")))
  qtl_type <- vec["QTL_type"]
  trait_type <- vec["Trait_type"]
  sumstat_path <- vec["Summary_statistic_path"]

  cat("[I]: Working on", trait_type, "summary statistics ...", "\n")
  save_to <- file.path(base_save_to, paste0("candidate_snp.", qtl_type, ".tsv"))
  if (!file.exists(save_to)) {
    qtl_sstab <- fread(cmd = paste("(head -1", sumstat_path, "; grep -Fw", top_snp_id, sumstat_path, ")"))
    fwrite(qtl_sstab, save_to)
  } else {
    qtl_sstab <- fread(save_to)
  }

  tryCatch({
    plot_tab <- dplyr::select(qtl_sstab, SNP, trait = gene, beta, t_stat = `t-stat`, p_value = `p-value`) %>%
      dplyr::filter(p_value < 0.05) %>%
      tidyr::separate(col = SNP, into = c("pos_id", "rs_id"), sep = ";") %>%
      dplyr::mutate(p_value_adj = p.adjust(p_value)) %>%
      dplyr::mutate(se = beta / t_stat, upper = beta + 1.96 * se, lower = beta - 1.96 * se) %>%
      dplyr::mutate(trait = forcats::fct_reorder(trait, -beta)) %>%
      dplyr::mutate(`FDR < 0.05` = ifelse(p_value_adj < 0.05, "Yes", "No") %>% factor(levels = c("Yes", "No")))

    if (nrow(plot_tab) > 15) { plot_tab <- dplyr::slice_min(plot_tab, p_value, n = 15) }
    p_qtl <- ggplot(data = plot_tab) +
      geom_pointrange(aes(x = beta, y = trait, color = `FDR < 0.05`, xmin = lower, xmax = upper)) +
      geom_vline(xintercept = 0, linetype = "dotted") +
      scale_color_npg() +
      labs(x = "Beta", y = NULL) +
      facet_grid(rs_id ~ ., scales = "free_y", space = "free_y") +
      theme_classic() +
      theme(legend.position = "top")

    plot_height <- nrow(plot_tab) * 0.2 + 1
    save_to <- file.path(base_save_to, paste0("cmv_cancidate_snp.", qtl_type, ".effect_size.pdf"))
    ggsave(save_to, plot = p_qtl, width = 3.5, height = plot_height)
  }, error = function(e) cat(e$message, "\n"))

  invisible(NA)
})


# Boxplots to show QTL effects of the top SNP associated with CMV+/-
snp_dosage_tab <- file.path(proj_dir, "outputs/associations/cmv_infection/pQTL/rs7180928.txt") %>%
  fread() %>%
  tidyr::pivot_longer(-SNP, names_to = "DonorID", values_to = "Dosage") %>%
  tidyr::separate(SNP, into = c("Chrom", "Pos", "Ref", "Alt", "RsID"), sep = "[:;]") %>%
  dplyr::mutate(Genotype = dplyr::case_when(Dosage <= 0.5 ~ paste0(Ref, Ref), Dosage >= 1.5 ~ paste0(Alt, Alt), TRUE ~ paste0(Ref, Alt)))

# Cytokine of IL22 M. tb
cytokine_tab <- file.path(proj_dir, "outputs/associations/cmv_infection/pQTL/cytokines_filtered.tsv") %>%
  fread() %>%
  dplyr::filter(stringr::str_detect(Cytokine, "il22_mtb")) %>%
  tidyr::pivot_longer(-Cytokine, names_to = "DonorID", values_to = "Cytokine_level")
cytokine_plot_tab <- snp_dosage_tab %>% dplyr::inner_join(cytokine_tab, by = "DonorID") %>% dplyr::filter(!is.na(Cytokine_level), !is.na(Genotype))
cytokine_boxpoot <- ggplot(cytokine_plot_tab) +
  geom_violin(aes(x = Genotype, y = Cytokine_level + 1)) +
  geom_boxplot(aes(x = Genotype, y = Cytokine_level + 1), width = 0.3) +
  geom_point(aes(x = Genotype, y = Cytokine_level, color = Genotype), position = position_jitter(width = 0.075), alpha = 0.5) +
  scale_color_npg() +
  scale_y_log10() +
  labs(x = "Genotype", y = "Cytokine levels") +
  theme_classic() +
  theme(legend.position = "top")
cytokine_boxpoot_saveto <- file.path(proj_dir, "outputs/associations/cmv_infection/gwas/function_analysis/cytokines_boxplot.pdf")
ggsave(cytokine_boxpoot_saveto, cytokine_boxpoot, width = 3, height = 4)


protein_tab <- file.path(proj_dir, "outputs/associations/cmv_infection/pQTL/proteins_filtered.tsv") %>%
  fread() %>%
  dplyr::filter(Protein == "KIR2DS4_Oncology_II") %>%
  tidyr::pivot_longer(-Protein, names_to = "DonorID", values_to = "Protein_level")
protein_plot_tab <- snp_dosage_tab %>% dplyr::inner_join(protein_tab, by = "DonorID") %>% dplyr::filter(!is.na(Genotype), !is.na(Protein_level))
protein_boxpoot <- ggplot(protein_plot_tab) +
  geom_violin(aes(x = Genotype, y = Protein_level)) +
  geom_boxplot(aes(x = Genotype, y = Protein_level), width = 0.3) +
  geom_point(aes(x = Genotype, y = Protein_level, color = Genotype), position = position_jitter(width = 0.075), alpha = 0.5) +
  scale_color_npg() +
  labs(x = "Genotype", y = "Protein adundance") +
  theme_classic() +
  theme(legend.position = "top")
protein_boxpoot_saveto <- file.path(proj_dir, "outputs/associations/cmv_infection/gwas/function_analysis/protein_boxplot.pdf")
ggsave(protein_boxpoot_saveto, protein_boxpoot, width = 3, height = 4)
