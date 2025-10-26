#!/usr/bin/env Rscript

# general
library(tidyverse)
library(data.table)

# plots
library(ggplot2)
library(ComplexHeatmap)
library(circlize)


#mofa <- readRDS("/vol/projects/CIIM/2000HIV/cQTL/mofa/out/model_corrected_scaled.rds")
# Extract MOFA factor scores
#factors_df <- as.data.frame(get_factors(mofa, factors = factors)) %>%
#  tibble::rownames_to_column(var = "Record.Id")
# Clean up factor column names by removing 'single_group.' prefix
#colnames(factors_df) <- sub("^single_group\\.", "", colnames(factors_df))

### Instead of above, getting factor values later from the clin_df file from Javi

outdir <- "/vol/projects/nvanunen/analysis/2000HIV_CMV/HIV_comorb/mofa_cor"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
setwd(outdir)

# ---- Correlate CMV serology measures with MOFA factor scores ----
# Load phenotype data
pheno <- fread("/vol/projects/BIIM/2000HIV/Phenotype/Phenotype_2000HIV_all_01.tsv") %>%
  dplyr::select(Record.Id, CMV_IgG_Serology, CMV_IgG_IU.mL, CD8_LATEST) %>%
  dplyr::rename(CMV_cat = CMV_IgG_Serology, CMV_IgG = CMV_IgG_IU.mL) %>%
  dplyr::mutate( # Clean CMV_IgG column: remove leading "<" and coerce to numeric (non-numeric to NA)
    CMV_IgG = as.numeric(sub("^<", "", CMV_IgG))
  )

clin_df <- fread("/vol/projects/CIIM/2000HIV/cQTL/mofa/out/clin_cor/clin_df_corrected_scaled.csv") %>% 
  dplyr::rename(Record.Id = "V1") %>% select(Record.Id, "MH_TR_RESP_D%COPD", matches("Factor"))

# Merge phenotype and factors
merge_df <- pheno %>% dplyr::inner_join(clin_df, by = "Record.Id")
dim(merge_df) # 996 26

# remove samples with missing CMV data
merge_df <- merge_df %>% filter(!is.na(CMV_cat))
dim(merge_df) # 994 26

# set CMV_IgG to NA for people who don't have CMV, so we only do the continuous analysis with CMV+
merge_df <- merge_df %>% mutate(CMV_IgG = ifelse(CMV_cat == 0, NA, CMV_IgG))
merge_df <- merge_df %>% mutate(log_CMV_IgG = log10(CMV_IgG))
dim(merge_df) # 994 27

# count amount of people with CMV+ and CMV-
merge_df %>% group_by(CMV_cat) %>% summarise(n = n()) # 63 vs 931

 # Compute correlations and linear regression for each factor
corr_results <- data.frame(
  Factor              = character(),
  b_CMV               = numeric(), p_CMV               = numeric(),
  b_CMV_adjCD8        = numeric(), p_CMV_adjCD8        = numeric(),
  b_CMV_IgG           = numeric(), p_CMV_IgG           = numeric(),
  b_CMV_IgG_adjCD8    = numeric(), p_CMV_IgG_adjCD8    = numeric(),
  stringsAsFactors = FALSE
)

cmv_cat <- merge_df$CMV_cat
cmv_igg <- merge_df$CMV_IgG

factors <- grep("Factor", names(merge_df), value = TRUE)

fac <- "Factor1"
for (fac in factors) {
  x <- merge_df[[fac]]
  # Spearman correlation for categorical (CMV_serology encoded as 0/1)
  # Linear regression for categorical serostatus
  lm_cat <- lm(x ~ CMV_cat, data = merge_df)
  beta_cat <- coef(lm_cat)["CMV_cat"]
  p_lm_cat <- summary(lm_cat)$coefficients["CMV_cat", "Pr(>|t|)"]

  # CMV (0/1) effect *after* adjusting for CD8
  lm_cat_adj <- lm(x ~ CMV_cat + CD8_LATEST, data = merge_df)
  beta_cat_adj <- coef(lm_cat_adj)["CMV_cat"]
  p_cat_adj    <- summary(lm_cat_adj)$coefficients["CMV_cat", "Pr(>|t|)"]

  # --- CMV IgG analysis restricted to CMV‑positive donors only ---
  pos_idx <- which(!is.na(merge_df$CMV_IgG) & !is.na(x))  # CMV_IgG is NA when CMV- so this already selects CMV+
  # Linear regression for continuous IgG (log‑transformed) in positives
  tmp_df <- data.frame(log_CMV_IgG = merge_df$log_CMV_IgG[pos_idx],
                       CD8_LATEST = merge_df$CD8_LATEST[pos_idx],
                       fac_score   = x[pos_idx])
  lm_IgG   <- lm(log_CMV_IgG ~ fac_score, data = tmp_df)
  beta_IgG <- coef(lm_IgG)["fac_score"]
  p_lm_IgG <- summary(lm_IgG)$coefficients["fac_score", "Pr(>|t|)"]

  # IgG effect after adjusting for CD
  lm_IgG_adj <- lm(fac_score ~ log_CMV_IgG + CD8_LATEST, data = tmp_df)
  beta_IgG_adj <- coef(lm_IgG_adj)["log_CMV_IgG"]
  p_IgG_adj    <- summary(lm_IgG_adj)$coefficients["log_CMV_IgG", "Pr(>|t|)"]

  corr_results[nrow(corr_results) + 1, ] <- list(
    fac,
    beta_cat,       p_lm_cat,
    beta_cat_adj,   p_cat_adj,
    beta_IgG,        p_lm_IgG,
    beta_IgG_adj,    p_IgG_adj
  )
}

 # Adjust p-values for multiple testing (BH)
corr_results <- corr_results %>%
  mutate(
    padj_CMV            = p.adjust(p_CMV,            method = "BH"),
    padj_CMV_adjCD8     = p.adjust(p_CMV_adjCD8,     method = "BH"),
    padj_CMV_IgG        = p.adjust(p_CMV_IgG,        method = "BH"),
    padj_CMV_IgG_adjCD8 = p.adjust(p_CMV_IgG_adjCD8, method = "BH")
  )
write.table(corr_results, file.path(outdir, "mofa_CMV_factor_correlations.tsv"), row.names=FALSE, sep = "\t")
message("Saved adjusted MOFA-CMV correlation results to ", file.path(outdir, "mofa_CMV_factor_correlations.tsv"))

# for downstream analysis
signif_cat_factors <- corr_results %>% filter(padj_CMV < 0.05) %>% pull(Factor)
signif_igg_factors <- corr_results %>% filter(padj_CMV_IgG < 0.05) %>% pull(Factor)

# ---- Downsampling analysis ----
nsim <- 100
set.seed(123)
downsample_results <- data.frame(
  Factor = character(),
  Simulation = integer(),
  Beta = numeric(),
  Pval = numeric(),
  stringsAsFactors = FALSE
)

cmv_cat <- merge_df$CMV_cat
for (fac in factors) {
  x <- merge_df[[fac]]
  pos_idx <- which(cmv_cat == 1)
  neg_idx <- which(cmv_cat == 0)
  n_min <- min(length(pos_idx), length(neg_idx))
  for (s in 1:nsim) {
    pos_idx_sub <- sample(pos_idx, n_min, replace = FALSE)
    neg_idx_sub <- sample(neg_idx, n_min, replace = FALSE)
    combined_idx <- c(pos_idx_sub, neg_idx_sub)
    x_sub <- x[combined_idx]
    IgG_sub <- cmv_cat[combined_idx]
    lm_ds <- lm(IgG_sub ~ x_sub)
    beta_ds <- coef(lm_ds)["x_sub"]
    p_ds    <- summary(lm_ds)$coefficients["x_sub", "Pr(>|t|)"]
    downsample_results[nrow(downsample_results) + 1, ] <- list(
      fac, s, beta_ds, p_ds
    )
  }
}

# Summarize proportion of significant p-values per factor
downsample_summary <- downsample_results %>%
  group_by(Factor) %>%
  summarise(prop_significant = mean(Pval < 0.05)) %>%
  ungroup()

write.table(downsample_results, file.path(outdir, "mofa_CMV_factor_downsampling_results.tsv"), row.names = FALSE, sep = "\t")
write.table(downsample_summary, file.path(outdir, "mofa_CMV_factor_downsampling_summary.tsv"), row.names = FALSE, sep = "\t")
message("Saved downsampling analysis results to ", file.path(outdir, "mofa_CMV_factor_downsampling_results.tsv"))
message("Saved downsampling summary to ", file.path(outdir, "mofa_CMV_factor_downsampling_summary.tsv"))


# ---- CD8 effect before vs after adjusting for CMV (linear models only) ----
after_cmv <- data.frame(
  Factor       = character(),
  CD8_b  = numeric(), CD8_b_CMVadj = numeric(), CD8_b_CMVnadj = numeric(), 
  CD8_p     = numeric(), CD8_p_CMVadj    = numeric(), CD8_p_CMVnadj = numeric(),
  stringsAsFactors = FALSE
)

for (fac in factors) {
  cc  <- complete.cases(merge_df[[fac]], merge_df$CD8_LATEST, cmv_cat)
  df  <- merge_df[cc, ]

  # CD8 effect without CMV
  lm_before  <- lm(as.formula(paste(fac, "~ CD8_LATEST")), data = df)
  beta_before <- coef(lm_before)["CD8_LATEST"]
  p_before    <- summary(lm_before)$coefficients["CD8_LATEST", "Pr(>|t|)"]

  # CD8 effect after adjusting for CMV
  lm_after   <- lm(as.formula(paste(fac, "~ CD8_LATEST + CMV_cat")), data = df)
  beta_after <- coef(lm_after)["CD8_LATEST"]
  p_after    <- summary(lm_after)$coefficients["CD8_LATEST", "Pr(>|t|)"]

  # CD8 effect after adjusting for CMV IgG
  lm_after_igg   <- lm(as.formula(paste(fac, "~ CD8_LATEST + CMV_IgG")), data = df)
  beta_after_igg <- coef(lm_after_igg)["CD8_LATEST"]
  p_after_igg    <- summary(lm_after_igg)$coefficients["CD8_LATEST", "Pr(>|t|)"]

  after_cmv[nrow(after_cmv) + 1, ] <- list(
    fac, beta_before, beta_after, beta_after_igg, p_before, p_after, p_after_igg
  )
}

write.table(after_cmv, file.path(outdir, "mofa_CD8_lm_before_after_CMV.tsv"),
            row.names = FALSE, sep = "\t")
message("Saved CD8 linear models before/after CMV adjustment to ",
        file.path(outdir, "mofa_CD8_lm_before_after_CMV.tsv"))

#print(sig_all4)
#"Factor9"  "Factor11" "Factor15"

# Boxplots for CMV_cat
for (fac in signif_cat_factors) {
  # Use only non-missing CMV_cat for plotting
  plot_df_cat <- merge_df %>% filter(!is.na(CMV_cat))
  # Create a factor with labels for plotting
  plot_df_cat <- plot_df_cat %>%
    dplyr::mutate(CMV_label = factor(CMV_cat, levels = c(0, 1), labels = c("CMV-", "CMV+")))
  # Compute t-test for CMV_cat
  grp0 <- plot_df_cat %>% filter(CMV_cat == 0) %>% pull(fac)
  grp1 <- plot_df_cat %>% filter(CMV_cat == 1) %>% pull(fac)
  ttest <- t.test(grp1, grp0)
  p_ttest <- signif(ttest$p.value, 2)

  # get downsampling significance proportion for this factor
  ds_prop <- downsample_summary %>% filter(Factor == fac) %>% pull(prop_significant)
  ds_label <- paste0("Downsample sig: ", signif(ds_prop * 100, 2), "%")

  p_cat <- ggplot(plot_df_cat, aes(x = CMV_label, y = .data[[fac]], fill = CMV_label)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA) +
    geom_jitter(aes(color = CMV_label), width = 0.2, alpha = 0.6) +
    scale_fill_manual(
      values = c("CMV-" = "#e1c741", "CMV+" = "#a970a0"),
      name = "CMV"
    ) +
    scale_color_manual(
      values = c("CMV-" = "#e1c741", "CMV+" = "#a970a0"),
      guide = FALSE
    ) +
    annotate("text", x = 1.5, y = max(plot_df_cat[[fac]], na.rm = TRUE) * 1.05,
             label = paste0("p = ", p_ttest), vjust = 0) +
    annotate("text", x = 1.5, y = max(plot_df_cat[[fac]], na.rm = TRUE) * 0.9,
             label = ds_label, vjust = 0) +
    labs(
      title = paste(fac, "by CMV serostatus"),
      x = "CMV",
      y = fac
    ) +
    theme_classic() +
    # hide legend
    theme(legend.position = "none")
  ggsave(file.path(outdir, paste0("plot_", fac, "_by_CMVcat.pdf")),
         plot = p_cat, width = 3, height = 4)
  message(file.path(outdir, paste0("plot_", fac, "_by_CMVcat.pdf")))
}

# Boxplots for residuals after regressing out CD8_LATEST
for (fac in signif_cat_factors) {
  # compute residuals
  cc <- complete.cases(merge_df[[fac]], merge_df$CD8_LATEST, cmv_cat)
  resid_df <- lm(as.formula(paste(fac, "~ CD8_LATEST")), data = merge_df[cc, ], na.action = na.exclude)
  resids <- resid(resid_df)
  plot_df_resid <- merge_df[cc, ] %>% 
    mutate(resid = resids,
           CMV_label = factor(CMV_cat, levels = c(0,1), labels = c("CMV-","CMV+")))
  # compute t-test on residuals
  grp0 <- plot_df_resid %>% filter(CMV_cat == 0) %>% pull(resid)
  grp1 <- plot_df_resid %>% filter(CMV_cat == 1) %>% pull(resid)
  ttest_resid <- t.test(grp1, grp0)
  p_resid <- signif(ttest_resid$p.value, 2)
  p_resid_plot <- ggplot(plot_df_resid, aes(x = CMV_label, y = resid, fill = CMV_label)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA) +
    geom_jitter(aes(color = CMV_label), width = 0.2, alpha = 0.6) +
    scale_fill_manual(values = c("CMV-"="#e1c741","CMV+"="#a970a0"), guide=FALSE) +
    scale_color_manual(values = c("CMV-"="#e1c741","CMV+"="#a970a0"), guide=FALSE) +
    annotate("text", x = 1.5, y = max(plot_df_resid$resid, na.rm=TRUE)*1.05,
             label = paste0("p = ", p_resid), vjust = 0) +
    labs(title = paste(fac, "residuals by CMV (adjusted CD8)"),
         x = "CMV", y = paste0(fac, "_resid")) +
    theme_classic() +
    theme(legend.position="none")
  ggsave(file.path(outdir, paste0("plot_", fac, "_resid_by_CMVcat.pdf")),
         plot = p_resid_plot, width = 4, height = 4)
  message(file.path(outdir, paste0("plot_", fac, "_resid_by_CMVcat.pdf")))
}
message("Saved residual boxplots for CMV by factor adjusted for CD8.")

# Scatter plots for CMV_IgG
for (fac in signif_igg_factors) {
  plot_df_IgG <- merge_df %>% filter(cmv_cat == 1 & !is.na(CMV_IgG))
  # Linear model for log_CMV_IgG vs. factor score
  lm_mod <- lm(as.formula(paste(fac, "~ log_CMV_IgG")), data = plot_df_IgG)
  lm_sum <- summary(lm_mod)
  p_lm <- signif(lm_sum$coefficients["log_CMV_IgG", "Pr(>|t|)"], 2)
  r2 <- signif(lm_sum$r.squared, 2)

  p_IgG <- ggplot(plot_df_IgG, aes_string(x = "log_CMV_IgG", y = fac)) +
    geom_point(alpha = 0.6, colour = "#a970a0") +
    geom_smooth(method = "lm", se = TRUE, colour = "#704369", linewidth = 1.1) +
    annotate(
      "text",
      x = -Inf, y = max(plot_df_IgG[[fac]], na.rm = TRUE) * 0.95,
      label = paste0("p = ", p_lm, "\nR² = ", r2),
      hjust = -0.3, vjust = 1
    ) +
    labs(title = paste(fac, "vs log10(CMV IgG)"),
         x = "log10(CMV IgG)", y = fac) +
    theme_classic()
  ggsave(file.path(outdir, paste0("plot_", fac, "_vs_CMV_IgG.pdf")),
         plot = p_IgG, width = 4.5, height = 4)
  message(file.path(outdir, paste0("plot_", fac, "_vs_CMV_IgG.pdf")))
}


 # 1. Give every data-frame unique, self-describing column names
corr_tbl <- corr_results #%>%
  #select(Factor,
  #       b_CMV,           p_CMV,
  #       b_CMV_adjCD8,    p_CMV_adjCD8,
  #       b_CMV_IgG,       p_CMV_IgG,
  #       b_CMV_IgG_adjCD8, p_CMV_IgG_adjCD8)

ds_tbl     <- downsample_summary %>% 
  rename(downsample_pct_signif = prop_significant)


cmv_adj_tbl <- after_cmv 

# 2. Merge everything by Factor
master_tbl <- list(corr_tbl, ds_tbl, cmv_adj_tbl) |>
  reduce(full_join, by = "Factor") |>
  mutate(Factor_num = as.integer(sub("Factor", "", Factor))) |>
  arrange(Factor_num) |>
  select(-Factor_num)


# 3. Format numbers: round to 2 decimals if |x| ≥ 0.01, otherwise scientific (2 sig figs)
master_tbl_fmt <- master_tbl %>%
  mutate(across(
    where(is.numeric),
    ~ ifelse(
        is.na(.),
        NA,
        ifelse(abs(.) < 0.01 & . != 0,
               formatC(., format = "e", digits = 2),
               formatC(round(., 2), format = "f", digits = 2))
      )
  ))

write.table(master_tbl_fmt,
            file.path(outdir, "mofa_factor_master_summary.tsv"),
            row.names = FALSE, sep = "\t", quote = FALSE)
file.path(outdir, "mofa_factor_master_summary.tsv")

 # ---- Directional p-value table (sign of beta determines sign of p) ----
#   start from the numeric master_tbl to keep beta signs numeric
master_tbl_dirp <- master_tbl %>%
  mutate(
    p_CMV_dir            = ifelse(b_CMV           >= 0, p_CMV,           -p_CMV),
    p_CMV_adjCD8_dir     = ifelse(b_CMV_adjCD8    >= 0, p_CMV_adjCD8,    -p_CMV_adjCD8),
    p_CMV_IgG_dir        = ifelse(b_CMV_IgG       >= 0, p_CMV_IgG,       -p_CMV_IgG),
    p_CMV_IgG_adjCD8_dir = ifelse(b_CMV_IgG_adjCD8>= 0, p_CMV_IgG_adjCD8,-p_CMV_IgG_adjCD8),
    cd8adj_p_before_dir  = ifelse(CD8_b           >= 0, CD8_p,           -CD8_p),
    cd8adj_p_after_dir   = ifelse(CD8_b_CMVadj    >= 0, CD8_p_CMVadj,    -CD8_p_CMVadj),
    cmvadj_p_before_dir  = ifelse(CD8_b           >= 0, CD8_p,           -CD8_p),
    cmvadj_p_after_dir   = ifelse(CD8_b_CMVnadj   >= 0, CD8_p_CMVnadj,   -CD8_p_CMVnadj)
  ) %>%
  select(Factor,
         p_CMV_dir, p_CMV_adjCD8_dir, p_CMV_IgG_dir, p_CMV_IgG_adjCD8_dir,
         cd8adj_p_before_dir, cd8adj_p_after_dir,
         cmvadj_p_before_dir, cmvadj_p_after_dir,
         downsample_pct_signif) %>% # keep proportion column as-is
  mutate(across(
    where(is.numeric),
    ~ ifelse(
        is.na(.),
        NA,
        ifelse(abs(.) < 0.01 & . != 0,
               formatC(., format = "e", digits = 2),
               formatC(round(., 2), format = "f", digits = 2))
      )
  ))# %>%
  #dplyr::filter(!is.na(downsample_pct_signif))

# Save directional p-value table (no further rounding, scientific as needed)
write.table(master_tbl_dirp,
            file.path(outdir, "mofa_factor_master_summary_dirp.tsv"),
            row.names = FALSE, sep = "\t", quote = FALSE)

message("Saved directional p-value summary to ",
        file.path(outdir, "mofa_factor_master_summary_dirp.tsv"))


## HEATMAP ## 

master_tbl_dirp_forplot <- master_tbl_dirp %>%                                # ensure correct factor order
  mutate(Factor_num = as.integer(sub("Factor", "", Factor))) %>% 
  arrange(Factor_num) %>% 
  select(-Factor_num)

heat_cols <- c("p_CMV_dir", "p_CMV_adjCD8_dir",
               "cd8adj_p_before_dir", "cd8adj_p_after_dir",
               "cmvadj_p_before_dir", "cmvadj_p_after_dir",
               "p_CMV_IgG_dir", "p_CMV_IgG_adjCD8_dir")

heat_mat <- master_tbl_dirp_forplot %>% 
  select(all_of(heat_cols)) %>% 
  mutate(across(everything(), as.numeric)) %>% 
  as.matrix()

rownames(heat_mat) <- master_tbl_dirp_forplot$Factor
colnames(heat_mat) <- c("CMV_cat",
                        "CMV_cat_adjCD8",
                        "CD8_before", "CD8_after",
                        "CMVadj_before", "CMVadj_after",
                        "CMV_IgG",
                        "CMV_IgG_adjCD8")

# Signed –log10(p)
heat_vals <- -log10(abs(heat_mat)) * sign(heat_mat)

# Bar-plot values (NA → 0)
ds_vals <- as.numeric(master_tbl_dirp_forplot$downsample_pct_signif)
ds_vals[is.na(ds_vals)] <- 0

# ----------------------------------
# 2. Build annotation & palette  ----
# ----------------------------------
max_abs <- ceiling(max(abs(heat_vals), na.rm = TRUE))
col_fun <- colorRamp2(c(-max_abs, 0, max_abs), c("steelblue", "white", "firebrick"))

right_anno <- rowAnnotation(
  Downsample = anno_barplot(ds_vals,
                            gp = gpar(fill = "grey40"),
                            border = FALSE,
                            axis_param = list(side = "top")),
  annotation_name_side = "top"
)

# ----------------------------------
# 3. Draw and save  -----------------
# ----------------------------------
ht <- Heatmap(heat_vals,
              name = "-log10(p)·sign",
              col = col_fun,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              right_annotation = right_anno,
              na_col = "grey90",
              heatmap_legend_param = list(
                at = c(-max_abs, 0, max_abs),
                labels = c(paste0("-", max_abs), "0", paste0("+", max_abs))
              ))

pdf(file.path(outdir, "heatmap_directional_pvalues_with_downsample.pdf"), width = 8, height = 8)
draw(ht)
dev.off()
file.path(outdir, "heatmap_directional_pvalues_with_downsample.pdf")




# ---- Association of CD8_LATEST with CMV_cat and CMV_IgG ----
cmv_cat <- merge_df$CMV_cat  # reuse 0/1 coding
cmv_igg <- merge_df$CMV_IgG  # reuse 0/1 coding

# CD8 vs CMV serostatus
lm_cd8_cat <- lm(CD8_LATEST ~ cmv_cat, data = merge_df)
beta_cat_cd8 <- coef(lm_cd8_cat)["cmv_cat"]
p_cat_cd8    <- summary(lm_cd8_cat)$coefficients["cmv_cat", "Pr(>|t|)"]

# CD8 vs continuous CMV IgG (CMV+ only)
IgG_df <- merge_df %>% filter(!is.na(CMV_IgG) & !is.na(CD8_LATEST))
lm_cd8_IgG <- lm(CD8_LATEST ~ log_CMV_IgG, data = IgG_df)
beta_IgG_cd8 <- coef(lm_cd8_IgG)["log_CMV_IgG"]
p_IgG_cd8    <- summary(lm_cd8_IgG)$coefficients["log_CMV_IgG", "Pr(>|t|)"]

cd8_cmv_tbl <- data.frame(
  Measure = c("CMV_cat", "log_CMV_IgG"),
  Beta    = c(beta_cat_cd8, beta_IgG_cd8),
  Pval    = c(p_cat_cd8,    p_IgG_cd8)
)

write.table(cd8_cmv_tbl,
            file.path(outdir, "CD8_vs_CMV_measures.tsv"),
            row.names = FALSE, sep = "\t")

message("Saved CD8 associations with CMV_cat and log(CMV_IgG) to ",
        file.path(outdir, "CD8_vs_CMV_measures.tsv"))


# ---- Plots: CD8_LATEST vs CMV measures ----

# 1) Boxplot of CD8 by CMV serostatus (0/1)
plot_df_cat <- merge_df %>%
  filter(!is.na(CMV_cat)) %>%
  mutate(CMV_label = factor(CMV_cat, levels = c(0, 1), labels = c("CMV-", "CMV+")))

# p‑value (two‑sided t‑test)
grp0 <- plot_df_cat %>% filter(CMV_cat == 0) %>% pull(CD8_LATEST)
grp1 <- plot_df_cat %>% filter(CMV_cat == 1) %>% pull(CD8_LATEST)
p_val_cat <- signif(t.test(grp1, grp0)$p.value, 2)

p_cd8_cat <- ggplot(plot_df_cat, aes(x = CMV_label, y = CD8_LATEST, fill = CMV_label)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, aes(color = CMV_label)) +
  scale_fill_manual(values = c("CMV-" = "#e1c741", "CMV+" = "#a970a0"), guide = FALSE) +
  scale_color_manual(values = c("CMV-" = "#e1c741", "CMV+" = "#a970a0"), guide = FALSE) +
  annotate("text", x = 1.5, y = max(plot_df_cat$CD8_LATEST, na.rm = TRUE) * 1.05,
           label = paste0("p = ", p_val_cat), vjust = 0) +
  labs(title = "CD8 counts by CMV serostatus", x = "CMV", y = "CD8_LATEST") +
  theme_classic()

ggsave(file.path(outdir, "plot_CD8_by_CMVcat.pdf"), plot = p_cd8_cat, width = 4, height = 4)
message(file.path(outdir, "plot_CD8_by_CMVcat.pdf"))

# 2) Scatter plot CD8 vs log10(CMV IgG) for CMV+ donors
plot_df_IgG <- merge_df %>% filter(!is.na(CMV_IgG))

lm_cd8_IgG <- lm(CD8_LATEST ~ log_CMV_IgG, data = plot_df_IgG)
p_val <- signif(summary(lm_cd8_IgG)$coefficients["log_CMV_IgG", "Pr(>|t|)"], 2)
r2_val <- signif(summary(lm_cd8_IgG)$r.squared, 2)

p_cd8_IgG <- ggplot(plot_df_IgG, aes(x = log_CMV_IgG, y = CD8_LATEST)) +
  geom_point(alpha = 0.6, colour = "#377eb8") +
  geom_smooth(method = "lm", se = TRUE, colour = "#e41a1c", linewidth = 1.1) +
  annotate("text",
           x = min(plot_df_IgG$log_CMV_IgG, na.rm = TRUE),
           y = max(plot_df_IgG$CD8_LATEST, na.rm = TRUE) * 0.95,
           hjust = 0, vjust = 1,
           label = paste0("p = ", p_val, "\nR² = ", r2_val)) +
  labs(title = "CD8 counts vs log10(CMV IgG)", x = "log10(CMV IgG)", y = "CD8_LATEST") +
  theme_classic()

ggsave(file.path(outdir, "plot_CD8_vs_CMV_IgG.pdf"),
       plot = p_cd8_IgG, width = 5, height = 4)
message(file.path(outdir, "plot_CD8_vs_CMV_IgG.pdf"))
