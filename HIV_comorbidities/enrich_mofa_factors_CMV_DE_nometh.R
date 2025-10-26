#!/usr/bin/env Rscript

# general
library(tidyverse)
library(data.table)

# plots
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

#BiocManager::install("MOFA2")
library(MOFA2)
#library(reticulate)
#use_condaenv("mofa_env", required = TRUE)
#py_config() # Check if Python is detected

# gene annot
library(AnnotationDbi)
library(org.Hs.eg.db)

outdir <- "/vol/projects/nvanunen/analysis/2000HIV_CMV/HIV_comorb/mofa_enrich_nometh"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
setwd(outdir)

mofa <- readRDS("/vol/projects/CIIM/2000HIV/cQTL/mofa/out/model_corrected_scaled.rds")

# ----------------------------------------
# Load differential DE lists from Excel
library(readxl)
xls_path <- "/vol/projects/nvanunen/analysis/2000HIV_CMV/HIV_comorb/202409_suppTables.xlsx"
# List all sheets (run once to inspect and adapt sheet names)
sheets <- excel_sheets(xls_path)
print(sheets)  # adjust sheet names below based on actual names


available_views <- setdiff(views_names(mofa), "meth")   # drop methylation view

# Nice labels for plots
view_label_map <- c(
  cyt   = "Cytokines",
  gex   = "Gene expression",
  metab = "Metabolites",
  prot  = "Proteins"
)
available_views
#[1] "cyt"   "gex"   "metab" "prot" 

# Read DE results from discovery sheets across omics
de_cyt_df     <- read_excel(xls_path, sheet = "S2_cytokines_discovery", skip = 2) %>% filter(!is.na(padj) & padj < 0.05) %>% dplyr::rename(feature = "...1")
de_gex_df     <- read_excel(xls_path, sheet = "S6_geneExpressions_discovery", skip = 2) %>% filter(!is.na(padj) & padj < 0.05) %>% dplyr::rename(feature = "gene", effectSize = "log2FoldChange")
de_metab_df   <- read_excel(xls_path, sheet = "S11_metaboliteLevels_discovery", skip = 2) %>% filter(!is.na(padj) & padj < 0.05) %>% dplyr::rename(feature = "...1")
de_prot_df    <- read_excel(xls_path, sheet = "S9_proteinLevels_discovery", skip = 2) %>% filter(!is.na(padj) & padj < 0.05) %>% dplyr::rename(feature = "...1")

# check if each discovery has the expected columns
expected_cols <- c("feature", "padj", "effectSize")
check_cols <- function(df, name) {
  missing_cols <- setdiff(expected_cols, colnames(df))
  if (length(missing_cols) > 0) {
    warning(paste("Missing columns in", name, ":", paste(missing_cols, collapse = ", ")))
  }
}
check_cols(de_cyt_df, "Cytokines")
check_cols(de_gex_df, "Gene Expression")
check_cols(de_metab_df, "Metabolites")
check_cols(de_prot_df, "Proteins")



# Enrichment test function (hypergeometric)
enrich_test <- function(feature_list, de_list, universe) {
  # restrict DE list to universe
  de_univ <- intersect(de_list, universe)
  k <- sum(feature_list %in% de_univ)
  n <- length(de_univ)
  N <- length(feature_list)
  M <- length(universe)
  # hypergeom p-value P(X >= k)
  phyper(k - 1, n, M - n, N, lower.tail = FALSE)
}

# print count of DE features in each discovery
cat(
  "DE Cytokines:", nrow(de_cyt_df), # 10
  "| DE Genes:", nrow(de_gex_df), # 1442 
  "| DE Metabolites:", nrow(de_metab_df), # 5
  "| DE Proteins:", nrow(de_prot_df),  # 38
  "\n" 
)

# Only test these factor-direction combinations
factor_dir_map <- c(Factor6="pos", Factor8="neg", Factor9="neg", Factor11="pos", Factor20="pos")

factor_labels <- c(
  Factor6 = "Plaque",
  Factor8 = "Hypertension & Myocardial Infarction",
  Factor9 = "CD8 latest",
  Factor11 = "COPD", # & CD8 latest",
  Factor20 = "Rapid Progressors" # & neg with CD8 latest
)

role_map <- c(pos = "risk", neg = "protective")

view_name <- "metab"
de_df <- de_metab_df 

view_name <- "gex"
de_df <- de_gex_df 
perc_thresholds <- seq(10, 90, by=10)
factor_dir <- factor_dir_map

# Plot enrichment per factor for each view
plot_factor_enrichment <- function(view_name, de_df, factor_dir, outdir, perc_thresholds, return_df=FALSE) {
  message(view_name)
  # Retrieve MOFA weights and universe
  weights <- get_weights(mofa, views = view_name, factors = paste0("Factor", 1:mofa@dimensions$K))[[view_name]]
  feats_universe <- rownames(weights)
  if (view_name == "gex") {
    orig_names <- feats_universe
    symbol_map <- mapIds(org.Hs.eg.db, keys = orig_names, column="SYMBOL", keytype="ENSEMBL", multiVals="first")
    is_ens <- grepl("^ENSG", orig_names); mapped <- !is.na(symbol_map)
    keep <- (!is_ens) | (is_ens & mapped)
    weights <- weights[keep, , drop=FALSE]
    rownames(weights) <- ifelse(is_ens[keep], symbol_map[keep], rownames(weights))
    feats_universe <- rownames(weights)
  } else if (view_name == "metab") {
    metabo_idx <- fread("/vol/projects/CIIM/2000HIV/Metabolites/metabo_index.tsv")
    map_formula <- setNames(metabo_idx$formula, metabo_idx$name)
    orig_metab <- feats_universe
    valid <- orig_metab %in% names(map_formula)
    weights <- weights[valid, , drop=FALSE]
    rownames(weights) <- map_formula[orig_metab[valid]]
    feats_universe <- rownames(weights)
  } else if (view_name == "prot") {
    clean_names <- gsub("_[^_]+(?:_II)?$", "", feats_universe)
    rownames(weights) <- clean_names
    feats_universe <- clean_names
  }

  # DE lists
  de_pos <- de_df %>% dplyr::filter(effectSize > 0) %>% pull(feature)
  de_neg <- de_df %>% dplyr::filter(effectSize < 0) %>% pull(feature)
  de_all <- unique(c(de_pos, de_neg))

  df_plot_list <- list()

  # Loop over all MOFA factors, defaulting to factor_dir_map when available
  #all_factors <- paste0("Factor", 1:mofa@dimensions$K)
  all_factors <- names(factor_dir_map)
  for (f in all_factors) {
    # Determine direction for this factor; default to "pos" if not specified
    direction <- if (f %in% names(factor_dir)) factor_dir[[f]] else "pos"
    message(f)
    # Define sign-specific sorted weight lists per factor
    universe_pos <- rownames(weights)[weights[, f] > 0]
    universe_neg <- rownames(weights)[weights[, f] < 0]
    w_pos <- universe_pos[order(weights[universe_pos, f], decreasing = TRUE)]
    w_neg <- universe_neg[order(weights[universe_neg, f], decreasing = FALSE)]
    df_list <- list()
    #pct <- perc_thresholds[[1]]
    for (pct in perc_thresholds) {
      N_pos <- max(1, floor(length(universe_pos) * pct / 100))
      N_neg <- max(1, floor(length(universe_neg) * pct / 100))
      feats_pos <- head(w_pos, N_pos)
      feats_neg <- head(w_neg, N_neg)
      combos <- list(
        posW_posDE = list(feats=feats_pos, de=de_pos), # risk factors with CMV+
        posW_negDE = list(feats=feats_pos, de=de_neg), # risk factors with CMV-
        negW_posDE = list(feats=feats_neg, de=de_pos), # prot factors with CMV+
        negW_negDE = list(feats=feats_neg, de=de_neg)  # prot factors with CMV-
      )
      #test <- "posW_posDE"  
      for (test in names(combos)) { # posW_posDE etc
        universe_test <- if (grepl("^posW", test)) universe_pos else universe_neg
        N <- if (grepl("^posW", test)) N_pos else N_neg

        L <- combos[[test]]
        k <- sum(L$feats %in% intersect(L$feats, L$de))
        perc <- k / N * 100 # % intersection of DEG list and top % feature list

        # choose universe based on weight direction
        pval <- enrich_test(L$feats, L$de, universe_test)
        df_list[[length(df_list)+1]] <- data.frame(
          Pct = pct,
          Test = test,
          Perc = perc,
          Pval = pval,
          Tested = N,
          DE_count = length(L$de)
        )
      }
    }
    df_plot <- do.call(rbind, df_list)
    df_plot$Pct <- factor(df_plot$Pct, levels=perc_thresholds)
    df_plot <- df_plot %>%
      dplyr::mutate(
        WeightRole = role_map[ifelse(direction == "pos", "pos", "neg")]
      )
    df_plot_list[[f]] <- df_plot

    if (FALSE) {
      # Only show points with P < 0.05
      #df_plot_sig <- df_plot %>% dplyr::filter(Pval < 0.05)
      df_plot_sig <- df_plot

      # Dot plot per factor (only significant points)
      pfile <- file.path(outdir, paste0("dotplot_", view_name, "_", f, ".pdf"))
      p <- ggplot(df_plot_sig, aes(x=Pct, y=Test, size=Perc, color=-log10(Pval))) +
        geom_point() +
        scale_color_viridis_c() +
        labs(title=paste(view_name, f, "enrichment"), x="Pct", y="Test") +
        theme_bw() +
        theme(axis.text.x=element_text(angle=45, hjust=1))
      ggsave(filename=pfile, plot=p, width=8, height=6)
      message(paste("Dot plot for", view_name, f, "saved to", pfile))

      # Prepare heatmap matrix: use unique Pval per Pct and Test
      df_pval <- df_plot %>%
        dplyr::select(Pct, Test, Pval) %>%
        dplyr::group_by(Pct, Test) %>%
        dplyr::slice(1) %>%
        dplyr::ungroup()
      mat <- df_pval %>%
        tidyr::pivot_wider(
          names_from = Test,
          values_from = Pval
        ) %>%
        tibble::column_to_rownames("Pct") %>%
        as.matrix()
      mat <- -log10(mat)
      mat[is.infinite(mat)] <- NA
      hfile <- file.path(outdir, paste0("heatmap_", view_name, "_", f, ".pdf"))
      pdf(hfile, width=6, height=6)
      ht <- Heatmap(mat, name="-log10(pval)", cluster_rows=FALSE, cluster_columns=FALSE,
              column_title=paste(view_name, f), row_title="Pct")
      draw(ht)
      dev.off()
      message(paste("Heatmap for", view_name, f, "saved to", hfile))
    }
  }
  # Write combined enrichment results for this view
  df_all <- do.call(rbind, lapply(names(df_plot_list), function(f) {
    df_plot_list[[f]] %>% dplyr::mutate(Factor = f)
  }))
  # After combining df_all and before writing CSV, adjust p-values
  df_all <- df_all %>%
    dplyr::group_by(Factor) %>%
    dplyr::mutate(padj = p.adjust(Pval, method = "BH")) %>%
    dplyr::ungroup()
  csv_path <- file.path(outdir, paste0("factor_enrichment_", view_name, ".csv"))
  write.csv(df_all, csv_path, row.names = FALSE)
  message(paste("Combined CSV for", view_name, "saved to", csv_path))
  if (return_df) return(df_plot_list)
}


for (view in available_views) {
  de_df <- switch(view,
                  cyt   = de_cyt_df,
                  gex   = de_gex_df,
                  metab = de_metab_df,
                  prot  = de_prot_df)
  # Use percentage thresholds
  perc_thresholds <- seq(5, 95, by=10)
  plot_factor_enrichment(view, de_df, factor_dir_map, outdir, perc_thresholds)
}


# Summarize significant enrichments (adjusted P < 0.05) across views for each factor
sig_enrich <- list()
for (view in available_views) {
  df_all <- read.csv(file.path(outdir, paste0("factor_enrichment_", view, ".csv"))) %>%
    dplyr::mutate(View = view) %>%
    dplyr::filter(padj < 0.05)
  sig_enrich[[view]] <- df_all
}
sig_df <- dplyr::bind_rows(sig_enrich)

# Initialize list to collect filtered results without CMV+/- tests
df_f_list <- list()

## Dotplots per factor, combining all views
# Loop over all factors present in sig_df
all_factors <- unique(sig_df$Factor)
for (f in all_factors) {
  df_f <- sig_df %>% filter(Factor == f)
  if (nrow(df_f) == 0) next

  # Annotate weight and DE directions
  df_f <- df_f %>%
    dplyr::mutate(
      View = dplyr::recode(View, !!!view_label_map),
      WeightDir = ifelse(grepl("^posW", Test), "pos", "neg"),
      DEdir = sub(".*_", "", Test),
      WeightRole = {
        dir_map <- if (f %in% names(factor_dir_map)) factor_dir_map[f] else "pos"
        ifelse(WeightDir == dir_map, "risk", "protective")
      },
      DEstatus = dplyr::case_when(
        DEdir == "posDE" ~ "CMV+",
        DEdir == "negDE" ~ "CMV-"
      ),
      TestLabel = paste(View, DEstatus, sep = "_")
    )

  # Store df_f for supplementary table
  df_f_list[[f]] <- df_f %>% dplyr::mutate(Factor = f)
  if (nrow(df_f) > 0) {
    # Set subtitle direction mapping, default to "pos" if not mapped
    subtitle_dir <- if (f %in% names(factor_dir_map)) factor_dir_map[f] else "pos"
    p2 <- ggplot(df_f, aes(x = Pct, y = TestLabel, size = Perc, color = -log10(padj))) +
      geom_point() +
      facet_grid(WeightRole ~ ., scales="free_y", space="free_y") +
      scale_color_viridis_c(name = "-log10(padj)") +
      scale_size_continuous(name = "% enriched") +
      labs(#title = paste("Enrichment for", f, " - associated with ", factor_labels[f]),
           subtitle = paste("Enrichment for ", f , " - associated with", factor_labels[f]),
           x = "Tested top % contributing features to factor", 
           y = ""
           ) +
      theme_bw() +
      # hide bg squares
      theme(
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())
    ggsave(
      file.path(outdir, paste0("summary_dotplot_", f, ".pdf")),
      plot = p2, width = 6, height = 4
    )
    message(paste("Summary dotplot for", f, "saved to", file.path(outdir, paste0("summary_dotplot_", f, ".pdf"))))
  }
}

# Write combined supplementary table
supp_df <- do.call(rbind, df_f_list)
write.csv(supp_df, file.path(outdir, "supplementary_enrichment.csv"), row.names = FALSE)
message("Supplementary table saved to ", file.path(outdir, "supplementary_enrichment.csv"))




all_enrich <- list()
for (view in available_views) {
  df_all <- read.csv(file.path(outdir, paste0("factor_enrichment_", view, ".csv"))) %>%
    dplyr::mutate(View = view)
  all_enrich[[view]] <- df_all
}
all_df <- dplyr::bind_rows(all_enrich)

all_factors <- unique(all_df$Factor)
for (f in all_factors) {
  df_f <- all_df %>% filter(Factor == f)
  if (nrow(df_f) == 0) next

  # Annotate weight and DE directions
  df_f <- df_f %>%
    dplyr::mutate(
      View = dplyr::recode(View, !!!view_label_map),
      WeightDir = ifelse(grepl("^posW", Test), "pos", "neg"),
      DEdir = sub(".*_", "", Test),
      WeightRole = {
        dir_map <- if (f %in% names(factor_dir_map)) factor_dir_map[f] else "pos"
        ifelse(WeightDir == dir_map, "risk", "protective")
      },
      DEstatus = dplyr::case_when(
        DEdir == "posDE" ~ "CMV+",
        DEdir == "negDE" ~ "CMV-"
      ),
      TestLabel = paste(View, DEstatus, sep = "_")
    )

  # Store df_f for supplementary table
  df_f_list[[f]] <- df_f %>% dplyr::mutate(Factor = f)
  if (nrow(df_f) > 0) {

    # ensure both facet rows exist
    df_f <- df_f %>%
      mutate(
        WeightRole = factor(WeightRole, levels = c("protective", "risk"))
      )

    this_sig_df  <- df_f %>% filter(padj < 0.05)
    nonsig_df <- df_f %>% filter(padj >= 0.05)

    subtitle_dir <- if (f %in% names(factor_dir_map)) factor_dir_map[f] else "pos"

    p2 <- ggplot(NULL, aes(x = Pct, y = TestLabel, size = Perc)) +
      # non‑significant dots in light grey
      geom_point(data = nonsig_df,
                 color = "grey80", alpha = 0.7) +
      # significant dots coloured by –log10(padj)
      geom_point(data = this_sig_df,
                 aes(color = -log10(padj))) +
      facet_grid(WeightRole ~ ., scales = "free_y", space = "free_y", drop = FALSE) +
      scale_color_viridis_c(name = "-log10(padj)") +
      scale_size_continuous(name = "% enriched") +
      labs(
        subtitle = paste("Enrichment for Factor", f, "- associated with", factor_labels[f]),
        x = "Tested top % contributing features to factor",
        y = ""
      ) +
      theme_bw() +
      theme(
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )

    ggsave(
      file.path(outdir, paste0("summary_dotplot_ALL_", f, ".pdf")),
      plot = p2, width = 6, height = 4
    )
    message(paste("Summary dotplot for", f, "saved to", file.path(outdir, paste0("summary_dotplot_ALL_", f, ".pdf"))))
  }
}



all_enrich <- list()
for (view in available_views) {
  df_all <- read.csv(file.path(outdir, paste0("factor_enrichment_", view, ".csv"))) %>%
    dplyr::mutate(View = view)
  all_enrich[[view]] <- df_all
}
all_df <- dplyr::bind_rows(all_enrich)

all_factors <- unique(all_df$Factor)
for (f in all_factors) {
  df_f <- all_df %>% filter(Factor == f)
  if (nrow(df_f) == 0) next

  # Annotate weight and DE directions
  df_f <- df_f %>%
    dplyr::mutate(
      View = dplyr::recode(View, !!!view_label_map),
      WeightDir = ifelse(grepl("^posW", Test), "pos", "neg"),
      DEdir = sub(".*_", "", Test),
      WeightRole = {
        dir_map <- if (f %in% names(factor_dir_map)) factor_dir_map[f] else "pos"
        ifelse(WeightDir == dir_map, "risk", "protective")
      },
      DEstatus = dplyr::case_when(
        DEdir == "posDE" ~ "CMV+",
        DEdir == "negDE" ~ "CMV-"
      ),
      TestLabel = paste(View, DEstatus, sep = "_")
    )

  # Store df_f for supplementary table
  df_f_list[[f]] <- df_f %>% dplyr::mutate(Factor = f)
  if (nrow(df_f) > 0) {

    # ensure both facet rows exist
    df_f <- df_f %>%
      mutate(
        WeightRole = factor(WeightRole, levels = c("protective", "risk"))
      )

    non_signif <- df_f %>% filter(Pval >= 0.05)
    nom_signif  <- df_f %>% filter(Pval < 0.05 & padj >= 0.05)
    fdr_signif <- df_f %>% filter(padj < 0.05)

    subtitle_dir <- if (f %in% names(factor_dir_map)) factor_dir_map[f] else "pos"

    # Make DEstatus an ordered factor so "CMV+" appears left of "CMV‑"
    df_f <- df_f %>%
      mutate(DEstatus = factor(DEstatus, levels = c("CMV+", "CMV-")))
    p2 <- ggplot(NULL, aes(x = Pct, y = View, size = Perc, color = -log10(padj))) +
      geom_point(data = non_signif, color = "grey80", alpha = 0.7, shape = 16) + # non‑significant circle in light grey
      #geom_point(data = nom_signif, shape = 1) + # nominal significant circles
      geom_point(data = fdr_signif, shape = 16) + # fdr significant triangles
      facet_grid(WeightRole ~ DEstatus, scales = "free_y", space = "free", drop = FALSE) +
      scale_color_viridis_c(name = "-log10(padj)", option = "viridis") +
      scale_size_continuous(name = "% enriched") +
      labs(
        subtitle = paste("Enrichment for", f, "- associated with", factor_labels[f]),
        x = "Tested top percentage contributing features to factor",
        y = ""
      ) +
      theme_bw() +
      theme(
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )

    ggsave(
      file.path(outdir, paste0("summary_dotplot_ALL_", f, "_SPLIT.pdf")),
      plot = p2, width = 8, height = 4
    )
    message(paste("Summary dotplot for", f, "saved to", file.path(outdir, paste0("summary_dotplot_ALL_", f, "_SPLIT.pdf"))))
  }
}





# ---- Aggregate strongest enrichment per Factor/Test ----

# Read combined enrichment results for each view
all_results <- lapply(available_views, function(view) {
  read.csv(file.path(outdir, paste0("factor_enrichment_", view, ".csv"))) %>%
    mutate(View = view)
}) %>% bind_rows()

# Combine Test and View into TestLabel if not already
all_results <- all_results %>%
  mutate(TestLabel = paste(View, Test, sep = "_"))

# For each Factor and TestLabel, take the minimum adjusted p-value
best_results <- all_results %>%
  group_by(Factor, TestLabel) %>%
  summarize(min_padj = min(padj, na.rm = TRUE), .groups = "drop")

# ---- Faceted heatmap: easier interpretation ----

# Split TestLabel into View and Test columns
library(tidyr)
library(dplyr)
best_results <- best_results %>%
  separate(TestLabel, into = c("View", "Test"), sep = "_", remove = FALSE)

# Set View and Test orderings
best_results$View <- factor(best_results$View,
                            levels = view_label_map[available_views])
# Specify Test order if desired
test_levels <- c("negW_negDE", "negW_posDE", "posW_negDE", "posW_posDE")
best_results$Test <- factor(best_results$Test, levels = test_levels)

# Rebuild heat_df from updated best_results (not used for plotting but for consistency)
heat_df <- best_results %>%
  mutate(LogP = -log10(min_padj)) %>%
  arrange(View, Test) %>%
  pivot_wider(id_cols = Factor, names_from = TestLabel, values_from = LogP) %>%
  column_to_rownames("Factor")

# Faceted tile plot
p_heat <- ggplot(best_results, aes(x = Test, y = Factor, fill = -log10(min_padj))) +
  geom_tile(color = "white") +
  facet_wrap(~ View, ncol = 1, scales = "free_x", strip.position = "top") +
  scale_fill_gradientn(
    name = "-log10(min padj)",
    colours = c("#384E78", "#5874DC", "#6AAB9C", "#FA9284", "#E06C78"),
    na.value = "grey90"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(),
    strip.background = element_rect(fill = "grey95", color = NA)
  )

# Save new faceted plot
ggsave(file.path(outdir, "summary_strongest_enrichment_heatmap_faceted.pdf"),
       p_heat, width = 8, height = 6)

message(file.path(outdir, "summary_strongest_enrichment_heatmap_faceted.pdf"))




if (FALSE) {


  # ---- Alternative visualizations ----

  # 1. Small-multiples line charts: percent overlap trends by View & DEstatus
  df_all <- supp_df
  df_all$DEstatus <- ifelse(grepl("posDE", df_all$Test), "CMV+", "CMV-")
  for (f in unique(df_all$Factor)) {
    p_line <- ggplot(df_all %>% filter(Factor==f),
                    aes(x=as.numeric(as.character(Pct)), y=Perc, color=View, linetype=DEstatus)) +
      geom_line() + geom_point() +
      facet_wrap(~ View+WeightRole, nrow=1) +
      labs(title=paste("Overlap trends for", f),
          x="Top percentile of features", y="% overlap") +
      theme_bw()
    ggsave(file.path(outdir, paste0("alt1_line_", f, ".pdf")),
          plot=p_line, width=8, height=3)
  }

  # 2. Heatmap of –log10(pval) with percent-overlap overlay
  library(ComplexHeatmap)
  library(circlize)
  for (f in unique(df_all$Factor)) {
    df_mat <- df_all %>% filter(Factor==f) %>%
      mutate(LogP=-log10(Pval)) %>%
      dplyr::select(Pct, TestLabel, LogP, Perc) %>%
      tidyr::pivot_wider(names_from=Pct, values_from=LogP) 
    mat <- as.matrix(df_mat[,-1])
    rownames(mat) <- df_mat$TestLabel
    ht <- Heatmap(mat, name="-log10(pval)", cell_fun = function(j, i, x, y, w, h, col) {
      grid.text(sprintf("%.1f%%", df_all$Perc[df_all$Factor==f &
          df_all$TestLabel==rownames(mat)[i] &
          df_all$Pct==colnames(mat)[j]]),
        x, y, gp=gpar(fontsize=8))
    },
    cluster_rows=FALSE, cluster_columns=FALSE,
    column_title=paste("Pct for", f))
    pdf(file.path(outdir, paste0("alt2_heatmap_", f, ".pdf")), width=6, height=6)
    draw(ht)
    dev.off()
  }

  # 3. Paired bar-and-dot infographic at best percentile
  for (f in unique(df_all$Factor)) {
    df_best <- df_all %>% filter(Factor==f) %>%
      group_by(TestLabel) %>%
      slice_min(order_by=Pval, with_ties=FALSE) %>%
      ungroup()
    p_bar <- ggplot(df_best, aes(x=TestLabel, y=Perc, fill=DEstatus)) +
      geom_bar(stat="identity") +
      geom_point(aes(size=-log10(Pval)), color="black") +
      coord_flip() +
      labs(title=paste("Best enrichment for", f),
          x="", y="% overlap") +
      theme_bw()
    ggsave(file.path(outdir, paste0("alt3_bar_", f, ".pdf")),
          plot=p_bar, width=6, height=4)
  }

}
