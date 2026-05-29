# =============================================================================
# Author      : Dr. rer. nat. Guillermo G. Torres <guigotoe@gmail.com>
# Affiliation : Institute of Clinical Molecular Biology (IKMB), Kiel University
# Project     : MILENA — Microbiome of Children (Pastoral Kids Cohort)
# Repository  : https://github.com/cobinet/Microbiome_Children_Milena_Project
# Script      : AssocAnalysis.R
# Version     : 1.0.0
# Date        : 2025
#
# Description :
#   Differential abundance and multi-contrast association analysis linking
#   functional microbiome profiles (KO, EC, MetaCyc pathways) to parasitic
#   infection status (Blastocystis hominis [BH], Giardia lamblia [GL],
#   Dientamoeba fragilis [DF]) in a pediatric cohort (6–59 months).
#
#   The analysis framework:
#     1. Loads CLR-transformed functional feature matrices + cleaned metadata
#     2. Residualizes each feature block against age-group and sex covariates
#        using limma::lmFit to remove demographic confounding
#     3. Constructs a 3-level infection outcome Y per pathogen target:
#        {None = no infections; <target> = single target infection only;
#         Coinfection = target + ≥1 other concurrent pathogen}
#     4. Fits multi-contrast linear models (limma eBayes) per feature block
#        with BH-adjusted p-values across three contrasts per pathogen
#     5. Exports publication-ready HTML tables (gt) and lollipop / volcano
#        plots (ggplot2) for each pathogen × contrast combination
#
# Statistical method :
#   Empirical Bayes moderated t-tests (limma) on CLR-scale data.
#   FDR correction via Benjamini-Hochberg (BH). Threshold: q ≤ 0.10.
#
# Inputs (all under DATA_DIR — NOT included in this repository; see data/README.md):
#   preprocessing_output/asv_table_clr_transformed.tsv
#     rows = ASV_IDs, cols = PC_code sample identifiers, values = CLR-transformed counts
#   preprocessing_output/KO_Functions_clr_transformed.tsv
#     rows = KEGG Orthology IDs (K#####), cols = PC_codes, values = CLR
#   preprocessing_output/EC_Functions_clr_transformed.tsv
#     rows = EC numbers (ec:#.#.#.#), cols = PC_codes, values = CLR
#   preprocessing_output/MetaCyc_Pathways_clr_transformed.tsv
#     rows = MetaCyc pathway IDs, cols = PC_codes, values = CLR
#   analysis_ready_data/asv_taxonomy_clean.csv
#     cols = ASV_ID, Kingdom, Phylum, Class, Order, Family, Genus, Species
#   analysis_ready_data/metadata_clean.csv
#     cols = PC_code, AgeM, AgeGroup, Sexchildren, KGcode,
#            Dietary_pattern, Nutritionalstatus, HEI, WAZ, HAZ, BAZ,
#            infectioncoinfection, BH, GL, DF
#
# Outputs (written to OUT_DIR):
#   {TARGET}_ALL_features_all_contrasts.csv   — full limma result table per pathogen
#   {TARGET}_publication_table.html           — top hits formatted for publication
#   {TARGET}_lollipop_{contrast}.png          — lollipop plots of top 20 features
#   volcano_{BLOCK}_{contrast}.png            — volcano plots per block × contrast
#   infect_coinfection_top_features_best_by_block.csv — combined best-hit summary
#   infect_coinfection_{BLOCK}_all_contrasts.csv      — per-block full tables
#
# Usage :
#   # Adjust DATA_DIR and OUT_DIR at the top of this script, then:
#   Rscript scripts/AssocAnalysis.R
#
# Dependencies :
#   R ≥ 4.3; see R_environment.txt for full sessionInfo()
#   CRAN: gt, data.table, dplyr, tidyr, stringr, readr, purrr, forcats,
#         ggplot2, gridExtra, RColorBrewer, scales, tibble
#   Bioconductor: limma, vegan, phyloseq, zCompositions, ANCOMBC, Maaslin2
#   CRAN: fgsea, plotly, mixOmics
# =============================================================================


# ── 0. Bootstrap Package Installer ───────────────────────────────────────────
# Utility: installs missing packages from CRAN/Bioc before loading them.
packages <- function(requirements, quiet = FALSE) {
  has <- requirements %in% rownames(installed.packages())
  if (any(!has)) {
    message("Installing missing packages: ", paste(requirements[!has], collapse = ", "))
    setRepositories(ind = c(1:7))
    r <- getOption("repos")
    r["CRAN"] <- "https://cran.uni-muenster.de/"
    install.packages(requirements[!has], repos = r)
  }
  for (r in requirements) {
    if (quiet) suppressMessages(require(r, character.only = TRUE))
    else message(r, ": ", suppressMessages(require(r, character.only = TRUE)))
  }
}

packages(
  c("gt", "data.table", "dplyr", "tidyr", "stringr", "readr", "purrr",
    "forcats", "ggplot2", "gridExtra", "RColorBrewer", "scales", "tibble",
    "limma", "vegan", "phyloseq", "zCompositions", "ANCOMBC", "Maaslin2",
    "fgsea", "plotly", "mixOmics"),
  quiet = TRUE
)


# ── 1. User-configurable Paths ────────────────────────────────────────────────
# Change DATA_DIR to the folder containing your preprocessed data outputs
# (produced by funAnalysis_metagPastoKids.ipynb) and OUT_DIR to where
# you want results written.

DATA_DIR   <- file.path("data", "microbiome_analysis_output")   # <-- adjust if needed
QC_DIR     <- file.path(DATA_DIR, "analysis_ready_data")
PRE_DIR    <- file.path(DATA_DIR, "preprocessing_output")
OUT_DIR    <- file.path(DATA_DIR, "infectionsFunctions_result")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)


# ── 2. Load CLR-Transformed Feature Matrices + Metadata ──────────────────────
# All matrices are already CLR-transformed (output of the preprocessing pipeline
# in funAnalysis_metagPastoKids.ipynb). Rows = features, cols = samples (PC_codes).

message("Loading data...")
asv  <- read_tsv(file.path(PRE_DIR, "asv_table_clr_transformed.tsv"),
                 show_col_types = FALSE)
tax  <- read_csv(file.path(QC_DIR,  "asv_taxonomy_clean.csv"),
                 show_col_types = FALSE)
ko   <- read_tsv(file.path(PRE_DIR, "KO_Functions_clr_transformed.tsv"),
                 show_col_types = FALSE)
ec   <- read_tsv(file.path(PRE_DIR, "EC_Functions_clr_transformed.tsv"),
                 show_col_types = FALSE)
mc   <- read_tsv(file.path(PRE_DIR, "MetaCyc_Pathways_clr_transformed.tsv"),
                 show_col_types = FALSE)
meta <- read_csv(file.path(QC_DIR,  "metadata_clean.csv"),
                 show_col_types = FALSE)


# ── 3. Metadata Factor Preparation ───────────────────────────────────────────
# Cast categorical variables to factors with explicit level ordering.
# AgeGroup follows developmental order (6–11m < 12–23m < 24–59m).

meta <- meta %>%
  mutate(
    Dietary_pattern   = factor(Dietary_pattern),
    Nutritionalstatus = factor(Nutritionalstatus),
    AgeGroup = factor(AgeGroup, levels = c("6-11m", "12-23m", "24-59m"))
  )

meta_df <- meta %>%
  mutate(
    Sexchildren = factor(Sexchildren),
    KGcode      = factor(KGcode),
    AgeM        = as.numeric(AgeM),
    AgeGroup    = factor(AgeGroup, levels = c("6-11m", "12-23m", "24-59m"))
  ) %>%
  column_to_rownames("PC_code")


# ── 4. Identify Common Samples Across All Blocks ─────────────────────────────
# Restrict analysis to samples present in ALL feature blocks AND the metadata.
# This ensures row alignment between feature matrices and the outcome vector Y.

sample_ids <- Reduce(
  intersect,
  list(colnames(ko)[-c(1, 2)],
       colnames(mc)[-c(1, 2)],
       colnames(asv)[-c(1, 2)],
       meta$PC_code)
)
message("Samples in all blocks: ", length(sample_ids))


# ── 5. Residualize Feature Blocks Against Demographic Covariates ──────────────
# Remove variation attributable to AgeGroup and Sexchildren (nuisance confounders)
# via limma::lmFit residuals. The residualized CLR values represent infection-
# associated functional variation independent of demographic structure.
# Z-scoring per feature is optional but aids interpretability across blocks.

design_resid <- model.matrix(~ AgeGroup + Sexchildren, data = meta_df)

resid_block <- function(M, design) {
  fit <- lmFit(M, design)
  res <- residuals(fit)
  t(scale(t(res)))  # z-score per feature row
}

# Subset each block to common samples and transpose to features × samples
KO_res  <- as.matrix(ko  %>% select(`function`, all_of(sample_ids)) %>%
                       column_to_rownames("function")  %>% select(all_of(sample_ids)))
MC_res  <- as.matrix(mc  %>% select(`pathway`,  all_of(sample_ids)) %>%
                       column_to_rownames("pathway")   %>% select(all_of(sample_ids)))
EC_res  <- as.matrix(ec  %>% select(`function`, all_of(sample_ids)) %>%
                       column_to_rownames("function")  %>% select(all_of(sample_ids)))
ASV_res <- as.matrix(asv %>% select(`ASV_ID`,   all_of(sample_ids)) %>%
                       column_to_rownames("ASV_ID")    %>% select(all_of(sample_ids)))

# Column-transposed versions (samples × features) kept for mixOmics
KO_block  <- t(KO_res)
MC_block  <- t(MC_res)
EC_block  <- t(EC_res)
ASV_block <- t(ASV_res)


# ── 6. Build Annotation Look-up Tables ───────────────────────────────────────
# Robustly detect ID and description column names for each functional database.

ko_id_col <- names(ko)[1]   # usually "function"
ec_id_col <- names(ec)[1]
mc_id_col <- names(mc)[1]   # usually "pathway"

ko_annot <- ko %>% select(id = all_of(ko_id_col), description = any_of("description")) %>% distinct()
ec_annot <- ec %>% select(id = all_of(ec_id_col), description = any_of("description")) %>% distinct()
mc_annot <- mc %>% select(id = all_of(mc_id_col), description = any_of("description")) %>% distinct()


# ── 7. Infection Outcome Construction ────────────────────────────────────────
# For each target pathogen (BH / GL / DF), classify each child as:
#   "None"         — no infection by any of the three parasites
#   "<target>"     — infected only by the target pathogen (single infection)
#   "Coinfection"  — infected by target AND ≥1 additional pathogen
# Samples positive for non-target parasites only are excluded (NA → dropped).

build_Y_infection <- function(meta_df,
                              target   = "BH",
                              inf_cols = c("BH", "GL", "DF")) {
  df      <- meta_df %>% mutate(across(all_of(inf_cols), ~ replace_na(as.numeric(.), 0)))
  others  <- setdiff(inf_cols, target)
  any_all   <- rowSums(df[, inf_cols,  drop = FALSE]) > 0
  target_on <- df[[target]] == 1
  other_sum <- rowSums(df[, others,   drop = FALSE])

  Y <- ifelse(
    !any_all,                        "None",
    ifelse(target_on & other_sum == 0, target,
    ifelse(target_on & other_sum >= 1, "Coinfection", NA))
  )

  keep <- !is.na(Y)
  list(
    Y    = factor(Y[keep], levels = c("None", target, "Coinfection")),
    keep = keep
  )
}


# ── 8. Multi-Contrast limma Differential Abundance Wrapper ───────────────────
# Fits a one-way ANOVA-style design (~ 0 + Y) and extracts all pairwise
# contrasts between infection groups. Returns full topTable results enriched
# with group means, SEs, 95% CIs, and optional functional annotations.

run_multiclass_limma2 <- function(M_res,
                                  Y,
                                  annot_df   = NULL,
                                  block_name = "KO",
                                  desc_col   = "description",
                                  target) {
  stopifnot(is.matrix(M_res), ncol(M_res) == length(Y))
  Y <- droplevels(Y)

  tabY <- table(Y)
  if (any(tabY < 3))
    message(block_name, ": small class sizes → results may be unstable: ",
            paste(names(tabY), tabY, collapse = ", "))

  # One-hot design matrix (no intercept) for clean contrast specification
  design <- model.matrix(~ 0 + Y)
  fit0   <- lmFit(M_res, design)

  # Build only contrasts whose group levels actually exist in this subset
  coefs        <- colnames(design)
  target_coef  <- paste0("Y", target)
  contr_list   <- list()

  if (all(c(target_coef, "YNone") %in% coefs))
    contr_list[[paste0(target, "_vs_None")]] <- makeContrasts(
      contrasts = sprintf("%s - YNone", target_coef), levels = design)

  if (all(c("YCoinfection", "YNone") %in% coefs))
    contr_list$Coinf_vs_None <- makeContrasts(
      contrasts = "YCoinfection - YNone", levels = design)

  if (all(c("YCoinfection", target_coef) %in% coefs))
    contr_list[[paste0("Coinf_vs_", target)]] <- makeContrasts(
      contrasts = sprintf("YCoinfection - %s", target_coef), levels = design)

  if (length(contr_list) == 0)
    stop("No valid contrasts — check infection class counts for target: ", target)

  contr <- do.call(cbind, contr_list)
  fit   <- eBayes(contrasts.fit(fit0, contr))

  # Gather per-contrast topTable results
  all_tbl <- map_dfr(seq_len(ncol(contr)), function(i) {
    topTable(fit, coef = i, number = Inf, sort.by = "P") %>%
      rownames_to_column("feature") %>%
      as_tibble() %>%
      rename(p = P.Value, q = adj.P.Val) %>%
      mutate(contrast = colnames(contr)[i])
  })

  # Group means (CLR) — dynamically named to the target pathogen
  grp_means <- sapply(levels(Y), function(g) {
    rowMeans(M_res[, Y == g, drop = FALSE], na.rm = TRUE)
  }) %>%
    as.data.frame() %>%
    rownames_to_column("feature") %>%
    {
      nm <- names(.)
      nm[nm == "None"]        <- "mean_None"
      nm[nm == target]        <- paste0("mean_", target)
      nm[nm == "Coinfection"] <- "mean_Coinf"
      names(.) <- nm; .
    }

  out <- all_tbl %>% left_join(grp_means, by = "feature")

  # Per-group SE for group means
  sd_none  <- if ("None"        %in% levels(Y)) apply(M_res[, Y == "None",        drop = FALSE], 1, sd, na.rm = TRUE) else NULL
  sd_tgt   <- if (target        %in% levels(Y)) apply(M_res[, Y == target,        drop = FALSE], 1, sd, na.rm = TRUE) else NULL
  sd_coinf <- if ("Coinfection" %in% levels(Y)) apply(M_res[, Y == "Coinfection", drop = FALSE], 1, sd, na.rm = TRUE) else NULL

  n_none  <- as.numeric(tabY["None"])
  n_tgt   <- as.numeric(tabY[target])
  n_coinf <- as.numeric(tabY["Coinfection"])

  out <- out %>%
    mutate(
      logFC_se      = logFC / t,
      mean_None_se  = if (!is.null(sd_none))  sd_none[feature]  / sqrt(n_none)  else NA_real_,
      tmp_target_se = if (!is.null(sd_tgt))   sd_tgt[feature]   / sqrt(n_tgt)   else NA_real_,
      mean_Coinf_se = if (!is.null(sd_coinf)) sd_coinf[feature] / sqrt(n_coinf) else NA_real_
    ) %>%
    rename(!!paste0("mean_", target, "_se") := tmp_target_se)

  # 95% CIs for logFC and group means
  crit               <- qt(0.975, df = df.residual(fit))
  target_mean_col    <- paste0("mean_", target)
  target_mean_se_col <- paste0("mean_", target, "_se")

  out <- out %>%
    mutate(
      logFC_lower = logFC - crit * (logFC / t),
      logFC_upper = logFC + crit * (logFC / t),
      mean_None_lower  = mean_None - crit * mean_None_se,
      mean_None_upper  = mean_None + crit * mean_None_se,
      !!paste0("mean_", target, "_lower") :=
        if (target %in% levels(Y)) (!!sym(target_mean_col)) - crit * (!!sym(target_mean_se_col)) else NA_real_,
      !!paste0("mean_", target, "_upper") :=
        if (target %in% levels(Y)) (!!sym(target_mean_col)) + crit * (!!sym(target_mean_se_col)) else NA_real_,
      mean_Coinf_lower = if ("Coinfection" %in% levels(Y)) mean_Coinf - crit * mean_Coinf_se else NA_real_,
      mean_Coinf_upper = if ("Coinfection" %in% levels(Y)) mean_Coinf + crit * mean_Coinf_se else NA_real_
    )

  if (!is.null(annot_df))
    out <- out %>% left_join(annot_df, by = c("feature" = "id"))

  out %>% mutate(block = block_name)
}


# ── 9. Publication Table Generator ───────────────────────────────────────────
# Formats significant hits (q ≤ 0.10, logFC > 0.15) as a gt HTML table.

fmt_q <- function(x) {
  case_when(
    x < 1e-4 ~ formatC(x, format = "e", digits = 2),
    x < 0.01 ~ sprintf("%.3f", x),
    TRUE     ~ sprintf("%.2f", x)
  )
}

make_pub_table <- function(df, contrast_labels,
                           top_n_per_block = 25,
                           file_html = "pub.html",
                           title = "Top features") {
  tab <- df %>%
    filter(contrast %in% contrast_labels) %>%
    group_by(block, contrast) %>%
    arrange(q, desc(abs(logFC)), .by_group = TRUE) %>%
    slice_head(n = top_n_per_block) %>%
    ungroup() %>%
    mutate(
      Direction = if_else(logFC > 0, "↑ left group higher", "↓ left group lower"),
      `q (FDR)` = fmt_q(q)
    ) %>%
    transmute(
      Block = block, Contrast = contrast, Feature = feature,
      Description = coalesce(description, NA_character_),
      `CLR logFC` = round(logFC, 2), `q (FDR)`,
      `Mean None`   = round(mean_None,   2),
      `Mean Single` = round(mean_Single, 2),
      `Mean Coinf`  = round(mean_Coinf,  2),
      Direction
    )

  gt_tab <- tab %>%
    gt(groupname_col = "Block") %>%
    tab_header(
      title    = md(paste0("**", title, "**")),
      subtitle = md("Effect sizes on CLR scale; FDR by BH")
    ) %>%
    cols_label(
      Feature = "ID", Description = "Description",
      `CLR logFC` = "logFC", `q (FDR)` = "q",
      `Mean None` = "Mean (None)", `Mean Single` = "Mean (Single)",
      `Mean Coinf` = "Mean (Coinf)"
    ) %>%
    tab_spanner(label = "Group means (CLR)",
                columns = c(`Mean None`, `Mean Single`, `Mean Coinf`)) %>%
    opt_row_striping() %>%
    opt_table_outline()

  gtsave(gt_tab, file_html)
  invisible(tab)
}


# ── 10. Lollipop Plot Generator ───────────────────────────────────────────────
# Horizontal lollipop chart of top N features per block, coloured by direction
# and shaped by FDR significance. Saved as high-resolution PNG.

make_lollipop <- function(df, contrast_label,
                          file_png = "lollipop.png", n = 20) {
  dat <- df %>%
    filter(contrast == contrast_label) %>%
    arrange(q, desc(abs(logFC))) %>%
    group_by(block) %>%
    slice_head(n = n) %>%
    ungroup() %>%
    mutate(
      label = if_else(!is.na(description),
                      paste0(feature, " — ", description), feature),
      label = stringr::str_trunc(label, 80),
      sig   = q <= 0.10,
      dir   = if_else(logFC > 0, "↑ left group", "↓ left group"),
      label = forcats::fct_reorder(label, logFC)
    )

  p <- ggplot(dat, aes(x = label, y = logFC, color = dir)) +
    geom_segment(aes(xend = label, y = 0, yend = logFC), linewidth = 0.4) +
    geom_point(aes(shape = sig), size = 2) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), name = "FDR ≤ 0.10") +
    facet_wrap(~ block, scales = "free_y", ncol = 1) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = contrast_label, x = NULL, y = "CLR logFC (left vs right)") +
    theme_bw(base_size = 11) +
    theme(strip.background = element_rect(fill = "grey95"),
          strip.text       = element_text(face = "bold"),
          legend.position  = "top") +
    coord_flip()

  ggsave(file_png, p, width = 7.5, height = 10, dpi = 300)
}


# ── 11. Per-Pathogen Driver Function ─────────────────────────────────────────
# Orchestrates the full analysis for a single target pathogen:
#   build outcome Y → subset matrices → run limma → combine → export results.

analyze_one_infection <- function(target     = "BH",
                                  KO_res, EC_res, MC_res,
                                  meta_df,
                                  out_path   = OUT_DIR,
                                  out_prefix = "BH") {
  by   <- build_Y_infection(meta_df, target = target)
  Y    <- by$Y
  keep <- by$keep

  if (length(Y) < 10) {
    message("[", target, "] Too few usable samples after filtering; skipping.")
    return(invisible(NULL))
  }
  message("[", target, "] Classes: ", paste(names(table(Y)), table(Y), collapse = " | "))

  # Restrict feature matrices to samples retained after NA exclusion
  KO_sub <- KO_res[, keep, drop = FALSE]
  EC_sub <- EC_res[, keep, drop = FALSE]
  MC_sub <- MC_res[, keep, drop = FALSE]

  # Differential abundance per functional block
  res_KO <- run_multiclass_limma2(KO_sub, Y, annot_df = ko_annot, block_name = "KO",      target = target)
  res_EC <- run_multiclass_limma2(EC_sub, Y, annot_df = ec_annot, block_name = "EC",      target = target)
  res_MC <- run_multiclass_limma2(MC_sub, Y, annot_df = mc_annot, block_name = "MetaCyc", target = target)

  all_res <- bind_rows(res_KO, res_EC, res_MC) %>%
    mutate(contrast = str_replace_all(contrast, "Y", ""))

  # Derive "<target>_vs_Coinf" by sign-flipping "Coinf_vs_<target>"
  coinf_contr <- paste0("Coinf_vs_", target)
  if (coinf_contr %in% unique(all_res$contrast)) {
    target_vs_coinf <- all_res %>%
      filter(contrast == coinf_contr) %>%
      mutate(contrast = paste0(target, "_vs_Coinf"),
             logFC = -logFC, t = -t)
    all_res <- bind_rows(all_res, target_vs_coinf)
  }

  readr::write_csv(all_res, file.path(out_path, paste0(target, "_ALL_features_all_contrasts.csv")))

  # Publication table: significant hits only
  contrats_pub <- c(paste0(target, " - None"), paste0("Coinfection - ", target))
  make_pub_table(
    df = all_res %>%
      filter(q <= 0.10, logFC > 0.15) %>%
      mutate(mean_Single = .data[[paste0("mean_", target)]]),
    contrast_labels = intersect(contrats_pub, unique(all_res$contrast)),
    file_html = file.path(out_path, paste0(target, "_publication_table.html")),
    title = paste0("Top functional features for ", target,
                   " (", target, " vs None; ", target, " vs Coinfection)")
  )

  # Lollipop plots for the two primary contrasts
  if (paste0(target, " - None") %in% unique(all_res$contrast))
    make_lollipop(all_res, paste0(target, " - None"),
                  file_png = file.path(out_path, paste0(target, "_lollipop_", target, "_vs_None.png")),
                  n = 20)

  if (paste0("Coinfection - ", target) %in% unique(all_res$contrast))
    make_lollipop(all_res, paste0("Coinfection - ", target),
                  file_png = file.path(out_path, paste0(target, "_lollipop_", target, "_vs_Coinf.png")),
                  n = 20)

  message("[", target, "] Complete — results written to: ", out_path)
}


# ── 12. Execute Analysis for All Three Pathogens ─────────────────────────────
# BH = Blastocystis hominis | GL = Giardia lamblia | DF = Dientamoeba fragilis

message("\n=== Starting per-pathogen differential abundance analysis ===\n")

BH_out <- analyze_one_infection(target = "BH", KO_res, EC_res, MC_res, meta_df)
GL_out <- analyze_one_infection(target = "GL", KO_res, EC_res, MC_res, meta_df)
DF_out <- analyze_one_infection(target = "DF", KO_res, EC_res, MC_res, meta_df)

message("\n=== All analyses complete. Results in: ", OUT_DIR, " ===\n")
