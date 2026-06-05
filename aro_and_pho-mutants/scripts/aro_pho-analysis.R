library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(pracma)
library(multcomp)

##user configurations:
EXPERIMENT  <- "aro"   # "aro" or "pho"
PLOT_DIR    <- "./aro_and_pho-mutants/plots/"
DATA_DIR    <- "./aro_and_pho-mutants/data/"

## ============================================================
## EXPERIMENT DEFINITIONS
## ============================================================
## Each named list entry defines one plot panel.
## Fields:
##   conditions   - strain_condition values to include
##   cols         - named colour vector
##   factor_order - factor level order for x-axis / legend
##   baseline_ref - the strain_condition used as AUC denominator
##   curve_out    - filename for the growth curve tiff
##   auc_out      - filename for the AUC bar chart tiff
##   auc_w/auc_h  - width/height of the AUC plot (inches)
##   dunnett_ref  - reference level for Dunnett's test
##                  (if NULL the first factor level is used)

ARO_PANELS <- list(
  main = list(
    conditions   = c("WT -Phe", "WT +8mM", "aro10 +8mM", "aro4 +8mM"),
    cols         = c("WT -Phe" = "gray50", "WT +8mM" = "#0072BD",
                     "aro4 +8mM" = "#D95319", "aro10 +8mM" = "#77AC30"),
    factor_order = c("WT -Phe", "WT +8mM", "aro10 +8mM", "aro4 +8mM"),
    baseline_ref = "WT -Phe",
    curve_out    = "aro4_curves.tiff",
    auc_out      = "aro4_auc.tiff",
    auc_w = 5, auc_h = 3.2,
    dunnett_ref  = "WT +8mM"   # second Dunnett run; first uses default (WT -Phe)
  ),

  neg_phe_controls = list(
    conditions   = c("WT -Phe", "aro4 -Phe", "aro10 -Phe"),
    cols         = c("WT -Phe" = "gray50", "aro4 -Phe" = "#D95319",
                     "aro10 -Phe" = "#77AC30"),
    factor_order = c("WT -Phe", "aro4 -Phe", "aro10 -Phe"),
    baseline_ref = "WT -Phe",
    curve_out    = "aro4_curves (neg_phe controls).tiff",
    auc_out      = "aro4_auc (neg_phe controls).tiff",
    auc_w = 4.3, auc_h = 2.8,
    dunnett_ref  = NULL
  ),

  aro4_mutant_control = list(
    conditions   = c("aro4 -Phe", "aro4 +8mM"),
    cols         = c("aro4 -Phe" = "gray50", "aro4 +8mM" = "#D95319"),
    factor_order = c("aro4 -Phe", "aro4 +8mM"),
    baseline_ref = "aro4 -Phe",
    curve_out    = "aro4_curves (aro4 controls).tiff",
    auc_out      = "aro4_auc (aro4 controls).tiff",
    auc_w = 4, auc_h = 2.5,
    dunnett_ref  = NULL
  ),

  aro10_mutant_control = list(
    conditions   = c("aro10 -Phe", "aro10 +8mM"),
    cols         = c("aro10 -Phe" = "gray50", "aro10 +8mM" = "#77AC30"),
    factor_order = c("aro10 -Phe", "aro10 +8mM"),
    baseline_ref = "aro10 -Phe",
    curve_out    = "aro4_curves (aro10 controls).tiff",
    auc_out      = "aro4_auc (aro10 controls).tiff",
    auc_w = 4, auc_h = 2.5,
    dunnett_ref  = NULL
  ),

  moderate_conc_aro10 = list(
    conditions   = c("WT +0.9mM", "aro10 +0.9mM", "aro10 +8mM"),
    cols         = c("WT +0.9mM" = "gray50", "aro10 +0.9mM" = "#0072BD",
                     "aro10 +8mM" = "#77AC30"),
    factor_order = c("WT +0.9mM", "aro10 +0.9mM", "aro10 +8mM"),
    baseline_ref = "WT +0.9mM",
    curve_out    = "aro10_curves_moderateconc.tiff",
    auc_out      = NULL,   # shared AUC plot below
    auc_w = NULL, auc_h = NULL,
    dunnett_ref  = NULL
  ),

  moderate_conc_aro4 = list(
    conditions   = c("WT +0.9mM", "aro4 +0.9mM", "aro4 +8mM"),
    cols         = c("WT +0.9mM" = "gray50", "aro4 +0.9mM" = "#D95319",
                     "aro4 +8mM" = "#EDB120"),
    factor_order = c("WT +0.9mM", "aro4 +0.9mM", "aro4 +8mM"),
    baseline_ref = "WT +0.9mM",
    curve_out    = "aro4_curves_moderateconc.tiff",
    auc_out      = NULL,
    auc_w = NULL, auc_h = NULL,
    dunnett_ref  = NULL
  )
)

PHO_PANELS <- list(
  main = list(
    conditions   = c("WT -Phe", "WT +Phe", "pho84 +Phe", "pho81 +Phe", "vtc4 +Phe"),
    cols         = c("WT -Phe" = "gray50", "WT +Phe" = "#0072BD",
                     "pho84 +Phe" = "#D95319", "pho81 +Phe" = "#77AC30",
                     "vtc4 +Phe" = "#EDB120"),
    factor_order = c("WT -Phe", "WT +Phe", "pho84 +Phe", "pho81 +Phe", "vtc4 +Phe"),
    baseline_ref = "WT -Phe",
    curve_out    = "pho_curves.tiff",
    auc_out      = "pho_auc.tiff",
    auc_w = 5, auc_h = 3.2,
    dunnett_ref  = "WT +Phe"
  ),

  neg_phe_controls = list(
    conditions   = c("WT -Phe", "pho84 -Phe", "pho81 -Phe", "vtc4 -Phe"),
    cols         = c("WT -Phe" = "gray50", "pho84 -Phe" = "#D95319",
                     "pho81 -Phe" = "#77AC30", "vtc4 -Phe" = "#EDB120"),
    factor_order = c("WT -Phe", "pho84 -Phe", "pho81 -Phe", "vtc4 -Phe"),
    baseline_ref = "WT -Phe",
    curve_out    = "pho_curves (neg_phe controls).tiff",
    auc_out      = "pho_auc (neg_phe controls).tiff",
    auc_w = 5.3, auc_h = 3.2,
    dunnett_ref  = NULL
  ),

  pho84_mutant_control = list(
    conditions   = c("pho84 -Phe", "pho84 +Phe"),
    cols         = c("pho84 -Phe" = "gray50", "pho84 +Phe" = "#D95319"),
    factor_order = c("pho84 -Phe", "pho84 +Phe"),
    baseline_ref = "pho84 -Phe",
    curve_out    = "pho_curves (pho84 controls).tiff",
    auc_out      = "pho_auc (pho84 controls).tiff",
    auc_w = 4.3, auc_h = 2.8,
    dunnett_ref  = NULL
  ),

  pho81_mutant_control = list(
    conditions   = c("pho81 -Phe", "pho81 +Phe"),
    cols         = c("pho81 -Phe" = "gray50", "pho81 +Phe" = "#77AC30"),
    factor_order = c("pho81 -Phe", "pho81 +Phe"),
    baseline_ref = "pho81 -Phe",
    curve_out    = "pho_curves (pho81 controls).tiff",
    auc_out      = "pho_auc (pho81 controls).tiff",
    auc_w = 4.3, auc_h = 2.8,
    dunnett_ref  = NULL
  ),

  vtc4_mutant_control = list(
    conditions   = c("vtc4 -Phe", "vtc4 +Phe"),
    cols         = c("vtc4 -Phe" = "gray50", "vtc4 +Phe" = "#EDB120"),
    factor_order = c("vtc4 -Phe", "vtc4 +Phe"),
    baseline_ref = "vtc4 -Phe",
    curve_out    = "pho_curves (vtc4 controls).tiff",
    auc_out      = "pho_auc (vtc4 controls).tiff",
    auc_w = 4.3, auc_h = 2.8,
    dunnett_ref  = NULL
  )
)

##load and clean data:
load_aro_data <- function(data_dir) {
  data <- read.csv(file.path(data_dir, "aro_expts_demo.csv"))
  colnames(data)[1] <- "Sample"
  data <- data %>%
    mutate(
      Sample = str_trim(Sample),
      Sample = str_replace_all(Sample, "\\s+", " "),
      Sample = str_replace_all(Sample, "\\+8 mM", "+8mM"),
      Sample = str_replace_all(Sample, "\\+0\\.9 mM", "+0.9mM"),
      Sample = str_replace_all(Sample, "- Phe", "-Phe")
    ) %>%
    mutate(
      Genotype = case_when(
        str_detect(Sample, "^WT")    ~ "WT",
        str_detect(Sample, "^aro4")  ~ "aro4",
        str_detect(Sample, "^aro9")  ~ "aro9",
        str_detect(Sample, "^aro10") ~ "aro10"
      ),
      Treatment = case_when(
        str_detect(Sample, "-Phe") ~ "-Phe",
        str_detect(Sample, "8mM")  ~ "+8mM",
        str_detect(Sample, "0.9")  ~ "+0.9mM"
      )
    )
}

load_pho_data <- function(data_dir) {
  data <- read.csv(file.path(data_dir, "pho_expts_demo.csv"))
  colnames(data)[1] <- "Sample"
  data$Sample[data$Sample == "WT colony #2"] <- "WT"
  data$Sample[data$Sample == "WT colony #1"] <- "WT"
  data <- data %>%
    filter(Sample %in% c("WT", "WT\nPhe 8mM",
                         "\u0394pho81", "\u0394pho81\nPhe 8mM",
                         "\u0394pho84", "\u0394pho84\nPhe 8mM",
                         "\u0394vtc4",  "\u0394vtc4\nPhe 8mM")) %>%
    mutate(
      Sample = str_trim(Sample),
      Sample = str_replace_all(Sample, "\nPhe 8M",   "+Phe"),
      Sample = str_replace_all(Sample, "\nPhe 8mM",  "+Phe"),
      Sample = str_replace_all(Sample, "\nPhe 8 mM", "+Phe")
    ) %>%
    mutate(
      Genotype = case_when(
        str_detect(Sample, "^WT")   ~ "WT",
        str_detect(Sample, "pho81") ~ "pho81",
        str_detect(Sample, "pho84") ~ "pho84",
        str_detect(Sample, "vtc4")  ~ "vtc4"
      ),
      Treatment = case_when(
        str_detect(Sample, "\\+Phe") ~ "+Phe",
        TRUE ~ "-Phe"
      )
    )
}

to_long <- function(data_clean, max_time = 30) {
  data_clean %>%
    pivot_longer(cols = starts_with("X"), names_to = "Time", values_to = "Value") %>%
    mutate(
      Time  = as.numeric(str_remove(Time, "X")),
      Value = as.numeric(Value)
    ) %>%
    filter(Time <= max_time)
}

##plot functions:
plot_curves <- function(data_long, panel, out_path) {
  summary_df <- data_long %>%
    group_by(Time, Genotype, Treatment) %>%
    summarise(mean_value = mean(Value, na.rm = TRUE),
              sd_value   = sd(Value,   na.rm = TRUE),
              .groups = "drop") %>%
    mutate(strain_condition = paste(Genotype, Treatment)) %>%
    filter(strain_condition %in% panel$conditions) %>%
    mutate(strain_condition = factor(strain_condition, levels = panel$factor_order))

  p <- ggplot(summary_df,
              aes(x = Time, y = mean_value,
                  color = strain_condition, fill = strain_condition)) +
    geom_ribbon(aes(ymin = mean_value - sd_value,
                    ymax = mean_value + sd_value),
                alpha = 0.2, color = NA) +
    geom_line(linewidth = 1.2, alpha = 0.85) +
    scale_color_manual(values = panel$cols) +
    scale_fill_manual(values  = panel$cols) +
    theme_minimal(base_size = 14)

  ggsave(out_path, p, width = 6, height = 3.5, dpi = 600, compression = "lzw")
  invisible(p)
}

compute_auc <- function(data_long, panel) {
  ##keep only conditions in this panel
  all_conditions <- panel$conditions
  auc_df <- data_long %>%
    mutate(strain_condition = paste(Genotype, Treatment)) %>%
    filter(strain_condition %in% all_conditions) %>%
    group_by(Genotype, Treatment, Replicate) %>%
    summarise(AUC = trapz(Time, Value), .groups = "drop")

  ##baseline = the condition labelled as baseline_ref
  baseline_geno  <- str_extract(panel$baseline_ref, "^\\S+")
  baseline_treat <- str_remove(panel$baseline_ref, "^\\S+\\s")
  baseline_auc <- auc_df %>%
    filter(Genotype == baseline_geno, Treatment == baseline_treat) %>%
    dplyr::select(Replicate, AUC_baseline = AUC)

  auc_ratio <- auc_df %>%
    left_join(baseline_auc, by = "Replicate") %>%
    mutate(
      fold_change     = AUC / AUC_baseline,
      strain_condition = paste(Genotype, Treatment)
    ) %>%
    filter(strain_condition %in% all_conditions) %>%
    mutate(strain_condition = factor(strain_condition, levels = panel$factor_order))

  auc_ratio
}

plot_auc <- function(auc_ratio, panel, out_path) {
  summary_auc <- auc_ratio %>%
    group_by(strain_condition) %>%
    summarise(mean_fc = mean(fold_change, na.rm = TRUE),
              se_fc   = sd(fold_change,   na.rm = TRUE) / sqrt(n()),
              .groups = "drop") %>%
    mutate(strain_condition = factor(strain_condition, levels = panel$factor_order))

  p1 <- ggplot(summary_auc,
               aes(x = strain_condition, y = mean_fc, fill = strain_condition)) +
    geom_bar(stat = "identity", color = "black", alpha = 0.6) +
    geom_errorbar(aes(ymin = mean_fc - se_fc, ymax = mean_fc + se_fc),
                  width = 0.2, linewidth = 0.8) +
    scale_fill_manual(values = panel$cols) +
    theme_minimal(base_size = 14)

  ggsave(out_path, p1,
         width = panel$auc_w, height = panel$auc_h,
         dpi = 600, compression = "lzw")
  invisible(p1)
}

run_dunnett <- function(auc_ratio, ref_level = NULL) {
  ## if ref_level is provided, relevel so that becomes the reference
  df <- auc_ratio
  if (!is.null(ref_level)) {
    other_levels <- setdiff(levels(df$strain_condition), ref_level)
    df$strain_condition <- factor(df$strain_condition,
                                  levels = c(ref_level, other_levels))
  }
  aov_res     <- aov(fold_change ~ strain_condition, data = df)
  dunnett_res <- glht(aov_res, linfct = mcp(strain_condition = "Dunnett"))
  cat("\n--- ANOVA summary (ref:", levels(df$strain_condition)[1], ") ---\n")
  print(summary(aov_res))
  cat("\n--- Dunnett ---\n")
  print(summary(dunnett_res))
  invisible(dunnett_res)
}

##run full analysis:
run_analysis <- function(experiment = EXPERIMENT,
                         data_dir   = DATA_DIR,
                         plot_dir   = PLOT_DIR) {
  if (experiment == "aro") {
    raw        <- load_aro_data(data_dir)
    data_long  <- to_long(raw)
    panels     <- ARO_PANELS
  } else if (experiment == "pho") {
    raw        <- load_pho_data(data_dir)
    data_long  <- to_long(raw)
    panels     <- PHO_PANELS
  } else {
    stop("experiment must be 'aro' or 'pho'")
  }

  for (panel_name in names(panels)) {
    panel <- panels[[panel_name]]
    cat("\n========== Panel:", panel_name, "==========\n")

    ## growth curve
    plot_curves(data_long, panel,
                file.path(plot_dir, panel$curve_out))

    ## AUC
    auc_ratio <- compute_auc(data_long, panel)

    if (!is.null(panel$auc_out)) {
      plot_auc(auc_ratio, panel,
               file.path(plot_dir, panel$auc_out))
    }

    ##stats: primary Dunnett (default reference = first factor level)
    run_dunnett(auc_ratio, ref_level = NULL)

    ## secondary Dunnett if a second reference is specified
    if (!is.null(panel$dunnett_ref)) {
      run_dunnett(auc_ratio, ref_level = panel$dunnett_ref)
    }
  }

  ##aro moderate-conc combined AUC plot (both genotypes together)
  if (experiment == "aro") {
    cat("\n========== Combined moderate-conc AUC ==========\n")
    auc_mc10 <- compute_auc(data_long, panels$moderate_conc_aro10)
    auc_mc4  <- compute_auc(data_long, panels$moderate_conc_aro4)
    auc_comb <- bind_rows(auc_mc10, auc_mc4) %>%
      distinct(Genotype, Treatment, Replicate, .keep_all = TRUE)

    all_mc_conds  <- c("WT +0.9mM", "aro10 +0.9mM", "aro10 +8mM",
                       "aro4 +0.9mM", "aro4 +8mM")
    all_mc_cols   <- c("WT +0.9mM" = "gray50", "aro10 +0.9mM" = "#0072BD",
                       "aro10 +8mM" = "#77AC30", "aro4 +0.9mM" = "#D95319",
                       "aro4 +8mM" = "#EDB120")

    summary_mc <- auc_comb %>%
      filter(strain_condition %in% all_mc_conds) %>%
      mutate(strain_condition = factor(strain_condition, levels = all_mc_conds),
             group = gsub(" .*", "", as.character(strain_condition)),
             group = factor(group, levels = c("WT", "aro10", "aro4"))) %>%
      group_by(strain_condition, group) %>%
      summarise(mean_fc = mean(fold_change, na.rm = TRUE),
                se_fc   = sd(fold_change,   na.rm = TRUE) / sqrt(n()),
                .groups = "drop")

    p_mc <- ggplot(summary_mc,
                   aes(x = strain_condition, y = mean_fc,
                       fill = strain_condition, group = group)) +
      geom_bar(stat = "identity", color = "black",
               position = position_dodge(width = 0.75), alpha = 0.6) +
      geom_errorbar(aes(ymin = mean_fc - se_fc, ymax = mean_fc + se_fc),
                    width = 0.2, linewidth = 0.8) +
      scale_fill_manual(values = all_mc_cols) +
      theme_minimal(base_size = 14)

    ggsave(file.path(plot_dir, "moderateconc_auc.tiff"),
           p_mc, width = 5, height = 3.2, dpi = 600, compression = "lzw")

    ## Dunnett for combined moderate-conc
    auc_comb_f <- auc_comb %>%
      filter(strain_condition %in% all_mc_conds) %>%
      mutate(strain_condition = factor(strain_condition,
                                       levels = c("WT +0.9mM", all_mc_conds[-1])))
    run_dunnett(auc_comb_f, ref_level = NULL)
  }

  cat("\nDone. Plots saved to:", plot_dir, "\n")
}

run_analysis()

