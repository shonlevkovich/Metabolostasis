library(dplyr)
library(tidyr)
library(pracma)
library(broom)
library(ggplot2)
library(ggbreak)
library(gridExtra)

##import AUC data:
auc_results = read.csv('./synergy-analysis/data/auc.csv')
auc_summary = read.csv('./synergy-analysis/data/auc_summary.csv')

##calculate fractional inhibition:
auc_results <- auc_results %>%
  left_join(
    auc_results %>% 
      filter(Condition == "WT") %>% 
      select(Batch, AUC_control = AUC),
    by = "Batch"
  ) %>%
  mutate(Fractional_Inhibition = 1 - (AUC / AUC_control))

##define combinations of interest:
combos_list <- list(
  "Phe.4mM...Iso.4mM" = c("Phe.4mM", "Iso.4mM"),
  "Phe.4mM..Thr.4mM"  = c("Phe.4mM", "Thr.4mM"),
  "Thr.4mM...Iso.4mM" = c("Thr.4mM", "Iso.4mM")
)

##calculate ΔBliss per replicate
delta_bliss_df <- data.frame()
for (combo_name in names(combos_list)) {
  singles <- combos_list[[combo_name]]
  A <- singles[1]
  B <- singles[2]
  # Extract replicate FI values
  df_A <- auc_results %>%
    filter(Condition == A) %>%
    select(Batch, FI_A = Fractional_Inhibition)
  df_B <- auc_results %>%
    filter(Condition == B) %>%
    select(Batch, FI_B = Fractional_Inhibition)
  df_combo <- auc_results %>%
    filter(Condition == combo_name) %>%
    select(Batch, FI_combo = Fractional_Inhibition)
  # Merge by replicate
  merged <- df_A %>%
    left_join(df_B, by = "Batch") %>%
    left_join(df_combo, by = "Batch")
  # Expected Bliss
  merged <- merged %>%
    mutate(
      FI_bliss = FI_A + FI_B - (FI_A * FI_B),
      delta_bliss = FI_combo - FI_bliss,
      Combo = combo_name
    )
  delta_bliss_df <- bind_rows(delta_bliss_df, merged)
}

##t-tests:
synergy_stats <- delta_bliss_df %>%
  group_by(Combo) %>%
  summarise(
    mean_delta = mean(delta_bliss),
    t_test = list(t.test(delta_bliss, mu = 0)),
    .groups = "drop"
  ) %>%
  mutate(
    tidy = lapply(t_test, tidy)
  ) %>%
  unnest(tidy) %>%
  select(Combo, mean_delta, p.value, statistic, conf.low, conf.high)
synergy_stats$p_adj <- p.adjust(synergy_stats$p.value, method = "BH")

write.csv(synergy_stats, './synergy-analysis/data/synergy_stats.csv',row.names=FALSE)

##plot dataframe:
plot_data_all <- data.frame()
for (combo_name in names(combos_list)) {
  single_conditions <- combos_list[[combo_name]]
  # Single compounds
  singles <- auc_results %>%
    filter(Condition %in% single_conditions) %>%
    group_by(Condition) %>%
    summarise(
      mean_FI = mean(Fractional_Inhibition),
      SE_FI   = sd(Fractional_Inhibition) / sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(Type = "Single",
           Combo_Group = combo_name)
  # Observed combination
  combo_df <- auc_results %>%
    filter(Condition == combo_name) %>%
    summarise(
      mean_FI = mean(Fractional_Inhibition),
      SE_FI   = sd(Fractional_Inhibition) / sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(Condition = "Combo",
           Type = "Combo",
           Combo_Group = combo_name)
  # Combine singles + combo
  plot_data_all <- bind_rows(plot_data_all, singles, combo_df)
}
plot_data_all$Condition = factor(plot_data_all$Condition, levels=c('Phe.4mM','Iso.4mM','Thr.4mM','Combo'))

##compute expected Bliss inhibition per combo
##bliss = FI_A + FI_B - FI_A*FI_B
bliss_df <- data.frame()
for (combo_name in names(combos_list)) {
  conds <- combos_list[[combo_name]]
  FI_A <- plot_data_all %>% filter(Condition == conds[1]) %>% pull(mean_FI)
  FI_B <- plot_data_all %>% filter(Condition == conds[2]) %>% pull(mean_FI)
  bliss_df <- (rbind(bliss_df, data.frame(Combo_Group = combo_name, Bliss_mean = FI_A + FI_B - FI_A*FI_B)))
}

##plot:
colors <- c("Phe.4mM" = "#D95319","Iso.4mM" = "#0072BD","Thr.4mM" = "#EDB120","Combo"   = "gray50")
bliss_plot <- ggplot(plot_data_all, aes(x = Condition, y = mean_FI, fill = Condition)) +
  geom_bar(stat = "identity", width = 0.5, color = "black", alpha=.6) +
  geom_errorbar(aes(ymin = mean_FI, ymax = pmin(1, mean_FI + SE_FI)), width = 0.2) +
  geom_hline(data = bliss_df, aes(yintercept = Bliss_mean), linetype = "dotted", color = "black", size = 1) +
  scale_fill_manual(values = colors) +
  labs(y = "Fractional Inhibition (FI)", x = "", fill = "") +
  facet_wrap(~Combo_Group, scales = "free_x") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(size = 12), legend.position = "none") +
  scale_y_break(c(0.2, 0.7))   # break from 0.05 to 0.25

##save:
ggsave("./synergy-analysis/plots/bliss.tiff", bliss_plot, 
       width = 8.5, height = 3.5, dpi = 600, compression = "lzw")
