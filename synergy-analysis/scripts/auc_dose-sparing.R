library(dplyr)
library(tidyr)
library(pracma)
library(broom)
library(ggplot2)
library(gridExtra)

##import AUC data:
auc_results = read.csv('./synergy-analysis/data/auc.csv')
auc_summary = read.csv('./synergy-analysis/data/auc_summary.csv')

##calculate growth:
auc_results <- auc_results %>%
  left_join(
    auc_results %>% 
      filter(Condition == "WT") %>% 
      select(Batch, AUC_control = AUC),
    by = "Batch"
  ) %>%
  mutate(growth = (AUC / AUC_control))

##define dose-sparing comparisons:
dose_sparing_list <- list(
  "Phe+Ile" = list(
    singles = c("Phe.8mM", "Iso.8mM"),
    combo   = "Phe.4mM...Iso.4mM"
  ),
  "Phe+Thr" = list(
    singles = c("Phe.8mM", "Thr.8mM"),
    combo   = "Phe.4mM..Thr.4mM"
  ),
  "Thr+Ile" = list(
    singles = c("Thr.8mM", "Iso.8mM"),
    combo   = "Thr.4mM...Iso.4mM"
  )
)

##plotting dataframe:
plot_data_all <- data.frame()
for (grp in names(dose_sparing_list)) {
  singles <- dose_sparing_list[[grp]]$singles
  combo   <- dose_sparing_list[[grp]]$combo
  # Singles (8 mM)
  singles_df <- auc_results %>%
    filter(Condition %in% singles) %>%
    group_by(Condition) %>%
    summarise(
      mean_growth = mean(growth),
      SE_growth   = sd(growth)/sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(Group = grp)
  # Combo (4+4 mM)
  combo_df <- auc_results %>%
    filter(Condition == combo) %>%
    summarise(
      mean_growth = mean(growth),
      SE_growth   = sd(growth)/sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(
      Condition = "Combo (4+4 mM)",
      Group = grp
    )
  plot_data_all <- bind_rows(plot_data_all, singles_df, combo_df)
}
plot_data_all$Condition <- factor(
  plot_data_all$Condition,
  levels = c("Phe.8mM", "Iso.8mM", "Thr.8mM", "Combo (4+4 mM)")
)

##plot:
colors <- c("Phe.8mM" = "#D95319","Iso.8mM" = "#0072BD","Thr.8mM" = "#EDB120","Combo (4+4 mM)" = "gray50")
dose_sparing_plot <- ggplot(plot_data_all,
                            aes(x = Condition, y = mean_growth, fill = Condition)) +
  geom_bar(stat = "identity", width = 0.5,
           color = "black", alpha = 0.6) +
  geom_errorbar(aes(ymin = mean_growth,
                    ymax = pmin(1, mean_growth + SE_growth)),
                width = 0.2) +
  scale_fill_manual(values = colors) +
  labs(y = "Growth", x = "") +
  facet_wrap(~Group, scales = "free_x") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(size = 11),
        legend.position = "none")

##save:
ggsave("./synergy-analysis/plots/dose_sparing_Growth.tiff", dose_sparing_plot,
  width = 8.5,height = 3,dpi = 600,compression = "lzw")

##test:
stats_data <- data.frame()
for (grp in names(dose_sparing_list)) {
  singles <- dose_sparing_list[[grp]]$singles
  combo   <- dose_sparing_list[[grp]]$combo
  tmp <- auc_results %>%
    filter(Condition %in% c(singles, combo)) %>%
    mutate(
      Group = grp,
      Condition = case_when(
        Condition == combo ~ "Combo",
        TRUE ~ Condition
      )
    )
  stats_data <- bind_rows(stats_data, tmp)
}

stats_results <- data.frame()
for (grp in unique(stats_data$Group)) {
  df <- stats_data %>% filter(Group == grp)
  combo_vals <- df %>% filter(Condition == "Combo") %>% pull(growth)
  singles <- unique(df$Condition[df$Condition != "Combo"])
  for (s in singles) {
    single_vals <- df %>% filter(Condition == s) %>% pull(growth)
    test <- t.test(combo_vals, single_vals)
    stats_results <- bind_rows(
      stats_results,
      data.frame(
        Group = grp,
        Comparison = paste("Combo vs", s),
        p_value = test$p.value,
        mean_combo = mean(combo_vals),
        mean_single = mean(single_vals),
        sd_single = sd(single_vals)
      )
    )
  }
}
stats_results$p_adj <- p.adjust(stats_results$p_value, method = "BH")

write.csv(stats_results, './synergy-analysis/data/dosesparing_stats.csv',row.names=FALSE)
