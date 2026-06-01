library(dplyr)
library(tidyr)
library(ggplot2)
library(pracma) 

##import data:
curves = read.delim('./growth-curves/data/phe_curves.txt')
cell_conc = read.delim('./growth-curves/data/FEBS_cell-conc.txt')

##auc:
auc_data <- curves %>%
  group_by(Strain, Dose) %>%
  summarise(AUC = trapz(Time, OD), .groups = "drop") %>%
  group_by(Strain) %>%
  mutate(AUC_norm = AUC / AUC[Dose == 0] * 100) %>%
  ungroup()

##plot auc barplot:
auc_data$Strain = factor(auc_data$Strain, levels = c("WT", "aro4Δ"))
p = ggplot(auc_data, aes(x = factor(Dose), y = AUC_norm, fill = Strain, group=Strain)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.85), width = 0.7, alpha=0.6, color='black') +
  labs(
    x = "Extracellular Phe (mM)",
    y = "Growth (% of control)",
    title = "Dose-dependent growth in WT vs aro4Δ"
  ) +
  scale_fill_manual(values = c("WT" = "#0072BD", "aro4Δ" = "#D95319")) +
  theme_minimal()

##save:
ggsave("./growth-curves/plots/curve_auc.tiff", p, 
       width = 4.5, height = 3, dpi = 600, compression = "lzw")

##intracellular Phe (FEBS) mean ± SE; fold-change:
conc_summary <- cell_conc %>%
  group_by(Strain, Phe) %>%
  summarise(
    mean_mM = mean(mM.per.cell, na.rm = TRUE),
    se_mM   = sd(mM.per.cell, na.rm = TRUE)/sqrt(n()),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = Phe, values_from = c(mean_mM, se_mM)) %>%
  mutate(
    fold_change = mean_mM_Yes / mean_mM_No,
    log2FC = log2(fold_change),
    se_log2FC = se_mM_Yes / (mean_mM_Yes * log(2))
  )

##plot:
conc_summary$Strain = factor(conc_summary$Strain, levels = c("WT", "aro4"))
p1 = ggplot(conc_summary, aes(x = Strain, y = log2FC, fill = Strain)) +
  geom_bar(stat = "identity", width = 0.5, alpha = 0.6, color='black') +
  geom_errorbar(aes(ymin = log2FC - se_log2FC, ymax = log2FC + se_log2FC),
                width = 0.2) +
  labs(y = "Intracellular Phe log2FC") +
  scale_fill_manual(values = c("WT" = "#0072BD", "aro4" = "#D95319")) +
  theme_minimal()

##save:
ggsave("./growth-curves/plots/phe_conc.tiff", p1, 
       width = 3, height = 2.5, dpi = 600, compression = "lzw")

