library(dplyr)
library(igraph)
library(ggplot2)

## ============================================================
## PATHS
## ============================================================

OUT_DIR     <- './network/data/'
NETWORK_DIR <- './network/data/'
PLOT_DIR    <- './network/plots/'
PROPS_FILE  <- './property-analysis/data/AA_properties.csv'

##helper functions:
min_max_norm <- function(x, new_min = 0, new_max = 1) {
  (x - min(x, na.rm = TRUE)) /
    (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) *
    (new_max - new_min) + new_min
}

##build graph & run pagerank:
edgelist <- read.delim(file.path(NETWORK_DIR, 'edgelist_final.txt'))
g           <- graph_from_data_frame(edgelist, directed = TRUE)
E(g)$weight <- edgelist$norm_weight
ppr <- page_rank(g, algo = "prpack", directed = TRUE,
                 weights = E(g)$weight, damping = 0.85)
ppr_scores <- data.frame(node = names(ppr$vector), ppr_score = ppr$vector)
ppr_scores$norm_ppr  <- min_max_norm(ppr_scores$ppr_score)
cutoff               <- mean(ppr_scores$norm_ppr)

##save ppr scores:
ppr_scores = ppr_scores[,c("node","ppr_score","norm_ppr")]
write.table(ppr_scores, file.path(OUT_DIR, 'ppr_scores.txt'), sep = '\t', row.names = FALSE)

##group:
ppr_scores$ppr_group <- ifelse(ppr_scores$norm_ppr > cutoff, "Central", "Peripheral")

##plot influence density:
dens_list  <- lapply(split(ppr_scores$norm_ppr, ppr_scores$ppr_group), density)
mean_lines <- do.call(rbind, lapply(names(dens_list), function(grp) {
  dens     <- dens_list[[grp]]
  mean_val <- mean(ppr_scores$norm_ppr[ppr_scores$ppr_group == grp])
  data.frame(ppr_group = grp, mean = mean_val,
             ymax = approx(dens$x, dens$y, xout = mean_val)$y)
}))

p_density <- ggplot(ppr_scores, aes(x = norm_ppr, fill = ppr_group)) +
  geom_density(alpha = 0.6, color = "black", linewidth = 0.8) +
  geom_segment(data = mean_lines,
               aes(x = mean, xend = mean, y = 0, yend = ymax),
               color = "#555555", linetype = "dashed", linewidth = 0.8) +
  scale_fill_manual(values = c("Central" = "#0072BD", "Peripheral" = "#D95319")) +
  theme_minimal() +
  theme(legend.position = 'none') +
  labs(x = "Influence (Normalised)", y = "Density")

#ggsave(file.path(PLOT_DIR, 'influence_density.tiff'),
#       p_density, width = 4.4, height = 3.8, dpi = 600)

##plot ppr score against log10(failure threshold):
thresh <- read.csv(PROPS_FILE) %>%
  select(Name3let, Mean_Failure_Thresh) %>%
  distinct() %>%
  rename(node = Name3let)

temp <- left_join(ppr_scores, thresh, by = 'node') %>%
  filter(!is.na(Mean_Failure_Thresh)) %>%
  mutate(ppr_score_norm = min_max_norm(ppr_score))

cor.test(log10(temp$Mean_Failure_Thresh), temp$ppr_score_norm)

p_corr <- ggplot(temp, aes(x = ppr_score_norm, y = log10(Mean_Failure_Thresh))) +
  geom_point(color = "#0072BD", size = 6, alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "#0072BD", linewidth = 1) +
  theme_minimal(base_size = 14) +
  labs(x = "minmax(PPR score)", y = "log10(Failure Threshold)")

#ggsave(file.path(PLOT_DIR, 'corr_ppr-vs-failurethresh.tiff'),
#       p_corr, width = 5.4, height = 5, dpi = 600)

