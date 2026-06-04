library(dplyr)
library(tidyr)
library(pracma)

##import kinetics data:
##data in long format, with columns: Batch, Time, Condition, Growth
data_long = read.csv('./synergy-analysis/data/kinetics.csv')

##calculate AUC for each curve:
auc_results = data_long %>%
  group_by(Batch, Condition) %>%
  summarise(
    AUC = trapz(Time, Growth),
    .groups = "drop"
  )

write.csv(auc_results, './synergy-analysis/data/auc.csv',row.names=FALSE)

##calculate AUC sd:
auc_summary = auc_results %>%
  group_by(Condition) %>%
  summarise(
    mean_AUC = mean(AUC),
    sd_AUC = sd(AUC),
    .groups = "drop"
  )

write.csv(auc_summary, './synergy-analysis/data/auc_summary.csv',row.names=FALSE)

