library(ggplot2)

##import data:
data = read.delim('./growth-curves/data/example_growth-curve.txt')

##set colour palette:
color_list <- c("#0072BD", "#D95319", "#77AC30", "#EDB120")

##plot:
p = ggplot(data, aes(x = Time, y = Growth, color = factor(Concentration))) +
  geom_line(alpha = 0.6, linewidth = 1.4) +
  scale_color_manual(values = color_list, name = 'Concentration') +
  labs(x = "Time", y = "Cell Density") +
  theme_minimal() +
  theme(
    legend.position = c(0.05, 0.95),
    legend.justification = c(0, 1),
    legend.background = element_rect(
      fill = "white",    
      color = "black", 
      size = 0.2, 
      linetype = "solid"
    )
  )

##save:
#ggsave("./growth-curves/plots/example_growth-curve.tiff", p, 
#       width = 4.5, height = 3.5, dpi = 600, compression = "lzw")



