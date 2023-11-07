### SCRIPT: Density plots BMN

## 07.11.23 Laura Sudupe , git @lsudupe

###### Libraries
library("Seurat")
library("ggplot2")
library("RColorBrewer")

#Data---------------------------------
se <- readRDS("./objects/heterogeneity/se_hierarchical.rds")

metadata <- as.data.frame(se@meta.data)
metadata <- subset(metadata, !(name %in% c("M3_F_1C", "M3_fem_1C")))

# If necessary, convert the 'name' and 'MM_MIC' columns to appropriate data types
metadata$name <- as.factor(metadata[["name"]])
metadata$MM_MIC <- as.numeric(metadata[["MM_MIC"]])

# Create the density plot
p <- ggplot(metadata, aes(x = MM_MIC, fill = name)) +
  geom_density(alpha = 0.5) +
  labs(x = "MM-PC Value", y = "Density", title = "Density Plot of MM-PC per Sample") +
  theme_minimal()

# Display the plot
print(p)
# Save the plot to a PDF file
ggsave(filename = "./results/density_plot/density_plot.pdf", plot = p)

####separate
p <- ggplot(metadata, aes(x = MM_MIC, fill = name)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(0, 0.7), breaks = seq(0, 0.7, by = 0.1)) +
  labs(x = "MM-PC Value", y = "Density") +
  theme_minimal() +
  facet_wrap(~ name, scales = "fixed") # Changed to 'fixed' to have the same scales

# Save the plot to a PDF file
ggsave(filename = "./results/density_plot/density_plots_colored.pdf", plot = p, width = 10, height = 8)


print(p)


