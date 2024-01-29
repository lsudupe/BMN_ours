## SCRIPT: Mouse quantile

## 28.11.23 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(ggplot2)
library(STutility)
library(gridExtra)
library(dplyr)
library(tidyr)
library(cowplot) 
library(ggpubr)
library(RColorBrewer)
color <- rev(brewer.pal(11,"Spectral"))

mouse <- readRDS("./objects/sp/se_hierarchical_signatures.rds")
df <- mouse@meta.data


# Convert PC_new_UCell and MM_MIC to numeric if they are not
df$PC_new_UCell <- as.numeric(df$PC_new_UCell)
df$MM_MIC <- as.numeric(df$MM_MIC)

# Convert name to a factor if it's not already
df$name <- as.factor(df$name)

# Get unique names
unique_names <- unique(df$name)

# Create a list to store plots
plot_list_PC <- list()
plot_list_MM <- list()

# Loop over unique names to create histograms for PC_new_UCell
for (n in unique_names) {
  plot_list_PC[[n]] <- ggplot(df[df$name == n, ], aes(x = PC_new_UCell)) +
    geom_histogram(bins = 30, fill = "blue", color = "black") +
    ggtitle(paste("PC_new_UCell for", n)) +
    theme(plot.title = element_text(size = 8))
}

# Loop over unique names to create histograms for MM_MIC
for (n in unique_names) {
  plot_list_MM[[n]] <- ggplot(df[df$name == n, ], aes(x = MM_MIC)) +
    geom_histogram(bins = 30, fill = "red", color = "black") +
    ggtitle(paste("MM_MIC for", n)) +
    theme(plot.title = element_text(size = 8))  # Adjust title size
}

# Combine plots into a grid
grid_layout_PC <- do.call(grid.arrange, c(plot_list_PC, ncol = length(unique_names)))
grid_layout_MM <- do.call(grid.arrange, c(plot_list_MM, ncol = length(unique_names)))

# Optionally, you can save the plots to a file
ggsave("./results/new_grouping/PC_new_UCell_histograms.pdf", grid_layout_PC, width = 15, height = 4)
ggsave("./results/new_grouping/MM_MIC_histograms.pdf", grid_layout_MM, width = 15, height = 4)  


# Create a combined histogram for MM_MIC
histogram_MM_MIC <- ggplot(df, aes(x = MM_MIC, fill = name)) +
  geom_histogram(bins = 100, position = "identity", alpha = 0.6) +
  scale_fill_brewer(palette = "Set1") +  # Change color palette as needed
  ggtitle("Combined Histogram of MM_MIC") +
  theme_minimal()

# Create a combined histogram for PC_new_UCell
histogram_PC_new_UCell <- ggplot(df, aes(x = PC_new_UCell, fill = name)) +
  geom_histogram(bins = 100, position = "identity", alpha = 0.6) +
  scale_fill_brewer(palette = "Set1") +  # Change color palette as needed
  ggtitle("Combined Histogram of PC_new_UCell") +
  theme_minimal()

# Save the plots to PDF
ggsave("./results/new_grouping/Combined_MM_MIC_histogram.pdf", histogram_MM_MIC, width = 10, height = 6)
ggsave("./results/new_grouping/Combined_PC_new_UCell_histogram.pdf", histogram_PC_new_UCell, width = 10, height = 6)



###SEPARATE DATAs
# Calcular el percentil 0.75 en las muestras M3_F_1C y M3_fem_1C
percentile_90 <- quantile(df$PC_new_UCell[df$name %in% c("M3_F_1C", "M3_fem_1C")], probs = 0.90)
# Usar el valor del percentil 0.75 como umbral para categorizar todas las filas
df <- df %>%
  mutate(new_groups = ifelse(PC_new_UCell < percentile_90, "rest", "notknown"))


st <- mouse
st <- AddMetaData(st, df)
FeatureOverlay(st, features = c("new_groups"), sampleids = 1:6, ncols = 2, 
               pt.size = 0.7, cols = c(rest="blue", notknown="white"))
#st <- SubsetSTData(st, idents = c("M1_fem_1C", "M2_F_2B", "M8_F2_1C", "M9_F2_1C"))

##
df_filtered <- st@meta.data
#df_filtered <- df %>% filter(name %in% c("M3_F_1C", "M3_fem_1C"))

# Crear una lista para almacenar los histogramas
plot_list_PC <- list()
plot_list_Neutrophils <- list()
plot_list_Tcell_Exh_UCell<- list()

# Función para agregar cuantiles a un histograma
add_quantiles <- function(plot, data, variable) {
  quantiles <- quantile(data[[variable]], probs = c(0.25, 0.5, 0.75))
  plot + geom_vline(xintercept = quantiles[1], linetype = "dashed", color = "red") +
    geom_vline(xintercept = quantiles[2], linetype = "solid", color = "green") +
    geom_vline(xintercept = quantiles[3], linetype = "dotted", color = "blue")
}

# Bucle para crear histogramas con cuantiles para PC_new_UCell
for (n in unique(df_filtered$name)) {
  plot <- ggplot(df_filtered[df_filtered$name == n, ], aes(x = PC_new_UCell)) +
    geom_histogram(bins = 30, fill = "blue", color = "black") +
    ggtitle(paste("PC_new_UCell for", n))
  plot_list_PC[[n]] <- add_quantiles(plot, df_filtered[df_filtered$name == n, ], "PC_new_UCell")
}

# Bucle para crear histogramas con cuantiles para Neutrophils
for (n in unique(df_filtered$name)) {
  plot <- ggplot(df_filtered[df_filtered$name == n, ], aes(x = Neutrophils)) +
    geom_histogram(bins = 30, fill = "green", color = "black") +
    ggtitle(paste("Neutrophils for", n))
  plot_list_Neutrophils[[n]] <- add_quantiles(plot, df_filtered[df_filtered$name == n, ], "Neutrophils")
}

# Bucle para crear histogramas con cuantiles para Tcell_Exh_UCell
for (n in unique(df_filtered$name)) {
  plot <- ggplot(df_filtered[df_filtered$name == n, ], aes(x = Tcell_Exh_UCell)) +
    geom_histogram(bins = 30, fill = "purple", color = "black") +
    ggtitle(paste("Tcell_Exh_UCell for", n))
  plot_list_Tcell_Exh_UCell[[n]] <- add_quantiles(plot, df_filtered[df_filtered$name == n, ], "Tcell_Exh_UCell")
}

# Combinar y organizar los histogramas en una cuadrícula
grid_layout_PC <- do.call(grid.arrange, c(plot_list_PC, ncol = length(unique(df_filtered$name))))
grid_layout_Neutrophils <- do.call(grid.arrange, c(plot_list_Neutrophils, ncol = length(unique(df_filtered$name))))
grid_layout_Tcell_Exh_UCell <- do.call(grid.arrange, c(plot_list_Tcell_Exh_UCell, ncol = length(unique(df_filtered$name))))

# Guardar los histogramas en archivos PDF
ggsave("./results/new_grouping/PC_new_UCell_histograms_quantiles_2.pdf", grid_layout_PC, width = 15, height = 4)
ggsave("./results/new_grouping/Neutrophils_histograms_quantiles_2.pdf", grid_layout_Neutrophils, width = 15, height = 4)
ggsave("./results/new_grouping/Tcell_Exh_UCell_histograms_quantiles_2.pdf", grid_layout_Tcell_Exh_UCell, width = 15, height = 4)


## seprate Q3
# Calcular los cuantiles 0.25, 0.5, y 0.75 de PC_new_UCell, excluyendo 'rest' rows
df_filtered <- df_filtered %>%
  mutate(row_name = rownames(df_filtered))

# Calcular los cuantiles por cada muestra en 'name', excluyendo 'rest'
df_quantiles <- df_filtered %>%
  filter(new_groups != "rest") %>%
  group_by(name) %>%
  summarise(
    Q1 = quantile(PC_new_UCell, probs = 0.25),
    Q2 = quantile(PC_new_UCell, probs = 0.5),
    Q3 = quantile(PC_new_UCell, probs = 0.75)
  ) %>%
  ungroup() # Desagrupar para evitar problemas en la unión

# Unir los cuantiles con el dataframe original
df_filtered_q3 <- df_filtered %>%
  left_join(df_quantiles, by = "name") %>%
  mutate(new_groups = case_when(
    new_groups == "rest" ~ "rest",
    PC_new_UCell < Q1 ~ "Remote",
    PC_new_UCell < Q2 ~ "Border",
    TRUE ~ "Hot_spot"
  )) %>%
  select(-c(Q1, Q2, Q3))  # Elimina las columnas de cuantiles

# Restaurar los nombres de fila desde la columna guardada
rownames(df_filtered_q3) <- df_filtered_q3$row_name
df_filtered_q3 <- df_filtered_q3 %>% select(-row_name) 

# Asumiendo que st es un objeto Seurat y df_filtered_q3 es el metadato a agregar
st_q3 <- AddMetaData(st, df_filtered_q3)

##plot areas colorblind
#library(colorBlindness)
#selectedColors <- c(rest="#E5B17E", Remote="#C3E57E", Border="#7EC3E5", Hot_spot="#CC5151",notknown="white")

p <-FeatureOverlay(st_q3, features = c("new_groups"),  ncols = 2, sampleids = 1:6,
               pt.size =0.6, cols = c(rest = "#E5B17E",   # Azul oscuro
                                      Remote = "#C3E57E",    # Gris claro
                                      Border = "#7EC3E5", # Rojo vino
                                       Hot_spot = "#CC5151")) # Verde esmeralda

saveRDS(st_q3, "./objects/sp/st_q3_all.rds")

pdf(paste("./results/distribution/spatial_areas_distribution.pdf",sep=""))
print(p)
dev.off()

##Boxplots
meta_q3 <-  st_q3@meta.data
#meta_q3$new_groups <- factor(meta_q3$new_groups, levels = c("rest", "Low", "Surronding", "Nuclei"))


# Identify unique samples
unique_names <- unique(meta_q3$name)
colors_palette <- setNames(rainbow(length(unique_names)), unique_names)

# Correct the spelling to match the data and ensure proper factor levels
meta_q3$new_groups <- factor(meta_q3$new_groups, levels = c("rest", "Low", "Surronding", "Hot_spot"))

# Function to create combined boxplots
create_combined_boxplots <- function(df, variables, group_var, sample_var, colors) {
  plots <- list()
  for (var in variables) {
    p <- ggplot(df, aes_string(x = group_var, y = var, fill = sample_var)) +
      geom_boxplot(position = position_dodge(width = 0.8)) +
      scale_fill_manual(values = colors) +
      ggtitle(paste("Combined Boxplot of", var)) +
      theme_minimal()
    plots[[var]] <- p
  }
  return(plots)
}

# Create and save combined boxplots for meta_q3
combined_boxplots_meta_q3 <- create_combined_boxplots(meta_q3, variables_to_plot, "new_groups", "name", colors_palette)
pdf("./results/new_grouping/combined_boxplots_meta_q3.pdf", width = 11, height = 8.5)
do.call(grid.arrange, c(combined_boxplots_meta_q3, ncol = 1))
dev.off()

# Function to create combined boxplots NO SAMPLE SEPARATION
create_combined_boxplots <- function(df, variables, group_var,colors) {
  plots <- list()
  for (var in variables) {
    p <- ggplot(df, aes_string(x = group_var, y = var)) +
      geom_boxplot(position = position_dodge(width = 0.8)) +
      scale_fill_manual(values = colors) +
      ggtitle(paste("Combined Boxplot of", var)) +
      theme_minimal()
    plots[[var]] <- p
  }
  return(plots)
}

combined_boxplots_meta_q3 <- create_combined_boxplots(df_filtered, "Tcell_Exh_UCell" , "new_groups", colors_palette)
pdf("./results/distribution/combined_boxplots_meta_q3_Tcell_Exh_UCell_together.pdf", width = 7, height = 10)
do.call(grid.arrange, c(combined_boxplots_meta_q3, ncol = 1))
dev.off()

# Create and save combined boxplots for meta_q3
# Assuming your data frame is named df
df_filtered <- meta_q3 %>% 
  filter(meta_q3$new_groups != "rest")

a <- ggplot(df_filtered, aes(x = new_groups, y = Tcell_Exh_UCell, fill = new_groups)) +
  geom_boxplot() +
  theme_minimal() +
  ggtitle("Tcell Exh score by Group")
b <- ggplot(df_filtered, aes(x = new_groups, y = PC_new_UCell, fill = new_groups)) +
  geom_boxplot() +
  theme_minimal() +
  ggtitle("PC score by Group")
c <- ggplot(df_filtered, aes(x = new_groups, y = Teff_UCell, fill = new_groups)) +
  geom_boxplot() +
  theme_minimal() +
  ggtitle("Teff score by Group")

pdf("./results/distribution/boxplot_signatures.pdf", width = 15, height = 6)
grid.arrange(a, b, c, ncol = 3)
dev.off()

##CORRELACION Q3
# Función para crear gráficos de dispersión y calcular la correlación de Spearman
create_cor_plots <- function(df, x_var, y_var, group_var) {
  plots <- list()
  unique_groups <- unique(df[[group_var]])
  
  for (group in unique_groups) {
    p <- ggplot(df[df[[group_var]] == group, ], aes_string(x = x_var, y = y_var)) +
      geom_point(aes(color = name)) +
      theme_cowplot() +
      stat_cor(method = "spearman", cor.coef.name = "rho") +
      labs(title = paste(group, "-", x_var, "vs", y_var)) +
      theme(plot.title = element_text(size=10))
    plots[[group]] <- p
  }
  
  return(plots)
}

##no color by sample
create_cor_plots <- function(df, x_var, y_var, group_var) {
  plots <- list()
  unique_groups <- unique(df[[group_var]])
  
  for (group in unique_groups) {
    p <- ggplot(df[df[[group_var]] == group, ], aes_string(x = x_var, y = y_var)) +
      geom_point() +  # Removed the color aesthetic
      theme_cowplot() +
      stat_cor(method = "spearman", cor.coef.name = "rho") +
      labs(title = paste(group, "-", x_var, "vs", y_var)) +
      theme(plot.title = element_text(size=10))
    plots[[group]] <- p
  }
  
  return(plots)
}

# PC_new_UCell vs Tcell_Exh_UCell
plots_PC_TcellExh <- create_cor_plots(meta_q3, "PC_new_UCell", "Tcell_Exh_UCell", "new_groups")
combined_plots_PC_TcellExh <- plot_grid(plotlist = plots_PC_TcellExh, ncol = 4)

pdf("./results/distribution/correlation_meta_q3_pcvstex.pdf", width = 17, height = 6)
print(combined_plots_PC_TcellExh, ncol = 4)
dev.off()

# PC_new_UCell vs Teff_UCell
plots_PC_Teff <- create_cor_plots(meta_q3, "PC_new_UCell", "Teff_UCell", "new_groups")
combined_plots_PC_Teff <- plot_grid(plotlist = plots_PC_Teff, ncol = 4)

pdf("./results/distribution/correlation_meta_q3_pcvsteff.pdf", width = 17, height = 6)
print(combined_plots_PC_Teff, ncol = 4)
dev.off()

# Tcell_Exh_UCell vs Teff_UCell
plots_TcellExh_Teff <- create_cor_plots(meta_q3, "Tcell_Exh_UCell", "Teff_UCell", "new_groups")
combined_plots_TcellExh_Teff <- plot_grid(plotlist = plots_TcellExh_Teff, ncol = 1)

pdf("./results/new_grouping/correlation_meta_q3_texvsteff.pdf", width = 5, height = 10)
print(combined_plots_TcellExh_Teff, ncol = 4)
dev.off()



## tendencias
# Reestructurar los datos para incluir PC_new_UCell, Tcell_Exh_UCell y Teff_UCell
long_data <- meta_q3 %>%
  select(name, new_groups, PC_new_UCell, Tcell_Exh_UCell, Teff_UCell) %>%
  pivot_longer(cols = c(PC_new_UCell, Tcell_Exh_UCell, Teff_UCell), names_to = "variable", values_to = "value")

# Crear el gráfico
ggplot(long_data, aes(x = new_groups, y = value, color = variable, group = variable)) +
  geom_point(position = position_dodge(width = 0.2)) +
  geom_line(position = position_dodge(width = 0.2)) +
  facet_wrap(~ name) +  # Si quieres separar por muestra
  theme_minimal() +
  labs(title = "Comparación de PC_new_UCell, Tcell_Exh_UCell y Teff_UCell por new_groups",
       y = "Valor",
       color = "Variable")


# Reestructurar los datos para tener Tcell_Exh_UCell y Teff_UCell en formato largo
long_data <- meta_q3 %>%
  select(name, new_groups, Tcell_Exh_UCell, Teff_UCell, Neutrophils) %>%
  pivot_longer(cols = c(Tcell_Exh_UCell, Teff_UCell, Neutrophils), names_to = "variable", values_to = "value")

# Crear el gráfico
ggplot(long_data, aes(x = new_groups, y = value, color = variable, group = variable)) +
  geom_point(position = position_dodge(width = 0.2)) +
  geom_line(position = position_dodge(width = 0.2)) +
  facet_wrap(~ name) +  # Si quieres separar por muestra
  theme_minimal() +
  labs(title = "Comparación de Tcell_Exh_UCell y Teff_UCell por new_groups",
       y = "Valor",
       color = "Variable")

########UMAP NUCLEI

st_nuclei <- st_q3
Idents(object = st_nuclei) <- "new_groups"
st_nuclei <- SubsetSTData(st_nuclei, idents = c("Hot_spot"))

DefaultAssay(st_nuclei) <- "RNA"
st_nuclei <- NormalizeData(st_nuclei, verbose = FALSE)
st_nuclei <- FindVariableFeatures(st_nuclei, selection.method = "vst", nfeatures = 2000)
st_nuclei <- ScaleData(st_nuclei, verbose = FALSE)
st_nuclei <- RunPCA(st_nuclei, features = VariableFeatures(object = st_nuclei))
st_nuclei <- RunUMAP(st_nuclei, dims = 1:10)
DimPlot(st_nuclei, reduction = "umap",group.by ="name")

DefaultAssay(st_nuclei) <- "SCT"
st_nuclei <- SCTransform(st_nuclei, vars.to.regress = "name")
st_nuclei <- RunPCA(st_nuclei, features = VariableFeatures(object = st_nuclei))
st_nuclei <- FindNeighbors(st_nuclei, dims = 1:10)
st_nuclei <- FindClusters(st_nuclei, resolution = 0.5)
st_nuclei <- RunUMAP(st_nuclei, dims = 1:10, )
DimPlot(st_nuclei, reduction = "umap",group.by ="name")
DimPlot(st_nuclei, reduction = "umap",group.by ="seurat_clusters")

FeatureOverlay(st_nuclei, features = c("seurat_clusters"), 
               sampleids = 1:6, ncols = 2, 
               pt.size = 0.9)

FeaturePlot(st_nuclei, reduction = "umap", features = "Xbp1", cols = color,pt.size = 1.3)
FeatureOverlay(st_nuclei, features = c("Xbp1"), 
               sampleids = 1:6, ncols = 2, 
               pt.size = 1.3)

##divide by sample
saveRDS(st_q3, "./objects/sp/st_q3.rds")
st_nuclei <- st_q3
Idents(object = st_nuclei) <- "new_groups"
st_nuclei <- SubsetSTData(st_nuclei, idents = c("Hot_spot"))

Idents(object = st_nuclei) <- "name"
name <- unique(st_nuclei@meta.data[["name"]])
objects <- c()
plots_list <- c()
plot_spatial <- c()

for (i in name){
  a <- SubsetSTData(st_nuclei, idents = i)
  a <- RunPCA(a, features = VariableFeatures(object = st_nuclei))
  a <- FindNeighbors(a, dims = 1:10)
  a <- FindClusters(a, resolution = 0.3)
  a <- RunUMAP(a, dims = 1:10, )
  plots_list[[i]] <- DimPlot(a, reduction = "umap",group.by ="seurat_clusters") +
                  ggtitle(paste("nuclei sep for", i))
  plot_spatial[[i]] <- FeatureOverlay(a, features = c("seurat_clusters"),
                 pt.size = 1.3)+
                  ggtitle(paste("spatial sep for", i))
  objects[[length(objects) + 1]] <- a
}

names(objects) <- name
grid_umaps <- do.call(grid.arrange, c(plots_list, ncol = length(name)))
grid_spatial <- do.call(grid.arrange, c(plot_spatial, ncol = length(name)))

FeatureOverlay(mouse, features = c("seurat_clusters"), 
               sampleids = 1:6, ncols = 2, 
               pt.size = 0.7)

FeatureOverlay(st_q3, features = c("new_groups"), 
               sampleids = 1:6, ncols = 2, 
               pt.size = 0.7)


