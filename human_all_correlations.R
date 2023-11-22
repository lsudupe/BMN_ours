## SCRIPT: Human correlation analysis

## 21.11.23 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(ggplot2)
library(UCell)
library(STutility)
library(RColorBrewer)
library(cowplot)
library(ggpubr)
color <- rev(brewer.pal(11,"Spectral"))

## Data
healthy <- readRDS("./objects/sp/integrated/healthy_clean.rds")
mm <- readRDS("./objects/sp/human/human_combined.rds")

all <- c(healthy, mm)
## Signatures
source("./signatures.R")
sig <- list(new_PC_markers_human, Teff_Juanro_human, TcellAct_Juanro_human, genes38_human, cytotoxicity_human)
names(sig) <- c("new_PC_markers_human", "Teff_Juanro_human", "TcellAct_Juanro_human", "genes38_human", "cytotoxicity_human")

## Enrichment score
post <- c()
for (i in 1:length(all)) {
  a <- all[[i]]
  b <- names(all)[i]
  plots_list <- list()
  
  for (signature_name in 1:length(sig)) {
    d <- names(sig[signature_name])
    genes <- sig[[signature_name]]
    
    # Añade el UCellScore para la firma actual
    vector <- ScoreSignatures_UCell(a@assays[["RNA"]]@counts, features = list(genes))
    a@meta.data[[d]] <- as.vector(vector)
    
    # Aplicar subconjuntos si es necesario
    # ... (tu código para aplicar subconjuntos) ...
    
    # Crea el gráfico
    lim <- c(min(a@meta.data[[d]]), max(a@meta.data[[d]]))
    p <- FeatureOverlay(a, features = c(d), ncols = 1, pt.size = 1.5, value.scale = "all", cols = color)
    if (lim[1] != lim[2]) {
      p <- p + scale_fill_gradientn(colours = color, breaks = lim, labels = lim, limits = lim)
    }
    
    plots_list <- c(plots_list, list(p))
  }
  
  # Combinar todos los gráficos en un único PDF para este objeto Seurat
  file_name_all_plots <- paste0("./results/human/correlation/spatial/", b, "_all_plots.pdf")
  if (length(plots_list) > 0) {
    combined_plot <- ggarrange(plotlist = plots_list, ncol = 1, nrow = length(plots_list))
    ggsave(file_name_all_plots, combined_plot, width = 11, height = 8.5 * length(plots_list), limitsize = FALSE)
  } else {
    message("No plots were created for ", b)
  }
  
  post[[length(post) + 1]] <- a
}

names(post) <- c( "BM_human_H13","BM_human_H6","BM_human_H7", "BM_human_H9","BM_human_AP-B00182_",
"BM_human_AP-B02149_" ,"BM_human_AP-B08041_", "BM_human_AP-B08805" , "BM_B000943" , "BM_B01320",      
"BM_B02817" ,"BM_B10395")

# Correlations
for (i in 1:length(post)) {
  seurat <- post[[i]]
  b <- names(post)[i]
  plots_list <- list()
  
  combinations <- list(
    list("new_PC_markers_human", "Teff_Juanro_human", "PC vs Teffectoras"),
    list("new_PC_markers_human", "TcellAct_Juanro_human", "PC vs Tactivadas"),
    list("new_PC_markers_human", "genes38_human", "PC vs Texhausted"),
    list("new_PC_markers_human", "cytotoxicity_human", "PC vs cytotoxicidad"),
    list("genes38_human", "Teff_Juanro_human", "Texhausted vs Teffectoras"),
    list("genes38_human", "TcellAct_Juanro_human", "Texhausted vs Tactivadas"),
    list("genes38_human", "cytotoxicity_human", "Texhausted vs cytotoxicidad")
  )
  
  for (comb in combinations) {
    signature1 <- comb[[1]]
    signature2 <- comb[[2]]
    plot_title <- comb[[3]]
    
    if (signature1 %in% names(seurat@meta.data) && signature2 %in% names(seurat@meta.data)) {
      p <- ggplot(seurat@meta.data, aes(!!sym(signature1), !!sym(signature2))) +
        geom_point() +
        theme_cowplot() +
        stat_cor(method = "spearman", cor.coef.name = "r") +
        ggtitle(plot_title)
      
      plots_list <- c(plots_list, list(p))
    } else {
      message("Skipping combination ", signature1, " vs ", signature2, " for ", b)
    }
  }
  num_plots <- length(plots_list)
  height_per_plot <- 8.5  # Altura en pulgadas para cada gráfico en orientación horizontal
  width_per_plot <- 11    # Ancho en pulgadas para cada gráfico en orientación horizontal
  max_height <- 36        # Altura máxima en pulgadas para el PDF en orientación horizontal
  
  # Calcular la altura total manteniendo dentro del límite máximo
  total_height <- min(num_plots * height_per_plot, max_height)
  
  # Combinar todos los gráficos y guardarlos como PDF en orientación horizontal
  if (num_plots > 0) {
    combined_plot <- ggarrange(plotlist = plots_list, ncol = 1, nrow = num_plots)
    file_name_all_plots <- paste0("./results/human/correlation/plot/", b, "_all_plots.pdf")
    ggsave(file_name_all_plots, combined_plot, width = width_per_plot, height = total_height, limitsize = FALSE)
  } else {
    message("No plots were created for ", b)
  }
  
}


## Check which genes are in each object
resultados <- list()

for (nombre_objeto in names(all)) {
  genes_counts <- all[[nombre_objeto]]@assays[["RNA"]]@counts@Dimnames[[1]]
  
  for (nombre_firma in names(sig)) {
    genes_interes <- sig[[nombre_firma]]
    genes_ausentes <- genes_interes[!genes_interes %in% genes_counts]
    resultados[[nombre_objeto]][[nombre_firma]] <- genes_ausentes
  }
}

# Convertir los resultados a un formato de tabla
tabla_resultados <- do.call(rbind, lapply(resultados, function(x) {
  sapply(x, function(y) paste(y, collapse = ", "))
}))

# Poner nombres a las filas y columnas
rownames(tabla_resultados) <- names(resultados)
colnames(tabla_resultados) <- names(sig)

# Mostrar la tabla
print(tabla_resultados)



