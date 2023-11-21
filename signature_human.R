## SCRIPT: Extract signature per cluster and check in human BM project

## 03.10.23 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(ggplot2)
library(STutility)
library(RColorBrewer)
library(tidyverse)
library(UCell)
library(readxl)
library(dplyr)

#Data---------------------------------
se <- readRDS("./objects/heterogeneity/se_hierarchical.rds")
human <- readRDS("./objects/sp/human/human_combined.rds")

# split human
B08041 <- human[["BM_human_AP-B08041_"]]
B08805 <- human[["BM_human_AP-B08805"]]
B10395 <- human[["BM_B10395"]]

Idents(se) <- "clustering"
all_markers <- FindAllMarkers(se, min.pct = 0.25)

top_markers_per_cluster <- all_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 50, wt = -log10(p_val_adj))

split_data <- split(top_markers_per_cluster, top_markers_per_cluster$cluster)

# Assume you have gene lists for each human sample
genes_B08041 <- rownames(B08041)  # Extract gene names
genes_B08805 <- rownames(B08805)
genes_B10395 <- rownames(B10395)

# Find common genes across all three samples
common_genes <- Reduce(intersect, list(genes_B08041, genes_B08805, genes_B10395))

# Function to get top genes ensuring they are in the common gene list
get_top_genes <- function(cluster_data, common_genes, n = 25) {
  # Convert gene names in the cluster data to uppercase
  cluster_genes_upper <- toupper(cluster_data$gene)
  # Ensure that common_genes are also in uppercase
  common_genes_upper <- toupper(common_genes)
  # Check for common genes in uppercase
  common_cluster_genes <- cluster_genes_upper[cluster_genes_upper %in% common_genes_upper]
  if(length(common_cluster_genes) >= n) {
    return(head(common_cluster_genes, n))
  } else {
    # Select additional genes if fewer than n are common
    additional_genes <- setdiff(cluster_genes_upper, common_cluster_genes)
    total_genes <- c(common_cluster_genes, head(additional_genes, n - length(common_cluster_genes)))
    return(total_genes)
  }
}

# Your existing loop
for(i in 1:length(split_data)) {
  cluster_name <- unique(split_data[[i]]$cluster)[1]
  top_genes <- get_top_genes(split_data[[i]], common_genes)
  assign(paste0("cluster", cluster_name), top_genes)
}

lista <- c(B08041, B08805, B10395)
names(lista) <- c("B08041","B08805", "B10395")

###Enrichment score
color <- brewer.pal(11,"Spectral")
color <- rev(color)

post <- c()
# Define the cluster names
cluster_names <- paste0("cluster", 1:7)

for (i in 1:length(lista)){
  v <- "sss"
  a <- lista[[i]]
  b <- names(lista[i])
  
  for (cluster_name in cluster_names) {
    # Define the variable names dynamically based on the cluster_name
    signature_var <- paste0("signature_1_", cluster_name)
    
    file_name <- paste0("./results/human/cluster_enrich/new/25/", b, "_", cluster_name, "_pc.pdf", sep = "")
    
    # Add UCellScore for the current cluster
    vector <- ScoreSignatures_UCell(a@assays[["RNA"]]@counts, features = list(get(cluster_name)))
    a@meta.data[[signature_var]] <- as.vector(vector)
    
    #c <- c(min(a@meta.data[[signature_var]]), max(a@meta.data[[signature_var]]))
    # Create the plot
    p <- FeatureOverlay(a, features = c(signature_var), ncols = 1, pt.size = 1.1, 
                        value.scale = "all" ,cols = color) +
      scale_fill_gradientn(colours = color,
                           breaks = c(0.0, 0.5),
                           labels = c("Min", "Max"),
                           limits = c(0.0, 0.5))
    
    # Save the plot to a PDF file
    pdf(file_name)
    print(p)
    dev.off()
  }
  
  #########
  # Call the function for each comparison
  plot_ratio(a, 1, 7, color, "./results/human/cluster_enrich/25_ratio/")
  plot_ratio(a, 5, 7, color, "./results/human/cluster_enrich/25_ratio/")
  plot_ratio(a, 5, 6, color, "./results/human/cluster_enrich/25_ratio/")
  plot_ratio(a, 6, 7, color, "./results/human/cluster_enrich/25_ratio/")
  
  post[[length(post) + 1]] <- a
}

names(post) <- c("B08041","B08805", "B10395")

## Check Isa signatures
# Leer los nombres de las hojas del archivo Excel
nombres_hojas <- excel_sheets("./data/itziar_genes.xlsx")

# Leer los datos de cada hoja y seleccionar los top 100 genes por avg_log2FC
top_genes_por_hoja <- lapply(nombres_hojas, function(nombre_hoja) {
  dataset <- read_excel("./data/itziar_genes.xlsx", sheet = nombre_hoja)
  top_genes <- dataset %>%
    arrange(desc(abs(avg_log2FC))) %>%
    slice_head(n = 50) %>%
    select(gene) %>% # Seleccionamos solo la columna de genes
    pull() # Convertimos la columna de genes en un vector
})

# Asignar los nombres de las hojas a los vectores de genes
names(top_genes_por_hoja) <- nombres_hojas
listas_marcadores <- list(
  Neutrop = top_genes_por_hoja$Neutrop,
  Erythrobl = top_genes_por_hoja$Erythrobl,
  Tcells = top_genes_por_hoja$Tcells,
  Dendritic = top_genes_por_hoja$Dendritic,
  Bcells = top_genes_por_hoja$Bcells
)

#list2env(objects,envir=.GlobalEnv)
###Enrichment score
post <- c()
for(i in 1:length(human)) {
  a <- human[[i]]
  b <- names(human)[i]  # Note: I've changed this to ensure it gets the name of the ith object
  
  # Iterar sobre cada lista de marcadores
  for(lista_name in names(listas_marcadores)) {
    marcadores <- listas_marcadores[[lista_name]] # Directly use the list
    
    # Calcular UCell Score para la lista actual de marcadores
    vector <- ScoreSignatures_UCell(a@assays[["RNA"]]@counts, features = list(marcadores))
    signature_name <- paste0("signature_", lista_name)
    a@meta.data[[signature_name]] <- as.vector(vector)
    
    # Crear visualización
    color <- rev(brewer.pal(11,"Spectral"))
    
    p <- FeatureOverlay(a, features = c(signature_name), ncols = 1, pt.size = 1.5, 
                        value.scale = "all", cols = color)
    
    # Crear directorios si no existen
    dir.create(file.path("./results/human/mm_celltype/", lista_name), recursive = TRUE, showWarnings = FALSE)
    
    # Guardar la visualización en PDF
    pdf_name <- paste("./results/human/mm_celltype/", lista_name, paste0(b, "_", lista_name, ".pdf"), sep = "/")
    pdf(pdf_name)
    print(p)
    dev.off()
  }
  
  # Añadir el objeto Seurat modificado a la lista post
  post[[length(post) + 1]] <- a
}

names(post) <- c("BM_human_AP-B00182_", "BM_human_AP-B02149_", "BM_human_AP-B08041_", "BM_human_AP-B08805",
                  "BM_B000943", "BM_B01320", "BM_B02817", "BM_B10395")

###same range
###Enrichment score
for(i in 1:length(post)) {
  a <- post[[i]]
  b <- names(post)[i] 
  
  #bcell
  p <- FeatureOverlay(a, features = c("signature_Bcells"), ncols = 1, pt.size = 1.1, 
                      value.scale = "all" ,cols = color) +
    scale_fill_gradientn(colours = color,
                         breaks = c(0.0,0.2),
                         labels = c("Min", "Max"),
                         limits = c(0.0,0.2))
  
  # Save the plot to a PDF filef
  pdf(paste("./results/human/healthy/samerange/", b,"_bcell.pdf",sep=""))
  print(p)
  dev.off()
  
  #dentritic
  p <- FeatureOverlay(a, features = c("signature_Dendritic"), ncols = 1, pt.size = 1.1, 
                      value.scale = "all" ,cols = color) +
    scale_fill_gradientn(colours = color,
                         breaks = c(0.0,0.2),
                         labels = c("Min", "Max"),
                         limits = c(0.0,0.2))
  
  # Save the plot to a PDF filef
  pdf(paste("./results/human/healthy/samerange/", b,"_dentritic.pdf",sep=""))
  print(p)
  dev.off()
  
  #Erythrobl
  p <- FeatureOverlay(a, features = c("signature_Erythrobl"), ncols = 1, pt.size = 1.1, 
                      value.scale = "all" ,cols = color) +
    scale_fill_gradientn(colours = color,
                         breaks = c(0.0,0.5),
                         labels = c("Min", "Max"),
                         limits = c(0.0,0.5))
  
  # Save the plot to a PDF filef
  pdf(paste("./results/human/healthy/samerange/", b,"_Erythrobl.pdf",sep=""))
  print(p)
  dev.off()
  
  #Neutrop
  p <- FeatureOverlay(a, features = c("signature_Neutrop"), ncols = 1, pt.size = 1.1, 
                      value.scale = "all" ,cols = color) +
    scale_fill_gradientn(colours = color,
                         breaks = c(0.0,0.30),
                         labels = c("Min", "Max"),
                         limits = c(0.0,0.30))
  
  # Save the plot to a PDF filef
  pdf(paste("./results/human/healthy/samerange/", b,"_Neutrop.pdf",sep=""))
  print(p)
  dev.off()
  
  #Tcells
  p <- FeatureOverlay(a, features = c("signature_Tcells"), ncols = 1, pt.size = 1.1, 
                      value.scale = "all" ,cols = color) +
    scale_fill_gradientn(colours = color,
                         breaks = c(0.0,0.31),
                         labels = c("Min", "Max"),
                         limits = c(0.0,0.31))
  
  # Save the plot to a PDF filef
  pdf(paste("./results/human/healthy/samerange/", b,"_Tcells.pdf",sep=""))
  print(p)
  dev.off()
  
}
