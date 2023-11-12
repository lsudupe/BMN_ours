## SCRIPT: Human sample Healthy enrichment scores

## 08.11.23 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(ggplot2)
library(Signac)
library(RColorBrewer)
library(STutility)
library(UCell)
library(readr)
library(dplyr)
library(gridExtra)
library(tibble)
library(readxl)


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
  MSC = top_genes_por_hoja$MSC,
  EC = top_genes_por_hoja$EC,
  Pericytes = top_genes_por_hoja$Pericytes,
  Megakar = top_genes_por_hoja$Megakar,
  Erythrobl = top_genes_por_hoja$Erythrobl,
  Tcells = top_genes_por_hoja$Tcells,
  Dendritic = top_genes_por_hoja$Dendritic,
  Bcells = top_genes_por_hoja$Bcells
)

# Crear variables en el entorno global para cada hoja
list2env(top_genes_por_hoja, envir = .GlobalEnv)

#Variables---------------------------------
DIR_ROOT <- file.path(getwd())
DIR_DATA <- file.path(DIR_ROOT, "data/data/")

# samples
all <- dir(path = DIR_DATA)
#samples <- samples[ samples %in% c("BM_human_AP-B08805")]

samples <- list.files(path = DIR_DATA, recursive = TRUE, pattern = "filtered",full.names = TRUE)
imgs <- list.files(path = DIR_DATA, pattern = "tissue_hires_image.png", recursive = TRUE, full.names = TRUE)
spotfiles <- list.files(path = DIR_DATA, pattern = "tissue_positions_list.csv", recursive = TRUE, full.names = TRUE)
json <- list.files(path = DIR_DATA, pattern = "scalefactors_json.json", recursive = TRUE, full.names = TRUE)
name <- dir(path = DIR_DATA)

infoTable <- data.frame(samples = samples, imgs = imgs, spotfiles = spotfiles, json = json, 
                        name = name)


infoTable <- subset(infoTable, name != "M1_fem_1C" &
                      name != "M1_tib_1A" & 
                      name != "M2_F_2B" &
                      name != "M3_F_1C" &
                      name != "M3_fem_1C" &
                      name != "M3_tib_2A" & 
                      name != "M8_F2_1C" & 
                      name !="M9_F2_1C" &
                      name !="BM_B000943" &
                      name !="BM_B01320" &
                      name !="BM_B02817" &
                      name !="BM_B10395" &
                      name !="BM_human_AP-B00182_" &
                      name !="BM_human_AP-B02149_" &
                      name !="BM_human_AP-B08041_" &
                      name !="BM_human_AP-B08805" )

se <- InputFromTable(infoTable,
                     platform =  "Visium")

st.object <- GetStaffli(se)
st.object
se <- LoadImages(se, time.resolve = FALSE)
saveRDS(se, "./objects/sp/integrated/se_human_healthy.rds")
se <- readRDS("./objects/sp/integrated/se_human_healthy.rds")

##divide by sample
Idents(object = se) <- "name"
name <- unique(se@meta.data[["name"]])
objects <- c()

for (i in name){
  a <- SubsetSTData(se, idents = i)
  a <- SCTransform(a)
  ## add object to list
  objects[[length(objects) + 1]] <- a
  
}
names(objects) <- name
## separate list
list2env(objects,envir=.GlobalEnv)

se <- SCTransform(se)

objects <- c(BM_human_H13, BM_human_H6, BM_human_H7, BM_human_H9)
names(objects) <- c("BM_human_H13", "BM_human_H6", "BM_human_H7", "BM_human_H9")

###Enrichment score
post <- c()
for(i in 1:length(objects)) {
  a <- objects[[i]]
  b <- names(objects)[i]  # Note: I've changed this to ensure it gets the name of the ith object
  
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
    dir.create(file.path("./results/human/healthy", lista_name), recursive = TRUE, showWarnings = FALSE)
    
    # Guardar la visualización en PDF
    pdf_name <- paste("./results/human/healthy", lista_name, paste0(b, "_", lista_name, ".pdf"), sep = "/")
    pdf(pdf_name)
    print(p)
    dev.off()
  }
  
  # Añadir el objeto Seurat modificado a la lista post
  post[[length(post) + 1]] <- a
}

names(post) <- c("BM_human_H13", "BM_human_H6", "BM_human_H7", "BM_human_H9")

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
                         breaks = c(0.0,0.4),
                         labels = c("Min", "Max"),
                         limits = c(0.0,0.4))
  
  # Save the plot to a PDF filef
  pdf(paste("./results/human/healthy/samerange/", b,"_dentritic.pdf",sep=""))
  print(p)
  dev.off()
  
  #EC
  p <- FeatureOverlay(a, features = c("signature_EC"), ncols = 1, pt.size = 1.1, 
                      value.scale = "all" ,cols = color) +
    scale_fill_gradientn(colours = color,
                         breaks = c(0.0,0.3),
                         labels = c("Min", "Max"),
                         limits = c(0.0,0.3))
  
  # Save the plot to a PDF filef
  pdf(paste("./results/human/healthy/samerange/", b,"_EC.pdf",sep=""))
  print(p)
  dev.off()
  
  #Erythrobl
  p <- FeatureOverlay(a, features = c("signature_Erythrobl"), ncols = 1, pt.size = 1.1, 
                      value.scale = "all" ,cols = color) +
    scale_fill_gradientn(colours = color,
                         breaks = c(0.0,0.4),
                         labels = c("Min", "Max"),
                         limits = c(0.0,0.4))
  
  # Save the plot to a PDF filef
  pdf(paste("./results/human/healthy/samerange/", b,"_Erythrobl.pdf",sep=""))
  print(p)
  dev.off()
  
  #Megakar
  p <- FeatureOverlay(a, features = c("signature_Megakar"), ncols = 1, pt.size = 1.1, 
                      value.scale = "all" ,cols = color) +
    scale_fill_gradientn(colours = color,
                         breaks = c(0.0,0.2),
                         labels = c("Min", "Max"),
                         limits = c(0.0,0.2))
  
  # Save the plot to a PDF filef
  pdf(paste("./results/human/healthy/samerange/", b,"_Megakar.pdf",sep=""))
  print(p)
  dev.off()
  
  #MSC
  p <- FeatureOverlay(a, features = c("signature_MSC"), ncols = 1, pt.size = 1.1, 
                      value.scale = "all" ,cols = color) +
    scale_fill_gradientn(colours = color,
                         breaks = c(0.0,0.2),
                         labels = c("Min", "Max"),
                         limits = c(0.0,0.2))
  
  # Save the plot to a PDF filef
  pdf(paste("./results/human/healthy/samerange/", b,"_MSC.pdf",sep=""))
  print(p)
  dev.off()
  
  #Neutrop
  p <- FeatureOverlay(a, features = c("signature_Neutrop"), ncols = 1, pt.size = 1.1, 
                      value.scale = "all" ,cols = color) +
    scale_fill_gradientn(colours = color,
                         breaks = c(0.0,0.25),
                         labels = c("Min", "Max"),
                         limits = c(0.0,0.25))
  
  # Save the plot to a PDF filef
  pdf(paste("./results/human/healthy/samerange/", b,"_Neutrop.pdf",sep=""))
  print(p)
  dev.off()
  
  #Pericytes
  p <- FeatureOverlay(a, features = c("signature_Pericytes"), ncols = 1, pt.size = 1.1, 
                      value.scale = "all" ,cols = color) +
    scale_fill_gradientn(colours = color,
                         breaks = c(0.0,0.35),
                         labels = c("Min", "Max"),
                         limits = c(0.0,0.35))
  
  # Save the plot to a PDF filef
  pdf(paste("./results/human/healthy/samerange/", b,"_Pericytes.pdf",sep=""))
  print(p)
  dev.off()
  
  #Tcells
  p <- FeatureOverlay(a, features = c("signature_Tcells"), ncols = 1, pt.size = 1.1, 
                      value.scale = "all" ,cols = color) +
    scale_fill_gradientn(colours = color,
                         breaks = c(0.0,0.2),
                         labels = c("Min", "Max"),
                         limits = c(0.0,0.2))
  
  # Save the plot to a PDF filef
  pdf(paste("./results/human/healthy/samerange/", b,"_Tcells.pdf",sep=""))
  print(p)
  dev.off()

}


###plot specific genes
# Read the Excel files
ec_data <- read_excel("./data/EC_mergedNORIBO_markers_ages.xlsx")
msc_data <- read_excel("./data/MSC_NEW_Age_05FCmarkers.xlsx")

# Define a function to filter and get top 5 genes
get_top_genes <- function(data) {
  data %>%
    filter(p_val_adj < 0.005, avg_log2FC > 0) %>%
    arrange(desc(avg_log2FC)) %>%
    slice_head(n = 10) %>%
    select(gene)
}

# Apply the function to each dataset
top_genes_ec <- get_top_genes(ec_data)
top_genes_msc <- get_top_genes(msc_data)

# Print the results
print("Top 5 genes for EC:")
print(top_genes_ec)

print("Top 5 genes for MSC:")
print(top_genes_msc)

for(i in 1:length(post)) {
  a <- post[[i]]
  b <- names(post)[i]
  
  # MSC
  for(gene in c("PLCG2", "SOX4", "PPP1R10", "RETREG1", "FAM43A")) {
    p <- FeatureOverlay(a, features = c(gene), ncols = 1, pt.size = 1.5, 
                        value.scale = "all", cols = color)
    
    # Guardar el plot
    pdf_name <- paste0("./human/healthy/genes/", b, "_", gene, "_MSC.pdf")
    pdf(pdf_name)
    print(p)
    dev.off()
  }
  
  # EC
  for(gene in c("PLCG2","MAFB", "SLC2A3")) {
    p <- FeatureOverlay(a, features = c(gene), ncols = 1, pt.size = 1.5, 
                        value.scale = "all", cols = color)
    
    # Guardar el plot
    pdf_name <- paste0("./human/healthy/genes/", b, "_", gene, "_EC.pdf")
    pdf(pdf_name)
    print(p)
    dev.off()
  }
}




