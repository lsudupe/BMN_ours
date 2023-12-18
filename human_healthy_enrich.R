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
color <- rev(brewer.pal(11,"Spectral"))

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

## Read list of objects after cleanning
objects <- readRDS("./objects/sp/integrated/healthy_clean.rds")

## separate list
#list2env(objects,envir=.GlobalEnv)
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
                         breaks = c(0.0,0.2),
                         labels = c("Min", "Max"),
                         limits = c(0.0,0.2))
  
  # Save the plot to a PDF filef
  pdf(paste("./results/human/healthy/samerange/", b,"_dentritic.pdf",sep=""))
  print(p)
  dev.off()
  
  #EC
  p <- FeatureOverlay(a, features = c("signature_EC"), ncols = 1, pt.size = 1.1, 
                      value.scale = "all" ,cols = color) +
    scale_fill_gradientn(colours = color,
                         breaks = c(0.0,0.25),
                         labels = c("Min", "Max"),
                         limits = c(0.0,0.25))
  
  # Save the plot to a PDF filef
  pdf(paste("./results/human/healthy/samerange/", b,"_EC.pdf",sep=""))
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
                         breaks = c(0.0,0.45),
                         labels = c("Min", "Max"),
                         limits = c(0.0,0.45))
  
  # Save the plot to a PDF filef
  pdf(paste("./results/human/healthy/samerange/", b,"_MSC.pdf",sep=""))
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
                         breaks = c(0.0,0.31),
                         labels = c("Min", "Max"),
                         limits = c(0.0,0.31))
  
  # Save the plot to a PDF filef
  pdf(paste("./results/human/healthy/samerange/", b,"_Tcells.pdf",sep=""))
  print(p)
  dev.off()

}


###plot specific genes
# Read the Excel files
ec_data <- read_excel("./data/EC_mergedNORIBO_markers_ages.xlsx", sheet = 1)
msc_data <- read_excel("./data/MSC_NEW_Age_05FCmarkers.xlsx", sheet = 1)

# Define a function to filter and get top 5 genes
get_top_genes <- function(data) {
  data %>%
    #filter(p_val_adj < 0.005, avg_log2FC > 0) %>%
   # arrange(desc(avg_log2FC)) %>%
   # slice_head(n = 10) %>%
    select(gene)
}

# Apply the function to each dataset
top_genes_ec <- get_top_genes(ec_data)
top_genes_msc <- get_top_genes(msc_data)

ec_genes <- pull(top_genes_ec, gene) 
msc_genes <- pull(top_genes_msc, gene)
# Print the results
print("Top 5 genes for EC:")
print(top_genes_ec)

print("Top 5 genes for MSC:")
print(top_genes_msc)

for(i in 1:length(objects)) {
  a <- objects[[i]]
  b <- names(objects)[i]
  
  # MSC
  for(gene in msc_genes) {
    tryCatch({
      p <- FeatureOverlay(a, features = c(gene), ncols = 1, pt.size = 1.5, 
                          value.scale = "all", cols = color)
      
      # Save the plot
      pdf_name <- paste0("./results/human/healthy/genes/young/", b, "_", gene, "_MSC.pdf")
      pdf(pdf_name)
      print(p)
      dev.off()
    }, error = function(e) {
      message(paste("Skipping:", gene, "in MSC - not found in dataset."))
    })
  }
  
  # EC
  for(gene in ec_genes) {
    tryCatch({
      p <- FeatureOverlay(a, features = c(gene), ncols = 1, pt.size = 1.5, 
                          value.scale = "all", cols = color)
      
      # Save the plot
      pdf_name <- paste0("./results/human/healthy/genes/young/", b, "_", gene, "_EC.pdf")
      pdf(pdf_name)
      print(p)
      dev.off()
    }, error = function(e) {
      message(paste("Skipping:", gene, "in EC - not found in dataset."))
    })
  }
}


#####CANONICAL CLEAN 

ec <- c("CD34", "CDH5", "EMCN", "KDR", "PECAM1", "IL6ST", "VWF", 
        "TEK", "CD9", "LY6E", "VCAM1", "CAVIN2", "ICAM2", "ADGRF5", 
        "ADGRL4", "LY6C1", "TIE2", "ARHGAP31", "CLDN5", "COTL1", "CTNNBIP1",
        "CXCL12", "EFNA1", "ESAM", "ETS2", "FAM167B", "FAM212A", "FCER1G", "GPIHBP1",
        "HSPA1A", "IFGBP7", "ITGB1", "KLF2", "KLF6", "LALBA", "LIMS2", "PDK1", "PLIN2",
        "PLVAP", "PREX2", "RPL13A", "SPARCL1", "TCN2", "THRSP", "TINAGL1", "TM4SF1", "TRP53I11", "TXNIP", "VIM", "YES1", "ICAM1")
msc <- c("LEPR", "CXCL12", "VCAM1", "ANGPT1", "KITLG", "ENG", "LPL", "CP", "COL1A1", 
         "PDGFRB", "HP", "PTH1R", "NT5E", "THY1", "NES", "PDGFRA", "GREM1")
oln <- c("COL1A1", "SPP1", "CLEC11A", "COL3A1", "COL5A2", "COL5A1", "COL12A1", "COL16A1", 
         "COL6A3", "COL11A1", "MMP13", "BGLAP", "RUNX2", "SP7")
adipocytes <- c("SOX9", "ACAN", "CILP")
pericytes <- c("RGS5", "ACTA2", "TAGLN", "MYL9", "TPM2", "NOTCH3", "HIGD1B", "GM13889", "SNCG",
"MUSTN1", "STEAP4", "RASL11A", "MYH11", "RRAD", "TINAGL1", "PTP4A3", "CYGB", "FILIP1L", "CALD1", 
"MYLK", "IL6", "NRIP2", "FST", "PCP4L1", "GUCY1A3", "H2-M9", "CSRP1", "THY1", "DES", "SPARCL1",
"PDGFRB", "OLFR558", "ATP1B2", "ASPN", "OLFML2A", "RGS4", "TPM1", "AOC3", "COL18A1", "CXCL1", 
"PPP1R14A", "EBF1", "COL5A3", "MND1", "NR2F2", "RARRES2", "COL4A1", "PPP1R12A", "GJA4", "LMOD1", 
"BCAM", "MEG3", "CSRP2", "LGALS1", "GUCY1B3", "ZAK", "FLNA", "ADAMTS4", "PDGFA", "PLN", "DTX3", 
"VSTM4", "ART3", "KCNK3", "TUSC5", "COL4A2", "HSPB1", "ELN", "NGF", "RGS16", "CYSTM1", "NRARP", 
"GJC1", "MYO1B", "RASL12", "VTN", "MYL6", "CEBPB", "MAP1LC3A", "EPAS1", "ITGA7", "43712", "KCNJ8",
"LBH", "PRSS23", "PHLDA1", "LHFP", "PDE3A", "PDLIM3", "RCAN2", "C1S1", "TBX2", "NR4A1", "NDUFA4L2",
"KCNE4", "TBX3OS1", "MIR143HG", "EDNRA", "RBPMS2", "BCAS3", "MGST3", "NAB1", "ABCC9", "UBA2", "HEYL",
"JAG1", "TSPAN12", "GADD45B", "MPRIP", "HSPB2", "AGTR1A", "MCAM", "RASGRP2", "ADAP2OS", 
"1500009L16RIK", "GM13861", "GM37800", "FAM162A", "PGF", "ENPEP", "IFITM1", "ANO1", "PDE1A",
"GSTM1", "CRIP1", "ID3", "S1PR3", "COX4I2", "ATF3", "DNAJB4", "NTRK2", "CRYAB", "MEF2C", 
"TESC", "GCNT2", "ADAMTS1", "RBPMS", "ANXA1", "TNFRSF21", "NDRG2", "LIN7A", "ZFHX3", "IRF1",
"OAZ2", "LPP", "CASQ2", "MRVI1", "CCK", "EMID1", "RGS7BP", "ADAP2", "ADCY6", "ATP2A2", 
"VASP", "RHOJ", "DLC1", "MYL1", "ARHGEF17", "NUDT4", "SERPINI1", "RASD1", "PRRX1", "GPC6", 
"ACTN1", "PDE5A", "PARM1", "TOB2", "KLF9", "PTEN", "SLC38A11", "CPE", "NPY1R", "PPP1R15A",
"GNG11", "VCL", "ACTB", "TGFB1I1", "TSPAN15", "PTRF", "RAMP1", "GNB4", "USP2", "ID4", "FRY",
"PPP1CB", "LIMD1", "MEOX2", "FILIP1", "WTIP", "PPP1R12B", "ARHGDIB", "SEPW1", "S100A11",
"TMEM51", "HSPA1A", "BTG2", "PLCE1", "RTN4RL1", "MAT2A", "FUS", "HRCT1", "FOXS1", "IGFBP7", "ZEB2", "TNS1", "REM1", "TM4SF1", "PLEKHG2", "43715", "EBF2", "ROCK1", "ISYNA1", "MAPRE2", "KLHDC8B", "ITGB1", "JUNB", "PRKAR1A", "COX8A", "PDE4D", "LIMS1", "SPRY4", "PCDH19", "GALNT16", "CRTC3", "LPL", "GEM", "EPHX3", "LAMB2", "SYNE2", "UTRN", "COG7", "PIP5K1B", "POSTN", "ADAMTS15", "ATPIF1", "MSRB1", "PIK3R1", "ECH1", "STAC", "CPNE2", "MALAT1", "ARHGEF25", "CAV1", "IFT43", "MAST4", "SEMA5A", "ACTN4", "SORBS2", "PHLDA3", "NR2F1", "EPB41L1", "SNRK", "TRPV2", "ATP1A2", "FXYD1", "PALLD", "PDLIM1", "PHLDB2", "PLA2R1", "NBEAL1", "KCNAB1", "FERMT2", "ZC2HC1A")




canonical <- list(ec, msc, oln,adipocytes, pericytes)
names(canonical) <- c("ec", "msc", "oln", "adipocytes", "pericytes")

###Enrichment score
post <- c()
for(i in 1:length(objects)) {
  a <- objects[[i]]
  b <- names(objects)[i]  # Note: I've changed this to ensure it gets the name of the ith object
  
  # Iterar sobre cada lista de marcadores
  for(lista_name in names(canonical)) {
    marcadores <- canonical[[lista_name]] # Directly use the list
    
    # Calcular UCell Score para la lista actual de marcadores
    vector <- ScoreSignatures_UCell(a@assays[["RNA"]]@counts, features = list(marcadores))
    signature_name <- paste0("signature_", lista_name)
    a@meta.data[[signature_name]] <- as.vector(vector)
    
    # Crear visualización
    color <- rev(brewer.pal(11,"Spectral"))
    
    p <- FeatureOverlay(a, features = c(signature_name), ncols = 1, pt.size = 1.5, 
                        value.scale = "all", cols = color)
    
    # Crear directorios si no existen
    dir.create(file.path("./results/human/healthy/clean_canonical", lista_name), recursive = TRUE, showWarnings = FALSE)
    
    # Guardar la visualización en PDF
    pdf_name <- paste("./results/human/healthy/clean_canonical", lista_name, paste0(b, "_", lista_name, ".pdf"), sep = "/")
    pdf(pdf_name)
    print(p)
    dev.off()
  }
  
  # Añadir el objeto Seurat modificado a la lista post
  post[[length(post) + 1]] <- a
}

names(post) <- c("BM_human_H13", "BM_human_H6", "BM_human_H7", "BM_human_H9")

seurat_objects <- post

# Features to plot
features <- c("signature_ec", "signature_msc", "signature_oln", "signature_adipocytes", "signature_pericytes")

# Prepare data frame for plotting
plot_data <- data.frame()
for (name in names(seurat_objects)) {
  obj <- seurat_objects[[name]]
  for (feature in features) {
    temp_df <- data.frame(sample = name, 
                          feature = feature, 
                          value = obj@meta.data[[feature]])
    plot_data <- rbind(plot_data, temp_df)
  }
}

# Create violin plots for each feature
plots <- lapply(features, function(feat) {
  ggplot(plot_data[plot_data$feature == feat, ], aes(x = sample, y = value, fill = sample)) +
    geom_violin(trim = FALSE) +
    scale_fill_brewer(palette = "Set1") +
    scale_y_continuous(limits = c(0, NA)) + # Establece el límite inferior de y en 0
    theme_minimal() +
    ggtitle(paste("Gráfico de Violín de", feat))
})

# Arrange and save plots to PDF
pdf("./results/human/healthy/violin_plots.pdf", width = 12, height = 8)
do.call(grid.arrange, c(plots, ncol = 2))
dev.off()

