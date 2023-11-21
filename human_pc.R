## SCRIPT: Human sample pc utility evaluation

## 24.05.23 Laura Sudupe , git @lsudupe

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

#pc_markers
pc_mm_markers <- c("SHH", "DHH", "IHH", "PTCH1", "PTCH2", "SMO", "SUFU", "GLI1", 
                   "GLI2", "GLI3", "CD19", "CD44", "CXCR4", "KLF4", "CD28", "CD33", "CD27")

human_tcell_exhausted <- c("PDCD1", "CTLA4", "TNFRSF9", "HAVCR2", "TOX", "TIGIT", "WARS", "RSAD2",
                           "MCM7", "MX1", "NDFIP2", "ENOSF1", "CCDC141", "STMN1", "TTN", "FASLG",
                           "MCM5", "NAB1", "PHLDA1", "MCM3", "PCNA", "GAPDH", "OASL", "IFI44L",
                           "TBC1D4", "SLC43A3", "PAM", "CCL3", "ACP5", "OAS3", "CD38", "TNFSF10",
                           "GBP2", "KIF20B", "CTSB")

neutro <- c("GNB2L1", "ATP5L", "C14ORF2", "ATP5E", "ATP5G2", "TMSB4X", "ATP5I", "MPO", "UQCR11.1",
            "RPS4Y1", "C19ORF43", "TCEB2", "GLTSCR2", "GPX1", "ATP5D", "SOX4", "ATP5G3", "ATP5J",
            "ATP5O", "HLA-DRA", "C6ORF48", "ATP5J2", "USMG5", "PRTN3", "SHFM1", "ATP5B", "DEK",
            "ATP5A1", "MATR3", "TMEM66", "C11ORF31", "LAPTM5", "AZU1", "RPSAP58", "SF3B14",
            "HLA-DPB1", "ELANE", "HMGB2", "HLA-DRB1", "ALDOA", "NME1-NME2", "AIF1", "SEPT7",
            "RP11-620J15.3", "C14ORF166", "ATP5G1", "NHP2L1", "CD99", "ATPIF1", "LSMD1")


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
                      name !="M9_F2_1C")

se <- InputFromTable(infoTable,
                     platform =  "Visium")

st.object <- GetStaffli(se)
st.object
se <- LoadImages(se, time.resolve = FALSE)
saveRDS(se, "./objects/sp/integrated/se_human.rds")
se <- readRDS("./objects/sp/integrated/se_human.rds")

##divide by sample
Idents(object = se) <- "name"
name <- unique(se@meta.data[["name"]])
objects <- c()

for (i in name){
  a <- SubsetSTData(se, idents = i)
  a <- SCTransform(a)
  a <- NormalizeData(a)
  a <- ScaleData(a)
  ## add object to list
  objects[[length(objects) + 1]] <- a
  
}
names(objects) <- name
## separate list
list2env(objects,envir=.GlobalEnv)

#####change some of the objects
df2 <- `BM_human_AP-B08041_`@meta.data
df <-  read_csv("./data/data/BM_human_AP-B08041_/spatial/on_tissue.csv")
# Replace NA in the 'on_tissue' column with 'off'
df$on_tissue[is.na(df$on_tissue)] <- "off"
# Convert rownames to a column 
df <- df %>%
  rownames_to_column(var = "index")
df$index <- df$Barcode
# Add "_7" to each index in df
df$index <- paste(df$index, "_7", sep = "")
# Suppose df2 is your second data frame
# Convert rownames to a column in df2
df2 <- df2 %>%
  rownames_to_column(var = "index")
# Add 'on_tissue' column to df2 based on index match and in the order of df2
df2 <- df2 %>%
  left_join(df[, c("index", "on_tissue")], by = "index") %>%
  filter(!is.na(on_tissue))
# In case you have NAs in df2 for 'on_tissue', replace them with 'off'
df2$on_tissue[is.na(df2$on_tissue)] <- "off"
df2 <- column_to_rownames(df2, var = "index")
`BM_human_AP-B08041_`@meta.data <- df2
#subset
Idents(object = `BM_human_AP-B08041_`) <- `BM_human_AP-B08041_`@meta.data[["on_tissue"]]
`BM_human_AP-B08041_` <- SubsetSTData(object = `BM_human_AP-B08041_`, ident = c("on"))

#####fin

prueba <- c(`BM_human_AP-B00182_`, `BM_human_AP-B02149_`,`BM_human_AP-B08041_`,`BM_human_AP-B08805`,
            `BM_B000943`,`BM_B01320`,`BM_B02817`,`BM_B10395`, se)
names(prueba) <- c("BM_human_AP-B00182_", "BM_human_AP-B02149_", "BM_human_AP-B08041_", "BM_human_AP-B08805",
                   "BM_B000943", "BM_B01320", "BM_B02817", "BM_B10395", "se")
lista <- prueba
lista$se <- NULL

#save
saveRDS(lista,"./objects/sp/human/human_combined.rds")
prueba <- readRDS("./objects/sp/human/human_combined.rds")

###Enrichment score
post <- c()
for (i in 1:length(prueba)){
  a <- prueba[[i]]
  b <- names(prueba[i])
  
  ## Add UCellScore PC
  vector<- ScoreSignatures_UCell(a@assays[["RNA"]]@counts, features = list(pc_mm_markers))
  a@meta.data[["signature_1_pc"]] <- as.vector(vector)
  
  color <- brewer.pal(11,"Spectral")
  color <- rev(color)
  
  p <- FeatureOverlay(a, features = c("signature_1_pc"), ncols = 1, pt.size = 1.5, 
                      value.scale = "all" ,cols = color)
  
  pdf(paste("./results/human/pc_enrich/", b,"_pc.pdf",sep=""))
  print(p)
  dev.off()
  
  ## Add UCellScore neutro
  vector<- ScoreSignatures_UCell(a@assays[["RNA"]]@counts, features = list(neutro))
  a@meta.data[["signature_1_neutro"]] <- as.vector(vector)
  
  p <- FeatureOverlay(a, features = c("signature_1_neutro"), ncols = 1, pt.size = 1.5, 
                      value.scale = "all" ,cols = color)
  
  pdf(paste("./results/human/neutro_enrich/", b,"_pc.pdf",sep=""))
  print(p)
  dev.off()
  
  ## Add UCellScore Tcell
  vector<- ScoreSignatures_UCell(a@assays[["RNA"]]@counts, features = list(human_tcell_exhausted))
  a@meta.data[["signature_1_tcell_ex"]] <- as.vector(vector)
  
  p <- FeatureOverlay(a, features = c("signature_1_tcell_ex"), ncols = 1, pt.size = 1.5, 
                      value.scale = "all" ,cols = color)
  
  pdf(paste("./results/human/tcell_ex_enrich/", b,"_tcell_ex.pdf",sep=""))
  print(p)
  dev.off()
  
  post[[length(post) + 1]] <- a
}

names(post) <- c("BM_human_AP-B00182_", "BM_human_AP-B02149_", "BM_human_AP-B08041_", "BM_human_AP-B08805",
                   "BM_B000943", "BM_B01320", "BM_B02817", "BM_B10395")

saveRDS(post,"./objects/sp/human/human_combined_enriched.rds")

se <- post[["se"]]

color <- brewer.pal(11,"Spectral")
color <- rev(color)

#pc
a <- FeatureOverlay(se, features = c("signature_1_pc"),sampleids = 1, pt.size = 0.5) +
  scale_fill_gradientn(colours = color,
                      breaks = c(0.0,0.25),
                      labels = c("Min", "Max"),
                      limits = c(0.0,0.25))
b <- FeatureOverlay(se, features = c("signature_1_pc"),sampleids = 2, pt.size = 0.5) +
  scale_fill_gradientn(colours = color,
                      breaks = c(0.0,0.25),
                      labels = c("Min", "Max"),
                      limits = c(0.0,0.25))
c <- FeatureOverlay(se, features = c("signature_1_pc"),sampleids = 3, pt.size = 0.5) +
  scale_fill_gradientn(colours = color,
                      breaks = c(0.0,0.25),
                      labels = c("Min", "Max"),
                      limits = c(0.0,0.25))
d <- FeatureOverlay(se, features = c("signature_1_pc"),sampleids = 4, pt.size = 0.5) +
  scale_fill_gradientn(colours = color,
                      breaks = c(0.0,0.25),
                      labels = c("Min", "Max"),
                      limits = c(0.0,0.25))
f <- FeatureOverlay(se, features = c("signature_1_pc"),sampleids = 5, pt.size = 0.5) +
  scale_fill_gradientn(colours = color,
                      breaks = c(0.0,0.25),
                      labels = c("Min", "Max"),
                      limits = c(0.0,0.25))
g <- FeatureOverlay(se, features = c("signature_1_pc"),sampleids = 6, pt.size = 0.5) +
  scale_fill_gradientn(colours = color,
                      breaks = c(0.0,0.25),
                      labels = c("Min", "Max"),
                      limits = c(0.0,0.25))
h <- FeatureOverlay(se, features = c("signature_1_pc"),sampleids = 7, pt.size = 0.5) +
  scale_fill_gradientn(colours = color,
                      breaks = c(0.0,0.25),
                      labels = c("Min", "Max"),
                      limits = c(0.0,0.25))
i <- FeatureOverlay(se, features = c("signature_1_pc"),sampleids = 8, pt.size = 0.5) +
  scale_fill_gradientn(colours = color,
                      breaks = c(0.0,0.25),
                      labels = c("Min", "Max"),
                      limits = c(0.0,0.25))

pdf(paste("./results/human/pc_enrich/se_pc_1.pdf",sep=""))
print(grid.arrange(a,b,c,d, ncol=2))
dev.off()
pdf(paste("./results/human/pc_enrich/se_pc_2.pdf",sep=""))
print(grid.arrange(f,g,h,i, ncol=2))
dev.off()

#neutro
a <- FeatureOverlay(se, features = c("signature_1_neutro"),sampleids = 1, pt.size = 0.5) +
  scale_fill_gradientn(colours = color,
                       breaks = c(0.0,0.2),
                       labels = c("Min", "Max"),
                       limits = c(0.0,0.2))
b <- FeatureOverlay(se, features = c("signature_1_neutro"),sampleids = 2, pt.size = 0.5) +
  scale_fill_gradientn(colours = color,
                       breaks = c(0.0,0.2),
                       labels = c("Min", "Max"),
                       limits = c(0.0,0.2))
c <- FeatureOverlay(se, features = c("signature_1_neutro"),sampleids = 3, pt.size = 0.5) +
  scale_fill_gradientn(colours = color,
                       breaks = c(0.0,0.2),
                       labels = c("Min", "Max"),
                       limits = c(0.0,0.2))
d <- FeatureOverlay(se, features = c("signature_1_neutro"),sampleids = 4, pt.size = 0.5) +
  scale_fill_gradientn(colours = color,
                       breaks = c(0.0,0.2),
                       labels = c("Min", "Max"),
                       limits = c(0.0,0.2))
f <- FeatureOverlay(se, features = c("signature_1_neutro"),sampleids = 5, pt.size = 0.5) +
  scale_fill_gradientn(colours = color,
                       breaks = c(0.0,0.2),
                       labels = c("Min", "Max"),
                       limits = c(0.0,0.2))
g <- FeatureOverlay(se, features = c("signature_1_neutro"),sampleids = 6, pt.size = 0.5) +
  scale_fill_gradientn(colours = color,
                       breaks = c(0.0,0.2),
                       labels = c("Min", "Max"),
                       limits = c(0.0,0.25))
h <- FeatureOverlay(se, features = c("signature_1_neutro"),sampleids = 7, pt.size = 0.5) +
  scale_fill_gradientn(colours = color,
                       breaks = c(0.0,0.2),
                       labels = c("Min", "Max"),
                       limits = c(0.0,0.2))
i <- FeatureOverlay(se, features = c("signature_1_neutro"),sampleids = 8, pt.size = 0.5) +
  scale_fill_gradientn(colours = color,
                       breaks = c(0.0,0.2),
                       labels = c("Min", "Max"),
                       limits = c(0.0,0.2))

pdf(paste("./results/human/neutro_enrich/se_neutro_1.pdf",sep=""))
print(grid.arrange(a,b,c,d, ncol=2))
dev.off()
pdf(paste("./results/human/neutro_enrich/se_neutro_2.pdf",sep=""))
print(grid.arrange(f,g,h,i, ncol=2))
dev.off()


#tcell
a <- FeatureOverlay(se, features = c("signature_1_tcell_ex"),sampleids = 1, pt.size = 0.5) +
  scale_fill_gradientn(colours = color,
                       breaks = c(0.0,0.2),
                       labels = c("Min", "Max"),
                       limits = c(0.0,0.2))
b <- FeatureOverlay(se, features = c("signature_1_tcell_ex"),sampleids = 2, pt.size = 0.5) +
  scale_fill_gradientn(colours = color,
                       breaks = c(0.0,0.2),
                       labels = c("Min", "Max"),
                       limits = c(0.0,0.2))
c <- FeatureOverlay(se, features = c("signature_1_tcell_ex"),sampleids = 3, pt.size = 0.5) +
  scale_fill_gradientn(colours = color,
                       breaks = c(0.0,0.2),
                       labels = c("Min", "Max"),
                       limits = c(0.0,0.2))
d <- FeatureOverlay(se, features = c("signature_1_tcell_ex"),sampleids = 4, pt.size = 0.5) +
  scale_fill_gradientn(colours = color,
                       breaks = c(0.0,0.2),
                       labels = c("Min", "Max"),
                       limits = c(0.0,0.2))
f <- FeatureOverlay(se, features = c("signature_1_tcell_ex"),sampleids = 5, pt.size = 0.5) +
  scale_fill_gradientn(colours = color,
                       breaks = c(0.0,0.2),
                       labels = c("Min", "Max"),
                       limits = c(0.0,0.2))
g <- FeatureOverlay(se, features = c("signature_1_tcell_ex"),sampleids = 6, pt.size = 0.5) +
  scale_fill_gradientn(colours = color,
                       breaks = c(0.0,0.2),
                       labels = c("Min", "Max"),
                       limits = c(0.0,0.25))
h <- FeatureOverlay(se, features = c("signature_1_tcell_ex"),sampleids = 7, pt.size = 0.5) +
  scale_fill_gradientn(colours = color,
                       breaks = c(0.0,0.2),
                       labels = c("Min", "Max"),
                       limits = c(0.0,0.2))
i <- FeatureOverlay(se, features = c("signature_1_tcell_ex"),sampleids = 8, pt.size = 0.5) +
  scale_fill_gradientn(colours = color,
                       breaks = c(0.0,0.2),
                       labels = c("Min", "Max"),
                       limits = c(0.0,0.2))

pdf(paste("./results/human/tcell_ex_enrich/se_tex_1.pdf",sep=""))
print(grid.arrange(a,b,c,d, ncol=2))
dev.off()
pdf(paste("./results/human/tcell_ex_enrich/se_tex_2.pdf",sep=""))
print(grid.arrange(f,g,h,i, ncol=2))
dev.off()



###Check interesting markers Cd81, Xbp1, Cd44, Mmp9
genes <- c("CD81", "XBP1", "CD44", "MMP9", "TNFRS17")

for (i in genes){
  a <- i
  
  color <- brewer.pal(11,"Spectral")
  color <- rev(color)
  
  p1 <- FeatureOverlay(se, features = i,sampleids = 1:4, ncols = 2,pt.size = 0.5) +
    scale_fill_gradientn(colours = color,
                         breaks = c(0.0,40),
                         labels = c("Min", "Max"),
                         limits = c(0.0,40))
  
  p2 <- FeatureOverlay(se, features = i,sampleids = 5:8, ncols = 2, pt.size = 0.5) + 
    scale_fill_gradientn(colours = color,
                         breaks = c(0.0,40),
                         labels = c("Min", "Max"),
                         limits = c(0.0,40))
  
  pdf(paste("./results/human/individual/", i,"_1.pdf",sep=""))
  print(p1)
  dev.off()
  
  pdf(paste("./results/human/individual/", i,"_2.pdf",sep=""))
  print(p2)
  dev.off()
  
  p3 <- FeatureOverlay(`BM_human_AP-B08041_`, features = i,pt.size = 1.3) + 
    scale_fill_gradientn(colours = color,
                         breaks = c(0.0,40),
                         labels = c("Min", "Max"),
                         limits = c(0.0,40))
  pdf(paste("./results/human/individual/", i,"_B08041.pdf",sep=""))
  print(p3)
  dev.off()
}

VlnPlot(object = se, features = c('CD81', 'XBP1', 'CD44', 'MMP9', "TNFRS17"))

VlnPlot(object = se, features = c("BCM", "BCMA", "CD269", "TNFRSF13A", "TNFRSF17"))

pdf(paste("./results/human/individual/", i,"_B08041.pdf",sep=""))
print(p3)
dev.off()

#correlation plot
# Create a scatter plot
ggplot(se@meta.data, aes(x = signature_1_neutro, y = signature_1_pc)) +
  geom_point() +
  labs(title = "PC vs neutro",
       x = "Signature Neutro",
       y = "Signature PC") +
  theme_minimal()

# Create a scatter plot
ggplot(se@meta.data, aes(x = se@meta.data[["signature_1_neutro"]], y = se@meta.data[["signature_1_pc"]] , color = se@meta.data[["name"]]))  +
  geom_point() +
  labs(title = "PC vs Texhausted",
       x = "Signature Texhausted",
       y = "Signature PC") +
  theme_minimal()

## Individual scatterplots
## separate list
list2env(post,envir=.GlobalEnv)

scatter <- c(`BM_human_AP-B00182_`, `BM_human_AP-B02149_`,`BM_human_AP-B08041_`,`BM_human_AP-B08805`,
            `BM_B000943`,`BM_B01320`,`BM_B02817`,`BM_B10395`)
names(scatter) <- c("BM_human_AP-B00182_", "BM_human_AP-B02149_", "BM_human_AP-B08041_", "BM_human_AP-B08805",
                   "BM_B000943", "BM_B01320", "BM_B02817", "BM_B10395")

###EScatter and correlation
for (i in 1:length(scatter)){
  a <- scatter[[i]]
  b <- names(scatter[i])
  
  # Create a scatter plot
  pdf(paste("./results/human/scatter/", b,".pdf",sep=""))
  print(ggplot(a@meta.data, aes(x = signature_1_tcell_ex, y = signature_1_pc))  +
    geom_point() +
    labs(title = paste("PC vs Texhausted", b, sep=" "),
         x = "Signature Texhausted",
         y = "Signature PC") +
    theme_minimal())
  dev.off()
  
  # Extract enrichment scores for two genes (replace 'Gene1' and 'Gene2' with your genes of interest)
  gene1_scores <- a@meta.data[["signature_1_pc"]]
  gene2_scores <- a@meta.data[["signature_1_tcell_ex"]]
  
  # Compute Spearman's rank correlation
  correlation_result <- cor(gene1_scores, gene2_scores, method="spearman")
  
  pdf(paste("./results/human/correlation/", b,"TexhaustedvsPC.pdf",sep=""))
  print(correlation_result)
  dev.off()
  
}

se <- NormalizeData(se)

FeatureOverlay(se, features = c('XBP1', 'CD44', 'MMP9'),sampleids = 1, ncols = 2,pt.size = 0.5)
se <- SCTransform(se, assay = "RNA", verbose = FALSE)
VlnPlot(object = se, features = c('XBP1', 'CD44', 'MMP9'), group.by = "name")


BM_B10395 <- prueba[["BM_B10395"]]
FeatureOverlay(BM_B10395, features = c('XBP1', 'CD44', 'MMP9'), ncols = 2,pt.size = 0.5)

###EScatter and correlation
for (i in 1:length(prueba)){
  a <- prueba[[i]]
  b <- names(prueba[i])
  a <- NormalizeData(a)
  
  # spatial
  pdf(paste("./results/human/new_pc_genes/", b,"_spatial.pdf",sep=""))
  print(FeatureOverlay(a, features = c('XBP1', 'CD44', 'MMP9'), ncols = 2,pt.size = 0.5))
  dev.off()

  # violin
  pdf(paste("./results/human/new_pc_genes/", b,"_violin.pdf",sep=""))
  print(VlnPlot(object = a, features = c('XBP1', 'CD44', 'MMP9')))
  dev.off()
}


# Iterate over each Seurat object in your list
for (i in seq_along(prueba)) {
  # Extract the current Seurat object
  se <- prueba[[i]]
  # Normalize the Seurat object
  aaaa <- NormalizeData(se)
  # Create a name for the PDF file based on the name of the Seurat object within the list
  pdf_name <- paste0("./results/human/new_pc_genes/", names(prueba)[i], ".pdf")
  # Open a PDF device to save the plot
  pdf(pdf_name)
  # Create and print the violin plot to the PDF device
  print(VlnPlot(object = aaaa, features = c('XBP1', 'CD44', 'MMP9'), assay = "RNA"))
  # Close the PDF device
  dev.off()
}


