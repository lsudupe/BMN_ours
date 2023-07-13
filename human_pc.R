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
library(tibble)

# Read data, human sample
#Variables---------------------------------
pc_mm <- readRDS("./data/single-cell/human/pc_mm/MM_multiome_list.rds")

# Extract the gene list from the Seurat object
seurat_genes <- rownames(pc_mm[["MM3"]]@assays$RNA@counts)

#pc_markers
pc_mm_markers <- c("SHH", "DHH", "IHH", "PTCH1", "PTCH2", "SMO", "SUFU", "GLI1", 
                   "GLI2", "GLI3", "CD19", "CD44", "CXCR4", "KLF4", "CD28", "CD33", "CD27")

#CD45,CD56,CD117
# Check which genes_to_check are in seurat_genes
genes_present <- pc_mm_markers %in% seurat_genes

# Print out the results
print(genes_present)

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
            `BM_B000943`,`BM_B01320`,`BM_B02817`,`BM_B10395`)
names(prueba) <- c("BM_human_AP-B00182_", "BM_human_AP-B02149_", "BM_human_AP-B08041_", "BM_human_AP-B08805",
                   "BM_B000943", "BM_B01320", "BM_B02817", "BM_B10395")


###Enrichment score
for (i in 1:length(prueba)){
  a <- prueba[[i]]
  b <- names(prueba[i])
  
  ## Add UCellScore
  vector<- ScoreSignatures_UCell(a@assays[["RNA"]]@counts, features = list(pc_mm_markers))
  a@meta.data[["signature_1_pc"]] <- as.vector(vector)
  
  color <- brewer.pal(11,"Spectral")
  color <- rev(color)
  
  p <- FeatureOverlay(a, features = c("signature_1_pc"), ncols = 1, pt.size = 1.5, 
                      value.scale = "all" ,cols = color)
  
  pdf(paste("./results/human/pc_enrich/", b,"_pc.pdf",sep=""))
  print(p)
  dev.off()
  
}
