## SCRIPT: check pc markers in human BM project

## 31.05.23 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(ggplot2)
library(STutility)
library(dplyr)

#Data---------------------------------
pc_mm_markers <- c("SHH", "DHH", "IHH", "PTCH1", "PTCH2", "SMO", "SUFU", "GLI1", 
                   "GLI2", "GLI3", "CD19", "CD44", "CXCR4", "KLF4", "CD28", "CD33", "CD27")



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

prueba <- c(`BM_human_AP-B00182_`, `BM_human_AP-B02149_`,`BM_human_AP-B08041_`,`BM_human_AP-B08805`)
names(prueba) <- c("BM_human_AP-B00182_", "BM_human_AP-B02149_", "BM_human_AP-B08041_", "BM_human_AP-B08805")

library(RColorBrewer)
color <- brewer.pal(11,"Spectral")
color <- rev(color)

for (i in 1:length(prueba)){
  a <- prueba[[i]]
  b <- names(prueba[i])
  
  for (g in pc_mm_markers){
    p <- FeatureOverlay(a, features = g, pt.size = 1.8, 
                        value.scale = "all" ,cols = color)
    
    pdf(paste("./results/human/pc_genes/", g,"_", b, ".pdf",sep=""))
    print(p)
    dev.off()
  }
}

##prueba
for (i in pc_mm_markers){
  
  p <- FeatureOverlay(se, features = i, sampleids = 1:4,
                      pt.size = 0.4, 
                      value.scale = "all" ,cols = color)
    
  pdf(paste("./results/human/prueba/", i, ".pdf",sep=""))
  print(p)
  dev.off()
  
}




