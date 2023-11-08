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


###markers
endothelial_markers <- c("CD34", "CDH5", "EMCN", "KDR", "PECAM1", "IL6ST", "VWF", 
                         "TEK", "CD9", "LY6E", "VCAM1", "CAVIN2", "ICAM2", "ADGRF5", 
                         "ADGRL4", "LY6C1", "TIE2", "ARHGAP31", "CLDN5", "COTL1", "CTNNBIP1", 
                         "CXCL12", "EFNA1", "ESAM", "ETS2", "FAM167B", "FAM212A", "FCER1G", 
                         "GPIHBP1", "HSPA1A", "IFGBP7", "ITGB1", "KLF2", "KLF6", "LALBA", 
                         "LIMS2", "PDK1", "PLIN2", "PLVAP", "PREX2", "RPL13A", "SPARCL1", 
                         "TCN2", "THRSP", "TINAGL1", "TM4SF1", "TRP53I11", "TXNIP", "VIM", 
                         "YES1", "ICAM1")

msc_markers <- c("LEPR", "CXCL12", "VCAM1", "ANGPT1", "KITLG", "ENG", "LPL", 
                 "CP", "COL1A1", "PDGFRB", "HP", "PTH1R", "NT5E", "THY1", "NES", 
                 "PDGFRA", "GREM1")

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

###Enrichment score
post <- c()
for (i in 1:length(objects)){
  a <- objects[[i]]
  b <- names(objects[i])
  
  ## Add UCellScore endo
  vector<- ScoreSignatures_UCell(a@assays[["RNA"]]@counts, features = list(endothelial_markers))
  a@meta.data[["signature_1_endo"]] <- as.vector(vector)
  
  color <- brewer.pal(11,"Spectral")
  color <- rev(color)
  
  p <- FeatureOverlay(a, features = c("signature_1_endo"), ncols = 1, pt.size = 1.5, 
                      value.scale = "all" ,cols = color)
  
  pdf(paste("./results/human/healthy/endo_enrich/", b,"_endo.pdf",sep=""))
  print(p)
  dev.off()
  
  ## Add UCellScore msc
  vector<- ScoreSignatures_UCell(a@assays[["RNA"]]@counts, features = list(msc_markers))
  a@meta.data[["signature_1_msc"]] <- as.vector(vector)
  
  p <- FeatureOverlay(a, features = c("signature_1_msc"), ncols = 1, pt.size = 1.5, 
                      value.scale = "all" ,cols = color)
  
  pdf(paste("./results/human/healthy/msc_enrich/", b,"_msc.pdf",sep=""))
  print(p)
  dev.off()
  
  post[[length(post) + 1]] <- a
}

names(post) <- c("BM_human_H13", "BM_human_H6", "BM_human_H7", "BM_human_H9")








