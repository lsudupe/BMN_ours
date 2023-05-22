## SCRIPT: PC cluster differentiation visium data BM project

## 18.05.23 Laura Sudupe , git @lsudupe

## Libraries
library(Seurat)
library(tidyverse)
library(STutility)
library(UCell)

# Data
se <- readRDS("./objects/heterogeneity/se_hierarchical.rds")

Idents(object = se) <- "name"
#M8
M8 <- SubsetSTData(se, ident = "M8_F2_1C")
meta <- M8@meta.data
rownames(meta) <- sub("_5", "", rownames(meta))
M8@meta.data <- meta
meta <- meta[, ncol(meta), drop = FALSE]

write.csv(meta, "./data/areas/m8_clustering.csv", row.names=TRUE)
M8_s <- ManualAnnotation(M8)
meta <- M8_s@meta.data
meta <- meta %>%
  mutate(pc_clusters = paste0("M8_", labels, "_cluster", clustering))
M8_s@meta.data <- meta
saveRDS(M8_s, "./objects/pc_clusters/M8_s.rds")
M8_s <- readRDS("./objects/pc_clusters/M8_s.rds")

#M2
M2 <- SubsetSTData(se, ident = "M2_F_2B")
meta <- M2@meta.data
rownames(meta) <- sub("_5", "", rownames(meta))
M2@meta.data <- meta
meta <- meta[, ncol(meta), drop = FALSE]

write.csv(meta, "./data/areas/m2_clustering.csv", row.names=TRUE)
M2_s <- ManualAnnotation(M2)
meta <- M2_s@meta.data
meta <- meta %>%
  mutate(pc_clusters = paste0("M2_", labels, "_cluster", clustering))
M2_s@meta.data <- meta
saveRDS(M2_s, "./objects/pc_clusters/M2_s.rds")
M2_s <- readRDS("./objects/pc_clusters/M2_s.rds")

#M9
M9 <- SubsetSTData(se, ident = "M9_F2_1C")
meta <- M9@meta.data
rownames(meta) <- sub("_5", "", rownames(meta))
M9@meta.data <- meta
meta <- meta[, ncol(meta), drop = FALSE]

write.csv(meta, "./data/areas/m9_clustering.csv", row.names=TRUE)
M9_s <- ManualAnnotation(M9)
meta <- M9_s@meta.data
meta <- meta %>%
  mutate(pc_clusters = paste0("M9_", labels, "_cluster", clustering))
M9_s@meta.data <- meta
saveRDS(M9_s, "./objects/pc_clusters/M9_s.rds")
M9_s <- readRDS("./objects/pc_clusters/M9_s.rds")

#M1
M1 <- SubsetSTData(se, ident = "M1_fem_1C")
meta <- M1@meta.data
rownames(meta) <- sub("_5", "", rownames(meta))
M1@meta.data <- meta
meta <- meta[, ncol(meta), drop = FALSE]

write.csv(meta, "./data/areas/m1_clustering.csv", row.names=TRUE)
#M1_s <- ManualAnnotation(M1)
meta <- M1_s@meta.data
meta <- meta %>%
  mutate(pc_clusters = paste0("M1_", labels, "_cluster", clustering))
M1_s@meta.data <- meta
saveRDS(M1_s, "./objects/pc_clusters/M1_s.rds")
M1_s <- readRDS("./objects/pc_clusters/M1_s.rds")

###combine the data
se_s <- MergeSTData(M1_s, y = c(M2_s, M8_s, M9_s), 
                  add.spot.ids = c("M1_s", "M2_s", "M8_s", "M9_s"), project = "BM")

saveRDS(se_s, "./objects/pc_clusters/combined_s.rds")

################ CHECK DORMANT SIGNATURE
objects <- c(M1_s,M2_s,M8_s,M9_s)
names(objects) <- c("M1_s","M2_s","M8_s","M9_s")
#Check dormant enrichment score
dormant <- c("C1qa", "Aif1" ,"Axl" , "Mpeg1", "H2-Eb1",
             "Nr1h3" ,"Sdc3" ,"Fcgr3" ,"Oas1a", "Ifit2", "Ly6a",
             "Ly86" ,"Bcl2a1b", "Cd19" , "Hebp1", "Sirpa" ,"Slfn2" ,"Tnfaip2" ,
             "Vpreb3", "Oasl2","Slc43a2", "Ifi44", "Gng10" ,"Cxcl16" ,"Ptprcap" ,
             "Anxa2" ,"Rgs2" ,"Tmed3" ,"Igll1", "Hpgd", "Glipr1", "Cd4" ,"Cd84" ,"Gbp2", "AB124611", "Slc44a2" ,
             "Samd9l" ,"Oas1g" ,"Fcgr1", "Pla2g15", "Tifa" ,"Pmp22", "Abcc3" ,"S100a10")

library(RColorBrewer)
color <- brewer.pal(11,"Spectral")
color <- rev(color)


new <- c()
for (i in 1:length(objects)){
  a <- objects[[i]]
  b <- names(objects[i])
  
  ## Add UCellScore
  vector<- ScoreSignatures_UCell(a@assays[["RNA"]]@counts, features = list(dormant))
  a@meta.data[["signature_1_dormant"]] <- as.vector(vector)
  
  meta <- a@meta.data
  
  lm <- lm(meta$signature_1_dormant ~ meta$MM_MIC, data =meta)
  residuals <- lm$residuals
  a@meta.data[["residuals_dormant_ucellscore"]] <- residuals
  
  ##plots
  p <- FeatureOverlay(a, features = c("signature_1_dormant"), pt.size = 1.8, 
                      value.scale = "all" ,cols = color)
  
  pdf(paste("./results/pc_clusters/", b,"_spatial_Ucell.pdf",sep=""))
  print(p)
  dev.off()
  
  p <- FeatureOverlay(a, features = c("residuals_dormant_ucellscore"), pt.size = 1.8, 
                      value.scale = "all" ,cols = color)
  
  pdf(paste("./results/pc_clusters/", b,"_spatial_Ucell_residuals.pdf",sep=""))
  print(p)
  dev.off()
  
  ## add object to list
  new[[length(new) + 1]] <- a
}

names(new) <- c("M1_s","M2_s","M8_s","M9_s")
## separate list
list2env(new,envir=.GlobalEnv)

###plots
###merge data
se_s <- MergeSTData(M1_s, y = c(M2_s, M8_s, M9_s), 
                    add.spot.ids = c("M1_s", "M2_s", "M8_s", "M9_s"), project = "BM")

saveRDS(se_s, "./objects/pc_clusters/combined_s_dormant.rds")
se_s <- readRDS("./objects/pc_clusters/combined_s_dormant.rds")

##subset 
Idents(object = se_s) <- se_s@meta.data[["clustering"]]
subset <- SubsetSTData(se_s, idents = c("5", "6", "7"))

se_s@meta.data[["pc_clusters"]] <- as.factor(se_s@meta.data[["pc_clusters"]])

meta <- se_s@meta.data
meta_s <- subset(meta, clustering %in% c("5", "6", "7"))
meta_s <- meta_s[meta_s$pc_clusters != "M2_bone_marrow_cluster7", ]
se_s@meta.data <- meta_s

se_s@meta.data[["pc_clusters"]] <- as.factor(se_s@meta.data[["pc_clusters"]])

##spatial plots
pdf(file.path("./results/pc_clusters/clusters.pdf"))
FeatureOverlay(subset, features = "pc_clusters", ncols = 2,pt.size = 0.7)
dev.off()


p <- FeatureOverlay(se_s, features = c("residuals_dormant_ucellscore"), pt.size = 1.8, 
                    value.scale = "all" ,cols = color)

pdf(paste("./results/pc_clusters/", b,"_spatial_Ucell_residuals.pdf",sep=""))
print(p)
dev.off()

#set colors
n = 7
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols = gg_color_hue(n)
cluster_color_map <- setNames(cols, unique(se_s@meta.data[["clustering"]]))

