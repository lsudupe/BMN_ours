## SCRIPT: Plotting script Bone Marrow project

## 21.06.23 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(STutility)
library(RColorBrewer)
library(UCell)

# Data
se <- readRDS("./objects/heterogeneity/se_hierarchical.rds")

b <- SetIdent(se, value = se@meta.data[["clustering"]])
b <- SubsetSTData(b, idents = c("5", "6","7"))
b@meta.data[["hotspots"]] <- b@active.ident

b <- SetIdent(b, value = b@meta.data[["name"]])
subset <- SubsetSTData(b, idents = c("M1_fem_1C", "M2_F_2B", "M8_F2_1C", "M9_F2_1C"))

#spatial
##divide by sample
Idents(object = subset) <- "name"
name <- unique(subset@meta.data[["name"]])

objects <- c()
subset <- SetIdent(subset, value = subset@meta.data[["name"]])

for (i in name){
  a <- SubsetSTData(subset, idents = i)
  ## add object to list
  objects[[length(objects) + 1]] <- a
}

names(objects) <- name
## separate list
list2env(objects,envir=.GlobalEnv)

###create new objects 
new_list <- list()
for (i in 1:length(objects)){
  a <- objects[[i]]
  a <- SetIdent(a, value = a@meta.data[["hotspots"]])
  unique_hotspots <- unique(a@meta.data[["hotspots"]])
  b <- names(objects[i])
  
  ##dormant
  ## Add UCellScore
  vector<- ScoreSignatures_UCell(a@assays[["RNA"]]@counts, features = list(dormant))
  a@meta.data[["signature_1_dormant"]] <- as.vector(vector)
  meta <- a@meta.data
  lm <- lm(meta$signature_1_dormant ~ meta$MM_MIC, data =meta)
  residuals <- lm$residuals
  a@meta.data[["residuals_dormant"]] <- residuals
  
  vector<- ScoreSignatures_UCell(a@assays[["RNA"]]@counts, features = list(genes38))
  a@meta.data[["signature_1_genes38"]] <- as.vector(vector)
  meta <- a@meta.data
  lm <- lm(meta$signature_1_genes38 ~ meta$MM_MIC, data =meta)
  residuals <- lm$residuals
  a@meta.data[["residuals_genes38"]] <- residuals
  
  ##g2m
  ## Add UCellScore
  vector<- ScoreSignatures_UCell(a@assays[["RNA"]]@counts, features = list(g2m))
  a@meta.data[["signature_1_g2m"]] <- as.vector(vector)
  meta <- a@meta.data
  lm <- lm(meta$signature_1_g2m ~ meta$MM_MIC, data =meta)
  residuals <- lm$residuals
  a@meta.data[["residuals_g2m"]] <- residuals
  
  ##s
  ## Add UCellScore
  vector<- ScoreSignatures_UCell(a@assays[["RNA"]]@counts, features = list(s))
  a@meta.data[["signature_1_s"]] <- as.vector(vector)
  meta <- a@meta.data
  lm <- lm(meta$signature_1_s ~ meta$MM_MIC, data =meta)
  residuals <- lm$residuals
  a@meta.data[["residuals_s"]] <- residuals

  
  p <- VlnPlot(a, features = c("Mmp9"), group.by = "clustering")
  pdf(paste("./results/plotting_bio/mmp9/", b,"_violin_Mmp9.pdf",sep=""))
  print(p)
  dev.off()
  
  p <- VlnPlot(a, features = c("Cd44"), group.by = "clustering")
  pdf(paste("./results/plotting_bio/cd44/", b,"_violin_Cd44.pdf",sep=""))
  print(p)
  dev.off()
  
  p <- VlnPlot(a, features = c("residuals_dormant"), group.by = "clustering")
  pdf(paste("./results/plotting_bio/dormant/", b,"_violin_dormant.pdf",sep=""))
  print(p)
  dev.off()
  
  p <- VlnPlot(a, features = c("residuals_genes38"), group.by = "clustering")
  pdf(paste("./results/plotting_bio/38genes/", b,"_violin_genes38.pdf",sep=""))
  print(p)
  dev.off()
  
  p <- VlnPlot(a, features = c("residuals_g2m"), group.by = "clustering")
  pdf(paste("./results/plotting_bio/proliferation/g2m/", b,"_violin_g2m.pdf",sep=""))
  print(p)
  dev.off()
  
  p <- VlnPlot(a, features = c("residuals_s"), group.by = "clustering")
  pdf(paste("./results/plotting_bio/proliferation/s/", b,"_violin_s.pdf",sep=""))
  print(p)
  dev.off()
  
  
  # Iterate over each unique hotspo
  for (j in unique_hotspots){
    # Subset the data for the current hotspot
    b <- SubsetSTData(a, idents = j)
    # Name the new object
    name <- paste(names(objects)[i], j, sep = "_")
    # Add it to the new list
    new_list[[name]] <- b
  }
}

##create a list by sample, for scaling purpouses

M1 <- c(new_list[["M1_fem_1C_5"]], new_list[["M1_fem_1C_6"]], new_list[["M1_fem_1C_7"]])
M2 <- c(new_list[["M2_F_2B_5"]], new_list[["M2_F_2B_6"]], new_list[["M2_F_2B_7"]])
M8 <- c(new_list[["M8_F2_1C_5"]], new_list[["M8_F2_1C_7"]], new_list[["M8_F2_1C_6"]])
M9 <- c(new_list[["M9_F2_1C_6"]], )

#######variables to check
Cd44 <- "Cd44"
Mmp9 <- "Mmp9"
genes38 <- c("Tnfrsf17", "Sdc1", "Cd79a", "Plpp5", 
           "Fcrl5", "Slamf7", "Cd38", "Cd40", "Tnfrsf13b", 
           "Cadm1", "Ms4a1", "Cd79b", "Ednrb", "Cav1", "Kcnn3", 
           "Gpr160", "Lsr", "Perp", "Fcrl2", "Gprc5d", "Itga8", 
           "Il5ra", "Parm1", "Lsamp", "Dpep1", "Tnfrsf13c", "Cd19", 
           "Slc1a4", "P2rx5", "Cd22", "Fecer2", "Adam28", "Cd24", 
           "Fcrl1", "Lax1", "Hvcn1", "Ccr10", "Dcc")
g2m <- c("Gtse1", "Ube2c", "Ccnb2", "Mki67", "Aurka", 
         "Tpx2", "Ndc80", "Anln", "Rangap1", "Cks2", "Cdca2", 
         "Ckap5", "Cenpa", "Tubb4b", "Cdca8", "Gas2l3", "Ckap2l", 
         "Kif20b", "G2e3", "Ctcf", "Hmmr", "Hmgb2", "Top2a", 
         "Kif11", "Nusap1", "Kif23", "Dlgap5", "Ckap2", "Ttk", 
         "Nek2", "Tacc3", "Ect2", "Cdca3", "Anp32e", "Ncapd2", 
         "Aurkb", "Bub1", "Cdc25c", "Lbr", "Smc4", "Cenpf",
         "Cdk1", "Cenpe", "Hjurp", "Cks1b", "Birc5", "Psrc1", 
         "Cdc20", "Nuf2", "Kif2c", "Cbx5")
s <- c("Chaf1b", "Pola1", "Pcna", "Tyms",
       "Cdc6", "Rad51", "Cdc45", "Mcm5", "Ubr7", 
       "E2f8", "Wdr76", "Hells", "Clspn", "Gmnn", "Mcm4",
       "Usp1", "Casp8ap2", "Msh2", "Mcm6", "Cdca7",
       "Fen1", "Dscc1", "Uhrf1", "Rrm1", "Slbp",
       "Rfc2", "Gins2", "Rpa2", "Prim1", "Mcm2", 
       "Dtl", "Tipin", "Brip1", "Nasp", "Blm", 
       "Rad51ap1", "Ccne2", "Exo1", "Rrm2", "Ung")
dormant <- c("C1qa", "Aif1" ,"Axl" , "Mpeg1", "H2-Eb1",
             "Nr1h3" ,"Sdc3" ,"Fcgr3" ,"Oas1a", "Ifit2", "Ly6a",
             "Ly86" ,"Bcl2a1b", "Cd19" , "Hebp1", "Sirpa" ,"Slfn2" ,"Tnfaip2" ,
             "Vpreb3", "Oasl2","Slc43a2", "Ifi44", "Gng10" ,"Cxcl16" ,"Ptprcap" ,
             "Anxa2" ,"Rgs2" ,"Tmed3" ,"Igll1", "Hpgd", "Glipr1", "Cd4" ,"Cd84" ,"Gbp2", "AB124611", "Slc44a2" ,
             "Samd9l" ,"Oas1g" ,"Fcgr1", "Pla2g15", "Tifa" ,"Pmp22", "Abcc3" ,"S100a10")


#######variables to check FIN

#####create a plotting loop
color <- brewer.pal(11,"Spectral")
color <- rev(color)

for (i in 1:length(new_list)){
  a <- new_list[[i]]
  b <- names(new_list[i])
  
  #Mmp9
  # Check if 'Mmp9' is a feature of the object
  if (!"Mmp9" %in% rownames(a@assays[["RNA"]]@counts)) {
    next
  }
  
  p <- FeatureOverlay(a, features = c("Mmp9"), pt.size = 1.8, 
                      min.cutoff = (0),
                      max.cutoff = (5),
                      value.scale = "all" ,cols = color)
  
  pdf(paste("./results/plotting_bio/mmp9/", b,"_Mmp9.pdf",sep=""))
  print(p)
  dev.off()
  
  #Cd44
  p <- FeatureOverlay(a, features = c(Cd44), pt.size = 1.8, 
                      value.scale = "all" ,cols = color)
  
  pdf(paste("./results/plotting_bio/cd44/", b,"_Cd44.pdf",sep=""))
  print(p)
  dev.off()
  
  ##dormant
  ## Add UCellScore
  vector<- ScoreSignatures_UCell(a@assays[["RNA"]]@counts, features = list(dormant))
  a@meta.data[["signature_1_dormant"]] <- as.vector(vector)
  meta <- a@meta.data
  lm <- lm(meta$signature_1_dormant ~ meta$MM_MIC, data =meta)
  residuals <- lm$residuals
  a@meta.data[["residuals_dormant"]] <- residuals
  ##plots
  p <- FeatureOverlay(a, features = c("residuals_dormant"), pt.size = 1.8, 
                      value.scale = "all" ,cols = color)
  pdf(paste("./results/plotting_bio/dormant/", b,"_spatial_dormant.pdf",sep=""))
  print(p)
  dev.off()
  
  ##genes38
  ## Add UCellScore
  vector<- ScoreSignatures_UCell(a@assays[["RNA"]]@counts, features = list(genes38))
  a@meta.data[["signature_1_genes38"]] <- as.vector(vector)
  meta <- a@meta.data
  lm <- lm(meta$signature_1_genes38 ~ meta$MM_MIC, data =meta)
  residuals <- lm$residuals
  a@meta.data[["residuals_genes38"]] <- residuals
  ##plots
  p <- FeatureOverlay(a, features = c("residuals_genes38"), pt.size = 1.8, 
                      value.scale = "all" ,cols = color)
  pdf(paste("./results/plotting_bio/38genes/", b,"_spatial_genes38.pdf",sep=""))
  print(p)
  dev.off()
  
  ##g2m
  ## Add UCellScore
  vector<- ScoreSignatures_UCell(a@assays[["RNA"]]@counts, features = list(g2m))
  a@meta.data[["signature_1_g2m"]] <- as.vector(vector)
  meta <- a@meta.data
  lm <- lm(meta$signature_1_g2m ~ meta$MM_MIC, data =meta)
  residuals <- lm$residuals
  a@meta.data[["residuals_g2m"]] <- residuals
  ##plots
  p <- FeatureOverlay(a, features = c("residuals_g2m"), pt.size = 1.8, 
                      value.scale = "all" ,cols = color)
  pdf(paste("./results/plotting_bio/proliferation/g2m/", b,"_spatial_g2m.pdf",sep=""))
  print(p)
  dev.off()
  
  ##s
  ## Add UCellScore
  vector<- ScoreSignatures_UCell(a@assays[["RNA"]]@counts, features = list(s))
  a@meta.data[["signature_1_s"]] <- as.vector(vector)
  meta <- a@meta.data
  lm <- lm(meta$signature_1_s ~ meta$MM_MIC, data =meta)
  residuals <- lm$residuals
  a@meta.data[["residuals_s"]] <- residuals
  ##plots
  p <- FeatureOverlay(a, features = c("residuals_s"), pt.size = 1.8, 
                      value.scale = "all" ,cols = color)
  pdf(paste("./results/plotting_bio/proliferation/s/", b,"_spatial_s.pdf",sep=""))
  print(p)
  dev.off()
  
  
}


for (i in 1:length(new_list)){
  a <- new_list[[i]]
  b <- names(objects[i])
  
  
}



FeatureOverlay(M8_F2_1C, features = c("Mmp9"), pt.size = 1.8, 
               value.scale = "all" ,cols = color)

p1 <- ST.FeaturePlot(new_list[["M2_F_2B_5"]], features = c("Mmp9"), ncol = 2, show.sb = FALSE, palette = "Spectral")
FeatureOverlay(new_list[["M2_F_2B_5"]], 
               features = c("Mmp9"), 
               cols = c("lightgray", "mistyrose", "red", "dark red", "black"), 
               type = "raw",
               ncols = 2)
p1 <- ST.FeaturePlot(new_list[["M2_F_2B_5"]], features = c("Mmp9"), ncol = 2, show.sb = FALSE, palette = "Spectral")

p2 <- ST.FeaturePlot(new_list[["M2_F_2B_5"]], features = "hotspots", split.labels = T, pt.size = 2) & 
  theme(plot.title = element_blank(), strip.text = element_blank())

p1 +p2
cowplot::plot_grid(p1, p2, ncol = 2)

prueba <- new_list[["M2_F_2B_5"]]
prueba@active.ident <- prueba@meta.data[["hotspots"]]
a <- SetIdent(a, value = a@meta.data[["hotspots"]])

VlnPlot(a, features = c("Mmp9"), group.by = "hotspots")









