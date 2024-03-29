## SCRIPT: Spatial object clustering to check heterogeneity and 
## analize strange cluster BM project

## 23.12.22 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(dplyr)
library(ggplot2) 
library(tidyverse)
library(base)
library(CARD)
library(gtools)
library(scatterpie)
library(ggcorrplot)
source(file = "./card.plot2.R")


## Spatial object
M1_fem_1C <- readRDS("./objects/card/M1_fem_1C_subgroup.rds")
#M1_tib_1A <- readRDS("./objects/card/M1_tib_1A_subgroup.rds")
M3_fem_1C <- readRDS("./objects/card/M3_fem_1C_subgroup.rds")
#M3_tib_2A <- readRDS("./objects/card/M3_tib_2A_subgroup.rds")

######Prepare data
#single cell
single_cell_bonemarrow <- readRDS("./objects/heterogeneity/single_cell_bonemarrow.rds")
sub_list <- levels(single_cell_bonemarrow@meta.data[["ident"]])
single_cell_bonemarrow_counts <- single_cell_bonemarrow@assays[["RNA"]]@counts
single_cell_bonemarrow_meta <- single_cell_bonemarrow@meta.data


##Merge them
combined <- merge(M1_fem_1C, y = c( M3_fem_1C), 
                  add.cell.ids = c("M1_fem_1C", "M3_fem_1C"), project = "BM")

##Merge them
#combined <- merge(M1_tib_1A, y = c( M3_tib_2A), 
#                  add.cell.ids = c("M1_tib_1A", "M3_tib_2A"), project = "BM")

##Merge them
#combined <- merge(M1_tib_1A, y = c( M3_tib_2A, M1_fem_1C, M3_fem_1C), 
#                  add.cell.ids = c("M1_tib_1A", "M3_tib_2A", "M1_fem_1C", "M3_fem_1C"), project = "BM")

x <- combined
x.image <- x@images
x@images[["M1_fem_1C"]]<- NULL
x@images[["M3_fem_1C"]]<- NULL

#x@images[["M1_tib_1A"]]<- NULL
#x@images[["M3_tib_2A"]]<- NULL

######seurat#####
list <- SplitObject(x, split.by = "type")
list <- lapply(X = list, FUN = SCTransform, assay="Spatial")
features <- SelectIntegrationFeatures(object.list = list, nfeatures = 3000)
list <- PrepSCTIntegration(object.list = list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = list, normalization.method = "SCT",
                                  anchor.features = features)
x <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

x <- RunPCA(x,npcs = 20, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:20) %>%
  FindClusters(resolution=0.8) 

x@images <- x.image

saveRDS(x, "./objects/heterogeneity/femur_integrated.rds")
x <- readRDS("./objects/heterogeneity/femur_integrated.rds")
samples <- c(x)
names(samples) <- c("femur")
#names(samples) <- c("tibia")
#names(samples) <- c("all")

for (i in 1:length(samples)){
  a <- samples[[i]]
  #umap separate in cell type
  #pdf(file.path("./results/heterogeneity/femur/",filename = paste0("umap_clusters_",names(samples[i]),".pdf")))
  #pdf(file.path("./results/heterogeneity/tibia/",filename = paste0("umap_clusters_",names(samples[i]),".pdf")))
  pdf(file.path("./results/heterogeneity/all/",filename = paste0("umap_clusters_",names(samples[i]),".pdf")))
  print(DimPlot(a, group.by = c("seurat_clusters"), label = T) + ggtitle("cell type"))
  dev.off()
  #umap separate in samples
  #pdf(file.path("./results/heterogeneity/femur/",filename = paste0("umap_samples_",names(samples[i]),".pdf")))
  #pdf(file.path("./results/heterogeneity/tibia/",filename = paste0("umap_samples_",names(samples[i]),".pdf")))
  pdf(file.path("./results/heterogeneity/all/",filename = paste0("umap_samples_",names(samples[i]),".pdf")))
  print(DimPlot(a, group.by = c("orig.ident"), label = T) + ggtitle("sample"))
  dev.off()
  #umap separate in samples
  #pdf(file.path("./results/clusters/integration/",filename = paste0("umap_area_",names(samples[i]),".pdf")))
  #print(DimPlot(a, group.by = c("area"),repel=TRUE, label = T) + ggtitle("sample"))
  #dev.off()
  #umap separate in disease
  #pdf(file.path("./results/heterogeneity/femur/",filename = paste0("umap_type_",names(samples[i]),".pdf")))
  #pdf(file.path("./results/heterogeneity/tibia/",filename = paste0("umap_type_",names(samples[i]),".pdf")))
  pdf(file.path("./results/heterogeneity/all/",filename = paste0("umap_type_",names(samples[i]),".pdf")))
  print(DimPlot(a, group.by = c("type"), label = T) + ggtitle("type"))
  dev.off()
  a <- SetIdent(a, value = a@meta.data[["seurat_clusters"]])
  #spatial umap 
  #pdf(file.path("./results/heterogeneity/femur/",filename = paste0("cluster_spatial_",names(samples[i]),".pdf")))
  #pdf(file.path("./results/heterogeneity/tibia/",filename = paste0("cluster_spatial_",names(samples[i]),".pdf")))
  pdf(file.path("./results/heterogeneity/all/",filename = paste0("cluster_spatial_",names(samples[i]),".pdf")))
  print(SpatialDimPlot(a, combine = FALSE,label.size = 1.5, label = T, crop = TRUE, pt.size.factor = 10))
  dev.off()
  #a <- SetIdent(a, value = a@meta.data[["area"]])
  #spatial umap 
  #pdf(file.path("./results/clusters/integration/",filename = paste0("area_spatial_",names(samples[i]),".pdf")))
  #print(SpatialDimPlot(a, combine = FALSE,label.size = 1.5, label = T))
  #dev.off()
  
}

Seurat::Idents(object = x) <- x@meta.data[["type"]]
x.image <- x@images
x@images[["M1_fem_1C"]]<- NULL
x@images[["M3_fem_1C"]]<- NULL
samples <- SplitObject(x, split.by = "type")
names(samples) <- c("M1_fem_1C", "M3_fem_1C")

  
for (i in 1:length(samples)){
  ####Deco
  a <- samples[[i]]
  a_count <- a@assays[["Spatial"]]@data
  a_location <- x.image[[names(samples[i])]]@coordinates
  a_location <- a_location[,2:3]
  colnames(a_location) <- c("x", "y")
  v<- ("jojoj")
  
  ###################CARD object creation
  CARD_obj_a = createCARDObject(
    sc_count = single_cell_bonemarrow_counts,
    sc_meta = single_cell_bonemarrow_meta,
    spatial_count = a_count,
    spatial_location = a_location,
    ct.varname = "ident",
    ct.select = unique(single_cell_bonemarrow_meta$ident),
    sample.varname = "metadata....experiment..",
    minCountGene = 100,
    minCountSpot = 5)
  ###################Deco
  CARD_a = CARD_deconvolution(CARD_object = CARD_obj_a)
  CARD_a_proportions <- CARD_a@Proportion_CARD
  ###################Plots
  p1 <- CARD.visualize.pie(proportion = CARD_a@Proportion_CARD,
                           spatial_location = CARD_a@spatial_location)
  pdf(paste("./results/CARD/heterogeneity/", names(samples[i]),"_1.pdf",sep=""))
  print(p1)
  dev.off()
  ## select the cell type that we are interested
  ct.visualize = sub_list
  ## visualize the spatial distribution of the cell type proportion
  p2 <- CARD.visualize.prop.2(
    proportion = CARD_a@Proportion_CARD,        
    spatial_location = CARD_a@spatial_location, 
    ct.visualize = ct.visualize,                 ### selected cell types to visualize
    colors = c("lightblue","lightyellow","red"), ### if not provide, we will use the default colors
    NumCols = 8)                                 ### number of columns in the figure panel
  pdf(paste("./results/CARD/heterogeneity/", names(samples[i]),"_2.pdf",sep=""))
  print(p2)
  dev.off()
  ## correlation
  p3 <- CARD.visualize.Cor(CARD_a@Proportion_CARD,colors = NULL) # if not provide, we will use the default colors
  pdf(paste("./results/CARD/heterogeneity/", names(samples[i]),"_3.pdf",sep=""))
  print(p3)
  dev.off()
  ## save object
  saveRDS(a,file = paste0("./objects/card/heterogeneity/",names(samples[i]),"_subgroup.rds"))
  ## save card results
  saveRDS(CARD_a,file = paste0("./objects/card/heterogeneity/",names(samples[i]),"_CARD_obj_subgroup.rds"))
  
}


##################################CELL TYPE PROPORTIONS HETEROGENEITY

## Spatial object
M1_fem_1C <- readRDS("./objects/card/heterogeneity/M1_fem_1C_subgroup.rds")
M3_fem_1C <- readRDS("./objects/card/heterogeneity/M3_fem_1C_subgroup.rds")

##add spatial image again
M1_fem_1C@images[["M1_fem_1C"]] <- x.image[["M1_fem_1C"]]
M3_fem_1C@images[["M3_fem_1C"]] <- x.image[["M3_fem_1C"]]

saveRDS(M1_fem_1C, "./objects/card/heterogeneity/M1_fem_1C_subgroup.rds")
saveRDS(M3_fem_1C, "./objects/card/heterogeneity/M3_fem_1C_subgroup.rds")

spatial_list <- c(M1_fem_1C, M3_fem_1C)
names(spatial_list) <- c("M1_fem_1C","M3_fem_1C")

## Card object
M1_fem_1C_CARD_obj <- readRDS("./objects/card/heterogeneity/M1_fem_1C_CARD_obj_subgroup.rds")
M3_fem_1C_CARD_obj <- readRDS("./objects/card/heterogeneity/M3_fem_1C_CARD_obj_subgroup.rds")

card_list <- c(M1_fem_1C_CARD_obj, M3_fem_1C_CARD_obj)
names(card_list) <- c("M1_fem_1C","M3_fem_1C")

##both
list <- list(spatial_list, card_list)
names(list) <- c("spatial_list", "card_list")

##extract proportions dataframe per sample
for (i in 1:length(list)){
  a <- list[[1]][[i]]
  cc <- list[[1]][i]
  c <- names(cc)
  b <- list[[2]][[i]]
  ## 
  a_meta <- a@meta.data
  a_area <- a_meta$area
  a_proportions <- b@Proportion_CARD
  a_proportions <- as.data.frame(a_proportions)
  cell_type <- colnames(a_proportions)
  a_proportions["area"] <- as.vector(a_area)
  write.csv(a_proportions, file = paste0("./results/CARD/heterogeneity/",c,"_pro_subgroup.csv"))
}

##proportions
M1_fem_1C <- read.csv("./results/CARD/heterogeneity/M1_fem_1C_pro_subgroup.csv")
M3_fem_1C <- read.csv("./results/CARD/heterogeneity/M3_fem_1C_pro_subgroup.csv")

pro_list <- list(M1_fem_1C, M3_fem_1C)
names(pro_list) <- c("M1_fem_1C","M3_fem_1C")

z <- colnames(M1_fem_1C)
z <- z[2:14]

######Proportions#######################################################################
#https://stackoverflow.com/questions/46205479/looping-over-multiple-lists-with-base-r
for (i in 1:length(pro_list)){
  a <- pro_list[[i]]
  v <- pro_list[i]
  # AC
  #BM_value <- a[grepl("BM", a[,14]),]
  #BM_value$area <- NULL
  
  value <- as.vector(unique(a$area))
  lista <- list()
  
  for (i in value){
    #select cluster or group of interest rows
    value_1 <- a[grepl(i, a$area),]
    value_1$area <- NULL
    value_1$X <- NULL
    
    ###create list to add content
    proportions <- c()
    for (o in colnames(value_1)){ 
      proportions <- c(proportions, (sum(value_1[[o]])*100/nrow(value_1)))
    }
    name <- paste('area:',i,sep='')
    lista[[name]] <- proportions
    
  }

}

######create a df 
rows = colnames(value_1)
df = data.frame(matrix(nrow = length(colnames(value_1)), ncol = 0)) 
rownames(df) = colnames(value_1)

##add list values to df
for (i in 1:length(lista)){
  a <- lista[[i]]
  df[, ncol(df) + 1] <- a
  names(df)[ncol(df)] <- names(lista[i])
}

#df <- t(df)
df <- cbind(celltype = rownames(df), df)
rownames(df) <- 1:nrow(df)

library(tidyr)
library(ggplot2)
DF <- data.frame(group = c(df$celltype),
                 areaBM = c(df$`area:BM`))
DFtall <- DF %>% gather(key = Area, value = Value, areaBM)
DFtall

pdf(paste("./results/CARD/heterogeneity/porcentages.pdf"))
print(ggplot(DFtall, aes(group, Value, fill = Area)) + geom_col(position = "dodge") + theme(text = element_text(size = 6)) )
dev.off()


#########PLOT deconvolution results with images

## Card object
M1_fem_1C_CARD_obj <- readRDS("./objects/card/heterogeneity/M1_fem_1C_CARD_obj_subgroup.rds")
M3_fem_1C_CARD_obj <- readRDS("./objects/card/heterogeneity/M3_fem_1C_CARD_obj_subgroup.rds")

card_list <- c(M1_fem_1C_CARD_obj, M3_fem_1C_CARD_obj)
names(card_list) <- c("M1_fem_1C","M3_fem_1C")

## Spatial object
M1_fem_1C <- readRDS("./objects/card/heterogeneity/M1_fem_1C_subgroup.rds")
M3_fem_1C <- readRDS("./objects/card/heterogeneity/M3_fem_1C_subgroup.rds")

spatial_list <- c(M1_fem_1C, M3_fem_1C)
names(spatial_list) <- c("M1_fem_1C","M3_fem_1C")

list <- list(spatial_list, card_list)
names(list) <- c("spatial_list", "card_list")

##extract proportions and add them to spatial object
for (i in 1:length(list)){
  a <- list[[1]][[i]]
  cc <- list[[1]][i]
  c <- names(cc)
  b <- list[[2]][[i]]
  ## 
  xx_df <- as.data.frame(b@Proportion_CARD)
  cells <- colnames(xx_df)
  for (o in cells){ 
    i_name <- o
    a@meta.data[[i_name]] <- xx_df[,i_name]
    i_new <- gsub("/", ".", i_name)
    ## plot
    library(BuenColors)
    color <- jdb_palette("brewer_spectra")
    
    u <- c(min(a@meta.data[[i_name]]),max(a@meta.data[[i_name]]))
    #u <- c(0,1)
    label <- c("min", "max")
    
    p1 <- SpatialFeaturePlot(a, features = i_name, pt.size.factor = 10, combine=FALSE,alpha = 0.8)
    fix.p1 <- scale_fill_gradientn(colours=color,breaks=u, labels = label,limits =u)
    p2 <- lapply(p1, function (x) x + fix.p1)
    
    pdf(paste("./results/CARD/heterogeneity/images/",c, "_" ,i_new,"_.pdf",sep=""))
    print(CombinePlots(p2))
    dev.off()
  }
  ## save object
  saveRDS(a,file = paste0("./objects/card/heterogeneity/",c,"_subgroup_deco.rds"))
}
##################################AMAIA PART HETEROGENEITY

#ERYTHROBLAST. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7689499/
CCL3 #up
KLF1 #down
GATA1 #down
# Osteoblast. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4580250/
DKK1 # up
MIP-1
CCR1 #which is one pathway contributing to their decreased osteogenic capacity
CXCL12

femur <- readRDS("./objects/heterogeneity/femur_integrated.rds")



DefaultAssay(femur) <- "SCT"
DefaultAssay(femur) <- "integrated"
SpatialDimPlot(femur, features = c("seurat_clusters"))
FeaturePlot(femur, features = "Ccr1")
DimPlot(femur, group.by = c("seurat_clusters"), label = T)
SpatialDimPlot(femur, label = T, crop = TRUE, pt.size.factor = 10)
SpatialFeaturePlot(femur,  features = c("Cxcl12"), pt.size = 10)
SpatialFeaturePlot(femur, features = c("Klf1"), pt.size = 10)
SpatialFeaturePlot(femur, features = c("Ccl3"), pt.size = 10)

Seurat::Idents(object = femur) <- femur@meta.data[["seurat_clusters"]]
markers_seurat_area <- Seurat::FindAllMarkers(object = femur, 
                                              assay = "integrated",
                                              slot = "data",
                                              verbose = TRUE, 
                                              only.pos = TRUE)


# Determine differentiating markers for cluster 6
cluster_6 <- FindMarkers(femur,
                          ident.1 = 6,
                          ident.2 = c(0,1,2,3,4,5)) 

cluster_6_f <- subset(cluster_6, p_val_adj < 0.05 & 1 < avg_log2FC)

#Top5
top50 <- cluster_6_f %>%
  top_n(n = 50,
        wt = avg_log2FC)

write.csv(top50, "./results/DE/femur/cluster6_amaia.csv", row.names =TRUE)

pdf(paste("./results/CARD/heterogeneity/images/doheatmap.pdf",sep=""))
print(DoHeatmap(femur, features = rownames(top20)) + 
  theme(text = element_text(size = 4.5)))
dev.off()

genes <- rownames(top20)

pdf(paste("./results/CARD/heterogeneity/images/featurepot_4.pdf",sep=""))
print(FeaturePlot(femur, 
            reduction = "umap", 
            features = genes[15:18],
            label = TRUE, 
            order = TRUE,
            min.cutoff = 'q10',
            repel = TRUE
))
dev.off()

cell_types <- as.list(as.vector(unique(top20_area$cluster)))
names(cell_types) <- as.vector(unique(top20_area$cluster))

M1_fem_1C <- readRDS("./objects/card/heterogeneity/M1_fem_1C_subgroup.rds")
M3_fem_1C <- readRDS("./objects/card/heterogeneity/M3_fem_1C_subgroup.rds")


DefaultAssay(M3_fem_1C) <- "SCT"
DefaultAssay(M1_fem_1C) <- "SCT"

plasma_cells_core <- c("Irf4", "Sdc1", "Xbp1")
plasma_cells_MM <- c( "Cd33", "Cd81", "Fcrl5", "Icam1",  "Ncam1")
plasma_cells_MM_mouse <- c( "Cd33", "Cd81", "Fcrl5", "Icam1",  "Ncam1")

for (i in plasma_cells_core){
  library(BuenColors)
  color <- jdb_palette("brewer_spectra")
  pdf(paste("./results/CARD/heterogeneity/images/features/PC/",i,"spatial_M1_fem_1C_plasmacore.pdf",sep=""))
  print(SpatialFeaturePlot(M1_fem_1C, features = i, pt.size.factor = 10, combine = FALSE))
  dev.off()
}

CD138
pdf(paste("./results/CARD/heterogeneity/images/features/PC/Cd27_spatial_M1_fem_1C_plasmacore.pdf",sep=""))
print(SpatialFeaturePlot(M1_fem_1C, features = c("Cd27"), pt.size.factor = 10, combine = FALSE))
dev.off()

b <- c(0,3)
p1 <- SpatialFeaturePlot(M3_fem_1C, features = "Cd81", pt.size.factor = 10, combine = FALSE)
fix.p1 <- scale_fill_gradientn(colours=color,
                               breaks=b,
                               labels=c("Min","Max"),
                               limits = b)
p2 <- lapply(p1, function (x) x + fix.p1)

pdf(paste("./results/CARD/heterogeneity/images/features/PC/Cd81_spatial_M3_fem_1C.pdf",sep=""))
print(CombinePlots(p2))
dev.off()


#####
genes_MIC <- read.csv("./data/single-cell/PC/MM_MIC_genes.csv")
genes_MIC <- as.vector(genes_MIC$x)
genes <- intersect(genes_MIC, rownames(M3_fem_1C@assays[["SCT"]]@data))
genes_MM <- intersect(genes_MIC, rownames(M1_fem_1C@assays[["SCT"]]@data))

a <- list(genes_MM)

DefaultAssay(M1_fem_1C) <- "SCT"
DefaultAssay(M1_fem_1C) <- "integrated"
DefaultAssay(M1_fem_1C) <- "Spatial"
library(GSEABase)
genes <- list(as.vector(genes))
genes_set <- GeneSet(genes)

M3_fem_1C_module <- AddModuleScore(M3_fem_1C, features = a, name = c("MIC_"))
M1_fem_1C_module <- AddModuleScore(M1_fem_1C, features = a, name = c("MIC_"))

a <-FeaturePlot(M3_fem_1C_module, features = "MIC_1", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = color)
b <-FeaturePlot(M1_fem_1C_module, features = "MIC_1", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = color)

pdf("./results/addModuleScore/heterogeneity/umap_score_MM_PC.pdf")
print(a + b)
dev.off()

##spatial plot

b <- c(0,0.6)
p1 <- SpatialFeaturePlot(M1_fem_1C_module, features = "MIC_1", pt.size.factor = 10, combine = FALSE)
fix.p1 <- scale_fill_gradientn(colours=color,
                               breaks=b,
                               labels=c("Min","Max"),
                               limits = b)
p2 <- lapply(p1, function (x) x + fix.p1)

pdf(paste("./results/addModuleScore/heterogeneity/satial_M1_fem_1C.pdf",sep=""))
print(CombinePlots(p2))
dev.off()

bl <- colorRampPalette(c("navy","royalblue","lightskyblue"))(200)                      
re <- colorRampPalette(c("mistyrose", "red2","darkred"))(200)

for (i in clusters){
  p1 <- SpatialFeaturePlot(M1_fem_1C_module, features = c("MIC_1"), pt.size.factor = 10, combine = FALSE)
  fix.p1 <- scale_fill_gradientn(colours=c(bl,"white", re), na.value = "grey98",limits = c(0,0.6))
  p2 <- lapply(p1, function (x) x + fix.p1)
  
  pdf(paste(i,".pdf",sep=""))
  print(CombinePlots(p2))
  dev.off()
}





##################################AMAIA PART HETEROGENEITY
