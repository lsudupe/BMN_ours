## SCRIPT: Deconvolution using CARD areas subgroups BM project

## 01.11.22 Laura Sudupe , git @lsudupe
#https://github.com/YingMa0107/CARD/

#Libraries---------------------------------
library(Seurat)
library(ggplot2)
library(dplyr)
#devtools::install_github('YingMa0107/CARD')
library(CARD)
library(gtools)
library(scatterpie)
library(ggcorrplot)

source(file = "./card.plot2.R")


#Data---------------------------------
spatial <- readRDS("./objects/sp/second/combined_filtered.rds")
single_cell_bonemarrow <- readRDS("./objects/sc/single_cell_bonemarrow.rds")


######Prepare data
#single cell
sub_list <- levels(single_cell_bonemarrow@meta.data[["cell"]])

single_cell_bonemarrow_counts <- single_cell_bonemarrow@assays[["RNA"]]@counts
single_cell_bonemarrow_meta <- single_cell_bonemarrow@meta.data
#single_cell_bonemarrow_meta <- single_cell_bonemarrow_meta[, 6:7]

#spatial
x <- spatial
x.image <- x@images
x@images[["M1_tib_1A"]]<- NULL
x@images[["M1_fem_1C"]]<- NULL
x@images[["M3_tib_2A"]]<- NULL
x@images[["M3_fem_1C"]]<- NULL
list <- SplitObject(x, split.by = "orig.ident")

for (i in 1:length(list)){
  a <- list[[i]]
  a <- SCTransform(a, assay="Spatial")
  a <- RunPCA(a, npcs = 20, verbose = FALSE, assay="SCT")
  a <- RunUMAP(a, reduction = "pca", dims = 1:20, verbose = FALSE)
  a <- FindNeighbors(a, reduction = "pca", dims = 1:20)
  a <- FindClusters(a, resolution=0.5)
  a_count <- a@assays[["Spatial"]]@data
  a_location <- x.image[[names(list[i])]]@coordinates
  a_location <- a_location[,2:3]
  colnames(a_location) <- c("x", "y")
  ###################CARD object creation
  CARD_obj_a = createCARDObject(
    sc_count = single_cell_bonemarrow_counts,
    sc_meta = single_cell_bonemarrow_meta,
    spatial_count = a_count,
    spatial_location = a_location,
    ct.varname = "cell",
    ct.select = unique(single_cell_bonemarrow_meta$cell),
    sample.varname = "metadata....experiment..",
    minCountGene = 100,
    minCountSpot = 5)
  v<- ("jojoj")
  ###################Deco
  CARD_a = CARD_deconvolution(CARD_object = CARD_obj_a)
  CARD_a_proportions <- CARD_a@Proportion_CARD
  ###################Plots
  p1 <- CARD.visualize.pie(proportion = CARD_a@Proportion_CARD,
                           spatial_location = CARD_a@spatial_location)
  pdf(paste("./results/CARD/subgroups/", names(list[i]),"_1.pdf",sep=""))
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
    NumCols = 4)                                 ### number of columns in the figure panel
  pdf(paste("./results/CARD/subgroups/", names(list[i]),"_2.pdf",sep=""))
  print(p2)
  dev.off()
  ## correlation
  p3 <- CARD.visualize.Cor(CARD_a@Proportion_CARD,colors = NULL) # if not provide, we will use the default colors
  pdf(paste("./results/CARD/subgroups/", names(list[i]),"_3.pdf",sep=""))
  print(p3)
  dev.off()
  ## save object
  saveRDS(a,file = paste0("./objects/card/",names(list[i]),"_subgroup.rds"))
  ## save card results
  saveRDS(CARD_a,file = paste0("./objects/card/",names(list[i]),"_CARD_obj_subgroup.rds"))
}

## Card object
M1_fem_1C_CARD_obj <- readRDS("./objects/card/M1_fem_1C_CARD_obj_subgroup.rds")
M1_tib_1A_CARD_obj <- readRDS("./objects/card/M1_tib_1A_CARD_obj_subgroup.rds")
M3_fem_1C_CARD_obj <- readRDS("./objects/card/M3_fem_1C_CARD_obj_subgroup.rds")
M3_tib_2A_CARD_obj <- readRDS("./objects/card/M3_tib_2A_CARD_obj_subgroup.rds")

card_list <- c(M1_fem_1C_CARD_obj, M1_tib_1A_CARD_obj, M3_fem_1C_CARD_obj, M3_tib_2A_CARD_obj)
names(card_list) <- c("M1_fem_1C_CARD_obj","M1_tib_1A_CARD_obj","M3_fem_1C_CARD_obj","M3_tib_2A_CARD_obj")

card_list <- c(M3_fem_1C_CARD_obj, M3_tib_2A_CARD_obj)
names(card_list) <- c("M3_fem_1C_CARD_obj","M3_tib_2A_CARD_obj")

## Spatial object
M1_fem_1C <- readRDS("./objects/card/M1_fem_1C_subgroup.rds")
M1_tib_1A <- readRDS("./objects/card/M1_tib_1A_subgroup.rds")
M3_fem_1C <- readRDS("./objects/card/M3_fem_1C_subgroup.rds")
M3_tib_2A <- readRDS("./objects/card/M3_tib_2A_subgroup.rds")

spatial_list <- c(M1_fem_1C, M1_tib_1A, M3_fem_1C, M3_tib_2A)
names(spatial_list) <- c("M1_fem_1C","M1_tib_1A","M3_fem_1C","M3_tib_2A")

spatial_list <- c(M3_fem_1C, M3_tib_2A)
names(spatial_list) <- c("M3_fem_1C","M3_tib_2A")

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
  write.csv(a_proportions, file = paste0("./results/CARD/",c,"_pro_subgroup.csv"))
}

M1_fem_1C <- read.csv("./results/CARD/M1_tib_1A_pro_subgroup.csv")
M1_tib_1A <- read.csv("./results/CARD/M1_tib_1A_pro_subgroup.csv")
M3_fem_1C <- read.csv("./results/CARD/M3_fem_1C_pro_subgroup.csv")
M3_tib_2A <- read.csv("./results/CARD/M3_tib_2A_pro_subgroup.csv")

pro_list <- list(M1_fem_1C, M1_tib_1A, M3_fem_1C, M3_tib_2A)
names(pro_list) <- c("M1_fem_1C","M1_tib_1A","M3_fem_1C","M3_tib_2A")

z <- colnames(M1_fem_1C)
z <- z[2:33]

######Proportions#######################################################################
#https://stackoverflow.com/questions/46205479/looping-over-multiple-lists-with-base-r
for (i in 1:length(pro_list)){
  a <- pro_list[[i]]
  v <- pro_list[i]
  v <- names(v)
  # AC
  BM_value <- a[grepl("BM", a[,34]),]
  BM_value$area <- NULL
    BM_Ery.Mk.prog. <- (sum(BM_value$Ery.Mk.prog.))*100/nrow(BM_value)
    BM_Neutro.prog. <- (sum(BM_value$Neutro.prog.))*100/nrow(BM_value)
    BM_Mono.prog. <- (sum(BM_value$Mono.prog.))*100/nrow(BM_value)
    BM_Gran.Mono.prog. <- (sum(BM_value$Gran.Mono.prog.))*100/nrow(BM_value)
    BM_LMPPs <- (sum(BM_value$LMPPs))*100/nrow(BM_value)
    
    BM_large.pre.B. <- (sum(BM_value$large.pre.B.))*100/nrow(BM_value)
    BM_Mk.prog. <- (sum(BM_value$Mk.prog.))*100/nrow(BM_value)
    BM_Erythroblasts <- (sum(BM_value$Erythroblasts))*100/nrow(BM_value)
    BM_Eo.Baso.prog. <- (sum(BM_value$Eo.Baso.prog.))*100/nrow(BM_value)
    BM_Monocytes <- (sum(BM_value$Monocytes))*100/nrow(BM_value)
    
    BM_Ery.prog. <- (sum(BM_value$Ery.prog.))*100/nrow(BM_value)
    BM_pro.B <- (sum(BM_value$pro.B))*100/nrow(BM_value)
    BM_T.cells <- (sum(BM_value$T.cells))*100/nrow(BM_value)
    BM_Neutrophils <- (sum(BM_value$Neutrophils))*100/nrow(BM_value)
    BM_Adipo.CAR <- (sum(BM_value$Adipo.CAR))*100/nrow(BM_value)
    
    BM_Ng2..MSCs <- (sum(BM_value$Ng2..MSCs))*100/nrow(BM_value)
    BM_Osteoblasts <- (sum(BM_value$SOsteoblasts))*100/nrow(BM_value)
    BM_Schwann.cells <- (sum(BM_value$Schwann.cells))*100/nrow(BM_value)
    BM_Arteriolar.fibro. <- (sum(BM_value$Arteriolar.fibro.))*100/nrow(BM_value)
    BM_Sinusoidal.ECs <- (sum(BM_value$Sinusoidal.ECs))*100/nrow(BM_value)
    
    BM_Osteo.CAR <- (sum(BM_value$Osteo.CAR))*100/nrow(BM_value)
    BM_small.pre.B. <- (sum(BM_value$small.pre.B.))*100/nrow(BM_value)
    BM_Chondrocytes <- (sum(BM_value$Chondrocytes))*100/nrow(BM_value)
    BM_Endosteal.fibro. <- (sum(BM_value$Endosteal.fibro.))*100/nrow(BM_value)
    BM_Fibro.Chondro.p. <- (sum(BM_value$Fibro.Chondro.p..))*100/nrow(BM_value)
    
    BM_Stromal.fibro. <- (sum(BM_value$Stromal.fibro.))*100/nrow(BM_value)
    BM_Arteriolar.ECs <- (sum(BM_value$Arteriolar.ECs))*100/nrow(BM_value)
    BM_Myofibroblasts <- (sum(BM_value$Myofibroblasts))*100/nrow(BM_value)
    BM_Smooth.muscle <- (sum(BM_value$Smooth.muscle))*100/nrow(BM_value)
    BM_Dendritic.cells <- (sum(BM_value$Dendritic.cells))*100/nrow(BM_value)
    
    BM_NK.cells <- (sum(BM_value$NK.cells))*100/nrow(BM_value)
    BM_B.cell <- (sum(BM_value$B.cell))*100/nrow(BM_value)
    
    
  BM_proportions <- c(BM_Ery.Mk.prog., BM_Neutro.prog. , BM_Mono.prog., BM_Gran.Mono.prog., BM_LMPPs,
                      BM_large.pre.B., BM_Mk.prog., BM_Erythroblasts, BM_Eo.Baso.prog., BM_Monocytes,
                      BM_Ery.prog., BM_pro.B, BM_T.cells, BM_Neutrophils, BM_Adipo.CAR,
                      BM_Ng2..MSCs, BM_Osteoblasts, BM_Schwann.cells, BM_Arteriolar.fibro., BM_Sinusoidal.ECs,
                      BM_Osteo.CAR, BM_small.pre.B., BM_Chondrocytes, BM_Endosteal.fibro., BM_Fibro.Chondro.p.,
                      BM_Stromal.fibro., BM_Arteriolar.ECs, BM_Myofibroblasts, BM_Smooth.muscle, BM_Dendritic.cells,
                      BM_NK.cells, BM_B.cell)
  
  ## dataframe
  pro_df_ <- data.frame(BM_proportions)
  rownames(pro_df_) <- z
  write.csv(pro_df_, file = paste0("./results/CARD/subgroups/",v,".csv"))
}

M1_fem_1C <- read.csv("./results/CARD/subgroups/M1_fem_1C.csv" ,row.names = 1, header= TRUE)
M1_tib_1A <- read.csv("./results/CARD/subgroups/M1_tib_1A.csv",row.names = 1, header= TRUE)
M3_fem_1C <- read.csv("./results/CARD/subgroups/M3_fem_1C.csv",row.names = 1, header= TRUE)
M3_tib_2A <- read.csv("./results/CARD/subgroups/M3_tib_2A.csv",row.names = 1, header= TRUE)

p1<-ggplot(M1_fem_1C, aes(x=rownames(M1_fem_1C),y=BM_proportions ,fill=BM_proportions)) +
  ggtitle("M1_fem_1C") +
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p2<-ggplot(M1_tib_1A, aes(x=rownames(M1_tib_1A),y=BM_proportions ,fill=BM_proportions)) +
  ggtitle("M1_tib_1A") +
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p3<-ggplot(M3_fem_1C, aes(x=rownames(M3_fem_1C),y=BM_proportions ,fill=BM_proportions)) +
  ggtitle("M3_fem_1C") +
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p4<-ggplot(M3_tib_2A, aes(x=rownames(M3_tib_2A),y=BM_proportions ,fill=BM_proportions)) +
  ggtitle("M3_tib_2A") +
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


(p1 + p2 + p3 + p4)



