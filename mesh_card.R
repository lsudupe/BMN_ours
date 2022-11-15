## SCRIPT: Deconvolution using CARD with spatial grid areas subgroups BM project

## 14.11.22 Laura Sudupe , git @lsudupe
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

## Spatial object
M1_fem_1C <- readRDS("./objects/card/M1_fem_1C_subgroup.rds")
M1_tib_1A <- readRDS("./objects/card/M1_tib_1A_subgroup.rds")
M3_fem_1C <- readRDS("./objects/card/M3_fem_1C_subgroup.rds")
M3_tib_2A <- readRDS("./objects/card/M3_tib_2A_subgroup.rds")

## Card object
M1_fem_1C_CARD_obj <- readRDS("./objects/card/M1_fem_1C_CARD_obj_subgroup.rds")
M1_tib_1A_CARD_obj <- readRDS("./objects/card/M1_tib_1A_CARD_obj_subgroup.rds")
M3_fem_1C_CARD_obj <- readRDS("./objects/card/M3_fem_1C_CARD_obj_subgroup.rds")
M3_tib_2A_CARD_obj <- readRDS("./objects/card/M3_tib_2A_CARD_obj_subgroup.rds")

CARD_obj = CARD.imputation(M1_fem_1C_CARD_obj,NumGrids = 375,ineibor = 10,exclude = NULL)

## Visualize the newly grided spatial locations to see if the shape is correctly detected. If not, the user can provide the row names of the excluded spatial location data into the CARD.imputation function
location_imputation = cbind.data.frame(x=as.numeric(sapply(strsplit(rownames(CARD_obj@refined_prop),split="x"),"[",1)),
                                       y=as.numeric(sapply(strsplit(rownames(CARD_obj@refined_prop),split="x"),"[",2)))
rownames(location_imputation) = rownames(CARD_obj@refined_prop)
library(ggplot2)
p4 <- ggplot(location_imputation, 
             aes(x = x, y = y)) + geom_point(shape=22,color = "#7dc7f5")+
  theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        legend.position="bottom",
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(colour = "grey89", fill=NA, size=0.5))
print(p4)

sub_list <- as.vector(levels(single_cell_bonemarrow@meta.data[["cell"]]))


p5 <- CARD.visualize.prop.2(
  proportion = CARD_obj@refined_prop,                         
  spatial_location = location_imputation,            
  ct.visualize = sub_list,                    
  colors = c("lightblue","lightyellow","red"),    
  NumCols = 7)
pdf(file.path("./results/CARD/mesh/sub_groups/",filename = ("Smesh.pdf")))
print(p5)
dev.off()

#####extract refined proportions
#refined_prop <- CARD_obj@refined_prop
#img <- GetImage(M3_tib_2A)
#coords <- GetTissueCoordinates(object = M3_tib_2A[[img]])
vv <- as.data.frame(CARD_obj@algorithm_matrix[["Res"]][["V"]])
vvv <- colnames(vv)
M1_fem_1C@meta.data[["Sinusoidal_ECs"]] <- vv[,20]

pdf(file.path("./results/CARD/mesh/sub_groups/",filename = ("Sinusoidal_ECs.pdf")))
print(SpatialFeaturePlot(M1_fem_1C, features = c("Sinusoidal_ECs"), pt.size.factor = 7, alpha = 0.6))
dev.off()


p1 <- SpatialFeaturePlot(fibro.new, features = i, combine = FALSE)
fix.p1 <- scale_fill_gradientn(colours=c(bl,"white", re), na.value = "grey98",limits = b,breaks=b, labels = c("min", "max"))
#fix.p1 <- scale_fill_gradientn(colours=myPalette(100),breaks=b, labels = c("min", "max"),limits = b)
p2 <- lapply(p1, function (x) x + fix.p1)

pdf(file.path("./results/Aucell/B_genes/",filename = paste(i,"spatial.b.pdf",sep="")))
print(CombinePlots(p2))
dev.off()


