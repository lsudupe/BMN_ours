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


CARD_obj = CARD.imputation(M3_tib_2A_CARD_obj,NumGrids = 2000,ineibor = 10,exclude = NULL)

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

sub_list <- as.vector(levels(single_cell_bonemarrow@meta.data[["groups"]]))
 
p5 <- CARD.visualize.prop.2(
  proportion = CARD_obj@refined_prop,                         
  spatial_location = location_imputation,            
  ct.visualize = sub_list,                    
  colors = c("lightblue","lightyellow","red"),    
  NumCols = 7)          



print(p5)
