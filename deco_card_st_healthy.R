## SCRIPT: Deconvolution using CARD areas subgroups BM project ST object ONLY HEALTHY

## 28.03.23 Laura Sudupe , git @lsudupe
#https://github.com/YingMa0107/CARD/

#Libraries---------------------------------
library(Seurat)
library(ggplot2)
library(dplyr)
library(CARD)
library(gtools)
library(scatterpie)
library(ggcorrplot)
library(STutility)

source(file = "./card.plot2.R")


#Data---------------------------------
se <- readRDS("./objects/sp/integrated/se.rds")
single_cell_bonemarrow <- readRDS("./objects/heterogeneity/single_cell_bonemarrow.rds")
single_cell_bonemarrow@meta.data[["orig.ident"]] <- "ref"
#single_cell_bonemarrow_all <- readRDS("./objects/heterogeneity/single_cell_bonemarrow_all_groups.rds")
#single_cell_bonemarrow_all@meta.data[["orig.ident"]] <- "ref"



######Prepare data
#single cell
sub_list <- levels(single_cell_bonemarrow@meta.data[["ident"]])
#sub_list <- levels(single_cell_bonemarrow_all@meta.data[["cell"]])


single_cell_bonemarrow_counts <- single_cell_bonemarrow@assays[["RNA"]]@counts
single_cell_bonemarrow_meta <- single_cell_bonemarrow@meta.data
#single_cell_bonemarrow_meta <- single_cell_bonemarrow_meta[, 6:7]


#spatial
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

objects <- objects[3:4]


## card
for (i in 1:length(objects)){
  a <- objects[[i]]
  a_count <- a@assays[["SCT"]]@data
  image.info <- a@tools[["Staffli"]]@meta.data
  a_location <- image.info[,1:2]

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
  v<- ("jojij")
  CARD_a_proportions <- CARD_a@Proportion_CARD
  
  ###################Plots
  p1 <- CARD.visualize.pie(proportion = CARD_a@Proportion_CARD,
                           spatial_location = CARD_a@spatial_location)
  pdf(paste("./results/ST/card/healthy/all/", names(objects[i]),"_1.pdf",sep=""))
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
  
  pdf(paste("./results/ST/card/healthy/all/", names(objects[i]),"_2.pdf",sep=""))
  print(p2)
  dev.off()
  ## correlation
  p3 <- CARD.visualize.Cor(CARD_a@Proportion_CARD,colors = NULL) # if not provide, we will use the default colors
  pdf(paste("./results/ST/card/healthy/all/", names(objects[i]),"_3.pdf",sep=""))
  print(p3)
  dev.off()
  ## save object
  saveRDS(a,file = paste0("./objects/card/last/healthy/all/",names(objects[i]),"_subgroup_ST.rds"))
  ## save card results
  saveRDS(CARD_a,file = paste0("./objects/card/last/healthy/all/",names(objects[i]),"_CARD_obj_subgroup_ST.rds"))
}

## read card objects
card <- c()
DIR_ROOT <- file.path(getwd())  
DIR_DATA <- file.path(DIR_ROOT, "/objects/card/last/healthy/all/")
o <- c()

########ENRICHMENT PLOT IF INTERESTED
for (i in 1:length(objects)){
  a <- objects[[i]]
  b <- names(objects[i])
  # read data
  c <- readRDS(paste0("./objects/card/last/healthy/all/",names(objects[i]), "_CARD_obj_subgroup_ST.rds"))
  pro <- as.data.frame(c@Proportion_CARD)
  mm_ic <- pro$MM_MIC
  a@meta.data[["mm_ic"]] <- mm_ic
  o[[length(o) + 1]] <- a
  
  pdf(paste0("./results/ST/card/enrich/healthy/ll/",names(objects[i]),"_enrich_mm.pdf"))
  print(FeatureOverlay(a, features = "mm_ic",
                 cols = c("lightgray", "mistyrose", "red", "dark red", "black"), ncol = 1, pt.size = 1.4))
  dev.off()
}
########ENRICHMENT PLOT FIN

#######EXTRACT PROPORTIONS

o <- c()
for (i in 1:length(objects)){
  a <- objects[[i]]
  b <- names(objects[i])
  ## read data
  c <- readRDS(paste0("./objects/card/last/healthy/all/",names(objects[i]), "_CARD_obj_subgroup_ST.rds"))
  pro <- as.data.frame(c@Proportion_CARD)
  cell.types <- colnames(pro)
  ## plot to add proportions
    for (i in 1:length(cell.types)){
        d <- cell.types[i]
        e <- pro[[d]]
        a@meta.data[[d]] <- e
    }
  
  o[[length(o) + 1]] <- a
  #pdf(paste0("./results/ST/card/enrich/",names(objects[i]),"_enrich_mm.pdf"))
  #print(FeatureOverlay(a, features = "mm_ic",
                       #cols = c("lightgray", "mistyrose", "red", "dark red", "black"), ncol = 1, pt.size = 1.4))
  #dev.off()
}
#######

name <- c("M3_F_1C", "M3_fem_1C")
names(o) <- name
## separate list
list2env(o,envir=.GlobalEnv)

se_merge <- MergeSTData(M3_F_1C, y = c(M3_fem_1C), 
                         add.spot.ids = c("M3_F_1C", "M3_fem_1C"), project = "BM")

####Proportion analysis
pro_meta <- se_merge@meta.data
pro_meta_ <- pro_meta[12:17]

library(tidyverse)
df_long <- pro_meta_ %>%
  gather(key = "cell_type", value = "value")

pdf(file.path("./results/ST/card/healthy/all/",filename = "violin_cell_types_all_onlyhealthy_horizontal.pdf"))
ggplot(df_long, aes(x = cell_type, y = value, fill = cell_type)) +
  geom_violin(scale = "width", trim = FALSE, show.legend = FALSE) +
  stat_summary(fun.data = mean_sdl, color = "red", geom = "pointrange", position = position_dodge(0.9)) +
  stat_summary(fun = median, color = "blue", geom = "point", position = position_dodge(0.9)) +
  geom_boxplot(width = 0.1, outlier.color = "black", outlier.shape = 16, outlier.size = 1) +
  theme_minimal() +
  labs(x = "Cell Type", y = "Value", title = "Violin Plot of Cell Types")
dev.off()

pdf(file.path("./results/ST/card/healthy/all/",filename = "violin_cell_types_all_onlyhealthy.pdf"))
ggplot(df_long, aes(x = cell_type, y = value, fill = cell_type)) +
  geom_violin(scale = "width", trim = FALSE, show.legend = FALSE) +
  coord_flip() +
  stat_summary(fun.data = mean_sdl, color = "red", geom = "pointrange", position = position_dodge(0.9)) +
  stat_summary(fun = median, color = "blue", geom = "point", position = position_dodge(0.9)) +
  geom_boxplot(width = 0.1, outlier.color = "black", outlier.shape = 16, outlier.size = 1) +
  theme_minimal() +
  labs(y = "Cell Type", x = "Value", title = "Violin Plot of Cell Types")
dev.off()

####Proportion analysis FIN

###porcentage plots
se_subset <- SubsetSTData(se, idents = c("M3_F_1C","M3_fem_1C"))


se_subset@meta.data[["Tcell"]] <- se_merge@meta.data[["Tcell"]]
se_subset@meta.data[["Bcell"]] <- se_merge@meta.data[["Bcell"]]
se_subset@meta.data[["Erythroblasts"]] <- se_merge@meta.data[["Erythroblasts"]]
se_subset@meta.data[["Monocytes"]] <- se_merge@meta.data[["Monocytes"]]
se_subset@meta.data[["Neutrophils"]] <- se_merge@meta.data[["Neutrophils"]]
#se_subset@meta.data[["MSC"]] <- se_merge@meta.data[["MSC"]]
#se_subset@meta.data[["EC"]] <- se_merge@meta.data[["EC"]]
se_subset@meta.data[["DC"]] <- se_merge@meta.data[["DC"]]
#se_subset@meta.data[["NK"]] <- se_merge@meta.data[["NK"]]

saveRDS(se_subset, "./objects/sc/integrated/se_deco_healthy.rds")
se_subset <- readRDS("./objects/sc/integrated/se_deco_healthy.rds")

###pruebas plots

pdf(file.path("./results/ST/mm_ic_merge.pdf"))
FeatureOverlay(se, features = "MM_MIC",
               cols = c("lightgray", "mistyrose", "red", "dark red", "black"), 
               sampleids = 1:6, ncols = 2,pt.size = 0.7)
dev.off()

pdf(file.path("./results/ST/UMAP.pdf"))
DimPlot(se, group.by = "seurat_clusters", label = TRUE, label.size = 8, reduction = "umap")
dev.off()

p1 <- DimPlot(se, group.by = "name", reduction = "umap")
p2 <- DimPlot(se, group.by = "seurat_clusters", label = TRUE, label.size = 8, reduction = "umap")
p3 <- DimPlot(se, group.by = "condition", label = TRUE, label.size = 8, reduction = "umap")
pdf(file.path("./results/ST/UMAP.pdf"))
p3
dev.off()

pdf(file.path("./results/ST/mm_ic_merge_prueba.pdf"))
FeatureOverlay(se, features = "labels", 
               min.cutoff = 2,
               max.cutoff = 3,
               sampleids = 1:6, ncols = 2,pt.size = 0.7)
dev.off()









