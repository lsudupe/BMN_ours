## SCRIPT: Deconvolution using CARD areas subgroups BM project ST object

## 19.02.23 Laura Sudupe , git @lsudupe
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
library(STutility)

source(file = "./card.plot2.R")


#Data---------------------------------
spatial <- readRDS("./objects/sp/integrated/integrated.harmony.rds")
se <- readRDS("./objects/sp/integrated/se.rds")
single_cell_bonemarrow <- readRDS("./objects/sc/integrated/integrated_sc_harmony.rds")
single_cell_bonemarrow <- readRDS("./objects/sc/integrated/single_cell_bonemarrow_all_groups_harmony.rds")

######Prepare data
#single cell
sub_list <- levels(single_cell_bonemarrow@meta.data[["split"]])

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

## separate list
list2env(objects,envir=.GlobalEnv)


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
    ct.varname = "split",
    ct.select = unique(single_cell_bonemarrow_meta$split),
    sample.varname = "orig.ident",
    minCountGene = 100,
    minCountSpot = 5)
  
  ###################Deco
  CARD_a = CARD_deconvolution(CARD_object = CARD_obj_a)
  CARD_a_proportions <- CARD_a@Proportion_CARD
  v<- ("jojoj")
  ###################Plots
  #p1 <- CARD.visualize.pie(proportion = CARD_a@Proportion_CARD,
                           #spatial_location = CARD_a@spatial_location)
  #pdf(paste("./results/ST/card/all/", names(objects[i]),"_1.pdf",sep=""))
  #print(p1)
  #dev.off()
  
  ## select the cell type that we are interested
  ct.visualize = sub_list
  ## visualize the spatial distribution of the cell type proportion
  p2 <- CARD.visualize.prop.2(
    proportion = CARD_a@Proportion_CARD,        
    spatial_location = CARD_a@spatial_location, 
    ct.visualize = ct.visualize,                 ### selected cell types to visualize
    colors = c("lightblue","lightyellow","red"), ### if not provide, we will use the default colors
    NumCols = 4)                                 ### number of columns in the figure panel
  
  pdf(paste("./results/ST/card/all/", names(objects[i]),"_2.pdf",sep=""))
  print(p2)
  dev.off()
  ## correlation
  p3 <- CARD.visualize.Cor(CARD_a@Proportion_CARD,colors = NULL) # if not provide, we will use the default colors
  pdf(paste("./results/ST/card/all/", names(objects[i]),"_3.pdf",sep=""))
  print(p3)
  dev.off()
  ## save object
  saveRDS(a,file = paste0("./objects/card/last/all/",names(objects[i]),"_subgroup_ST.rds"))
  ## save card results
  saveRDS(CARD_a,file = paste0("./objects/card/last/all/",names(objects[i]),"_CARD_obj_subgroup_ST.rds"))
}

## read card objects
card <- c()
DIR_ROOT <- file.path(getwd())  
DIR_DATA <- file.path(DIR_ROOT, "/objects/card/last/all/")
o <- c()


for (i in 1:length(objects)){
  a <- objects[[i]]
  b <- names(objects[i])
  # read data
  c <- readRDS(paste0("./objects/card/last/all/",names(objects[i]), "_CARD_obj_subgroup_ST.rds"))
  pro <- as.data.frame(c@Proportion_CARD)
  mm_ic <- pro$MM_MIC
  a@meta.data[["mm_ic"]] <- mm_ic
  o[[length(o) + 1]] <- a
  
  pdf(paste0("./results/ST/card/enrich/all/",names(objects[i]),"_enrich_mm.pdf"))
  print(FeatureOverlay(a, features = "mm_ic",
                 cols = c("lightgray", "mistyrose", "red", "dark red", "black"), ncol = 1, pt.size = 1.4))
  dev.off()
}


#######EXTRACT PROPORTIONS

cell.types <- colnames(pro)
o <- c()
for (i in 1:length(objects)){
  a <- objects[[i]]
  b <- names(objects[i])
  ## read data
  c <- readRDS(paste0("./objects/card/last/all/",names(objects[i]), "_CARD_obj_subgroup_ST.rds"))
  pro <- as.data.frame(c@Proportion_CARD)
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

names(o) <- name
## separate list
list2env(o,envir=.GlobalEnv)

se_merge <- MergeSTData(M1_fem_1C, y = c(M2_F_2B, M3_F_1C,M3_fem_1C ,M8_F2_1C , M9_F2_1C), 
                         add.spot.ids = c("M1_fem_1C", "M2_F_2B", "M3_F_1C", "M3_fem_1C", "M8_F2_1C", "M9_F2_1C"), project = "BM")

####Proportion analysis
pro_meta <- se_merge@meta.data
pro_meta_ <- pro_meta[12:44]

library(tidyverse)
df_long <- pro_meta_ %>%
  gather(key = "cell_type", value = "value")

pdf(file.path("./results/ST/card/all/",filename = "violin_cell_types_all_onlyMM.pdf"))
ggplot(df_long, aes(x = cell_type, y = value, fill = cell_type)) +
  geom_violin(scale = "width", trim = FALSE, show.legend = FALSE) +
  coord_flip() +
  stat_summary(fun.data = mean_sdl, color = "red", geom = "pointrange", position = position_dodge(0.9)) +
  stat_summary(fun = median, color = "blue", geom = "point", position = position_dodge(0.9)) +
  geom_boxplot(width = 0.1, outlier.color = "black", outlier.shape = 16, outlier.size = 1) +
  theme_minimal() +
  labs(y = "Cell Type", x = "Value", title = "Violin Plot of Cell Types")
dev.off()


###porcentage plots
mm_porcentge <- (sum(pro$MM_MIC))*100/nrow(pro)
neutro_porcentge <- (sum(pro$Neutrophils))*100/nrow(pro)

meta_pro <- se_merge@meta.data

# Change violin plot line colors by groups
p1<-ggplot(meta_pro, aes(x=meta_pro$name, y=meta_pro$MM_MIC, color=meta_pro$condition)) +
  geom_violin(trim=FALSE) +
  labs(title="plasma cells % in different samples",x="samples", y = "porcentage")

p2<-ggplot(meta_pro, aes(x=meta_pro$name, y=meta_pro$MM_MIC, color=meta_pro$condition)) +
  geom_violin(trim=FALSE) +
  labs(title="MM_MIC % in different clusters",x="clusters", y = "porcentage")
p2<-ggplot(meta_pro, aes(x=meta_pro$name, y=meta_pro$Bcell, color=meta_pro$condition)) +
  geom_violin(trim=FALSE) +
  labs(title="Bcell % in different clusters",x="clusters", y = "porcentage")
p2<-ggplot(meta_pro, aes(x=meta_pro$name, y=meta_pro$Erythroblasts, color=meta_pro$condition)) +
  geom_violin(trim=FALSE) +
  labs(title="Erythroblasts % in different clusters",x="clusters", y = "porcentage")
p2<-ggplot(meta_pro, aes(x=meta_pro$name, y=meta_pro$Monocytes, color=meta_pro$condition)) +
  geom_violin(trim=FALSE) +
  labs(title="Monocytes % in different clusters",x="clusters", y = "porcentage")
p2<-ggplot(meta_pro, aes(x=meta_pro$name, y=meta_pro$Tcell, color=meta_pro$condition)) +
  geom_violin(trim=FALSE) +
  labs(title="Tcell % in different clusters",x="clusters", y = "porcentage")
p2<-ggplot(meta_pro, aes(x=meta_pro$name, y=meta_pro$Neutrophils, color=meta_pro$condition)) +
  geom_violin(trim=FALSE) +
  labs(title="Neutrophils % in different clusters",x="clusters", y = "porcentage")
p2<-ggplot(meta_pro, aes(x=meta_pro$name, y=meta_pro$MSC, color=meta_pro$condition)) +
  geom_violin(trim=FALSE) +
  labs(title="MSC % in different clusters",x="clusters", y = "porcentage")
p2<-ggplot(meta_pro, aes(x=meta_pro$name, y=meta_pro$EC, color=meta_pro$condition)) +
  geom_violin(trim=FALSE) +
  labs(title="EC % in different clusters",x="clusters", y = "porcentage")
p2<-ggplot(meta_pro, aes(x=meta_pro$name, y=meta_pro$DC, color=meta_pro$condition)) +
  geom_violin(trim=FALSE) +
  labs(title="DC % in different clusters",x="clusters", y = "porcentage")
p2<-ggplot(meta_pro, aes(x=meta_pro$name, y=meta_pro$NK, color=meta_pro$condition)) +
  geom_violin(trim=FALSE) +
  labs(title="NK % in different clusters",x="clusters", y = "porcentage")
p2

##


se@meta.data[["Tcell"]] <- se_merge@meta.data[["Tcell"]]
se@meta.data[["Bcell"]] <- se_merge@meta.data[["Bcell"]]
se@meta.data[["MM_MIC"]] <- se_merge@meta.data[["MM_MIC"]]
se@meta.data[["Erythroblasts"]] <- se_merge@meta.data[["Erythroblasts"]]
se@meta.data[["Monocytes"]] <- se_merge@meta.data[["Monocytes"]]
se@meta.data[["Neutrophils"]] <- se_merge@meta.data[["Neutrophils"]]
se@meta.data[["MSC"]] <- se_merge@meta.data[["MSC"]]
se@meta.data[["EC"]] <- se_merge@meta.data[["EC"]]
se@meta.data[["DC"]] <- se_merge@meta.data[["DC"]]
se@meta.data[["NK"]] <- se_merge@meta.data[["NK"]]


se <- readRDS("./objects/sc/integrated/se_deco.rds")
saveRDS(se, "./objects/sc/integrated/se_deco.rds")



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


###neighbours
library(data.table)
se_copy <- copy(se)
se_copy <- SetIdent(se_copy, value = "seurat_clusters")
se_copy <- RegionNeighbours(se, id = "1", keep.within.id = T,verbose = TRUE)
neigh_1 <- as.factor(se_copy@meta.data[["nbs_1"]])
se_copy <- RegionNeighbours(se, id = "2", keep.within.id = T,verbose = TRUE)
neigh_2 <- as.factor(se_copy@meta.data[["nbs_2"]])

se_copy@meta.data[["neigh_1"]] <- neigh_1
se_copy@meta.data[["neigh_2"]] <- neigh_2


pdf(file.path("./results/ST/neigh_2.pdf"))
FeatureOverlay(se_copy, features = "neigh_2", ncols = 2, sampleids = 1:6, cols = c("red", "lightgray"), pt.size = 0.7)
dev.off()

meta_pro <- se_copy@meta.data

# Change violin plot line colors by groups
ggplot(meta_pro, aes(x=meta_pro$name, y=meta_pro$MM_MIC, color=meta_pro$neigh_1)) +
  geom_violin(trim=FALSE) +
  labs(title="plasma cells % in neigh_1",x="samples", y = "porcentage")
ggplot(meta_pro, aes(x=meta_pro$name, y=meta_pro$MM_MIC, color=meta_pro$neigh_2)) +
  geom_violin(trim=FALSE) +
  labs(title="plasma cells % in neigh_2",x="samples", y = "porcentage")


##DE
library(magrittr)
library(dplyr)

se <- SetIdent(se_copy, value = "nbs_2")
nbs_2.markers <- FindMarkers(se, ident.1 = "2", ident.2 = "nbs_2")
nbs_2.markers <- FindMarkers(se, ident.1 = "NA", ident.2 = "nbs_2")
nbs_2.markers$gene <- rownames(nbs_2.markers)
se.subset <- SubsetSTData(se, expression = nbs_2 %in% c("2", "nbs_2"))
sorted.marks <- nbs_2.markers %>% arrange(-avg_log2FC) %>% top_n(n = 40, wt = abs(avg_log2FC))

pdf(file.path("./results/ST/doheatmap.pdf"))
DoHeatmap(se.subset, features = sorted.marks$gene, group.colors = c("red", "lightgray"), disp.min = -2, disp.max = 2)
dev.off()

pdf(file.path("./results/ST/Bpgm.pdf"))
FeatureOverlay(se.subset, features = c("Bpgm", "Retnlg", "Cd177", "Csf3r"), pt.size = 0.7,  
               sampleids = 1:2,
               ncols = 3, cols = c("darkblue", "cyan", "yellow", "red", "darkred"))
dev.off()

## with pseudobulk approach
Seurat::Idents(object = se) <- se@meta.data[["neigh_2"]]

a <- AggregateExpression(se, 
                         group.by = c("neigh_2", "name" ),
                         assays = 'SCT',
                         slot = "data",
                         return.seurat = FALSE)

cts <- a$SCT
cts <- as.data.frame(cts)






