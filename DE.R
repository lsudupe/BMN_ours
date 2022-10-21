## SCRIPT: DE analysis cell type BM project

## 21.10.22 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(ggplot2)
library(patchwork)
library(cowplot)


#Data---------------------------------
integrated_seurat <- readRDS("./objects/sp/integrated/second/integrated.seurat_type.rds")

###Markers

x <- integrated_seurat
x <- SetIdent(x, value = x@meta.data[["seurat_clusters"]])
markers_seurat_area <- Seurat::FindAllMarkers(object = x, 
                                              assay = "integrated",
                                              slot = "data",
                                              verbose = TRUE, 
                                              only.pos = TRUE)

#Filter
markers_area <- subset(markers_seurat_area, p_val_adj < 0.05 & 0.5 < avg_log2FC)

#Top5
top10_area <- markers_area %>%
  group_by(cluster) %>%
  top_n(n = 10,
        wt = avg_log2FC)

pdf(file.path("./results/DE/second/",filename = "heatmap.pdf"))
DoHeatmap(integrated_seurat, features = top10_area$gene) + 
  theme(text = element_text(size = 6.5))
dev.off()

#Top5
top3 <- markers_area %>%
  group_by(cluster) %>%
  top_n(n = 3,
        wt = avg_log2FC)


cell_types <- as.list(as.vector(unique(top3$cluster)))
names(cell_types) <- as.vector(unique(top3$cluster))

for (i in 1:length(cell_types)){
  a <- cell_types[[i]]
  top3$cluster <- as.character(top3$cluster)
  b <- top3[top3$cluster %in% c(a), ]
  b <- as.vector(b$gene)
  cell_types[[i]] <- b
}

##################################3
#https://rpubs.com/DarrenVan/628853
## remove the x-aNxis text and tick 
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj, feature, pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}

## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1), 
          axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
} 

#########################

for (i in 1:length(cell_types)){
  a <- cell_types[[i]]
  c <- paste0(names(cell_types[i]), "cluster")
  
  #pdf(file.path("./results/DE/second/",filename = paste0("violin" ,c,".pdf")))
  #print(VlnPlot(integrated_seurat, features = a, group.by = "seurat_clusters", split.by = "type",ncol = 1))
  #dev.off()
  
  pdf(file.path("./results/DE/second/",filename = paste0("2violin" ,c,".pdf")))
  print(StackedVlnPlot(obj = integrated_seurat, features = a, split.by = "type"))
  dev.off()
  #spatial 
  pdf(file.path("./results/DE/second/",filename = paste0("spatial_seurat_",c,".pdf")))
  print(SpatialFeaturePlot(x,features=a,combine = FALSE))
  dev.off()
}

StackedVlnPlot(obj = pbmc, features = features)





