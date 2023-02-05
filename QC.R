## SCRIPT: QC of the seurat spatial objects BMN project

## 24.06.22 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)

## Read data
combined <- readRDS("./objects/sp/combined.rds")

## subset data, take out bone
Seurat::Idents(object = combined) <- combined@meta.data[["area"]]
combined <- subset(x = combined, idents = c("bone_marrow"))

##Visualization
# Visualize the number of spots counts per sample
pdf(file.path("./results/QC",filename = "Number of spot per sample.pdf"))
combined@meta.data%>% 
  ggplot(aes(x=orig.ident, fill=orig.ident)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NSpots")
dev.off()

# Visualize the distribution of genes detected per spot via boxplot
pdf(file.path("./results/QC",filename = "genes detected per spot boxplot.pdf"))
combined@meta.data %>% 
  ggplot(aes(x=orig.ident, y=log10(nFeature_Spatial), fill=orig.ident)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Ngenes vs Npots")
dev.off()

# Visualize the distribution of genes detected per spot via histogram
pdf(file.path("./results/QC",filename = "genes detected per spot histogram.pdf"))
combined@meta.data %>% 
  ggplot(aes(color=orig.ident, x=nFeature_Spatial, fill=orig.ident)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)
dev.off()

# Visualize the number UMIs/transcripts per cell
pdf(file.path("./results/QC",filename = "number UMIs sati transcripts per spot.pdf"))
combined@meta.data %>% 
  ggplot(aes(color=orig.ident, x=nCount_Spatial, fill=orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)
dev.off()


## spatial plots
library(BuenColors)
color <- jdb_palette("brewer_spectra")

feature.list <- c("nCount_Spatial", "nFeature_Spatial")
sample <- c(unique(combined$orig.ident))
samples <- c()
Seurat::Idents(object = combined) <- combined@meta.data[["orig.ident"]]
images <- combined@images

for (i in sample){
combined@images[[i]] <- NULL
}

for (i in sample){
  a <- subset(x = combined, idents = i)
  a@images[[i]] <- images[[i]]
  samples[[length(samples) + 1]] <- a
}
names(samples) <- sample

## separate list
list2env(prueba,envir=.GlobalEnv)

for (i in 1:length(samples)){
  a <- samples[[i]]
  b <- c(min(a@meta.data[["nFeature_Spatial"]]), max(a@meta.data[["nFeature_Spatial"]]))
  p1 <- SpatialFeaturePlot(a, features = c("nFeature_Spatial"), combine = FALSE, ncol = 1, pt.size.factor = 7)
  fix.p1 <- scale_fill_gradientn(colours=color,
                                 breaks=b,
                                 labels=c("Min","Max"),
                                 limits = b)
  p2 <- lapply(p1, function (x) x + fix.p1)
  
  pdf(paste("./results/QC/unique/",names(samples[i]),"_feature_spatial.pdf",sep=""))
  print(CombinePlots(p2))
  dev.off()
}

for (i in 1:length(samples)){
  a <- samples[[i]]
  b <- c(min(a@meta.data[["nCount_Spatial"]]), max(a@meta.data[["nCount_Spatial"]]))
  p1 <- SpatialFeaturePlot(a, features = c("nCount_Spatial"), combine = FALSE, ncol = 1, pt.size.factor = 7)
  fix.p1 <- scale_fill_gradientn(colours=color,
                                 breaks=b,
                                 labels=c("Min","Max"),
                                 limits = b)
  p2 <- lapply(p1, function (x) x + fix.p1)
  
  pdf(paste("./results/QC/unique/",names(samples[i]),"_count_spatial.pdf",sep=""))
  print(CombinePlots(p2))
  dev.off()
}

combined.meta <- combined@meta.data

###########count/feature correlation

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
pdf(paste("./results/QC/count_features_cor.pdf",sep=""))
combined.meta %>% 
  ggplot(aes(x=nCount_Spatial, y=nFeature_Spatial)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~orig.ident)
dev.off()

pdf(paste("./results/QC/features_violinplot.pdf",sep=""))
VlnPlot(object = combined, features = 'nFeature_Spatial', split.by = 'orig.ident')
dev.off()

pdf(paste("./results/QC/count_violinplot.pdf",sep=""))
VlnPlot(object = combined, features = 'nCount_Spatial', split.by = 'orig.ident')
dev.off()


## add image
combined@images <- images

# Filter out low quality cells using selected thresholds - these will change with experiment
filtered_combined <- subset(x = combined, 
                          subset= (nCount_Spatial >= 200 & nCount_Spatial <= 100000) & 
                            (nFeature_Spatial >= 150 & nFeature_Spatial <= 10000))

#####save combined
saveRDS(filtered_combined,"./objects/sp/combined_filtered.rds")
