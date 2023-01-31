## SCRIPT: QC of the spatial data BMN project

## 15.09.22 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)

#Data--------------------------------------
M1_tib_1A  <- readRDS("./objects/sp/M1_tib_1A.rds")
M1_fem_2B  <- readRDS("./objects/sp/M1_fem_2B.rds")
#M3_tib_1A  <- readRDS("./objects/sp/M3_tib_1A.rds")
M3_fem_1C  <- readRDS("./objects/sp/M3_fem_1C.rds")


M1_tib_1A@meta.data[["type"]] <- "MM"
M1_tib_1A@meta.data[["orig.ident"]] <- "M1_tib_1A" 

M1_fem_2B@meta.data[["type"]] <- "MM"
M1_fem_2B@meta.data[["orig.ident"]] <- "M1_fem_2B" 

M3_fem_1C@meta.data[["type"]] <- "control"
M3_fem_1C@meta.data[["orig.ident"]] <- "M3_fem_1C" 

M3_tib_1A@meta.data[["type"]] <- "control"
M3_tib_1A@meta.data[["orig.ident"]] <- "M3_tib_1A" 

##Merge them
combined <- merge(M1_tib_1A, y = c(M1_fem_2B, M3_fem_1C), 
                  add.cell.ids = c("M1_tib_1A", "M1_fem_2B", "M3_fem_1C"), project = "BM")

###clean
a <- combined
a@meta.data[["area"]] <- as.factor(a@meta.data[["area"]])
levels(a@meta.data[["area"]])[20] <- "NS"
a <- SetIdent(a, value = a@meta.data[["area"]])
combined <- a

saveRDS(combined, "./objects/sp/combined.rds")

###PLOTS###
######quality metrics per area
pdf(file.path("./results/QC/",filename = "Number of spot per area.pdf"))
combined@meta.data%>% 
  ggplot(aes(x=area, fill=orig.ident)) + 
  geom_bar(alpha=0.8) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NSpots")
dev.off()

pdf(file.path("./results/QC/",filename = "Number of spot per area.tissue.pdf"))
combined@meta.data%>% 
  ggplot(aes(x=orig.ident, fill=area)) + 
  geom_bar(alpha=0.8) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NSpots")
dev.off()

# Visualize the distribution of genes detected per spot via histogram
pdf(file.path("./results/QC/combined/",filename = "genes detected per spot histogram.pdf"))
combined@meta.data %>% 
  ggplot(aes(color=orig.ident, x=nFeature_Spatial, fill=orig.ident)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 100)
dev.off()

# Visualize the distribution of genes detected per spot via histogram
pdf(file.path("./results/QC/",filename = paste("genes boxplot per area.pdf",sep="")))
combined@meta.data%>% 
  ggplot(aes(x=nFeature_Spatial, y= area, fill=area)) + 
  geom_boxplot(alpha=0.8) +  
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlim(0, 1500) +
  geom_vline(xintercept = 100) +
  theme(plot.title = element_text(hjust=0.5, face="bold"))
dev.off()

# Visualize the distribution of UMI detected per spot via histogram
pdf(file.path("./results/QC/combined/",filename = "UMIs detected per spot histogram.pdf"))
combined@meta.data %>% 
  ggplot(aes(color=orig.ident, x=nCount_Spatial, fill=orig.ident)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() +
  geom_vline(xintercept = 200)
dev.off()

# Visualize the distribution of UMI detected per spot via histogram
pdf(file.path("./results/QC/",filename = paste("UMI boxplot per area.pdf",sep="")))
combined@meta.data%>% 
  ggplot(aes(x=nCount_Spatial, y= area, fill=area)) + 
  geom_boxplot(alpha=0.8) +  
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlim(0, 1500) +
  geom_vline(xintercept = 100) +
  theme(plot.title = element_text(hjust=0.5, face="bold"))
dev.off()

pdf(file.path("./results/QC/combined",filename = "features.pdf"))
FeatureScatter(combined, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial", group.by = "orig.ident")
dev.off()

##########SPATIAL plots
library(BuenColors)
nuria <- jdb_palette("brewer_spectra")
#solar <- jdb_palette("solar_extra")

color <- nuria
image.trans <- 0.5
spot.trans <- 1

##specific plots

a <- as.vector(combined@meta.data[["nFeature_Spatial"]])
a <- log(a)
hist(a)
combined@meta.data[["log_nFeature_Spatial"]] <- a
combined <- a

b <- c(2,9)
label <- c("min", "max")
p1 <- SpatialFeaturePlot(combined, features = "log_nFeature_Spatial",combine = FALSE)
fix.p1 <- scale_fill_gradientn(colours=color,breaks=b, labels = label,limits =b)
p2 <- lapply(p1, function (x) x + fix.p1)

pdf(paste("./results/QC/combined/log_nFeature_Spatial.pdf",sep=""))
CombinePlots(p2)
dev.off()

a <- as.vector(combined@meta.data[["nCount_Spatial"]])
a <- log(a)
combined@meta.data[["log_nCount_Spatial"]] <- a
hist(a)

b <- c(3,11)
label <- c("min", "max")
p1 <- SpatialFeaturePlot(combined, features = "log_nCount_Spatial",combine = FALSE, alpha=0.9)
fix.p1 <- scale_fill_gradientn(colours=color,breaks=b, labels = label,limits =b)
p2 <- lapply(p1, function (x) x + fix.p1)

pdf(paste("./results/QC/combined/log_nCount_Spatial.pdf",sep=""))
CombinePlots(p2)
dev.off()


# Filter out low quality cells using selected thresholds - these will change with experiment
combined <- subset(x = combined, idents = c("3","4","5","7","8","9","10"))
combined <- subset(x = combined, 
                            subset= (nCount_Spatial >= 200 & nCount_Spatial <= 25000) & 
                              (nFeature_Spatial >= 100 & nFeature_Spatial <= 5500))

saveRDS(combined, "./objects/sp/combined_filtered.rds")



