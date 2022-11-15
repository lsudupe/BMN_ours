## SCRIPT: Spatial object creation visium  2 data BM project

## 20.10.22 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(ggplot2)

#Variables---------------------------------
DIR_ROOT <- file.path(getwd())
DIR_DATA <- file.path(DIR_ROOT, "data/2.data/")

#Data--------------------------------------

samples <- dir(path = DIR_DATA)
# Create each individual Seurat object for every sample
lista <- c(samples)
names(lista) <- samples
prueba <-c()

for (i in lista){
  a <- Load10X_Spatial(data.dir = paste0(DIR_DATA, i),
                          filename = "filtered_feature_bc_matrix.h5",
                          assay = "Spatial",
                          slice = i,
                          filter.matrix = TRUE)
  ######add area data
  area_a <- read.csv(paste0(DIR_DATA, i, "/new_area.csv"))
  area_a <- as.vector(area_a$new_area)
  a@meta.data["area"] <- as.factor(area_a)
  ######subset data
  Seurat::Idents(object = a) <- a@meta.data[["area"]]
  a@meta.data[["area"]] <- a@active.ident
  ######add object to list
  prueba[[length(prueba) + 1]] <- a
}
names(prueba) <- samples

###modify image data
for (i in 1:length(prueba)) {
  #https://stackoverflow.com/questions/73131436/error-in-funleft-right-non-numeric-argument-to-binary-operator-when-runni
  a <- prueba[[i]]
  b <- names(prueba[i])
  a@images[[b]]@coordinates[["tissue"]] <- as.integer(a@images[[b]]@coordinates[["tissue"]])
  a@images[[b]]@coordinates[["row"]] <- as.integer(a@images[[b]]@coordinates[["row"]])
  a@images[[b]]@coordinates[["col"]] <- as.integer(a@images[[b]]@coordinates[["col"]])
  a@images[[b]]@coordinates[["imagerow"]] <- as.integer(a@images[[b]]@coordinates[["imagerow"]])
  a@images[[b]]@coordinates[["imagecol"]] <- as.integer(a@images[[b]]@coordinates[["imagecol"]])
  saveRDS(a,file = paste0("./objects/sp/second/",names(lista[i]),".rds"))
}


#####################QC#############################

DIR_sp <- file.path(DIR_ROOT, "./objects/sp/second/")

#Data--------------------------------------
objects <- dir(path = DIR_sp)
objects <- objects[! objects %in% c("combined_filtered.rds")]

lista <- c(objects)
names(lista) <- objects
prueba <-c()

for (i in lista){
  b <- i
  a <- readRDS(paste0(DIR_sp, i))
  prueba[[length(prueba) + 1]] <- a
  names(prueba) <- i
}
names(prueba) <- samples

###
names(prueba)

M1_tib_1A <- prueba[[2]]
M1_tib_1A@meta.data[["type"]] <- "MM"
M1_tib_1A@meta.data[["orig.ident"]] <- "M1_tib_1A" 
M1_tib_1A@images[["M1_tib_1A"]]@scale.factors$lowres = M1_tib_1A@images[["M1_tib_1A"]]@scale.factors$hires

M1_fem_1C <- prueba[[1]]
M1_fem_1C@meta.data[["type"]] <- "MM"
M1_fem_1C@meta.data[["orig.ident"]] <- "M1_fem_1C" 
M1_fem_1C@images[["M1_fem_1C"]]@scale.factors$lowres = M1_fem_1C@images[["M1_fem_1C"]]@scale.factors$hires

M3_tib_2A <- prueba[[4]]
M3_tib_2A@meta.data[["type"]] <- "control"
M3_tib_2A@meta.data[["orig.ident"]] <- "M3_tib_2A" 
M3_tib_2A@images[["M3_tib_2A"]]@scale.factors$lowres = M3_tib_2A@images[["M3_tib_2A"]]@scale.factors$hires

M3_fem_1C <- prueba[[3]]
M3_fem_1C@meta.data[["type"]] <- "control"
M3_fem_1C@meta.data[["orig.ident"]] <- "M3_fem_1C" 
M3_fem_1C@images[["M3_fem_1C"]]@scale.factors$lowres = M3_fem_1C@images[["M3_fem_1C"]]@scale.factors$hires


##Merge them
combined <- merge(M1_tib_1A, y = c(M1_fem_1C, M3_tib_2A, M1_tib_1A), 
                  add.cell.ids = c("M1_tib_1A", "M1_fem_1C", "M3_tib_2A", "M3_fem_1C"), project = "BM")

combined <- subset(x =combined, idents = c("Muscle", "Bone", "BM", "GP", "AC"))


###PLOTS###
# Visualize the distribution of genes detected per spot via histogram
pdf(file.path("./results/QC/combined.second/",filename = "genes detected per spot histogram.pdf"))
combined@meta.data %>% 
  ggplot(aes(color=orig.ident, x=nFeature_Spatial, fill=orig.ident)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 100)
dev.off()

# Visualize the distribution of UMI detected per spot via histogram
pdf(file.path("./results/QC/combined.second/",filename = "UMIs detected per spot histogram.pdf"))
combined@meta.data %>% 
  ggplot(aes(color=orig.ident, x=nCount_Spatial, fill=orig.ident)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() +
  geom_vline(xintercept = 200)
dev.off()

pdf(file.path("./results/QC/combined.second",filename = "features.pdf"))
FeatureScatter(combined, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial", group.by = "orig.ident")
dev.off()

pdf(paste("./results/QC/combined.second/nFeature_Spatial.pdf",sep=""))
SpatialFeaturePlot(combined, features = "nFeature_Spatial",combine = T, pt.size.factor = 10)
dev.off()

pdf(paste("./results/QC/combined.second/nCount_Spatial.pdf",sep=""))
SpatialFeaturePlot(combined, features = "nCount_Spatial",combine = T, pt.size.factor = 10)
dev.off()

##########SPATIAL plots
library(BuenColors)
nuria <- jdb_palette("brewer_spectra")
color <- nuria
image.trans <- 0.5
spot.trans <- 1

##specific plots

a <- as.vector(combined@meta.data[["nFeature_Spatial"]])
a <- log(a)
hist(a)
combined@meta.data[["log_nFeature_Spatial"]] <- a


b <- c(5,9)
label <- c("min", "max")
p1 <- SpatialFeaturePlot(combined, features = "log_nFeature_Spatial",combine = FALSE)
fix.p1 <- scale_fill_gradientn(colours=color,breaks=b, labels = label,limits =b)
p2 <- lapply(p1, function (x) x + fix.p1)

pdf(paste("./results/QC/combined.second/log_nFeature_Spatial.pdf",sep=""))
CombinePlots(p2)
dev.off()

a <- as.vector(combined@meta.data[["nCount_Spatial"]])
a <- log(a)
combined@meta.data[["log_nCount_Spatial"]] <- a

b <- c(4,11)
label <- c("min", "max")
p1 <- SpatialFeaturePlot(combined, features = "log_nCount_Spatial",combine = FALSE, alpha=0.9)
fix.p1 <- scale_fill_gradientn(colours=color,breaks=b, labels = label,limits =b)
p2 <- lapply(p1, function (x) x + fix.p1)

pdf(paste("./results/QC/combined.second/log_nCount_Spatial.pdf",sep=""))
CombinePlots(p2)
dev.off()


# Filter out low quality cells using selected thresholds - these will change with experiment
combined <- subset(x = combined, 
                   subset= (nCount_Spatial >= 200 & nCount_Spatial <= 40000) & 
                     (nFeature_Spatial >= 100 & nFeature_Spatial <= 6000))

saveRDS(combined, "./objects/sp/second/combined_filtered.rds")



