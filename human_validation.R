## SCRIPT: Human sample utility evaluation

## 27.02.23 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(ggplot2)

# Read data, mouse sample/human sample
#Variables---------------------------------
DIR_ROOT <- file.path(getwd())
DIR_DATA <- file.path(DIR_ROOT, "data/data/")

# samples
samples <- dir(path = DIR_DATA)
samples <- samples[ samples %in% c("M1_fem_1C","BM_human_AP-B08805")]

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
  prueba[[length(prueba) + 1]] <- a
}
names(prueba) <- samples

# separate samples
list2env(prueba,envir=.GlobalEnv)

# add area data to mouse and only get bone marrow area
## add area data
area_a <- read.csv(paste0(DIR_DATA, "M1_fem_1C/area.csv"))
area_a <- as.vector(area_a$area)
M1_fem_1C@meta.data["area"] <- as.factor(area_a)

# subset data
Seurat::Idents(object = M1_fem_1C) <- M1_fem_1C@meta.data[["area"]]
M1_fem_1C <- subset(x =M1_fem_1C, idents = c("bone_marrow"))


# add id
`BM_human_AP-B08805`@meta.data[["orig.ident"]] <- "human"
human <- `BM_human_AP-B08805`
M1_fem_1C@meta.data[["orig.ident"]] <- "mouse"

lista <- c(human, M1_fem_1C)
names(lista) <- c("human","M1_fem_1C")
prueba <-c()

for (i in 1:length(lista)){
  a <- lista[[i]]
  # Normalize the counts
  a <- NormalizeData(a)
  
  # agregate the expression in both
  a <- AggregateExpression(a, 
                             group.by = c("orig.ident"),
                             assays = 'Spatial',
                             slot = "data",
                             return.seurat = FALSE)
  
  # add to the list
  prueba[[length(prueba) + 1]] <- a
}
names(prueba) <- samples


# separate matrixes
cts <- prueba[["BM_human_AP-B08805"]]$Spatial
human <- as.data.frame(cts)
cts <- prueba[["M1_fem_1C"]]$Spatial
mouse <- as.data.frame(cts)

# order genes
mouse$all <- mouse[order(mouse$all,decreasing=TRUE), ]
human$all <- human[order(human$all,decreasing=TRUE), ]

# extract 100 top
top_mouse <- tolower(rownames(mouse)[1:500])
top_human <- toupper(rownames(human)[1:500])

write.table(top_human, file = "./human.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE,
            quote = FALSE)

# Read in the data
mouse_ortho <- scan("./results/biomart/mart_export.txt", what="", sep="\n")
mouse_ortho <- tolower(as.vector(t(as.matrix(mouse_ortho[-1]))))
a <- intersect(mouse_ortho, top_mouse)






