## SCRIPT: Human sample Healthy cleanning and QC

## 08.11.23 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(ggplot2)
library(Signac)
library(RColorBrewer)
library(STutility)
library(UCell)
library(readr)
library(dplyr)
library(gridExtra)
library(tibble)
library(readxl)
color <- rev(brewer.pal(11,"Spectral"))

## Read the data
se <- readRDS("./objects/sp/integrated/se_human_healthy.rds")

##divide by sample
Idents(object = se) <- "name"
name <- unique(se@meta.data[["name"]])
objects <- c()

for (i in name){
  a <- SubsetSTData(se, idents = i)
  a <- SCTransform(a)
  a <- NormalizeData(a)
  a <- ScaleData(a)
  ## add object to list
  objects[[length(objects) + 1]] <- a
  
}
names(objects) <- name
a <- objects
objects <- a
## separate list
#list2env(objects,envir=.GlobalEnv)

process_seurat_object <- function(seurat_obj, sample_name) {
  # Read the area.csv file
  area_file_path <- paste0("./data/data/", sample_name, "/area.csv")
  area_df <- read_csv(area_file_path, show_col_types = FALSE)
  
  # Replace NA in the 'area' column with 'bm'
  area_df <- area_df %>% mutate(area = ifelse(is.na(area), "bm", area))
  
  # Extract metadata and convert rownames to a column
  metadata_df <- seurat_obj@meta.data %>%
    rownames_to_column(var = "index")
  
  # Extract the sample number from the 'index' column
  metadata_df$sample_number <- as.numeric(sub(".*_", "", metadata_df$index))
  
  # Append the extracted sample number to the 'index' column in area_df
  area_df$index <- paste0(area_df$Barcode, "_", metadata_df$sample_number[1])
  
  # Merge the metadata with area data based on index
  updated_metadata_df <- metadata_df %>%
    left_join(area_df, by = "index")
  
  # Set the updated metadata back to the Seurat object
  rownames(updated_metadata_df) <- updated_metadata_df$index
  seurat_obj@meta.data <- updated_metadata_df
  
  return(seurat_obj)
}

# Apply the function to each Seurat object
for (sample_name in names(objects)) {
  objects[[sample_name]] <- process_seurat_object(objects[[sample_name]], sample_name)
}

for (i in seq_along(objects)) {
  seurat_obj <- objects[[i]]
  seurat_obj <- subset(seurat_obj, subset = area == "bm")
  # Update the object in the list
  objects[[i]] <- seurat_obj
}

## separate list
list2env(objects,envir=.GlobalEnv)

## Plots
a <- BM_human_H13
v <- "nCount_RNA"
c <- c(min(a@meta.data[v]), max(a@meta.data[v]))
FeatureOverlay(a, features = v, ncols = 1, pt.size = 1.7, slot = "count",
               value.scale = "all" ,cols = color) +
  scale_fill_gradientn(colours = color,
                       breaks = c,
                       labels = c,
                       limits = c)

## Combine meta
combined_metadata <- data.frame()

# Loop through each object and extract metadata
for (sample_name in names(objects)) {
  metadata <- objects[[sample_name]]@meta.data
  # Add a column to indicate the sample name
  metadata$Sample <- sample_name
  # Combine with the existing metadata
  combined_metadata <- rbind(combined_metadata, metadata)
}
# Remove 'index' and 'Barcode' columns if they exist
combined_metadata$index <- NULL
combined_metadata$Barcode <- NULL

# Violin plot for nCount_RNA
ggplot(combined_metadata, aes(x = Sample, y = nCount_RNA, fill = Sample)) +
  geom_violin(trim = TRUE) +
  theme_bw() +
  labs(title = "nCount_RNA Distribution", x = "Sample", y = "nCount_RNA") +
  scale_fill_brewer(palette = "Set1")

# Violin plot for nFeature_RNA
ggplot(combined_metadata, aes(x = Sample, y = nFeature_RNA, fill = Sample)) +
  geom_violin(trim = TRUE) +
  theme_bw() +
  labs(title = "nFeature_RNA Distribution", x = "Sample", y = "nFeature_RNA") +
  scale_fill_brewer(palette = "Set1")

# Save your list to a file
saveRDS(objects, file = "./objects/sp/integrated/healthy_clean.rds")

