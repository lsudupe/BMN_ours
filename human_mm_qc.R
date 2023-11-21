## SCRIPT: Human sample MM cleanning and QC

## 21.11.23 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(dplyr)
library(hrbrthemes)
color <- rev(brewer.pal(11,"Spectral"))

## Data
mm <- readRDS("./objects/sp/human/human_combined.rds")

## Spatial plots
for (i in 1:length(mm)){
  a <- mm[[i]]
  b <- names(mm[i])
  
  # count plot
  count <- "nCount_RNA"
  c <- c(min(a@meta.data[count]), max(a@meta.data[count]))
  p <- FeatureOverlay(a, features = count, ncols = 1, pt.size = 1.5, slot = "count",
                 value.scale = "all" ,cols = color) +
    scale_fill_gradientn(colours = color,
                         breaks = c,
                         labels = c,
                         limits = c)
  # Save the plot to a PDF file
  pdf(paste0("./results/human/QC/mm/", b, "_ncount.pdf", sep = ""))
  print(p)
  dev.off()
  
  # feature plot
  feature <- "nFeature_RNA"
  c <- c(min(a@meta.data[feature]), max(a@meta.data[feature]))
  p <- FeatureOverlay(a, features = feature, ncols = 1, pt.size = 1.5, slot = "count",
                 value.scale = "all" ,cols = color) +
    scale_fill_gradientn(colours = color,
                         breaks = c,
                         labels = c,
                         limits = c)
  # Save the plot to a PDF file
  pdf(paste0("./results/human/QC/mm/", b, "_nfeature.pdf", sep = ""))
  print(p)
  dev.off()
  
}


## Combine meta
combined_metadata <- data.frame()

# Identify common columns
common_cols <- Reduce(intersect, lapply(mm, function(x) names(x@meta.data)))

# Loop through each object and extract common metadata columns
for (sample_name in names(mm)) {
  metadata <- mm[[sample_name]]@meta.data[, common_cols, drop = FALSE]
  metadata$Sample <- sample_name
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
  scale_fill_brewer(palette = "Set1")+
  coord_cartesian(ylim = c(0, 5000))

# Violin plot for nFeature_RNA
# Prepare the data
plot_data <- combined_metadata %>%
  group_by(Sample) %>%
  mutate(num = n()) %>%
  ungroup() %>%
  mutate(myaxis = paste0(Sample, "\n", "n=", num))

# Create the plot
  ggplot(combined_metadata, aes(x = Sample, y = nFeature_RNA, fill = Sample)) +
    geom_violin(width = 1.4) +
    geom_boxplot(width = 0.1, color = "grey", alpha = 0.2) +
    scale_fill_viridis(discrete = TRUE) +
    coord_cartesian(ylim = c(0, 5000))
    theme_ipsum() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 11)
    ) +
    ggtitle("A Violin wrapping a boxplot") +
    xlab("")

## Check quartiles
for (i in 1:length(mm)){
  a <- mm[[i]]
  b <- names(mm[i])
  metadata <- a@meta.data
      
  # Calculate the quartiles for nFeature_RNA
  Q1 <- quantile(metadata$nFeature_RNA, 0.1)
  Q4 <- quantile(metadata$nFeature_RNA, 0.9)
      
  # Filter rows
  filtered_metadata <- metadata[metadata$nFeature_RNA >= Q1 & metadata$nFeature_RNA <= Q4, ]
      
  # Update the Seurat object with the filtered metadata
  a@meta.data <- filtered_metadata
      
  # Plot the feature plot
  feature <- "nFeature_RNA"
  c <- c(min(a@meta.data[feature]), max(a@meta.data[feature]))
  p <- FeatureOverlay(a, features = feature, ncols = 1, pt.size = 1.5, slot = "count",
                   value.scale = "all" ,cols = color) +
        scale_fill_gradientn(colours = color,
                             breaks = c,
                             labels = c,
                             limits = c)
      
  # Save the plot to a PDF file
  pdf_filename <- paste0("./results/human/QC/mm_quartiles/", b, "_nfeature.pdf")
   ggsave(pdf_filename, plot = p, device = "pdf")
}
    

