



library(colorBlindness)
library(cowplot)
library(ggplot2)
library(patchwork)

distr_groups <- readRDS("./objects/sp/st_q3.rds")

distr_groups@meta.data[["new_groups"]] <- as.factor(distr_groups@meta.data[["new_groups"]])

# Rename the levels
levels(distr_groups@meta.data[["new_groups"]])[levels(distr_groups@meta.data[["new_groups"]]) == "Surronding"] <- "Border_zone"
levels(distr_groups@meta.data[["new_groups"]])[levels(distr_groups@meta.data[["new_groups"]]) == "Low"] <- "Remote_zone"
levels(distr_groups@meta.data[["new_groups"]])[levels(distr_groups@meta.data[["new_groups"]]) == "Hot_spot"] <- "Hotspot"

# Check the updated unique values
unique(distr_groups@meta.data[["new_groups"]])

Idents(distr_groups) <- "new_groups"
distr_groups_subset <- subset(distr_groups, subset = new_groups != "rest")
distr_groups_subset@meta.data[["new_groups"]] <- droplevels(distr_groups_subset@meta.data[["new_groups"]])
unique(distr_groups_subset@meta.data[["new_groups"]])

selectedColors <- c(rest="#E5B17E", Remote_zone="#C3E57E", Border_zone="#7EC3E5", Hotspot="#CC5151",notknown="white")


# Create a new color mapping for 'new_groups'
# Assuming 'new_groups' is a factor with 4 levelsnr
color_mapping <- setNames(selectedColors, levels(M9@meta.data[["new_groups"]]))

# Plot using FeatureOverlay
a<- FeatureOverlay(distr_groups, features = "new_groups",sampleids = 1:4, ncols = 2,
                    pt.size = 1.3,cols = color_mapping)

pdf("./results/distribution/comunities_spatial.pdf", width = 12, height = 8)
print(a)
dev.off()

# Display the plot
print(p)

VlnPlot(distr_groups_subset, features = "Tnfrsf17", 
                    group.by = "new_groups", cols = selectedColors,
                    pt.size = 1) 


my_list <- c("Tnfrsf17", "Plpp5", "Sdc1", "Fcrl5", "Slamf7", 
             "Cd38", "Cd79B", "Cd79A", "Ms4A1", "Tnfrsf13B", 
             "Kcnn3", "CadM1", "Slc1A4", "Cav1", "Lsr", "Gpr160", 
             "Tnfrsf13C", "EdnrB", "P2Rx5", "Perp", "Cd22", "Fcer2",
             "Fcrl2", "Gprc5D", "Itga8", "Il5Ra", "Lsamp", "Adam28", 
             "Cd24", "Fcrl1", "Lax1", "Cd19", "Hvcn1", "Cd40", "Parm1", "Ccr10", "Dpep1", "Dcc")

for (element in my_list) {
  print(element, quote = FALSE)
}


mouse <- readRDS("./objects/sp/st_q3_all.rds")
distr_groups <- mouse


####spatial plots
distr_groups
genes <- c("Sdc1", "Cd81", "Slamf7", "Prdm1", "Cd38", "Tnfrsf17", "Xbp1", "Mzb1", "Cd52")
samples <- unique(distr_groups_subset@meta.data$name)
levels(distr_groups_subset@meta.data$new_groups) <- c("Remote_zone", "Border_zone", "Hotspot")
selectedColors <- c(rest="#E5B17E", Remote_zone="#C3E57E", Border_zone="#7EC3E5", Hotspot="#CC5151",notknown="white")


for (gen in genes) {
  # Verificar si el gen está en el objeto Seurat
  if(gen %in% rownames(distr_groups_subset@assays$SCT@data)) {
    # Calcula el valor máximo para el gen en todas las muestras
    max_expression <- max(distr_groups_subset@assays$SCT@data[gen, ])
    
    # Define el valor mínimo y máximo para la escala
    min_value <- 0
    max_value <- max_expression
    
    # Crear el gráfico
    p <- FeatureOverlay(distr_groups, features = c(gen), sampleids = 1:6, ncols = 2, pt.size = 1, 
                        value.scale = "all", cols = color) +
      scale_fill_gradientn(colours = color,
                           breaks = c(min_value, max_value / 2, max_value),
                           labels = c("Min", max_value / 2, "Max"),
                           limits = c(min_value, max_value))
    
    # Guardar el gráfico en un archivo PDF
    pdf(paste0("./results/genes_plots/all_areas/", gen, ".pdf"), width = 12, height = 10)
    print(p)
    dev.off()
    
    # Initialize an empty list to store plots for each sample
    sample_plots <- list()
    
    for (sample in samples) {
      # Filter the Seurat object for the current sample
      sample_data <- subset(distr_groups_subset, subset = name == sample)
      
      # Create a violin plot for the current gene and sample
      vln_plot <- VlnPlot(sample_data, features = gen, 
                          group.by = "new_groups", cols = selectedColors,
                          pt.size = 1) + NoLegend() + ggtitle(sample)
      
      # Add the plot to the list
      sample_plots[[sample]] <- vln_plot
    }
    
    # Combine all sample plots into one plot
    combined_plot <- wrap_plots(sample_plots, ncol = 1)
    
    # Save the combined plot to a single page in a PDF
    pdf(paste0("./results/genes_plots/all_areas/", gen, "_violin.pdf"), width = 8, height = 6 * length(samples))
    print(combined_plot)
    dev.off()
    
  } else {
    print(paste("El gen", gen, "no se encuentra en el objeto Seurat."))
  }
}

M1 <- readRDS("./objects/sp/st_s1_community.rds")
M2 <- readRDS("./objects/sp/st_s2_community.rds")
M8 <- readRDS("./objects/sp/st_s8_community.rds")
M9 <- readRDS("./objects/sp/st_s9_community.rds")

##communities
samples <- c(M1,M2,M8,M9)
names(samples) <- c("M1","M2","M8","M9")
samples <- list(M1 = M1, M2 = M2, M8 = M8, M9 = M9)

combined_samples <- Reduce(function(x, y) merge(x, y), samples)
combined_samples@meta.data$community <- factor(combined_samples@meta.data$community, 
                                               levels = c("Pg1", "Pg2", "Pg3", "Pg4", "Pg5", "Pg6", "Pg7", "Pg8", "Pg9", "Pg10", "Pg11", "Pg12"))



color_mapping <- c(Pg1="#000000", Pg2="#004949", Pg3="#009292", Pg4="#ff6db6", Pg5="#ffb6db",
                   Pg6="#490092", Pg7="#006ddb", Pg8="#b66dff", Pg9="#6db6ff", Pg10="#b6dbff",
                   Pg11="#920000", Pg12="#924900")

genes <- c("Xbp1", "Mzb1", "Cd74", "Hsp90b1", "S100a9", "S100a8", "B2m", "Tmsb4x", "Sec11c", "Tpt1",
              "Cybb", "Fth1", "Cd52", "Ighg2b", "Ighg3", "Cyba", "Ly6k", "S100a6")



# Function to get max expression of a gene from a Seurat object
get_max_expression <- function(seurat_object, gene_name) {
  if (gene_name %in% rownames(seurat_object@assays$SCT@data)) {
    return(max(seurat_object@assays$SCT@data[gene_name, ], na.rm = TRUE))
  } else {
    return(NA)
  }
}

# Initialize an empty list to store max expressions for each gene
max_expressions <- sapply(genes, function(gene) {
  max(sapply(samples, get_max_expression, gene_name = gene), na.rm = TRUE)
}, USE.NAMES = TRUE)

for (gene in genes) {
  plot_list <- list()
  
  max_value <- max_expressions[gene]
  
  # Check if max_value is NA or not a finite number
  if (is.na(max_value) || !is.finite(max_value)) {
    warning(paste("Max value for gene", gene, "is not finite. Skipping this gene."))
    next  # Skip this iteration of the loop
  }
  
  # Round the max_value up to the nearest 0.5 if it's not already a multiple of 0.5
  max_value <- ceiling(max_value * 2) / 2
  
  # Create a sequence of breaks from 0 to max_value, with steps of 0.5
  breaks <- seq(0, max_value, by = 1)
  
  # Create labels for each break
  labels <- as.character(breaks)
  
  for (sample_name in names(samples)) {
    seurat_object <- samples[[sample_name]]
    
    # Create the FeatureOverlay plot for the current gene and sample
    p <- FeatureOverlay(seurat_object, features = c(gene), sampleids = 1:4, ncols = 2, pt.size = 1.1, 
                        value.scale = "all", cols = color) +
      scale_fill_gradientn(colours = color,
                           breaks = breaks,
                           labels = labels,
                           limits = c(0, max_value))
    plot_list[[sample_name]] <- p
  }
  
  # Combine the plots into a single plot
  combined_plot <- plot_grid(plotlist = plot_list, ncol = 2)
  
  # Save the combined plot to a PDF
  pdf_filename <- paste0("./results/genes_plots/hotspot/", gene, "_community.pdf")
  pdf(pdf_filename, width = 17, height = 8)
  print(combined_plot)
  dev.off()
}

###BY SAMPLE
# Function to create a violin plot for a given gene across different samples
create_violin_plot <- function(gene, samples) {
  plot_list <- list()
  
  for (sample_name in names(samples)) {
    seurat_object <- samples[[sample_name]]
    vln_plot <- VlnPlot(seurat_object, features = gene,  cols = color_mapping,
                        group.by = "community", pt.size = 1) + NoLegend() + ggtitle(sample_name)
    plot_list[[sample_name]] <- vln_plot
  }
  
  combined_vln_plot <- wrap_plots(plot_list, ncol = 2)
  return(combined_vln_plot)
}

# Adding violin plots for each gene
for (gene in genes) {
  violin_plot <- create_violin_plot(gene, samples)
  
  # Save the violin plot to a PDF
  pdf_filename <- paste0("./results/genes_plots/hotspot/violin_", gene, "_community.pdf")
  pdf(pdf_filename, width = 12, height = 6)
  print(violin_plot)
  dev.off()
  
}

##ALL TOGHETER
samples <- lapply(names(samples), function(x) {
  samples[[x]][["sample"]] <- x
  return(samples[[x]])
})

color_mapping <- c(Pg1="#000000", Pg2="#004949", Pg3="#009292", Pg4="#ff6db6", Pg5="#ffb6db",
                   Pg6="#490092", Pg7="#006ddb", Pg8="#b66dff", Pg9="#6db6ff", Pg10="#b6dbff",
                   Pg11="#920000", Pg12="#924900")


# Create violin plots for each gene, separated by 'community'
for (gene in genes) {
  if (gene %in% rownames(combined_samples@assays$SCT@data)) {
    vln_plot <- VlnPlot(combined_samples, features = gene,  cols = color_mapping,
                        group.by = "community", pt.size = 1) + NoLegend()
    
    # Save the violin plot to a PDF
    pdf_filename <- paste0("./results/genes_plots/hotspot/ALLviolin_", gene, "_community.pdf")
    pdf(pdf_filename, width = 8, height = 4)
    print(vln_plot)
    dev.off()
  } else {
    warning(paste("Gene", gene, "not found in the dataset. Skipping."))
  }
}




