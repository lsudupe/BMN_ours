## SCRIPT: Heterogeneity BONE MARROW

## 05.12.23 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(Matrix)
library(sva)
library(genefilter)
library(RColorBrewer)
library(gridExtra)
library(ggplot2)
library(STutility)
library(ComplexHeatmap)
library(DESeq2)
library(enrichplot)
library(clusterProfiler)
library(org.Mm.eg.db)
library(colorBlindness)
color <- rev(brewer.pal(11,"Spectral"))
#source("./regress_out_function.R")

## Data
M1 <- readRDS("./objects/sp/st_s1_module.rds")
M2 <- readRDS("./objects/sp/st_s2_module.rds")
M8 <- readRDS("./objects/sp/st_s8_module.rds")
M9 <- readRDS("./objects/sp/st_s9_module.rds")

identifiers <- c("M1", "M2", "M8", "M9")
# Use MergeSTData to merge the objects
se.merged <- MergeSTData(
  x = M1,
  y = list(M2, M8, M9),
  add.spot.ids = identifiers,
  merge.data = TRUE
)

# new names
se.merged@meta.data[["community"]] <- as.factor(se.merged@meta.data[["community"]])
meta <- se.merged@meta.data
current_levels <- levels(meta$community)
new_levels_map <- setNames(
  c("Pg1", "Pg2", "Pg3", "Pg4", "Pg5", "Pg6", "Pg7", "Pg8", "Pg9", "Pg10", "Pg11", "Pg12"),
  current_levels
)
meta$community <- factor(meta$community, levels = current_levels, labels = new_levels_map[current_levels])

se.merged <- AddMetaData(se.merged,meta)  

#spatial plot
FeatureOverlay(se.merged, features = c("community"), 
               sampleids = 1:4, ncols = 2, 
               pt.size = 1.3)

## quitar inmunos
#immunoglobulin_genes <- grep("^Igh", rownames(se.merged), value = TRUE)
# Remove these genes from the matrix
#se.merged <- se.merged[!rownames(se.merged) %in% immunoglobulin_genes, ]

## Aggregate counts by community en este caso con el regress por name
Seurat::Idents(object = se.merged) <- se.merged@meta.data[["community"]]
cts <- AggregateExpression(se.merged, 
                           group.by = c("community"),
                           assays = 'RNA',
                           slot = "count",
                           return.seurat = FALSE)

aggregated_counts <- cts$RNA
head(aggregated_counts)

## regress out sample with Combat
batch<- c("S1", "S1", "S2", "S2","S2","S2","S2" ,"S8","S8","S8","S9","S9")
modules <- c("M1", "M1","M2","M2","M2","M2","M2","M8","M8","M9","M9","M9")
#corrected_matrix <- ComBat_seq(counts=aggregated_counts, batch = batch, group = NULL)

head(corrected_matrix)
corrected_matrix <- aggregated_counts

colData_df <- DataFrame(modules)
dds_matrix <- DESeqDataSetFromMatrix(countData = corrected_matrix, colData = colData_df, design = ~modules)
dds_matrix <- DESeq(dds_matrix)

## Extract normalize counts
final_mat <- counts(dds_matrix, normalized = T)
head(final_mat)

## Estandar error
se <- function(x) sd(x)/sqrt(length(x))
t_final_mat <- Matrix::t(final_mat)
se_mat <- apply(t_final_mat, 2, se)
se_mat <- se_mat[order(se_mat, decreasing = T)]
hist(se_mat, main="Histogram of se_mat", xlab="standard error values")
log_se_mat <- log1p(se_mat) # log1p is used to avoid log(0) which is undefined


## select top
# Calculate the 99th quantile on the log-transformed data
quantile_99 <- quantile(log_se_mat, probs = 0.99)
hist(log_se_mat, main="Log-transformed Histogram of se_mat", xlab="Log-transformed  standard error values")
abline(v=quantile_99, col="red", lwd=2, lty=2)
threshold_value <- quantile_99[1]
# To get the threshold value on the original scale of se_mat
original_threshold_value <- exp(threshold_value) - 1
selected_genes <- se_mat[se_mat > original_threshold_value]

## Ordenar por los genes
top <- names(sort(selected_genes, decreasing = TRUE))#[1:100]
final_mat_se <- final_mat[top, ] #pseudobulk ordered most var genes
dim(final_mat_se)

## ana heatmap 
matrix <- final_mat[top, ]
ordered_cols <- c("Pg1", "Pg2", "Pg3", "Pg4", "Pg5", "Pg6", "Pg7", "Pg8", "Pg9", "Pg10", "Pg11", "Pg12")
matrix_ordered <- matrix[, ordered_cols]

pheatmap(matrix_ordered, scale = "row", fontsize = 4)

## pca
pca<- prcomp(t(matrix),center = FALSE, scale. = FALSE) 
PC1_and_PC2<- data.frame(PC1=pca$x[,1], PC2= pca$x[,2])
PC1_and_PC2$sample <- c("M1","M1", "M2", "M2","M2","M2","M2", "M8","M8","M9","M9","M9")
head(PC1_and_PC2)

pdf("./results/pca_sample.pdf")
ggplot(PC1_and_PC2, aes(x=PC1, y=PC2)) + 
  geom_point(aes(color = sample)) +
  theme_bw(base_size = 14)
dev.off()

color_mapping <- c(Pg1="#000000", Pg2="#004949", Pg3="#009292", Pg4="#ff6db6", Pg5="#ffb6db",
                   Pg6="#490092", Pg7="#006ddb", Pg8="#b66dff", Pg9="#6db6ff", Pg10="#b6dbff",
                   Pg11="#920000", Pg12="#924900")

rownames = factor(c("Pg1", "Pg2", "Pg3", "Pg4", "Pg5", "Pg6", "Pg7", "Pg8", "Pg9", "Pg10", "Pg11", "Pg12"),
                  levels = c("Pg1", "Pg2", "Pg3", "Pg4", "Pg5", "Pg6", "Pg7", "Pg8", "Pg9", "Pg10", "Pg11", "Pg12"))

# Plot
pdf("./results/pca_pg.pdf")
ggplot(PC1_and_PC2, aes(x=PC1, y=PC2, color = rownames)) + 
  geom_point(size = 3) +
  scale_color_manual(values = color_mapping) +
  theme_bw(base_size = 14) +
  labs(color = 'Sample')
dev.off()
 
## PLOTS SPATIAL Y VIOLIN
top5 <- top[1:100]
top5

top_batcheffect_name <- top
top_NObatcheffect_name <- top

genesComun <- intersect(top_batcheffect_name, top_NObatcheffect_name)

x <- "Mzb1"
m1 <- FeatureOverlay(M1, features = c(x),pt.size = 1.3, col=color)
m2 <- FeatureOverlay(M2, features = c(x),pt.size = 1.3, col=color)
m3 <- FeatureOverlay(M8, features = c(x),pt.size = 1.3, col=color)
m4 <- FeatureOverlay(M9, features = c(x),pt.size = 1.3, col=color)

grid.arrange(m1, m2, m3, m4, ncol = 2)

"Xbp1" %in% top_NObatcheffect_name
"Mzb1" %in% top_NObatcheffect_name

VlnPlot(se.merged, features = x, group.by = "community", assay = "SCT")


## community spatial plot
# change names
M1@meta.data[["community"]] <- factor(M1@meta.data[["community"]],
                                      levels = c("S1_M1", "S1_M2"),
                                      labels = c("Pg1", "Pg2"))
M2@meta.data[["community"]] <- factor(M2@meta.data[["community"]],
                                      levels = c("S2_M1", "S2_M2", "S2_M3", "S2_M4", "S2_M5"),
                                      labels = c("Pg3", "Pg4", "Pg5", "Pg6","Pg7"))
M8@meta.data[["community"]] <- factor(M8@meta.data[["community"]],
                                      levels = c("S8_M1", "S8_M2"),
                                      labels = c("Pg8", "Pg9"))
M9@meta.data[["community"]] <- factor(M9@meta.data[["community"]],
                                      levels = c("S9_M1", "S9_M2", "S9_M3"),
                                      labels = c("Pg10", "Pg11", "Pg12"))

# Save data
saveRDS(M1,"./objects/sp/st_s1_community.rds")
saveRDS(M2,"./objects/sp/st_s2_community.rds")
saveRDS(M8,"./objects/sp/st_s8_community.rds")
saveRDS(M9, "./objects/sp/st_s9_community.rds")

all_communities <- c("Pg1", "Pg2", "Pg3", "Pg4", "Pg5", "Pg6", "Pg7", "Pg8", "Pg9", "Pg10", "Pg11", "Pg12")
num_communities <- length(unique(all_communities))
paletteMartinColors <- colorBlindness::paletteMartin
color_mapping <- setNames(paletteMartinColors, unique(all_communities))

p1 <- FeatureOverlay(M1, features = "community", pt.size = 1.3,cols = color_mapping)
p2 <- FeatureOverlay(M2, features = "community",pt.size = 1.3,cols = color_mapping)
p3 <- FeatureOverlay(M8, features = "community", pt.size = 1.3,cols = color_mapping)
p4 <- FeatureOverlay(M9, features = "community",pt.size = 1.3,cols = color_mapping)
cowplot::plot_grid(p1, p2, p3, p4, ncol = 2)

pdf("./results/distribution/pg_groups_spatial.pdf", width = 12, height = 8)
print(cowplot::plot_grid(p1, p2, p3, p4, ncol = 2))
dev.off()

FeatureOverlay(M9, features = "new_groups", pt.size = 1.3,cols = color_mapping)

## ORA
sGenes <- names(selected_genes)
entrez_ids <- bitr(sGenes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
head(entrez_ids)

ego <- enrichGO(gene          = entrez_ids$ENTREZID,
                OrgDb         = org.Mm.eg.db,
                ont           = "CC",  # Cellular Component
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)
goplot(ego)

barplot(ego, showCategory=20) 
mutate(ego, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore")

dp <- dotplot(ego, showCategory=30) + ggtitle("dotplot for ORA")
dp <- dp + 
  theme(axis.text.x = element_text(size = 10),  # Adjust x-axis text size
        axis.text.y = element_text(size = 8),  # Adjust y-axis text size
        plot.title = element_text(size = 10)) 
print(dp)

##save data
library(openxlsx)
# Convertir la lista en un data frame
genes_df <- data.frame(Gene = names(selected_genes), Expression = selected_genes)
# Guardar el data frame en un archivo de Excel
write.xlsx(genes_df, file = "./194_genes.xlsx")

##chequear 38 genes
genes38 <- c("Tnfrsf17", "Plpp5", "Sdc1", "Fcrl5", "Slamf7", 
            "Cd38", "Cd79B", "Cd79A", "Ms4A1", "Tnfrsf13B", 
            "Kcnn3", "CadM1", "Slc1A4", "Cav1", "Lsr", "Gpr160", 
            "Tnfrsf13C", "EdnrB", "P2Rx5", "Perp", "Cd22", "Fcer2",
            "Fcrl2", "Gprc5D", "Itga8", "Il5Ra", "Lsamp", "Adam28", 
            "Cd24", "Fcrl1", "Lax1", "Cd19", "Hvcn1", "Cd40", "Parm1", "Ccr10", "Dpep1", "Dcc")

genes38_orto <- c("Tnfrsf17", "Plpp5","Plpp4", "Sdc1", "Fcrl5", "Slamf7", 
             "Cd38", "Cd79b", "Cd79a", "Ms4a1", "Tnfrsf13b", 
             "Kcnn3","Kcnn1", "Kcnn2", "Cadm1", "Slc1a4", "Cav1", "Lsr", "Gpr160", 
             "Tnfrsf13c", "Ednrb", "P2rx5", "Perp", "Cd22", "Fcer2a",
             "Fcrl1","Fcrl5","Fcrl6", "Gprc5d", "Itga8","Itga2b","Itga5","Itgav", "Il5ra", "Lsamp", "Adam28", 
             "Cd24a", "Fcrl1", "Lax1", "Cd19", "Hvcn1", "Cd40", "Parm1", "Ccr10", "Dpep1", "Dcc", "Neo1")


count(genes38_orto %in% names(selected_genes))

"Ighm" %in% selected_genes
"Ighm" %in% names(selected_genes)
"Cd24" %in% names(selected_genes)







## pathways y genes a excell
# Asumiendo que ego es tu objeto y ya tiene los resultados
# Extrae los top 50 pathways
top_pathways <- ego@result[["ID"]][1:50]
top_pathway_descriptions <- ego@result[["Description"]][1:50]

# Preparar un data frame para almacenar los resultados
pathway_genes_df <- data.frame(PathwayID = character(), Description = character(), Genes = character(), stringsAsFactors = FALSE)

# Iterar sobre cada pathway y extraer los genes asociados
for (i in 1:length(top_pathways)) {
  pathway_id <- top_pathways[i]
  genes_in_pathway <- ego@geneSets[[pathway_id]]
  # Filtrar para asegurar que el conjunto de genes no esté vacío
  if (length(genes_in_pathway) > 0) {
    gene_symbols <- ego@gene2Symbol[genes_in_pathway]
    # Filtrar cualquier NA
    gene_symbols <- gene_symbols[!is.na(gene_symbols)]
    # Solo añadir a la lista si hay símbolos de genes disponibles
    if (length(gene_symbols) > 0) {
      pathway_genes_df <- rbind(pathway_genes_df, data.frame(PathwayID = pathway_id, Description = top_pathway_descriptions[i], Genes = paste(gene_symbols, collapse = ", ")))
    }
  }
}

# Escribir el data frame en un archivo Excel
write.xlsx(pathway_genes_df, file = "Top50_Pathways_and_Genes.xlsx")

