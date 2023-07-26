## SCRIPT: Tcell exhausted correlations BM project

## 26.07.23 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(ggplot2)
library(STutility)
library(dplyr)
library(UCell)
library(tidyr)


##individual
M1_s <- readRDS("./objects/pc_clusters/M1_s_dormant.rds")
M2_s <- readRDS("./objects/pc_clusters/M2_s_dormant.rds")
M8_s <- readRDS("./objects/pc_clusters/M8_s_dormant.rds")
M9_s <- readRDS("./objects/pc_clusters/M9_s_dormant.rds")

candidate <- c("Pdcd1", "Ctla4", "Tnfrsf9", "Havcr2", "Tox", "Tigit", "Wars", "Rsad2",
               "Mcm7", "Mx1", "Ndfip2", "Enosf1", "Ccdc141", "Stmn1", "Ttn", "Faslg",
               "Mcm5", "Nab1", "Phlda1", "Mcm3", "Pcna", "Gapdh", "Oasl", "Ifi44l",
               "Tbc1d4", "Slc43a3", "Pam", "Ccl3", "Acp5", "Oas3", "Cd38", "Tnfs10",
               "Gbp2", "Kif20b", "Ctsb")

#M1
Idents(object = M1_s) <- M1_s@meta.data[["pc_clusters"]]
M1_subset <- SubsetSTData(object = M1_s, ident = c("M1_group_1_cluster5","M1_group_1_cluster6", "M1_group_1_cluster7"))

#M2
Idents(object = M2_s) <- M2_s@meta.data[["pc_clusters"]]
M2_subset <- SubsetSTData(object = M2_s, ident = c("M2_group_1_cluster5","M2_group_1_cluster6", "M2_group_1_cluster7",
                                                   "M2_group_2_cluster5", "M2_group_2_cluster6", "M2_group_2_cluster7",
                                                   "M2_group_3_cluster5","M2_group_3_cluster6", "M2_group_3_cluster7",
                                                   "M2_group_4_cluster5","M2_group_4_cluster6", "M2_group_4_cluster7"
))

#M8
Idents(object = M8_s) <- M8_s@meta.data[["pc_clusters"]]
M8_subset <- SubsetSTData(object = M8_s, ident = c("M8_group_1_cluster5","M8_group_1_cluster6", "M8_group_1_cluster7",
                                                   "M8_group_2_cluster5", "M8_group_2_cluster6", "M8_group_2_cluster7"))


#M9
Idents(object = M9_s) <- M9_s@meta.data[["pc_clusters"]]
M9_subset <- SubsetSTData(object = M9_s, ident = c("M9_group_1_cluster6",
                                                   "M9_group_2_cluster5", "M9_group_2_cluster6",
                                                   "M9_group_3_cluster6"))


###merge data
se_merged <- MergeSTData(M1_s, y = c(M2_s, M8_s, M9_s), 
                         add.spot.ids = c("M1_s", "M2_s", "M8_s", "M9_s"), project = "BM")

# Assume 'seurat_obj' is your Seurat object
# Normalize data
seurat_obj <- NormalizeData(se_merged)

###########DAVID values AND EXCELL#########

df <- seurat_obj@meta.data

# Group the data by 'labels', 'name', 'clustering', and 'pc_clusters', then calculate median values
df_grouped <- df %>%
  select(name, labels, clustering, pc_clusters, Tcell:DC, signature_1_dormant, signature_1_candidate) %>%
  group_by(name, labels, clustering, pc_clusters) %>%
  summarise(across(Tcell:DC, median, na.rm = TRUE),
            across(signature_1_dormant, median, na.rm = TRUE),
            across(signature_1_candidate, median, na.rm = TRUE),
            .groups = "keep")

# Save 'pc_clusters' as last column
df_grouped <- df_grouped %>% relocate(pc_clusters, .after = last_col())

# Create two dataframes based on 'clustering' values
df_group1 <- df_grouped %>% filter(clustering %in% c(1,2,3,4))
df_group2 <- df_grouped %>% filter(clustering %in% c(5,6,7))

#####group them
# Create a new column 'sample_group' combining 'name' and 'labels'
df_group1$sample_group <- paste(df_group1$name, df_group1$labels, sep = "_")
df_group2$sample_group <- paste(df_group2$name, df_group2$labels, sep = "_")

# Group the data by 'sample_group' and calculate median values
df_group1_grouped <- df_group1 %>%
  group_by(sample_group) %>%
  summarise(across(Tcell:DC, median, na.rm = TRUE),
            across(signature_1_dormant, median, na.rm = TRUE),
            across(signature_1_candidate, median, na.rm = TRUE),
            .groups = "keep")

df_group2_grouped <- df_group2 %>%
  group_by(sample_group) %>%
  summarise(across(Tcell:DC, median, na.rm = TRUE),
            across(signature_1_dormant, median, na.rm = TRUE),
            across(signature_1_candidate, median, na.rm = TRUE),
            .groups = "keep")

# Load required packages
library(corrplot)
library(gplots)

# Determine the columns with non-zero standard deviation
#cols_to_keep1 <- apply(df_group1_grouped[,3:ncol(df_group1_grouped)], 2, sd) != 0
#cols_to_keep2 <- apply(df_group2_grouped[,3:ncol(df_group2_grouped)], 2, sd) != 0

# Subset the data frames
df_grouped1_grouped <- df_group1_grouped[, c(TRUE, TRUE, cols_to_keep1)]
df_grouped2_grouped <- df_group2_grouped[, c(TRUE, TRUE, cols_to_keep2)]

# Compute correlations for df_grouped1_grouped and df_grouped2_grouped
corr1 <- cor(df_grouped1_grouped[,3:ncol(df_grouped1_grouped)], use = "pairwise.complete.obs")
corr2 <- cor(df_grouped2_grouped[,3:ncol(df_grouped2_grouped)], use = "pairwise.complete.obs")

hc_order <- hclust(dist(corr1))$order

# Set graphical parameters to make the labels smaller
par(cex=0.6)

# Generate correlation plot for df_grouped1_grouped
pdf("./results/correlation/pc_groups/Correlation_df_group1_grouped_Tcell.pdf")
corrplot(corr1, method="color", col= colorRampPalette(c("blue", "white", "red"))(200),
         type="upper", addCoef.col = "black", tl.col="black", tl.srt=45, number.cex=0.70)
dev.off()

# Generate correlation plot for df_grouped2_grouped
pdf("./results/correlation/pc_groups/Correlation_df_group2_grouped_Tcell.pdf")
corrplot(corr2, method="color", col= colorRampPalette(c("blue", "white", "red"))(200),
         type="upper", addCoef.col = "black", tl.col="black", tl.srt=45, number.cex=0.70)
dev.off()
##cor fin
