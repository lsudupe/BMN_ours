## SCRIPT: PC cluster differentiation visium data BM project

## 18.05.23 Laura Sudupe , git @lsudupe

## Libraries
library(Seurat)
library(tidyverse)
library(STutility)
library(UCell)
library(RColorBrewer)

# Data
se <- readRDS("./objects/heterogeneity/se_hierarchical.rds")

Idents(object = se) <- "name"
#M8
M8 <- SubsetSTData(se, ident = "M8_F2_1C")
meta <- M8@meta.data
#rownames(meta) <- sub("_5", "", rownames(meta))
M8@meta.data <- meta
meta <- meta[, ncol(meta), drop = FALSE]

write.csv(meta, "./data/areas/m8_clustering.csv", row.names=TRUE)
M8_s <- ManualAnnotation(M8)
meta <- M8_s@meta.data
meta <- meta %>%
  mutate(pc_clusters = paste0("M8_", labels, "_cluster", clustering))
M8_s@meta.data <- meta
saveRDS(M8_s, "./objects/pc_clusters/M8_s.rds")
M8_s <- readRDS("./objects/pc_clusters/M8_s.rds")

#M2
M2 <- SubsetSTData(se, ident = "M2_F_2B")
meta <- M2@meta.data
#rownames(meta) <- sub("_5", "", rownames(meta))
M2@meta.data <- meta
meta <- meta[, ncol(meta), drop = FALSE]

write.csv(meta, "./data/areas/m2_clustering.csv", row.names=TRUE)
#M2_s <- ManualAnnotation(M2)
meta <- M2_s@meta.data
meta <- meta %>%
  mutate(pc_clusters = paste0("M2_", labels, "_cluster", clustering))
M2_s@meta.data <- meta
saveRDS(M2_s, "./objects/pc_clusters/M2_s.rds")
M2_s <- readRDS("./objects/pc_clusters/M2_s.rds")

#M9
M9 <- SubsetSTData(se, ident = "M9_F2_1C")
meta <- M9@meta.data
#rownames(meta) <- sub("_5", "", rownames(meta))
M9@meta.data <- meta
meta <- meta[, ncol(meta), drop = FALSE]

write.csv(meta, "./data/areas/m9_clustering.csv", row.names=TRUE)
#M9_s <- ManualAnnotation(M9)
meta <- M9_s@meta.data
meta <- meta %>%
  mutate(pc_clusters = paste0("M9_", labels, "_cluster", clustering))
M9_s@meta.data <- meta
saveRDS(M9_s, "./objects/pc_clusters/M9_s.rds")
M9_s <- readRDS("./objects/pc_clusters/M9_s.rds")

#M1
M1 <- SubsetSTData(se, ident = "M1_fem_1C")
meta <- M1@meta.data
#rownames(meta) <- sub("_5", "", rownames(meta))
M1@meta.data <- meta
meta <- meta[, ncol(meta), drop = FALSE]

write.csv(meta, "./data/areas/m1_clustering.csv", row.names=TRUE)
#M1_s <- ManualAnnotation(M1)
meta <- M1_s@meta.data
meta <- meta %>%
  mutate(pc_clusters = paste0("M1_", labels, "_cluster", clustering))
M1_s@meta.data <- meta
saveRDS(M1_s, "./objects/pc_clusters/M1_s.rds")
M1_s <- readRDS("./objects/pc_clusters/M1_s.rds")

###combine the data
se_s <- MergeSTData(M1_s, y = c(M2_s, M8_s, M9_s), 
                  add.spot.ids = c("M1_s", "M2_s", "M8_s", "M9_s"), project = "BM")

saveRDS(se_s, "./objects/pc_clusters/combined_s.rds")

################ CHECK DORMANT SIGNATURE
objects <- c(M1_s,M2_s,M8_s,M9_s)
names(objects) <- c("M1_s","M2_s","M8_s","M9_s")
#Check dormant enrichment score
dormant <- c("C1qa", "Aif1" ,"Axl" , "Mpeg1", "H2-Eb1",
             "Nr1h3" ,"Sdc3" ,"Fcgr3" ,"Oas1a", "Ifit2", "Ly6a",
             "Ly86" ,"Bcl2a1b", "Cd19" , "Hebp1", "Sirpa" ,"Slfn2" ,"Tnfaip2" ,
             "Vpreb3", "Oasl2","Slc43a2", "Ifi44", "Gng10" ,"Cxcl16" ,"Ptprcap" ,
             "Anxa2" ,"Rgs2" ,"Tmed3" ,"Igll1", "Hpgd", "Glipr1", "Cd4" ,"Cd84" ,"Gbp2", "AB124611", "Slc44a2" ,
             "Samd9l" ,"Oas1g" ,"Fcgr1", "Pla2g15", "Tifa" ,"Pmp22", "Abcc3" ,"S100a10")

candidate <- c("Pdcd1", "Ctla4", "Tnfrsf9", "Havcr2", "Tox", "Tigit", "Wars", "Rsad2",
               "Mcm7", "Mx1", "Ndfip2", "Enosf1", "Ccdc141", "Stmn1", "Ttn", "Faslg",
               "Mcm5", "Nab1", "Phlda1", "Mcm3", "Pcna", "Gapdh", "Oasl", "Ifi44l",
               "Tbc1d4", "Slc43a3", "Pam", "Ccl3", "Acp5", "Oas3", "Cd38", "Tnfs10",
               "Gbp2", "Kif20b", "Ctsb")

library(RColorBrewer)
color <- brewer.pal(12,"Spectral")
color <- rev(color)
color <- rainbow(12)


new <- c()
for (i in 1:length(objects)){
  a <- objects[[i]]
  b <- names(objects[i])
  
  ## Add UCellScore
  vector<- ScoreSignatures_UCell(a@assays[["RNA"]]@counts, features = list(dormant))
  a@meta.data[["signature_1_dormant"]] <- as.vector(vector)
  
  meta <- a@meta.data
  
  lm <- lm(meta$signature_1_dormant ~ meta$MM_MIC, data =meta)
  residuals <- lm$residuals
  a@meta.data[["residuals_dormant_ucellscore"]] <- residuals
  
  ##plots
  p <- FeatureOverlay(a, features = c("signature_1_dormant"), pt.size = 1.8, 
                      value.scale = "all" ,cols = color)
  
  pdf(paste("./results/pc_clusters/", b,"_spatial_Ucell.pdf",sep=""))
  print(p)
  dev.off()
  
  p <- FeatureOverlay(a, features = c("residuals_dormant_ucellscore"), pt.size = 1.8, 
                      value.scale = "all" ,cols = color)
  
  pdf(paste("./results/pc_clusters/", b,"_spatial_Ucell_residuals.pdf",sep=""))
  print(p)
  dev.off()
  
  ####add Tcell exhausted signature
  vector<- ScoreSignatures_UCell(a@assays[["RNA"]]@counts, features = list(candidate))
  a@meta.data[["signature_1_candidate"]] <- as.vector(vector)
  
  ## add object to list
  new[[length(new) + 1]] <- a
}

names(new) <- c("M1_s","M2_s","M8_s","M9_s")
## separate list
list2env(new,envir=.GlobalEnv)

#save the objects
saveRDS(M1_s, "./objects/pc_clusters/M1_s_dormant.rds")
saveRDS(M2_s,"./objects/pc_clusters/M2_s_dormant.rds")
saveRDS(M8_s, "./objects/pc_clusters/M8_s_dormant.rds")
saveRDS(M9_s, "./objects/pc_clusters/M9_s_dormant.rds")

M1_s <- readRDS("./objects/pc_clusters/M1_s_dormant.rds")
M2_s <- readRDS("./objects/pc_clusters/M2_s_dormant.rds")
M8_s <- readRDS("./objects/pc_clusters/M8_s_dormant.rds")
M9_s <- readRDS("./objects/pc_clusters/M9_s_dormant.rds")

###plots
###merge data
se_s <- MergeSTData(M1_s, y = c(M2_s, M8_s, M9_s), 
                    add.spot.ids = c("M1_s", "M2_s", "M8_s", "M9_s"), project = "BM")

saveRDS(se_s, "./objects/pc_clusters/combined_s_dormant.rds")
se_s <- readRDS("./objects/pc_clusters/combined_s_dormant.rds")

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


#######TWO LOOPS FOR TWO LISTS
normal <- c(M1_s,M2_s,M8_s,M9_s)
names(normal) <- c("M1_s","M2_s","M8_s","M9_s")

for (i in 1:length(normal)){
  a <- normal[[i]]
  b <- names(normal[i])
  
  color <- rainbow(30)
  p <- FeatureOverlay(a, features = c("pc_clusters"), ncols = 1, pt.size = 1.5,cols = color)
  
  pdf(paste("./results/pc_clusters/all/", b,"_clusters567.pdf",sep=""))
  print(p)
  dev.off()
  
  color <- brewer.pal(11,"Spectral")
  color <- rev(color)
  
  p <- FeatureOverlay(a, features = c("signature_1_dormant"), ncols = 1, pt.size = 1.5, 
                      value.scale = "all" ,cols = color)
  
  pdf(paste("./results/pc_clusters/all/", b,"_spatial_dormant.pdf",sep=""))
  print(p)
  dev.off()
  
  p <- FeatureOverlay(a, features = c("residuals_dormant_ucellscore"), ncols = 1, pt.size = 1.5, 
                      value.scale = "all" ,cols = color)
  
  pdf(paste("./results/pc_clusters/all/", b,"_spatial_dormant_residuals.pdf",sep=""))
  print(p)
  dev.off()
  
  ##meta plots dormant signature
  df <- a@meta.data
  df_selected <- df %>% select(name, pc_clusters, signature_1_dormant, labels)
  
  df_grouped <- df_selected %>%
    group_by(name, pc_clusters, labels) %>%
    summarise(avg_signature_1_dormant = mean(signature_1_dormant, na.rm = TRUE))
  
  # Barplot
  plot <- ggplot(df_grouped, aes(x = pc_clusters, y = avg_signature_1_dormant, fill = labels)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(title = "Average 'signature_1_dormant' with respect to 'group_clusters' for each 'name'",
         x = "Group Clusters", 
         y = "Average Signature 1 Dormant", 
         fill = "Label") +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1)) +  # Adjust text size and angle here
    facet_wrap(~name)  # Separate plots by 'name'
  
  pdf(paste("./results/pc_clusters/all/", b,"_dormant_clusters_barplot.pdf",sep=""))
  print(plot)
  dev.off()
  
  #Violinplot
  df_selected <- df %>% select(name, pc_clusters, signature_1_dormant, labels)
  
  # Plot
  plot <- ggplot(df_selected, aes(x = pc_clusters, y = signature_1_dormant, fill = labels)) +
    geom_violin(scale = "width", adjust = 1) +
    labs(title = "'Signature_1_dormant' distributions with respect to 'group_clusters' for each 'name'",
         x = "Group Clusters", 
         y = "Signature 1 Dormant", 
         fill = "Label") +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1)) +  # Adjust text size and angle here
    facet_wrap(~name)  # Separate plots by 'name'
  
  pdf(paste("./results/pc_clusters/all/", b,"_dormant_clusters_violinplot.pdf",sep=""))
  print(plot)
  dev.off()
  
  # Plot
  plot <- ggplot(df_selected, aes(x = pc_clusters, y = signature_1_dormant, fill = labels)) +
    geom_boxplot() +
    labs(title = "'Signature_1_dormant' distributions with respect to 'group_clusters' for each 'name'",
         x = "Group Clusters", 
         y = "Signature 1 Dormant", 
         fill = "Label") +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1)) +  # Adjust text size and angle here
    facet_wrap(~name)  # Separate plots by 'name'
  
  pdf(paste("./results/pc_clusters/all/", b,"_dormant_clusters_boxplot.pdf",sep=""))
  print(plot)
  dev.off()
  
  ##meta plots dormant RESIDUALS signature
  df <- a@meta.data
  df_selected <- df %>% select(name, pc_clusters, residuals_dormant_ucellscore, labels)
  
  df_grouped <- df_selected %>%
    group_by(name, pc_clusters, labels) %>%
    summarise(avg_signature_1_dormant = mean(residuals_dormant_ucellscore, na.rm = TRUE))
  
  # Barplot
  plot <- ggplot(df_grouped, aes(x = pc_clusters, y = avg_signature_1_dormant, fill = labels)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(title = "Average 'signature_1_dormant' with respect to 'group_clusters' for each 'name'",
         x = "Group Clusters", 
         y = "Average Signature 1 Dormant", 
         fill = "Label") +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1)) +  # Adjust text size and angle here
    facet_wrap(~name)  # Separate plots by 'name'
  
  pdf(paste("./results/pc_clusters/all/", b,"_dormant_residuals_clusters_barplot.pdf",sep=""))
  print(plot)
  dev.off()
  
  #Violinplot
  df_selected <- df %>% select(name, pc_clusters, residuals_dormant_ucellscore, labels)
  
  # Plot
  plot <- ggplot(df_selected, aes(x = pc_clusters, y = residuals_dormant_ucellscore, fill = labels)) +
    geom_violin(scale = "width", adjust = 1) +
    labs(title = "'Signature_1_dormant' distributions with respect to 'group_clusters' for each 'name'",
         x = "Group Clusters", 
         y = "Signature 1 Dormant", 
         fill = "Label") +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1)) +  # Adjust text size and angle here
    facet_wrap(~name)  # Separate plots by 'name'
  
  pdf(paste("./results/pc_clusters/all/", b,"_dormant_residuals_clusters_violinplot.pdf",sep=""))
  print(plot)
  dev.off()
  
  # Plot
  plot <- ggplot(df_selected, aes(x = pc_clusters, y = residuals_dormant_ucellscore, fill = labels)) +
    geom_boxplot() +
    labs(title = "'Signature_1_dormant' distributions with respect to 'group_clusters' for each 'name'",
         x = "Group Clusters", 
         y = "Signature 1 Dormant", 
         fill = "Label") +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1)) +  # Adjust text size and angle here
    facet_wrap(~name)  # Separate plots by 'name'
  
  pdf(paste("./results/pc_clusters/all/", b,"_dormant_residuals_clusters_boxplot.pdf",sep=""))
  print(plot)
  dev.off()
  
  
  #########Each cluster plots
  p <- ST.FeaturePlot(a, features = "pc_clusters",   ncol =3,
                      split.labels = T, pt.size = 0.3) & theme(plot.title = element_blank(), strip.text = element_blank())
  
  pdf(paste("./results/pc_clusters/all/", b,"_each_cluster.pdf",sep=""))
  print(p)
  dev.off()
  
}

subsets <- c(M1_subset,M2_subset, M8_subset, M9_subset)
names(subsets) <- c("M1_subset","M2_subset", "M8_subset", "M9_subset")

for (i in 1:length(subsets)){
  a <- subsets[[i]]
  b <- names(subsets[i])
  
  color <- rainbow(12)
  p <- FeatureOverlay(a, features = c("pc_clusters"), ncols = 1, pt.size = 1.5,cols = color)
  
  pdf(paste("./results/pc_clusters/", b,"_clusters567.pdf",sep=""))
  print(p)
  dev.off()
  
  color <- brewer.pal(11,"Spectral")
  color <- rev(color)
  
  p <- FeatureOverlay(a, features = c("signature_1_dormant"), ncols = 1, pt.size = 1.5, 
                      value.scale = "all" ,cols = color)
  
  pdf(paste("./results/pc_clusters/", b,"_spatial_dormant.pdf",sep=""))
  print(p)
  dev.off()
  
  p <- FeatureOverlay(a, features = c("residuals_dormant_ucellscore"), ncols = 1, pt.size = 1.5, 
                      value.scale = "all" ,cols = color)
  
  pdf(paste("./results/pc_clusters/", b,"_spatial_dormant_residuals.pdf",sep=""))
  print(p)
  dev.off()
  
  ##meta plots dormant signature
  df <- a@meta.data
  df_selected <- df %>% select(name, pc_clusters, signature_1_dormant, labels)
  
  df_grouped <- df_selected %>%
    group_by(name, pc_clusters, labels) %>%
    summarise(avg_signature_1_dormant = mean(signature_1_dormant, na.rm = TRUE))
  
  # Barplot
  plot <- ggplot(df_grouped, aes(x = pc_clusters, y = avg_signature_1_dormant, fill = labels)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(title = "Average 'signature_1_dormant' with respect to 'group_clusters' for each 'name'",
         x = "Group Clusters", 
         y = "Average Signature 1 Dormant", 
         fill = "Label") +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1)) +  # Adjust text size and angle here
    facet_wrap(~name)  # Separate plots by 'name'
  
  pdf(paste("./results/pc_clusters/", b,"_dormant_clusters_barplot.pdf",sep=""))
  print(plot)
  dev.off()
  
  #Violinplot
  df_selected <- df %>% select(name, pc_clusters, signature_1_dormant, labels)
  
  # Plot
  plot <- ggplot(df_selected, aes(x = pc_clusters, y = signature_1_dormant, fill = labels)) +
    geom_violin(scale = "width", adjust = 1) +
    labs(title = "'Signature_1_dormant' distributions with respect to 'group_clusters' for each 'name'",
         x = "Group Clusters", 
         y = "Signature 1 Dormant", 
         fill = "Label") +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1)) +  # Adjust text size and angle here
    facet_wrap(~name)  # Separate plots by 'name'
  
  pdf(paste("./results/pc_clusters/", b,"_dormant_clusters_violinplot.pdf",sep=""))
  print(plot)
  dev.off()
  
  # Plot
  plot <- ggplot(df_selected, aes(x = pc_clusters, y = signature_1_dormant, fill = labels)) +
    geom_boxplot() +
    labs(title = "'Signature_1_dormant' distributions with respect to 'group_clusters' for each 'name'",
         x = "Group Clusters", 
         y = "Signature 1 Dormant", 
         fill = "Label") +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1)) +  # Adjust text size and angle here
    facet_wrap(~name)  # Separate plots by 'name'
  
  pdf(paste("./results/pc_clusters/", b,"_dormant_clusters_boxplot.pdf",sep=""))
  print(plot)
  dev.off()
  
  ##meta plots dormant RESIDUALS signature
  df <- a@meta.data
  df_selected <- df %>% select(name, pc_clusters, residuals_dormant_ucellscore, labels)
  
  df_grouped <- df_selected %>%
    group_by(name, pc_clusters, labels) %>%
    summarise(avg_signature_1_dormant = mean(residuals_dormant_ucellscore, na.rm = TRUE))
  
  # Barplot
  plot <- ggplot(df_grouped, aes(x = pc_clusters, y = avg_signature_1_dormant, fill = labels)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(title = "Average 'signature_1_dormant' with respect to 'group_clusters' for each 'name'",
         x = "Group Clusters", 
         y = "Average Signature 1 Dormant", 
         fill = "Label") +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1)) +  # Adjust text size and angle here
    facet_wrap(~name)  # Separate plots by 'name'
  
  pdf(paste("./results/pc_clusters/", b,"_dormant_residuals_clusters_barplot.pdf",sep=""))
  print(plot)
  dev.off()
  
  #Violinplot
  df_selected <- df %>% select(name, pc_clusters, residuals_dormant_ucellscore, labels)
  
  # Plot
  plot <- ggplot(df_selected, aes(x = pc_clusters, y = residuals_dormant_ucellscore, fill = labels)) +
    geom_violin(scale = "width", adjust = 1) +
    labs(title = "'Signature_1_dormant' distributions with respect to 'group_clusters' for each 'name'",
         x = "Group Clusters", 
         y = "Signature 1 Dormant", 
         fill = "Label") +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1)) +  # Adjust text size and angle here
    facet_wrap(~name)  # Separate plots by 'name'
  
  pdf(paste("./results/pc_clusters/", b,"_dormant_residuals_clusters_violinplot.pdf",sep=""))
  print(plot)
  dev.off()
  
  # Plot
  plot <- ggplot(df_selected, aes(x = pc_clusters, y = residuals_dormant_ucellscore, fill = labels)) +
    geom_boxplot() +
    labs(title = "'Signature_1_dormant' distributions with respect to 'group_clusters' for each 'name'",
         x = "Group Clusters", 
         y = "Signature 1 Dormant", 
         fill = "Label") +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1)) +  # Adjust text size and angle here
    facet_wrap(~name)  # Separate plots by 'name'
  
  pdf(paste("./results/pc_clusters/", b,"_dormant_residuals_clusters_boxplot.pdf",sep=""))
  print(plot)
  dev.off()
  
  
  #########Each cluster plots
  p <- ST.FeaturePlot(a, features = "pc_clusters",  ncol =3,
                       split.labels = T, pt.size = 0.3) & theme(plot.title = element_blank(), strip.text = element_blank())
  
  pdf(paste("./results/pc_clusters/", b,"_each_cluster.pdf",sep=""))
  print(p)
  dev.off()
  
}


#######PLAY WITH THE METADATA FIN

###each cluster plot



#subset individual
M8_s@meta.data[["pc_clusters"]] <- as.factor(M8_s@meta.data[["pc_clusters"]])
Idents(object = M8_s) <- M8@meta.data[["pc_clusters"]]

pdf(file.path("./results/pc_clusters/M1_clusters.pdf"))
FeatureOverlay(M1_s, features = "pc_clusters", pt.size = 1.5)
dev.off()



##spatial plots
pdf(file.path("./results/pc_clusters/clusters.pdf"))
FeatureOverlay(subset, features = "seurat_clusters", ncols = 2,pt.size = 0.7)
dev.off()


p <- FeatureOverlay(subset, features = c("residuals_dormant_ucellscore"), pt.size = 1.8, 
                    value.scale = "all" ,cols = color)

pdf(paste("./results/pc_clusters/", b,"_spatial_Ucell_residuals.pdf",sep=""))
print(p)
dev.off()


