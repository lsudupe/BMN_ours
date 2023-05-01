## SCRIPT: lot canonical markers BM project

## 19.03.23 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(ggplot2)
library(STutility)
library(dplyr)
library(UCell)

#Data---------------------------------
objects <- readRDS("./objects/sp/regress_out/list_regressout_ST.rds")
all <- readRDS("./objects/heterogeneity/se_hierarchical.rds")

names <- names(objects)
#Check dormant enrichment score
###Add module score
library(RColorBrewer)
color <- brewer.pal(11,"Spectral")
color <- rev(color)

#set colors
n = 7
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols = gg_color_hue(n)
cluster_color_map <- setNames(cols, unique(all@meta.data[["clustering"]]))

dormant <- c("C1qa", "Aif1" ,"Axl" ,"II18bp", "Glul", "Mpeg1", "H2-Eb1",
             "Nr1h3" ,"Lilrb4" ,"Sdc3" ,"Fcgr3" ,"Oas1a", "Ifit2", "Ly6a",
             "Ly86" ,"Bcl2a1b", "Cd19" ,"Gm4955", "Hebp1", "Sirpa" ,"Slfn2" ,"Tnfaip2" ,
             "Vpreb3", "ligp1", "Oasl2","am213b" ,"Slc43a2", "Ifi44", "Gng10" ,"Cxcl16" ,"Ptprcap" ,
             "Anxa2" ,"Rgs2" ,"Tmed3" ,"Igll1", "Hpgd", "Glipr1", "Cd4" ,"Cd84" ,"Gbp2", "AB124611", "Slc44a2" ,
             "Samd9l" ,"Oas1g" ,"Fcgr1", "Pla2g15", "Tifa" ,"Pmp22", "Abcc3" ,"S100a10")


new <- c()
for (i in 1:length(objects)){
  a <- objects[[i]]
  b <- names(objects[i])
  
  a <- AddModuleScore(a,
                      genes.pool = a@assays[["RNA"]]@data,
                      features = list(dormant),
                      name="dormant")

  meta <- a@meta.data
  
  lm <- lm(meta$dormant1 ~ meta$MM_MIC, data =meta)
  residuals <- lm$residuals
  a@meta.data[["residuals_dormant_addmodulescore"]] <- residuals
  
  ##plots
  p <- FeatureOverlay(a, features = c("dormant1"), pt.size = 1.8, 
                      value.scale = "all" ,cols = color)
  
  pdf(paste("./results/ST/dormant/", b,"_spatial_modulescore.pdf",sep=""))
  print(p)
  dev.off()
  
  p <- a@meta.data%>% 
    ggplot(aes(x=dormant1, y= clustering, fill=clustering)) + 
    geom_boxplot(aes(fill=clustering)) +  
    scale_fill_manual(values =cluster_color_map ) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    #xlim(-0.5, 0.5) +
    theme(plot.title = element_text(hjust=0.5, face="bold")) 
  
  pdf(paste("./results/ST/dormant/", b,"_boxplot_modulescore.pdf",sep=""))
  print(p)
  dev.off()
  
  p <- FeatureOverlay(a, features = c("residuals_dormant_addmodulescore"), pt.size = 1.8, 
                      value.scale = "all" ,cols = color)
  
  pdf(paste("./results/ST/dormant/", b,"_spatial_modulescore_residuals.pdf",sep=""))
  print(p)
  dev.off()
  
  p <- a@meta.data%>% 
    ggplot(aes(x=residuals_dormant_addmodulescore, y= clustering, fill=clustering)) + 
    geom_boxplot(aes(fill=clustering)) +  
    scale_fill_manual(values =cluster_color_map ) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    #xlim(-0.5, 0.5) +
    theme(plot.title = element_text(hjust=0.5, face="bold")) 
  
  pdf(paste("./results/ST/dormant/", b,"_boxplot_modulescore_residuals.pdf",sep=""))
  print(p)
  dev.off()
  
  ## Add UCellScore
  
  vector<- ScoreSignatures_UCell(a@assays[["RNA"]]@data, features = list(dormant))
  a@meta.data[["signature_1_dormant"]] <- as.vector(vector)
  
  meta <- a@meta.data
  
  lm <- lm(meta$signature_1_dormant ~ meta$MM_MIC, data =meta)
  residuals <- lm$residuals
  a@meta.data[["residuals_dormant_ucellscore"]] <- residuals
  
  ##plots
  p <- FeatureOverlay(a, features = c("signature_1_dormant"), pt.size = 1.8, 
                      value.scale = "all" ,cols = color)
  
  pdf(paste("./results/ST/dormant/", b,"_spatial_Ucell.pdf",sep=""))
  print(p)
  dev.off()
  
  p <- FeatureOverlay(a, features = c("residuals_dormant_ucellscore"), pt.size = 1.8, 
                      value.scale = "all" ,cols = color)
  
  pdf(paste("./results/ST/dormant/", b,"_spatial_Ucell_residuals.pdf",sep=""))
  print(p)
  dev.off()
  
  p <- a@meta.data%>% 
    ggplot(aes(x=signature_1_dormant, y= clustering, fill=clustering)) + 
    geom_boxplot(aes(fill=clustering)) +  
    scale_fill_manual(values =cluster_color_map ) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    #xlim(-0.5, 0.5) +
    theme(plot.title = element_text(hjust=0.5, face="bold")) 
  
  pdf(paste("./results/ST/dormant/", b,"_boxplot_Ucell.pdf",sep=""))
  print(p)
  dev.off()
  
  p <- a@meta.data%>% 
    ggplot(aes(x=residuals_dormant_ucellscore, y= clustering, fill=clustering)) + 
    geom_boxplot(aes(fill=clustering)) +  
    scale_fill_manual(values =cluster_color_map ) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    #xlim(-0.5, 0.5) +
    theme(plot.title = element_text(hjust=0.5, face="bold")) 
  
  pdf(paste("./results/ST/dormant/", b,"_boxplot_Ucell_residuals.pdf",sep=""))
  print(p)
  dev.off()
  
  
  ## add object to list
  new[[length(new) + 1]] <- a
}

names(new) <- names


library(reshape2)
library(ggplot2)
melted <- reshape2::melt(u.scores)
colnames(melted) <- c("Cell","Signature","UCell_score")
p <- ggplot(melted, aes(x=Signature, y=UCell_score)) + 
  geom_violin(aes(fill=Signature), scale = "width") +
  geom_boxplot(width=0.1, outlier.size=0) +
  theme_bw() + theme(axis.text.x=element_blank())
p


##Plots
custom_theme <- theme(legend.position = c(0.45, 0.8), # Move color legend to top
                      legend.direction = "horizontal", # Flip legend
                      legend.text = element_text(angle = 30, hjust = 1), # rotate legend axis text
                      strip.text = element_blank(), # remove strip text
                      plot.title = element_blank(), # remove plot title
                      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm")) # remove plot margins

p <-FeatureOverlay(se, features = c("Retnlg"), pt.size = 0.7,  
               sampleids = 3:4, ncols = 1, cols = c("darkblue", "cyan", "yellow", "red", "darkred"))
p <-FeatureOverlay(se, features = c("Mmp8"), pt.size = 0.7,  
                   sampleids = 3, ncols = 2, cols = c("darkblue", "cyan", "yellow", "red", "darkred"))
p <- ST.FeaturePlot(se, features = c("Adpgk"),indices  = 4, pt.size = 1.5,ncol = 1, palette = "Spectral")

pdf(file.path("./results/ST/canonical/dentrAdpgk_healthy4_mouse.pdf"))
print(p)
dev.off()

FeatureOverlay(se, features = c("Retnlg","Mmp8"), ncols = 2, sampleids = 3:4, cols = c("Spectral"), pt.size = 0.7)
FeatureOverlay(se.subset, features = c("Retnlg","Mmp8"), pt.size = 0.7,  
               sampleids = 1:2, ncols = 3, cols = c("darkblue", "cyan", "yellow", "red", "darkred"))