## SCRIPT: Plot gradients results BM project

## 19.03.23 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(ggplot2)
library(STutility)
library(dplyr)
library(hrbrthemes)
library(extrafont)
library(ggthemes)
library(viridis)




#Data---------------------------------
se <- readRDS("./objects/sc/integrated/se_deco.rds")

library(RColorBrewer)
display.brewer.all()
color <- brewer.pal(11,"Spectral")
color <- rev(color)

p <-FeatureOverlay(se, features = c("MM_MIC"), sampleids = 1:6,pt.size = 0.7,ncol = 2 , 
                   value.scale = "all" ,cols = color)
p <-FeatureOverlay(se, features = c("MM_MIC"), sampleids = 1:6,pt.size = 0.7,ncol = 2 , 
                   value.scale = "all" ,cols = c("lightgray", "mistyrose", "red", "dark red", "black"))

p <- ST.FeaturePlot(se, features = c("MM_MIC"),indices  = 1:6, pt.size = 0.80,ncol = 2 , grid.ncol = 1, palette = "Spectral")
p <- ST.FeaturePlot(se, features = c("MM_MIC"),indices  = 1:6, pt.size = 0.80,ncol = 2 , grid.ncol = 1, cols = c("lightgray", "mistyrose", "red", "dark red", "black"))



pdf(file.path("./results/ST/gradient/gradient_rojo.pdf"))
print(p)
dev.off()

meta <- se@meta.data

## density plots
# Make the histogram
M1_fem_1C <- meta[grepl("M1_fem_1C", meta$name),]
M2_F_2B <- meta[grepl("M2_F_2B", meta$name),]
M3_F_1C <- meta[grepl("M3_F_1C", meta$name),]
M3_fem_1C <- meta[grepl("M3_fem_1C", meta$name),]
M8_F2_1C <- meta[grepl("M8_F2_1C", meta$name),]
M9_F2_1C <- meta[grepl("M9_F2_1C", meta$name),]

pdf(file.path("./results/ST/gradient/density_together.pdf"))
print(ggplot(data=meta, aes(x=MM_MIC, group=name, fill=name)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  theme_ipsum()+
  ylab("") +
  xlab("MM values (%)"))
dev.off()

# Using Small multiple
pdf(file.path("./results/ST/gradient/density_separate.pdf"))
ggplot(data=meta, aes(x=MM_MIC, group=name, fill=name)) +
  geom_density(adjust=1.5) +
  theme_ipsum() +
  facet_wrap(~name) +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    axis.ticks.x=element_blank()
  )
dev.off()

##add list values to df
for (i in 1:length(lista)){
  a <- lista[[i]]
  pdf(paste("./results/ST/gradient/", names(lista[i]),"_density.pdf",sep=""))
  a %>%
    ggplot( aes(x=MM_MIC)) +
    geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8) +
    ylab("") +
    xlab("MM values (%)")
  dev.off()
}

