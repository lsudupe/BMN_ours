## SCRIPT: Dormant genes individual plot  BM project

## 01.06.23 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(ggplot2)
library(STutility)
library(dplyr)
library(UCell)

#Data---------------------------------
objects <- readRDS("./objects/sp/regress_out/list_regressout_ST.rds")
all <- readRDS("./objects/heterogeneity/se_hierarchical.rds")

##individual
M1_s <- readRDS("./objects/pc_clusters/M1_s_dormant.rds")
M2_s <- readRDS("./objects/pc_clusters/M2_s_dormant.rds")
M8_s <- readRDS("./objects/pc_clusters/M8_s_dormant.rds")
M9_s <- readRDS("./objects/pc_clusters/M9_s_dormant.rds")


dormant <- c("C1qa", "Aif1" ,"Axl" ,"II18bp", "Glul", "Mpeg1", "H2-Eb1",
             "Nr1h3" ,"Lilrb4" ,"Sdc3" ,"Fcgr3" ,"Oas1a", "Ifit2", "Ly6a",
             "Ly86" ,"Bcl2a1b", "Cd19" ,"Gm4955", "Hebp1", "Sirpa" ,"Slfn2" ,"Tnfaip2" ,
             "Vpreb3", "ligp1", "Oasl2","am213b" ,"Slc43a2", "Ifi44", "Gng10" ,"Cxcl16" ,"Ptprcap" ,
             "Anxa2" ,"Rgs2" ,"Tmed3" ,"Igll1", "Hpgd", "Glipr1", "Cd4" ,"Cd84" ,"Gbp2", "AB124611", "Slc44a2" ,
             "Samd9l" ,"Oas1g" ,"Fcgr1", "Pla2g15", "Tifa" ,"Pmp22", "Abcc3" ,"S100a10")

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



subsets <- c(M1_subset,M2_subset, M8_subset, M9_subset)
names(subsets) <- c("M1_subset","M2_subset", "M8_subset", "M9_subset")
