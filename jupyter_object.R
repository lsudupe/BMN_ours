## SCRIPT: Save data to open in Jupyter, spatial data BMN project

## 18.09.22 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(dplyr)
library(ggplot2) 
library('magrittr')
library(tidyverse)
library(base)

#Data--------------------------------------
M1_tib_1A  <- readRDS("./objects/sp/M1_tib_1A_normalized.rds")
x <- M1_tib_1A

coordinates <- x@images[["M1_tib_1A"]]@coordinates
meta <- x@meta.data
counts <- t(as.data.frame(x@assays[["SCT"]]@counts))

##save the data
write.csv(coordinates , './coordinates.csv', row.names = TRUE)
write.csv(meta , './meta.csv', row.names = TRUE)
write.csv(counts , './counts.csv', row.names = TRUE)
..


