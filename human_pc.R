## SCRIPT: Human sample pc utility evaluation

## 24.05.23 Laura Sudupe , git @lsudupe

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

#pc_markers
pc_mm_markers <- c("SHH", "DHH", "IHH", "PTCH1", "PTCH2", "SMO", "SUFU", "GLI1", 
                   "GLI2", "GLI3", "CD19", "CD44", "CXCR4", "KLF4", "CD28", "CD33", "CD27")

human_tcell_exhausted <- c("Pdcd1", "Ctla4", "Tnfrsf9", "Havcr2", "Tox", "Tigit", "Wars", "Rsad2",
                           "Mcm7", "Mx1", "Ndfip2", "Enosf1", "Ccdc141", "Stmn1", "Ttn", "Faslg",
                            "Mcm5", "Nab1", "Phlda1", "Mcm3", "Pcna", "Gapdh", "Oasl", "Ifi44l",
                            "Tbc1d4", "Slc43a3", "Pam", "Ccl3", "Acp5", "Oas3", "Cd38", "Tnfsf10",
                             "Gbp2", "Kif20b", "Ctsb")

neutro <- c(c("Rpl11", "Nbpf10", "Hist2h2aa3", "Hist2h2ac", "Rps27", "Zc3h11a", "Rps27a", "Nbeal1", 
              "Rpl32", "Rpl34", "Matr3", "Rbm27", "Rps14", "Rack1", "Sox4", "Hist1h1c", "Hist1h4c", 
              "Hist1h1e", "Hist1h1d", "Hist1h2bk", "Hist1h2aj", "Hist1h2am", "Rps18", "Eef1a1", 
              "Npy", "Polr2j3", "C7orf55-luc7l2", "Rpl10", "Defa4", "Ms4a3", "Malat1", "Atp5mg", 
              "Rpl41", "Hnrnpa1l2", "Rnase3", "Ctsg", "Cox16", "B2m", "Rplp1", "Rps15a", "Aldoa",
              "Rpl13", "Eif4a1", "Nme1-nme2", "Mpo", "Serpinb10", "Atp5f1e", "Azu1", "Prtn3", "Elane"))

#Variables---------------------------------
DIR_ROOT <- file.path(getwd())
DIR_DATA <- file.path(DIR_ROOT, "data/data/")

# samples
all <- dir(path = DIR_DATA)
#samples <- samples[ samples %in% c("BM_human_AP-B08805")]

samples <- list.files(path = DIR_DATA, recursive = TRUE, pattern = "filtered",full.names = TRUE)
imgs <- list.files(path = DIR_DATA, pattern = "tissue_hires_image.png", recursive = TRUE, full.names = TRUE)
spotfiles <- list.files(path = DIR_DATA, pattern = "tissue_positions_list.csv", recursive = TRUE, full.names = TRUE)
json <- list.files(path = DIR_DATA, pattern = "scalefactors_json.json", recursive = TRUE, full.names = TRUE)
name <- dir(path = DIR_DATA)

infoTable <- data.frame(samples = samples, imgs = imgs, spotfiles = spotfiles, json = json, 
                        name = name)


infoTable <- subset(infoTable, name != "M1_fem_1C" &
                    name != "M1_tib_1A" & 
                      name != "M2_F_2B" &
                      name != "M3_F_1C" &
                      name != "M3_fem_1C" &
                      name != "M3_tib_2A" & 
                      name != "M8_F2_1C" & 
                      name !="M9_F2_1C")

se <- InputFromTable(infoTable,
                     platform =  "Visium")

st.object <- GetStaffli(se)
st.object
se <- LoadImages(se, time.resolve = FALSE)
saveRDS(se, "./objects/sp/integrated/se_human.rds")
se <- readRDS("./objects/sp/integrated/se_human.rds")

##divide by sample
Idents(object = se) <- "name"
name <- unique(se@meta.data[["name"]])
objects <- c()

for (i in name){
  a <- SubsetSTData(se, idents = i)
  a <- SCTransform(a)
  ## add object to list
  objects[[length(objects) + 1]] <- a
  
}
names(objects) <- name
## separate list
list2env(objects,envir=.GlobalEnv)

#####change some of the objects
df2 <- `BM_human_AP-B08041_`@meta.data
df <-  read_csv("./data/data/BM_human_AP-B08041_/spatial/on_tissue.csv")
# Replace NA in the 'on_tissue' column with 'off'
df$on_tissue[is.na(df$on_tissue)] <- "off"
# Convert rownames to a column 
df <- df %>%
  rownames_to_column(var = "index")
df$index <- df$Barcode
# Add "_7" to each index in df
df$index <- paste(df$index, "_7", sep = "")
# Suppose df2 is your second data frame
# Convert rownames to a column in df2
df2 <- df2 %>%
  rownames_to_column(var = "index")
# Add 'on_tissue' column to df2 based on index match and in the order of df2
df2 <- df2 %>%
  left_join(df[, c("index", "on_tissue")], by = "index") %>%
  filter(!is.na(on_tissue))
# In case you have NAs in df2 for 'on_tissue', replace them with 'off'
df2$on_tissue[is.na(df2$on_tissue)] <- "off"
df2 <- column_to_rownames(df2, var = "index")
`BM_human_AP-B08041_`@meta.data <- df2
#subset
Idents(object = `BM_human_AP-B08041_`) <- `BM_human_AP-B08041_`@meta.data[["on_tissue"]]
`BM_human_AP-B08041_` <- SubsetSTData(object = `BM_human_AP-B08041_`, ident = c("on"))

#####fin

prueba <- c(`BM_human_AP-B00182_`, `BM_human_AP-B02149_`,`BM_human_AP-B08041_`,`BM_human_AP-B08805`,
            `BM_B000943`,`BM_B01320`,`BM_B02817`,`BM_B10395`, se)
names(prueba) <- c("BM_human_AP-B00182_", "BM_human_AP-B02149_", "BM_human_AP-B08041_", "BM_human_AP-B08805",
                   "BM_B000943", "BM_B01320", "BM_B02817", "BM_B10395", "se")


###Enrichment score
post <- c()
for (i in 1:length(prueba)){
  a <- prueba[[i]]
  b <- names(prueba[i])
  
  ## Add UCellScore PC
  vector<- ScoreSignatures_UCell(a@assays[["RNA"]]@counts, features = list(pc_mm_markers))
  a@meta.data[["signature_1_pc"]] <- as.vector(vector)
  
  color <- brewer.pal(11,"Spectral")
  color <- rev(color)
  
  p <- FeatureOverlay(a, features = c("signature_1_pc"), ncols = 1, pt.size = 1.5, 
                      value.scale = "all" ,cols = color)
  
  pdf(paste("./results/human/pc_enrich/", b,"_pc.pdf",sep=""))
  print(p)
  dev.off()
  
  ## Add UCellScore neutro
  vector<- ScoreSignatures_UCell(a@assays[["RNA"]]@counts, features = list(neutro))
  a@meta.data[["signature_1_neutro"]] <- as.vector(vector)
  
  p <- FeatureOverlay(a, features = c("signature_1_neutro"), ncols = 1, pt.size = 1.5, 
                      value.scale = "all" ,cols = color)
  
  pdf(paste("./results/human/neutro_enrich/", b,"_pc.pdf",sep=""))
  print(p)
  dev.off()
  
  ## Add UCellScore Tcell
  vector<- ScoreSignatures_UCell(a@assays[["RNA"]]@counts, features = list(human_tcell_exhausted))
  a@meta.data[["signature_1_tcell_ex"]] <- as.vector(vector)
  
  p <- FeatureOverlay(a, features = c("signature_1_tcell_ex"), ncols = 1, pt.size = 1.5, 
                      value.scale = "all" ,cols = color)
  
  pdf(paste("./results/human/tcell_ex_enrich/", b,"_tcell_ex.pdf",sep=""))
  print(p)
  dev.off()
  
  post[[length(post) + 1]] <- a
}

names(post) <- c("BM_human_AP-B00182_", "BM_human_AP-B02149_", "BM_human_AP-B08041_", "BM_human_AP-B08805",
                   "BM_B000943", "BM_B01320", "BM_B02817", "BM_B10395", "se")
se <- post[["se"]]

color <- brewer.pal(11,"Spectral")
color <- rev(color)

#pc

a <- FeatureOverlay(se, features = c("signature_1_pc"),sampleids = 1, pt.size = 0.5) +
  scale_fill_gradientn(colours = color,
                      breaks = c(0.0,0.25),
                      labels = c("Min", "Max"),
                      limits = c(0.0,0.25))
b <- FeatureOverlay(se, features = c("signature_1_pc"),sampleids = 2, pt.size = 0.5) +
  scale_fill_gradientn(colours = color,
                      breaks = c(0.0,0.25),
                      labels = c("Min", "Max"),
                      limits = c(0.0,0.25))
c <- FeatureOverlay(se, features = c("signature_1_pc"),sampleids = 3, pt.size = 0.5) +
  scale_fill_gradientn(colours = color,
                      breaks = c(0.0,0.25),
                      labels = c("Min", "Max"),
                      limits = c(0.0,0.25))
d <- FeatureOverlay(se, features = c("signature_1_pc"),sampleids = 4, pt.size = 0.5) +
  scale_fill_gradientn(colours = color,
                      breaks = c(0.0,0.25),
                      labels = c("Min", "Max"),
                      limits = c(0.0,0.25))
f <- FeatureOverlay(se, features = c("signature_1_pc"),sampleids = 5, pt.size = 0.5) +
  scale_fill_gradientn(colours = color,
                      breaks = c(0.0,0.25),
                      labels = c("Min", "Max"),
                      limits = c(0.0,0.25))
g <- FeatureOverlay(se, features = c("signature_1_pc"),sampleids = 6, pt.size = 0.5) +
  scale_fill_gradientn(colours = color,
                      breaks = c(0.0,0.25),
                      labels = c("Min", "Max"),
                      limits = c(0.0,0.25))
h <- FeatureOverlay(se, features = c("signature_1_pc"),sampleids = 7, pt.size = 0.5) +
  scale_fill_gradientn(colours = color,
                      breaks = c(0.0,0.25),
                      labels = c("Min", "Max"),
                      limits = c(0.0,0.25))
i <- FeatureOverlay(se, features = c("signature_1_pc"),sampleids = 8, pt.size = 0.5) +
  scale_fill_gradientn(colours = color,
                      breaks = c(0.0,0.25),
                      labels = c("Min", "Max"),
                      limits = c(0.0,0.25))



pdf(paste("./results/human/pc_enrich/se_pc_1.pdf",sep=""))
print(grid.arrange(a,b,c,d, ncol=2))
dev.off()
pdf(paste("./results/human/pc_enrich/se_pc_2.pdf",sep=""))
print(grid.arrange(f,g,h,i, ncol=2))
dev.off()

#neutro
FeatureOverlay(se, features = c("signature_1_neutro"), pt.size = 1) +
  scale_fill_gradient(colours = color,
                      breaks = c(0.0,2.5),
                      labels = c("Min", "Max"),
                      limits = c(0.0,2.5))

pdf(paste("./results/human/neutro_enrich/se_neutro.pdf",sep=""))
print(p)
dev.off()

#tcell
FeatureOverlay(se, features = c("signature_1_tcell_ex"),sampleids = 1:8, pt.size = 1) +
  scale_fill_gradientn(colours = color,
                      breaks = c(0.0,2.5),
                      labels = c("Min", "Max"),
                      limits = c(0.0,2.5))

pdf(paste("./results/human/tcell_ex_enrich/se_tcell_ex.pdf",sep=""))
print(p)
dev.off()



















