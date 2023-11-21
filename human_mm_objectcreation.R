## SCRIPT: Human sample pc utility evaluation

## 24.05.23 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(ggplot2)
library(STutility)

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
  a <- NormalizeData(a)
  a <- ScaleData(a)
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
lista <- prueba
lista$se <- NULL

#save
saveRDS(lista,"./objects/sp/human/human_combined.rds")