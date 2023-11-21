## SCRIPT: Human correlation analysis

## 21.11.23 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(ggplot2)


## Data
healthy <- readRDS("./objects/sp/integrated/healthy_clean.rds")
mm <- readRDS("./objects/sp/human/human_combined.rds")

all <- c(healthy, mm)

## Signatures
Tcell35genes <- c("Pdcd1", "Ctla4", "Tnfrsf9", "Havcr2", "Tox", "Tigit", "Wars", "Rsad2",
                  "Mcm7", "Mx1", "Ndfip2", "Enosf1", "Ccdc141", "Stmn1", "Ttn", "Faslg",
                  "Mcm5", "Nab1", "Phlda1", "Mcm3", "Pcna", "Gapdh", "Oasl", "fi44l",
                  "Tbc1d4", "Slc43a3", "Pam", "Ccl3", "Acp5", "Oas3", "Cd38", "Tnfs10",
                  "Gbp2", "Kif20b", "Ctsb")
Tcell_exh_Isa_human <- c("HAVCR2","CTLA4","PDCD1","LAG3", "LAYN", "PD1", "TIGIT", "TIM3")
cytotoxicity_human <- c("NKG7", "CCL4", "CST7", "PRF1", "GZMA", "GZMB", "IFNG", "CCL3")
Teff_Juanro_human <- c("ALOX5AP","ANXA1","ARL4C","ATP6V1G1","B2M","BIN2","CD3D","CD48","CD52","CTSC","CTSW",
                       "CXCR6","CYBA","EMP3","EOMES","EVL","FTL","FXYD5","GBP5","GRAP2","GZMA","GZMB","HCST",
                       "IFITM2","IFNG","IL32","IL4R","IRF1","ISG15","ISG20","ITM2B","LGALS3","LST1","LTB",
                       "MALAT1","MT-CO2","NCR3","PDLIM2","PRF1","RABAC1","RPL27A","RPS27","S100A4","S100A6",
                       "SAMHD1","SDF2L1","SH3BGRL3","SLFN5","STK17B","TC2N","TMSB4X","TPT1","TRAC","XBP1")
TcellAct_Juanro_human <- c("IFNG","CCL4","CCL3","IL3","XCL1","CSF2","GZMB","FABP5","XCL2","LTA","LAG3",
                           "MIR155HG","TNFRSF4","TNFRSF9","PIM3","CIITA","CD74","HLA-DRB1","HLA-DRB5",
                           "HLA-DMB","HLA-DPB1","HLA-DQA2","HLA-DOA","HLA-DRA","CD25","CD69","CCL5","SP140",
                           "TIGIT","CD247","CCR5","NLRC3","PTPRCAP","STAT4","IL12RB1","TRAT1","SLA2","CXCR3",
                           "ZAP70","NKG7","SIRPG","ICOS","IL18R1","SLAMF7","PTPN7","ITK","CRTAM","GZMA",
                           "CD3E","GPR171","TARP","CD3D","LCK","SLAMF1","CD3G")

## Enrichment score
for (i in 1:length(all)){
  a <- all[[i]]
  b <- names(all[i])
  
  for (cluster_name in cluster_names) {
    # Define the variable names dynamically based on the cluster_name
    signature_var <- paste0("signature_1_", cluster_name)
    
    file_name <- paste0("./results/human/cluster_enrich/new/25/", b, "_", cluster_name, "_pc.pdf", sep = "")
    
    # Add UCellScore for the current cluster
    vector <- ScoreSignatures_UCell(a@assays[["RNA"]]@counts, features = list(get(cluster_name)))
    a@meta.data[[signature_var]] <- as.vector(vector)
    
    #c <- c(min(a@meta.data[[signature_var]]), max(a@meta.data[[signature_var]]))
    # Create the plot
    p <- FeatureOverlay(a, features = c(signature_var), ncols = 1, pt.size = 1.1, 
                        value.scale = "all" ,cols = color) +
      scale_fill_gradientn(colours = color,
                           breaks = c(0.0, 0.5),
                           labels = c("Min", "Max"),
                           limits = c(0.0, 0.5))
    
    # Save the plot to a PDF file
    pdf(file_name)
    print(p)
    dev.off()
    
  }
  
}




ggplot(subset(x = seurat@meta.data, name %in% c("M1_fem_1C","M2_F_2B","M8_F2_1C","M9_F2_1C")), 
       aes(x=PC_new_UCell , y=Tcell_Exh_UCell)) + geom_point(mapping = aes(color = name)) + 
  cowplot::theme_cowplot() + stat_cor(method = "spearman", cor.coef.name = "rho")


