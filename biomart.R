## SCRIPT: Biomart for BM project

## 07.02.23 Laura Sudupe , git @lsudupe

#Libraries---------------------------------

library(biomaRt)
mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
ig_genes<-getBM(attributes = c('external_gene_name','gene_biotype'), 
                mart = mart)
ig_genes<-ig_genes[grep(c("^IG"),ig_genes$gene_biotype),]


