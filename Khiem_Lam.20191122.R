# R script for manuscript

rm(list = ls())  #to clean up before loading object

library(plyr)
library(stringr)
library(reshape2)
library(ggplot2)
library(Seurat)

options(stringsAsFactors = FALSE)

options(scipen=999)

### read in data

# hashtag label
dir <- 'hashtag2'
files <- list.files(dir)
files <- files[str_detect(files,'xls$')]
files <- files[str_detect(files,scaf)]

hto <- lapply(files, function(x) {
    file <- file.path(dir, x)
    df <- read.delim(file=file,sep="\t",header=TRUE,as.is=TRUE)
})
names(hto) <- files

# 10X data
scaf <- 'SCAF927|SCAF929'
dir <- '../data/CS025333/02_PrimaryAnalysisOutput/00_FullCellrangerOutputs'
files <- list.files(dir)
files <- files[str_detect(files,scaf)]

my.list <- lapply(files, function(x) {
    subdir <- paste(dir,x,'outs/filtered_feature_bc_matrix',sep='/')
    dat <- Read10X(data.dir = subdir)
    return(dat[[1]])
})
names(my.list) <- files

for (i in 1:length(my.list)) {
    my.list[[i]] <- CreateSeuratObject(counts = my.list[[i]],
                    project = names(my.list)[i])
    my.list[[i]][["percent.mt"]] <- PercentageFeatureSet(my.list[[i]], pattern = "^mt-")
    hashtag <- hto[[i]]$MULTI_ID
    names(hashtag) <- hto[[i]]$Barcode
    my.list[[i]][["hto"]] <- hashtag
    my.list[[i]] <- subset(my.list[[i]],
                    subset = percent.mt < 15 & hto != 'multi' & hto != 'neg')
    my.list[[i]] <- NormalizeData(my.list[[i]], verbose = FALSE)
    my.list[[i]] <- FindVariableFeatures(my.list[[i]], selection.method = "vst",
        nfeatures = 5000, verbose = FALSE)
}

reference.list <- my.list
my.anchors <- FindIntegrationAnchors(object.list = reference.list, 
                                     anchor.features = 5000,dims = 1:50)

my.integrated <- IntegrateData(anchorset = my.anchors, dims = 1:50)

DefaultAssay(my.integrated) <- "integrated"

my.integrated <- ScaleData(my.integrated, verbose = FALSE)
my.integrated <- RunPCA(my.integrated, npcs = 50, verbose = FALSE)
my.integrated <- RunUMAP(my.integrated, reduction = "pca", dims = 1:50)
my.integrated <- RunTSNE(my.integrated, reduction = "pca", dims = 1:50)

my.integrated <- FindNeighbors(my.integrated, dims = 1:10)
my.integrated <- FindClusters(my.integrated, resolution = 0.5)


### extended data Fig 2A

file <- 'paper/extended.data.Fig.2A.umap.by.cluster.pdf'
pdf(file)
DimPlot(my.integrated, reduction = "umap",label = TRUE, group.by = c("seurat_clusters"))
dev.off()

