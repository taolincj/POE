###--------------------cluster analysis: which cell type is associated with trait?------------------
# source activate senic_env
library(Seurat)
library(data.table)
library(dplyr)
library(reshape2)
library(tibble)
library(tidyr)
library(ggplot2)
setwd("/data_group/xiehaibing/xiehaibing2/magma/ljbscrna")
count2seurat = function(genecount, resolution){
  jz = read.csv(genecount, sep = "\t")
  jz = jz[,-1]
  jz <- as(as.matrix(jz), "dgCMatrix")
  pbmc <- CreateSeuratObject(counts = jz, min.cells = 3, min.features = 200)
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
  pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 25000 & percent.mt < 5)
  pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
  pbmc <- NormalizeData(pbmc)
  pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(pbmc)
  pbmc <- ScaleData(pbmc, features = all.genes)
  pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
  pbmc <- FindNeighbors(pbmc, dims = 1:40)
  pbmc <- FindClusters(pbmc, resolution = resolution)
  pbmc <- RunUMAP(pbmc, dims = 1:40)
  return(pbmc)
}

rib = count2seurat("version11_expr_matrix_allRib22", resolution=1.55)
saveRDS(rib, "rib.rds")