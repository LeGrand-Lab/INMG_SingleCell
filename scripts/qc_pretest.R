#!/usr/bin/env Rscript
library(tidyverse)
library(Seurat)
library(scran)
library(scater)
library(SingleCellExperiment)

setwd("~/INMG_SingleCell/")
resu = "/results/intgratedD0/"
seu <- readRDS("rds/integratedD0/integrated_seu_fitsne.rds")
veccell <- readRDS("rds/integratedD0/integr_celltype_Vector.rds")
seu <- RenameIdents(seu,veccell)

# do not use as.SingleCellExperiment, instead create from scratch:
dgcmatrix <- seu@assays[["RNA"]]@counts
head(colnames(dgcmatrix)) 
head(rownames(dgcmatrix))
sce = SingleCellExperiment(assays=list(counts=dgcmatrix))

rowData(sce) <- DataFrame(genes_names = rownames(sce))

head(colData(sce))
if ((all(rownames(colData(sce))==rownames(seu@meta.data)))){
  colData(sce) = DataFrame(seu@meta.data)
}

rowData(sce)$expressed <- scater::nexprs(sce,byrow=TRUE)>0
per_cell <- perCellQCMetrics(sce[rowData(sce)$expressed,])
colData(sce) <- cbind(colData(sce),per_cell)
head(rowData(sce))
colData(sce)$keep_total <- scater::isOutlier(colData(sce)$sum,type = "lower", log=TRUE)
table(colData(sce)$keep_total) # TRUE are OUTLIERS
sce <- scater::addPerFeatureQC(sce)

head(rowData(sce))

plotdet <- scater::plotColData(sce,x="sum",y="detected")
pdf(paste0(resu,"plotColData_sumVSdetected.pdf"))
plotdet
dev.off()

print("Finding doublets")
sce <- computeSumFactors(sce)

sce <- logNormCounts(sce)

dbl_dens <- doubletCells(sce)
sce$doublet_score <- 0
sce$doublet_score <- log10(dbl_dens + 1)

pdf(paste0(resu, "tsne_doublets.pdf"))
plotTSNE(sce, colour_by="doublet_score")
dev.off()

pdf(paste0(resu, "feat_doublets.pdf"))
detfeat <- scater::plotColData(sce, x="sum",y="detected",colour_by="doublet_score")
dev.off()

saveRDS(sce, "mysce_doublets.rds")