#!/usr/bin/env Rscript
library(tidyverse)
library(Seurat)
library(scran)
library(scater)
library(SingleCellExperiment)
library(RColorBrewer)
library(cowplot)

setwd("~/INMG_SingleCell/")
resu = "results/integratedD0/DCscran/"
seu <- readRDS("rds/integratedD0/integrated_seu_fitsne.rds")
veccell <- readRDS("rds/integratedD0/integr_celltype_Vector.rds")
seu <- RenameIdents(seu,veccell)
seu@meta.data$celltype = seu@active.ident
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

print("Finding doublets")
sce <- computeSumFactors(sce) # by scran: scaling normalization (implements deconvolution)
sce <- logNormCounts(sce)
dbl_dens <- doubletCells(sce)
sce$doublet_score <- 0
sce$doublet_score <- log10(dbl_dens + 1)

sce <- runTSNE(sce, perplexity=150, PCA=T, num_threads=4)

f <- scater::plotColData(sce, x="sum",y="detected",colour_by="doublet_score")
t <- plotTSNE(sce, colour_by="doublet_score")
o <- plotTSNE(sce, colour_by="orig.ident")
g <- plotTSNE(sce, colour_by="celltype")
pdf("results/integratedD0/DCscran/tsne_sce.pdf", width=12)
plot_grid(f,t,o,g, ncol=2)
dev.off()

saveRDS(sce, "rds/doubletsD0/mysce_doublets.rds")

# sce <- readRDS("rds/doubletsD0/mysce_doublets.rds")
# transfer FItSNE:
head(sce@int_colData@listData[["reducedDims"]]@listData[["TSNE"]])
head(seu@reductions[["tsne"]]@cell.embeddings)
typeof(sce@int_colData@listData[["reducedDims"]]@listData[["TSNE"]])
typeof(seu@reductions[["tsne"]]@cell.embeddings)


# quqntiles:
head(as.numeric(as.character(seu@meta.data$doublet_score)))
x <- as.numeric(as.character(seu@meta.data$doublet_score))

quantile(x,0.995)
# 99.5% 
# 4.942606 
quantile(x,0.95)
# 95% 
# 3.583153 

bicho = as.data.frame(colData(sce))
bicho <-  bicho %>% mutate(doubinterval = case_when(
  doublet_score <= 3.5 ~ "1. <= 3.5",
  doublet_score > 3.5 & doublet_score < 5 ~ "2. 3.5 - 5",
  TRUE ~ "3. >= 5" # doublet_score > 4.5
))
colData(sce)$doubinterval = as.factor(bicho$doubinterval)

# visualization not very good:
a <- plotReducedDim(sce, dimred="FItSNE", colour_by = "celltype")
b <- plotReducedDim(sce, dimred="FItSNE", colour_by = "doulet_score")
c <- plotReducedDim(sce, dimred= "FItSNE", colour_by = "doubinterval")

thecols = c(rgb(0.9, 0.2, 0.3, 1),
            rgb(0.2, 0.5, 0.8, 0.5), 
            rgb(0.1, 0.7, 0.3, 0.1))

seu@meta.data$doublet_score = as.factor(colData(sce)$doublet_score)
seu@meta.data$doubinterval = as.factor(colData(sce)$doubinterval)
levels(seu@meta.data$doubinterval) = c( "<= 3.5", "3.5 - 5", ">= 5" )
b <- plotReducedDim(sce, dimred="FItSNE", colour_by = "doublet_score")
l <- qplot(colData(sce)$doublet_score)
m <- DimPlot(seu, group.by = "doubinterval", cols=rev(thecols))
n <- DimPlot(seu, group.by = "celltype", cols=
               colorRampPalette(brewer.pal(8,"Set2"))(length(levels(seu$celltype))),
             label=T, repel=T, size=2)+NoLegend()


pdf(paste0(resu,"scranMeth.pdf"), width=11)
plot_grid(m,n,l,b, ncol=2)
dev.off()