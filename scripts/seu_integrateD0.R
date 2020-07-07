#!/usr/bin/env Rscript
# INTEGRATED ANALYSIS of three datasets, in homeostasis
# following instructions from:
# https://satijalab.org/seurat/pancreas_integration_label_transfer.html
# but adapting to our datasets
# --
# JohaGL

library(dplyr)
library(Seurat)
library(ggplot2)
library(sctransform)
library(cowplot)
library(RColorBrewer)
source(file="~/INMG_SingleCell/scripts/functions_stock.R",local=T)

prloc="~/INMG_SingleCell/"
setwd(prloc)
resu="results/integratedD0/"
rdsdir = "rds/integratedD0/"
dir.create(resu,recursive=T)
dir.create(rdsdir, recursive=T)
nbDIM=13 
nft = 2500
resolution=0.5
print(paste0("Will use *",nbDIM,"* princ components (dims), as we already runned this script", 
     "accordingly to elbowplot and good balance identities-nbclusters"))
print(" ")
print(paste0("this integration uses subset nFeature_RNA >",nft, 
   "and mitochondrial features < 5 % ==> this is stringent"))
 
# open files
# --*1*--

list.d.seu <- sapply(c("data/DellOrsoD0/dorsowt1","data/DellOrsoD0/dorsowt2"),
  function(d){
  tmp <- Read10X(data.dir=d)
  seu <- CreateSeuratObject(tmp, project="DellOrso",min.cells = 3,min.features = 200)
  return(seu) })
dorso <- merge(list.d.seu[[1]], y=list.d.seu[[2]], 
               add.cells.ids=c("wt1","wt2"),project="DellOrso")
rm(list.d.seu)

# --*2*--

list.g.seu <- sapply(c("data/GiordaniD0/GSM3520458_20171018_uninjured_wt_filtered.csv",
                       "data/GiordaniD0/GSM3520459_20180917_uninjured_wt_filtered.csv"), 
                     function(g){
    g = read.csv(g, sep=",", header=TRUE, row.names=1)
  seu <- CreateSeuratObject(g, project="Giordani",min.cells = 3,min.features = 200)
  return(seu) })
gio <- merge(list.g.seu[[1]], y=list.g.seu[[2]], 
             add.cell.ids=c("wt1","wt2"), project = "Giordani")
rm(list.g.seu)

# --*3*--

dmizerol <- read.table(paste0("data/DeMicheliD0/rawdataD0.txt"), sep="\t",
                       header=T, row.names=1)
dmich <- CreateSeuratObject(dmizerol, project="DeMicheli", min.cells=3, min.features=200)
rm(dmizerol)
# ====================================================================================
print("adding metadata concerning origin of datasets")
print("dellOrso METADATA already in object")
gio@meta.data[["orig.ident"]] <- rep("Giordani", length( gio@meta.data[["orig.ident"]] ))
dmich@meta.data[["orig.ident"]] <- rep("DeMicheli", length( dmich@meta.data[["orig.ident"]] ))
# ====================================================================================

# integration
# ====================================================================================

print("integrating")
muscle.list <- c(dorso, gio, dmich)
muscle.list <- sapply(muscle.list, function(seu){
  seu[["percent.mt"]] <- PercentageFeatureSet(seu,pattern="^mt-")
  seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < nft & percent.mt < 5)
  seu <- NormalizeData(seu,verbose=FALSE)
  seu <- FindVariableFeatures(seu, selection.method="vst",nfeatures=nft, verbose=FALSE)
  return(seu)}
)
muscle.anchors <- FindIntegrationAnchors(object.list = muscle.list, dims=1:30)
muscle.integrated <- IntegrateData(anchorset = muscle.anchors, dims = 1:30)

# switch to integrated assay. The variable features of this assay are
# automatically set during IntegrateData
DefaultAssay(object = muscle.integrated) <- "integrated"

# =====================================================================================
# find clusters
# =====================================================================================
print("scaling, reducing dim (PCA), printing ElbowPlot into preprocessIntegrated.pdf")
muscle.integrated <- ScaleData(muscle.integrated, verbose = FALSE)
muscle.integrated <- RunPCA(muscle.integrated, npcs = 30, verbose = FALSE)
#visuals
pdf(paste0(resu,"preprocessIntegrated.pdf"))
ElbowPlot(muscle.integrated)
DimHeatmap(muscle.integrated, dims=1:30, cells = 500, balanced = TRUE)
dev.off()

muscle.integrated <- KNNplusFITSNE(muscle.integrated, nbDIM, resolution)
saveRDS(muscle.integrated,file=paste0(rdsdir,"integrated_seu_fitsne.rds"))

# =====================================================================================
# markers : both, positive and negative LFC
# =====================================================================================
muscle.integrated.markers <- FindAllMarkers(muscle.integrated, only.pos = FALSE, 
                                            min.pct = 0.25, logfc.threshold = 0.25) #time consuming
# save into two tables:
threshold = 0.5
positive <- muscle.integrated.markers %>% filter(avg_logFC > threshold)
negative <- muscle.integrated.markers %>% filter(avg_logFC < -threshold)

write.table(positive, paste0(resu,"ALLMARKERS_Pos_integratedD0.txt"), sep="\t")
write.table(negative, paste0(resu,"ALLMARKERS_Neg_integratedD0.txt"), sep="\t")

# example getting 8 most relevant genes Pos and Neg :
most.neg <- negative %>% group_by(cluster) %>% top_n(n=4, wt= -avg_logFC)
most.pos <- positive %>% group_by(cluster) %>% top_n(n=4, wt= avg_logFC)

# =====================================================================================
# preliminary plots into results
# =====================================================================================
pdf(paste0(resu,"FITSNE.pdf"))
DimPlot(muscle.integrated, reduction="tsne", label=T,
        cols=definecolors(muscle.integrated@active.ident)) + ggtitle("Integrated D0")
dev.off()
pos4top <- positive %>% group_by(cluster) %>% top_n(n=4, wt= avg_logFC)
pdf(paste0(resu,"HEATMAP.pdf"), width=13)
DoHeatmap(muscle.integrated, features = pos4top$gene, 
          group.colors=definecolors(muscle.integrated@active.ident)) 
dev.off()

print("finished")
