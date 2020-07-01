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
print(paste0("Will use *",nbDIM,"* princ components (dims), as we already runned this script", 
     "accordingly to elbowplot and good balance identities-nbclusters"))
 
# ====================================================================================
# open files
# --**--
dmizerol <- read.table(paste0("data/DeMicheliD0/rawdataD0.txt"), sep="\t",
                       header=T, row.names=1)
dmich <- CreateSeuratObject(dmizerol, project="DeMicheli", min.cells=3, min.features=200)
rm(dmizerol)
# --**--
dorso1.data <- Read10X(data.dir=paste0("data/DellOrsoD0/dorsowt1"))
dorso2.data <- Read10X(data.dir=paste0("data/DellOrsoD0/dorsowt2"))

dorso1 <- CreateSeuratObject(dorso1.data, project="DellOrso", min.cells=3, min.features=200)
dorso2 <- CreateSeuratObject(dorso2.data, project="DellOrso", min.cells=3, min.features=200)
rm(dorso1.data)
rm(dorso2.data)
dorso <- merge(dorso1,y=dorso2, add.cells.ids=c("wt1","wt2"),project="DellOrso")
rm(dorso1)
rm(dorso2)
# --**--
wt1 <- read.csv(paste0("data/GiordaniD0/GSM3520458_20171018_uninjured_wt_filtered.csv"), 
                sep=",", header=TRUE, row.names=1)
wt2 <- read.csv(paste0("data/GiordaniD0/GSM3520459_20180917_uninjured_wt_filtered.csv"),
                sep=",", header=TRUE, row.names=1)
gio1 <- CreateSeuratObject(wt1, project="Giordani", min.cells=3, min.features=200 )
gio1
gio2 <- CreateSeuratObject(wt2, project="Giordani", min.cells=3, min.features=200)
gio2
gio <- merge(gio1, y=gio2, add.cell.ids=c("wt1","wt2"), project = "Giordani")
gio
rm(wt1,wt2,gio1,gio2)
# ====================================================================================
print("adding metadata concerning origin of datasets")
print("dellOrso METADATA already in object")
gio@meta.data[["orig.ident"]] <- rep("Giordani", length( gio@meta.data[["orig.ident"]] ))
dmich@meta.data[["orig.ident"]] <- rep("DeMicheli", length( dmich@meta.data[["orig.ident"]] ))
# ====================================================================================
print("filter mitochondrial genes and low quality cells")
dmich[["percent.mt"]] <- PercentageFeatureSet(dmich, pattern= "^mt-")
dorso[["percent.mt"]] <- PercentageFeatureSet(dorso, pattern= "^mt-")
gio[["percent.mt"]] <- PercentageFeatureSet(gio, pattern= "^mt-")

dmich <- subset(dmich, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
dorso <- subset(dorso, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
gio <- subset(gio, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# ====================================================================================
# integration
print("integrating")
muscle.list <- c(dorso, gio, dmich)
# preprocessing: log-normalisation and variable features
for (i in 1:length(x=muscle.list)) {
  muscle.list[[i]] <- NormalizeData(muscle.list[[i]],verbose=FALSE)
  muscle.list[[i]] <- FindVariableFeatures(muscle.list[[i]], 
                                           selection.method="vst",nfeatures=2500, verbose=FALSE)
}
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

muscle.integrated <- KNNplusFITSNE(muscle.integrated, nbDIM, 0.5)
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
nb.clus = max(as.integer(levels(muscle.integrated@meta.data$seurat_clusters)))+1
pdf(paste0(resu,"FITSNE.pdf"))
DimPlot(muscle.integrated, reduction="tsne", label=T,
        cols=definecolors(muscle.integrated@active.ident)) + ggtitle("Integrated D0")
dev.off()
pos4top <- positive %>% group_by(cluster) %>% top_n(n=4, wt= avg_logFC)
pdf(paste0(resu,"HEATMAP.pdf"), width=13)
DoHeatmap(muscle.integrated, features = pos4top$gene, 
          group.colors=definecolors(muscle.integrated@active.ident)) 
dev.off()

