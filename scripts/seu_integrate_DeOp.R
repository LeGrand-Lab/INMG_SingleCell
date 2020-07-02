# INTEGRATED ANALYSIS  ** Oprescu and DeMicheli **
# Integrating Oprescu and DeMicheli D0 data. This is the main integration job, 
# reason: both on tibialis anterior, similar experimental design

# --
# JohaGL

library(dplyr)
library(tidyverse)
library(Seurat)
library(ggplot2)
library(sctransform)
library(cowplot)
library(RColorBrewer)
source(file="~/INMG_SingleCell/scripts/functions_stock.R",local=T)

prloc="~/INMG_SingleCell/"
resu="results/integratedD0_DeOp/"
rdsdir = "rds/integratedD0_DeOp/"
nbDIM=15
resol=0.4
projectstr= "Oprescu_D0"
threshold = 0.5 #threshold log2foldChange for saving markers tables
newfeat.intgr = 6000
maxUMIop = 35000
maxUMIde = 25000
permit = 10

setwd(prloc)
dir.create(resu,recursive=T)
dir.create(rdsdir, recursive=T)

print("Integrating Oprescu and DeMicheli D0 data. This is the main integration job, 
  reason: both on tibialis anterior m, similar experimental design")
print(paste0("Will use *",nbDIM,"* princ components (dims), accordingly to  ", 
             "elbowplot and good balance identities-nbclusters seen in independent runs"))

opre.data <- read.table("data/OprescuD0/oprescu_Noninjured_raw.txt", sep="\t",
                        header=T, row.names = 1)
opre <- CreateSeuratObject(opre.data, project=projectstr, min.cells=3, min.features=200)

dmich.data <- read.table("data/DeMicheliD0/rawdataD0.txt", sep="\t",
                    header=T, row.names=1)
dmich <- CreateSeuratObject(dmich.data, project="DeMicheli", min.cells=3,min.features = 200)

# ====================================================================================
print("adding metadata concerning origin of datasets")
opre@meta.data[["orig.ident"]] <- rep("Oprescu", length( opre@meta.data[["orig.ident"]] ))
dmich@meta.data[["orig.ident"]] <- rep("DeMicheli", length( dmich@meta.data[["orig.ident"]] ))
# ====================================================================================
print("filter mitochondrial genes and low quality cells")
opre [["percent.mt"]] <- PercentageFeatureSet(opre, pattern= "^mt-")
dmich[["percent.mt"]] <- PercentageFeatureSet(dmich, pattern= "^mt-")

print(paste0("doing data filtering, admit max nbUMI: ", maxUMIop, " oprescu,  ",
             maxUMIde," demicheli. For both, admited mitochondrial % ",permit))
opre <- subset(opre, subset=nFeature_RNA > 200 & nCount_RNA <= maxUMIop & percent.mt < permit)
dmich <- subset(dmich, subset=nFeature_RNA > 200 & nCount_RNA <= maxUMIde & percent.mt < permit)


# ====================================================================================
# integration

print(paste0("integrating, using nbfeatures= ",newfeat.integr))

muscle.list <- c(opre,dmich)
# preprocessing: log-normalisation and variable features
for (i in 1:length(x=muscle.list)) {
  muscle.list[[i]] <- NormalizeData(muscle.list[[i]],verbose=FALSE)
  muscle.list[[i]] <- FindVariableFeatures(muscle.list[[i]], 
                                           selection.method="vst",
                                           nfeatures= newfeat.intgr, verbose=FALSE)
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

muscle.integrated <- KNNplusFITSNE(muscle.integrated, nbDIM, resol)
saveRDS(muscle.integrated,file=paste0(rdsdir,"integrated_seu_fitsne.rds"))
# =====================================================================================
# also run UMAP (slot FITSNE kept available)
# =====================================================================================
muscle.integrated <- RunUMAP(muscle.integrated, umap.method="uwot", dims=1:nbDIM)
# =====================================================================================
# markers : both, positive and negative LFC
# =====================================================================================
muscle.integrated.markers <- FindAllMarkers(muscle.integrated, only.pos = FALSE, 
                                            min.pct = 0.25, logfc.threshold = 0.25) #time consuming
# save into two tables:

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

pdf(paste0(resu,"FITSNEandUMAP.pdf"))
DimPlot(muscle.integrated, reduction="tsne", label=T,
        cols=definecolors(muscle.integrated@active.ident)) + 
  ggtitle(str_replace(projectstr,"_",""))
DimPlot(muscle.integrated, reduction="umap", label=T,
        cols=definecolors(muscle.integrated@active.ident)) + 
  ggtitle(str_replace(projectstr,"_",""))
dev.off()

pos4top <- positive %>% group_by(cluster) %>% top_n(n=4, wt= avg_logFC)
pdf(paste0(resu,"HEATMAP.pdf"), width=13)
DoHeatmap(muscle.integrated, features = pos4top$gene, 
          group.colors=definecolors(muscle.integrated@active.ident)) 
dev.off()

print("finished")


# ====================================================================================
sink(paste0(resu,"sessionInfo.txt"), append=TRUE) 
sink(paste0(resu,"sessionInfo.txt"), append=TRUE, type="message")
sessionInfo()
sink()
sink(type="message")


