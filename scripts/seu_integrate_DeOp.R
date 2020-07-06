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
prefix = "DeOp"
resu="results/integratedD0_DeOp/"
rdsdir = "rds/integratedD0_DeOp/"
nbDIM=15
resol=0.4
threshold = 0.5 #threshold log2foldChange for saving markers tables
permit = 10
newfeat.intgr = 6000

dir.create(resu,recursive=T)
dir.create(rdsdir, recursive=T)

authors = c("DeMicheliD0","OprescuD0")
# define vector .rds files SAME ORDER and assign names:
fi.name = c("data/DeMicheliD0/rawdataD0.txt","data/OprescuD0/oprescu_Noninjured_raw.txt")
names(fi.name) = authors
maxUMI = c(25000,35000)
names(maxUMI) = authors

print("checking vector as you defined it, to run analysis")
print(rds.name)
setwd(prloc)

print("opening files, loading seu objects into list  and filtering")
muscle.list <- list()
muscle.list <- lapply(names(fi.name), function(auth){
  cgmat <- read.table(fi.name[auth], sep="\t", header=T, row.names = 1)
  seu <- CreateSeuratObject(cgmat, project=auth, min.cells=3,min.features = 200)
  return(seu)
})

muscle.list <- mapply( function(seu,auth){
  seu@meta.data[["orig.ident"]] <- rep(auth,length(seu@meta.data[["orig.ident"]]))
  seu[["percent.mt"]] <- PercentageFeatureSet(seu,pattern="^mt-")
  seu <- subset(seu, subset=nFeature_RNA>200 & 
                  nCount_RNA <=maxUMI[auth] & percent.mt < permit)
  return(seu)}, muscle.list, names(fi.name)
)

# ====================================================================================
# integration
# ====================================================================================
title <- ''; for (i in 1:length(authors)){title <- paste(title, authors[i])}
print(paste0("Integrating ",title," data. This is the main integration job, 
  reason: both on tibialis anterior m, similar experimental design"))
print(paste0("Will use *",nbDIM,"* princ components (dims), accordingly to  ", 
             "elbowplot and good balance identities-nbclusters seen in independent runs."))
print(paste0("   Using nbfeatures= ",newfeat.integr))


muscle.list <- lapply(muscle.list, function(seu){
  seu <- NormalizeData(seu,verbose=FALSE)
  seu <- FindVariableFeatures(seu, selection.method="vst",
           nfeatures= newfeat.intgr, verbose=FALSE)
  return(seu)
  })

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
pdf(paste0(resu,prefix,"_preprocessIntegrated.pdf"))
ElbowPlot(muscle.integrated)
DimHeatmap(muscle.integrated, dims=1:30, cells = 500, balanced = TRUE)
dev.off()

muscle.integrated <- KNNplusFITSNE(muscle.integrated, nbDIM, resol)

# =====================================================================================
# also run UMAP (slot FITSNE kept available)
# =====================================================================================
muscle.integrated <- RunUMAP(muscle.integrated, umap.method="uwot", dims=1:nbDIM)

saveRDS(muscle.integrated,file=paste0(rdsdir, prefix, "_integrated_seu_fitsne.rds"))
# =====================================================================================
# markers : both, positive and negative LFC
# =====================================================================================
muscle.integrated.markers <- FindAllMarkers(muscle.integrated, only.pos = FALSE, 
                                            min.pct = 0.25, logfc.threshold = 0.25) #time consuming
# save into two tables:

positive <- muscle.integrated.markers %>% filter(avg_logFC > threshold)
negative <- muscle.integrated.markers %>% filter(avg_logFC < -threshold)

write.table(positive, paste0(resu,prefix,"_ALLMARKERS_Pos_integratedD0.txt"), sep="\t")
write.table(negative, paste0(resu,prefix,"_ALLMARKERS_Neg_integratedD0.txt"), sep="\t")

# example getting 8 most relevant genes Pos and Neg :
most.neg <- negative %>% group_by(cluster) %>% top_n(n=4, wt= -avg_logFC)
most.pos <- positive %>% group_by(cluster) %>% top_n(n=4, wt= avg_logFC)

# =====================================================================================
# preliminary plots into results
# =====================================================================================

pdf(paste0(resu,prefix, "_FITSNEandUMAP.pdf"), width=13)
fit <- DimPlot(muscle.integrated, reduction="tsne", label=T,
        cols=definecolors(muscle.integrated@active.ident)) + 
  ggtitle(title) + theme(title=element_text(size=9))
uma <- DimPlot(muscle.integrated, reduction="umap", label=T,
        cols=definecolors(muscle.integrated@active.ident)) + 
  ggtitle(title) + theme(title=element_text(size=9))
plot_grid(fit,uma)
dev.off()

pos4top <- positive %>% group_by(cluster) %>% top_n(n=4, wt= avg_logFC)
pdf(paste0(resu,prefix,"_HEATMAP.pdf"), width=13)
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


