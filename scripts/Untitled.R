#!/usr/bin/env Rscript
# --
# JohaGL

library(Seurat)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(cowplot)
library(scales)
library(viridis)
library(Matrix)
library(reticulate)
library(monocle)

source("~/INMG_SingleCell/scripts/functions_stock.R", local=T)
prloc = "~/INMG_SingleCell/"
rdsdir= "rds/OprescuTimePoints/"
resu = "results/OprescuTimePoints/"
subdirmarkers = "tablesMarkersSubPops/"
dir.create(paste0(resu,subdirmarkers), recursive = T)
CELLTYPEcol = "celltype" # the name I give to this metadata column

setwd(prloc)

fapsTeno <- readRDS(paste0(rdsdir, "FAPS_Teno_seu.rds"))
fapsalone <- subset(fapsTeno, idents="Tenocytes", invert=T) #exclude tenocytes
saveRDS(fapsalone, paste0(rdsdir,"fapsalone.rds"))

filesrdslist = c("MuSC_SC_seu.rds", "FAPS_Teno_seu.rds", "fapsalone.rds")
names(filesrdslist) = c("musc","fapsTeno","fapsalone")
### running SEu for MUSC and FAPS

exitstatus <- sapply(names(filesrdslist), function(x){
  print("$$$$"); print(filesrdslist[[x]]);print("ùùùùùù")
  seu = readRDS(paste0(rdsdir,filesrdslist[[x]]))
  seu@meta.data$orig.ident <- factor(x=seu@meta.data$orig.ident,
   levels=c( "0.5 DPI", "2 DPI", "3.5 DPI", "5 DPI", "10 DPI", "21 DPI", "Noninjured"))
  
  #safely keep given celltypes
  seu@meta.data$Icelltype <- seu@meta.data[[CELLTYPEcol]]
  
  # determine sub-sub-clusters , do again seurat
  seu <- doSCT_PCA_UMAP(seu, 20)
  seu <- KNNplusFITSNE(seu, 20, 0.4)
  seu.markers <- FindAllMarkers(seu, only.pos=TRUE, min.pct=0.25, logfc.threschold= 0.25 )
  
  FourMainMarkers = seu.markers %>% group_by(cluster) %>% top_n(n=4, wt = avg_logFC)
  write.table(FourMainMarkers,paste0(resu,subdirmarkers, x, "oprescutable_main_markers.txt"), sep="\t", 
              col.names = T, row.names = F)
  
  # use top 2 markers for each calculated cluster
  twomain.markers <- seu.markers %>% group_by(cluster) %>% top_n(n=2, wt = avg_logFC)
  concatmarkers <- data.frame(cluster=unique(twomain.markers$cluster))
  tmp <- c()
  i <- 1
  base <- twomain.markers$gene
  while (i < length(base)){
    a = paste0(base[i],"_",base[i+1])
    tmp = c(tmp,a)
    i = i +2
  }
  tmp
  concatmarkers$concatmarkers = tmp
  
  write.table(concatmarkers,paste0(resu,subdirmarkers,x, "_checkMarkersAndClusters.txt"), sep="\t", 
              col.names = T, row.names = F)
  
  saveRDS(seu, paste0(rdsdir, x, "postSEU.rds"))
  return("ok")
}
)





