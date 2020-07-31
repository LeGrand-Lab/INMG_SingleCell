#!/usr/bin/env Rscript
# --
# JohaGL
library(ggplot2)
library(tidyverse)
library(Seurat)
library(scales)
library(viridis)
library(Matrix)
library(reticulate)
library(monocle)
library(patchwork)
library(cowplot)
library(RColorBrewer)

prloc = "~/INMG_SingleCell/"
datadir = "data/Oprescu/"
rdsdir= "rds/OprescuTimePoints/"
resu = "results/OprescuTimePoints/"
CELLTYPEcol = "celltype" # the name I give to this metadata column
REDUC = "tsne"  #if desired use "tsne" ==> is FIt-SNE indeed !!

setwd(prloc)
print("this analysis runs on fapsalone object")
seu <- readRDS(paste0(rdsdir,"fapsalonepostSEU.rds"))
a <- DimPlot(seu,reduction="tsne",group.by = "celltype", label = T)
b <- DimPlot(seu,reduction="tsne",group.by = "seurat_clusters", label=T)
a+b

correctednoms = c("Dlk1_FAPs", "Cfh_FAPs", #formerOSR1
                  "Sfrp2_FAPs", # Wisp1
                  "Cxcl14_FAPs", "Dpp4_FAPs",
                  "Fibroblasts", "ActivatedFAPs",
                  "ActivatedFAPs", "8: Cd74",
                  "Cxcl14_FAPs", "10: Myl1"
                  )
names(correctednoms) = seq(0,10)
seu <- RenameIdents(seu,correctednoms)

source("~/INMG_SingleCell/scripts/functions_stock.R", local=T)
cds <- domono(seu, "RNA") 
# Error in as.igraph.vs(graph, nodes) : Invalid vertex names !
saveRDS(cds,paste0(rdsdir,"fapsA_postMONOcds.rds")) #failed


# take away "8" '10 clusters, are not FAPS
seu <- subset(seu,ident=c("8: Cd74","10: Myl1"), invert=T)
cds <- domono(seu, "RNA") 
saveRDS(cds,paste0(rdsdir,"fapsB_postMONOcds.rds"))


