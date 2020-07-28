# Gives plots showing general landscape of oprescu's cells
# and saves  FAPS and MuSC+SM sub-populations into each .rds
# --
# JohaGL 

library(Seurat)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(cowplot)
library(inlmisc)
source("~/INMG_SingleCell/scripts/functions_stock.R", local=T)
prloc = "~/INMG_SingleCell/"
datadir = "data/Oprescu/"
rdsdir= "rds/OprescuTimePoints/"
resu = "results/OprescuTimePoints/"
library(viridis)
daysorder = c("0.5 DPI", "2 DPI", "3.5 DPI", "5 DPI", "10 DPI", "21 DPI", "Noninjured")
dayscols = viridis(length(daysorder), alpha=0.6)

CELLTYPEcol = "celltype" # the name I give to this metadata column

setwd(prloc)
print("passing celltype names from saved rds vector (produced by scripts/impute_oprescuFULL.R) to seurat object")
filtered.seu <- readRDS(paste0(rdsdir, "filtered_opreFITSNEUMAP.rds"))
levels(filtered.seu@active.ident) #[1] "0"  "1"  "2"  "3"  "4"  "5"  ... "25" "26"
new.cluster.ids <- readRDS(paste0(rdsdir,"opreFULL_transf_celltype_Vector"))
filtered.seu.imp <- RenameIdents(filtered.seu, new.cluster.ids) 
# this immediatly changes @active.ident:
cat(levels(filtered.seu.imp@active.ident), sep=" ") # M2 Osr1_FAPs ... Pericytes Cd8_Tcells
rm(filtered.seu)

print("passing imputed cell types into column $celltype")
filtered.seu.imp@meta.data[[CELLTYPEcol]] = filtered.seu.imp@active.ident
colscelltype = colorRampPalette(brewer.pal(8,"Set2"))(length(levels(filtered.seu.imp@meta.data[[CELLTYPEcol]])))

# creating plots

A1 <- DimPlot(filtered.seu.imp, group.by="orig.ident", order = daysorder, cols = dayscols, reduction = "umap") +
  theme(legend.position = "none")
A2 <- DimPlot(filtered.seu.imp, group.by="orig.ident", order = daysorder, cols = dayscols, reduction = "tsne") 
legendTP <- cowplot::get_legend(A2)
B1 <- DimPlot(filtered.seu.imp, group.by=CELLTYPEcol, label=T, repel=T, reduction = "umap",
        cols = colscelltype, label.size=3) + theme(legend.position = "none")
B2 <- DimPlot(filtered.seu.imp, group.by=CELLTYPEcol, label=T, repel=T, reduction = "tsne",
        cols = colscelltype , label.size=3) 
legendclusters <- cowplot::get_legend(B2)


musc.bcL <- WhichCells(filtered.seu.imp,idents="MuSC")
myonu.bcL <- WhichCells(filtered.seu.imp,idents="SM")
faps.bcL <- list()
faps.bcL <- lapply(c(c("Osr1_FAPs","Wisp1_FAPs", "Cxcl14_FAPs",  "Dpp4_FAPs",
                      "Dlk1_FAPs" , "Fibroblasts", "ActivatedFAPs" , "Tenocytes" )), 
                  function(x) WhichCells(filtered.seu.imp, idents=x)) 

mu_myoLab <- doDimPlotHighlight(filtered.seu.imp, list(musc.bcL,myonu.bcL), c("cadetblue","orange2"), 
                                "tsne","MuSC and Myonuclei(SM)" )
faps_onlyLab <- doDimPlotHighlight(filtered.seu.imp, faps.bcL, 
                                GetColors(length(faps.bcL)+2,scheme="iridescent", rev=T, alpha=0.7)[1:length(faps.bcL)], 
                                "tsne", "FAPs distinct sub-clusters" ) # last argument "" as must exist a title

# create temporary seu to plot clusters numbers:
tmpseu <- filtered.seu.imp; inversedclus = names(new.cluster.ids); names(inversedclus)=new.cluster.ids
tmpseu <- RenameIdents(tmpseu,inversedclus) #only inside function, no change is returned
musc.bc <-  WhichCells(tmpseu, idents=15)
myonu.bc <- WhichCells(tmpseu, idents=14)
faps.bc <- list()
faps.bc <- lapply(c(1,6,10,11,13,16,19,17), function(x) WhichCells(tmpseu, idents=x))
mu_myo <- doDimPlotHighlight(tmpseu, list(musc.bc,myonu.bc), c("cadetblue","orange2"), "tsne","."  )
faps_only <- doDimPlotHighlight(tmpseu, faps.bc, 
                                GetColors(length(faps.bc)+2,scheme="iridescent", rev=T, alpha=0.7)[1:length(faps.bc)], 
                                "tsne"," "  )

mu_myo <- doDimPlotHighlight(tmpseu, list(musc.bc,myonu.bc), c("cadetblue","orange2"), "tsne",
                             " "  )
rm(tmpseu)

##print(paste0("saving selected sub-populations .rds into"),rds)  ### pending extracting subpopulations: faps+tenocyts faps musc+sm

# PRINTS
# ===========================================================================
print(paste0("Printing general landscape plots into  ",resu, "plots_landscape.pdf"))
givenTtl = "Cellular populations along datapoints"
pdf(paste0(resu, "plots_landscapeNew.pdf"), width=13, height=14)
plot_grid(
  plot_grid( A1,A2+theme(legend.position = "none"), legendTP,
             B1,B2+theme(legend.position="none"), legendclusters, nrow= 2, 
             rel_widths = c(2,2,1.5), labels=c("A","B","","C","D","") ) ,
  
  plot_grid(mu_myoLab,faps_onlyLab,NULL,
            mu_myo,faps_only,NULL,nrow=2,
            rel_widths = c(2,2,1.5), labels=c("E","F","", "G","H","") ) ,
  
  nrow=2
)  + plot_annotation(givenTtl)
dev.off()
# 
# plot_grid(plot_grid(A1,A2+theme(legend.position = "none"), legendTP, nrow=1,rel_widths = c(2,2,1.5)),
#           plot_grid(B1,B2+theme(legend.position="none"), legendclusters, nrow= 1, rel_widths = c(2,2,1.5)), 
#           nrow=2, labels=c("A","B")) + plot_annotation("blah")

