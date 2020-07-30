# --
# JohaGL

library(Seurat)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(cowplot)
source("~/INMG_SingleCell/scripts/functions_stock.R", local=T)
prloc = "~/INMG_SingleCell/"
datadir = "data/Oprescu/"
rdsdir= "rds/OprescuTimePoints/"
resu = "results/OprescuTimePoints/"
CELLTYPEcol = "celltype" # the name I give to this metadata column

daysorder = c("0.5 DPI", "2 DPI", "3.5 DPI", "5 DPI", "10 DPI", "21 DPI", "Noninjured")
dayscols = viridis(length(daysorder), alpha=0.6)

setwd(prloc)
musc.sm.seu <- readRDS(paste0(rdsdir, "MuSC_SC_seu.rds"))
#fix levels factor DPI: 
musc.sm.seu@meta.data$orig.ident <- factor(x=musc.sm.seu@meta.data$orig.ident,
                                   levels=c( "0.5 DPI", "2 DPI", "3.5 DPI", "5 DPI", "10 DPI", "21 DPI", "Noninjured"))

dim(musc.sm.seu)
#take appart only satellite cells:
onlymusc = subset(musc.sm.seu, idents="MuSC")
# do violins (known markers): 


#######
# grem 1 genes across time:
myfeat_grem = c("Grem1", "Nog", "Bmp2", "Bmp4", "Bmp7", "Bmpr1a",
         "Bmpr1b", "Bmpr2", "Id1", "Id2", "Id3", "Id3") # repeated last one to make pair nb

violins <- VlnPlot(onlymusc, group.by = "orig.ident", 
                   ncol=length(myfeat_grem)/2, pt.size = 0,  #log=TRUE,
                   combine=FALSE, cols=rev(dayscols),
                   assay="SCT",
                   features=myfeat_grem) 

aesthetics = list(rep(0,length(violins)))
aesthetics[1:6] <- lapply(
  X = violins[1:6],
  FUN = function(p) p +theme(axis.title=element_blank(),
                             axis.text.x=element_blank(),
                             axis.text.y=element_text(size=8),
                             legend.position = "none"))

aesthetics[[7]] <- violins[[7]] + theme(axis.text=element_text(size=8),
                                        axis.title=element_text(size=8),
                                        legend.position = "none")
aesthetics[8:11] <- lapply(
  X = violins[8:11],
  FUN = function(p) p  +
    theme(axis.text.x=element_text(size=8),
          axis.text.y=element_text(size=8),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position = "none"))
aesthetics[[12]] <- cowplot::get_legend(violins[[12]])

vln_grem <- wrap_plots(aesthetics,ncol=length(myfeat_grem)/2) +
  plot_annotation(title= "Grem1 signaling across timepoints in MuSC cells")

#######
#### known genes
myfeat_vln = c("Myod1","Six1", "Myf5", "Vcam1","Pax7","Myh3" )  

violins <- VlnPlot(onlymusc, group.by = "orig.ident", 
                   ncol=length(myfeat_vln)/2, pt.size = 0,  #log=TRUE,
                   combine=FALSE, cols=rev(dayscols),
                   assay="SCT",
                   features=myfeat_vln) 

aesthetics = list(rep(0,length(violins)))
aesthetics[1:3] <- lapply(
  X = violins[1:3],
  FUN = function(p) p +theme(axis.title=element_blank(),
                             axis.text.x=element_blank(),
                             axis.text.y=element_text(size=7),
                             legend.position = "none"))

aesthetics[[4]] <- violins[[4]] + theme(axis.text=element_text(size=7),
                                        axis.title=element_text(size=7),
                                        legend.position = "none")
aesthetics[5:6] <- lapply(
  X = violins[5:6],
  FUN = function(p) p  +
    theme(axis.text.x=element_text(size=7),
          axis.text.y=element_text(size=7),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position = "none"))

vlnknown <- wrap_plots(aesthetics,ncol=length(myfeat_vln)/2) +
  plot_annotation(title= "Known markers across timepoints in MuSC cells")

##### PRINTING ####
pdf(paste0(resu, "MuSC_ViolinpltGrem1genes.pdf"), width=12)
vln_grem
dev.off()


