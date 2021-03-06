#!/usr/bin/env Rscript
# Visualize Grem1 pathway genes on fitsne maps, 
# also known marker genes Satellite cells
# and other genes "at demand"
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
REDUCTIONdef= "tsne"
CELLTYPEcol = "celltype" # the name I give to this metadata column

setwd(prloc)
filtered.seu <- readRDS(paste0(rdsdir, "filtered_opreFITSNEUMAP.rds"))
levels(filtered.seu@active.ident) #[1] "0"  "1"  "2"  "3"  "4"  "5"  ... "25" "26"
print("passing celltype names from saved rds vector (produced by scripts/impute_oprescuFULL.R) to seurat object")
new.cluster.ids <- readRDS(paste0(rdsdir,"opreFULL_transf_celltype_Vector"))
filtered.seu.imp <- RenameIdents(filtered.seu, new.cluster.ids)
filtered.seu.imp@meta.data[[CELLTYPEcol]] = filtered.seu.imp@active.ident

rm(filtered.seu)
#  to ease job, lets call "seu" the object 'filtered.seu.imp' and delete it
seu <- filtered.seu.imp; rm(filtered.seu.imp) 
#fix levels days injury:
seu@meta.data$orig.ident <- factor(x=seu@meta.data$orig.ident,
              levels=c( "0.5 DPI", "2 DPI", "3.5 DPI", "5 DPI", "10 DPI", "21 DPI", "Noninjured"))

wholeplot <- function(gene,seu){
  allmix <-  FeaturePlot(seu,features=gene, pt.size = 0.4, reduction = REDUCTIONdef)  + 
    ggtitle(paste(gene," : all timepoints")) +
    theme(axis.title = element_text(size=6),
          axis.text= element_text(size = 4),
          title=element_text(size=8),
          legend.position = "right",
          legend.text = element_text(size = 6))
  return(allmix) # one single plot
}
spliplot <- function(gene,seu){
  bigtmplist <- FeaturePlot(seu,features=gene, split.by = "orig.ident",
                            reduction = REDUCTIONdef,
                            ncol=length(levels(seu@meta.data$orig.ident)), 
                            pt.size = 0.4, combine=F)
  omg <- list()  ;  numplots= length(levels(seu@meta.data$orig.ident))
  for (i in 1:numplots){
    omg[[i]] <- bigtmplist[[i]] + 
      theme(text=element_text(size=7),axis.title = element_blank(),
             axis.text= element_blank(),legend.position = "none")  
    }
  return(omg)
}

# ¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬
# Grem1 signaling plots
# ¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬

# Grem1 
# ============================================================================

grem.A <- plot_grid(wholeplot("Grem1",seu), ncol=1)
grem.B <- plot_grid(plotlist=spliplot("Grem1",seu), nrow=2)

# plot_grid(grem.A,grem.B, ncol= 2, rel_widths = c(2,5))
# ================================================================================

# Noggin
# ================================================================================
nog.A <- plot_grid(wholeplot("Nog",seu), ncol=1)
nog.B <- plot_grid(plotlist=spliplot("Nog",seu), nrow=2)

# plot_grid(nog.A,nog.B, ncol= 2, rel_widths = c(2,5))

# ================================================================================

# Bmp2
# ================================================================================
bmp2.A <- plot_grid(wholeplot("Bmp2",seu),ncol=1)
bmp2.B <- plot_grid(plotlist=spliplot("Bmp2",seu), nrow=2)

# plot_grid(bmp2.A, bmp2.B, ncol=2, rel_widths = c(2,5))
# ================================================================================

# Bmp4
# ================================================================================
bmp4.A <- plot_grid(wholeplot("Bmp4",seu),ncol=1)
bmp4.B <- plot_grid(plotlist=spliplot("Bmp4",seu), nrow=2)

# plot_grid(bmp4.A, bmp4.B, ncol=2, rel_widths = c(2,5))
# ================================================================================
# Bmp7
# ================================================================================
bmp7.A <- plot_grid(wholeplot("Bmp7",seu),ncol=1)
bmp7.B <- plot_grid(plotlist=spliplot("Bmp7",seu), nrow=2)

# plot_grid(bmp7.A, bmp7.B, ncol=2, rel_widths = c(2,5))

# ================================================================================

# Bmpr1a
# ================================================================================
bmpr1a.A <- plot_grid(wholeplot("Bmpr1a",seu),ncol=1)
bmpr1a.B <- plot_grid(plotlist=spliplot("Bmpr1a",seu), nrow=2)

# plot_grid(bmpr1a.A , bmpr1a.B, ncol=2, rel_widths = c(2,5))

# Bmpr1b
# ================================================================================

# Bmpr1b
# ================================================================================
bmpr1b.A <- plot_grid(wholeplot("Bmpr1b",seu),ncol=1)
bmpr1b.B <- plot_grid(plotlist=spliplot("Bmpr1b",seu), nrow=2)

# plot_grid(bmpr1b.A , bmpr1b.B, ncol=2, rel_widths = c(2,5))

# ================================================================================
# Bmpr2
# ================================================================================
bmpr2.A <- plot_grid(wholeplot("Bmpr2",seu),ncol=1)
bmpr2.B <- plot_grid(plotlist=spliplot("Bmpr2",seu), nrow=2)

# plot_grid(bmpr2.A , bmpr2.B, ncol=2, rel_widths = c(2,5))
# ================================================================================

# Id1
# ================================================================================
Id1.A <- plot_grid(wholeplot("Id1",seu),ncol=1)
Id1.B <- plot_grid(plotlist=spliplot("Id1",seu), nrow=2)

# plot_grid(Id1.A, Id1.B, ncol=2, rel_widths = c(2,5))

# ================================================================================

# Id2
# ================================================================================
Id2.A <- plot_grid(wholeplot("Id2",seu),ncol=1)
Id2.B <- plot_grid(plotlist=spliplot("Id2",seu), nrow=2)

# plot_grid(Id2.A, Id2.B, ncol=2, rel_widths = c(2,5))

# ================================================================================

# Id3
# ================================================================================
Id3.A <- plot_grid(wholeplot("Id3",seu),ncol=1)
Id3.B <- plot_grid(plotlist=spliplot("Id3",seu), nrow=2)

# plot_grid(Id3.A, Id3.B, ncol=2, rel_widths = c(2,5))

# ================================================================================

pdf(paste0(resu,"WholeCells_genesGrem1signaling.pdf"), width=14, height= 24)
plot_grid(
  plot_grid(grem.A,grem.B, ncol= 2, rel_widths = c(2,5)),
  plot_grid(nog.A,nog.B, ncol= 2, rel_widths = c(2,5)),
  plot_grid(bmp2.A, bmp2.B, ncol=2, rel_widths = c(2,5)),
  plot_grid(bmp4.A, bmp4.B, ncol=2, rel_widths = c(2,5)),
  plot_grid(bmp7.A, bmp7.B, ncol=2, rel_widths = c(2,5)),
  plot_grid(NULL, NULL, ncol=2, rel_widths = c(2,5)),
  ncol=1) + plot_annotation("Grem1 signaling genes, 1 of 2 pages")
plot_grid(
  plot_grid(bmpr1a.A , bmpr1a.B, ncol=2, rel_widths = c(2,5)),
  plot_grid(bmpr1b.A , bmpr1b.B, ncol=2, rel_widths = c(2,5)),
  plot_grid(bmpr2.A , bmpr2.B, ncol=2, rel_widths = c(2,5)),
  plot_grid(Id1.A, Id1.B, ncol=2, rel_widths = c(2,5)),
  plot_grid(Id2.A, Id2.B, ncol=2, rel_widths = c(2,5)),
  plot_grid(Id3.A, Id3.B, ncol=2, rel_widths = c(2,5)), 
  ncol = 1
)  + plot_annotation("Grem1 signaling genes, 2 of 2 pages")
dev.off()

# ¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬
# known marker genes plots
# ¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬
print("printing also known genes across timepoints")
## Prkaa1 == AMPK
knowngenes = c("Pax7",  "Myf5", "Myod1", "Myh3", "Six1", "Prkaa1")
Lwhole = lapply(knowngenes, function(x) return(
  plot_grid(wholeplot(x, seu),ncol=1)))
Swhole = lapply(knowngenes, function(x) return(
  plot_grid(plotlist=spliplot(x, seu), nrow=2)
))
pdf(paste0(resu,"WholecellsKnowngenes.pdf"), width=14, height= 24)
plot_grid(
  plot_grid(Lwhole[[1]], Swhole[[1]], rel_widths=c(2,5)),
  plot_grid(Lwhole[[2]], Swhole[[2]], rel_widths=c(2,5)),
  plot_grid(Lwhole[[3]], Swhole[[3]], rel_widths=c(2,5)),
  plot_grid(Lwhole[[4]], Swhole[[4]], rel_widths=c(2,5)),
  plot_grid(Lwhole[[5]], Swhole[[5]], rel_widths=c(2,5)),
  plot_grid(Lwhole[[6]], Swhole[[6]], rel_widths=c(2,5)),
   ncol= 1) + plot_annotation("Known marker genes")
dev.off()

# ¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬
# Plots at demand
# ¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬

# Tek
# ================================================================================
Tek.A <- plot_grid(wholeplot("Tek",seu),ncol=1)
Tek.B <- plot_grid(plotlist=spliplot("Tek",seu), nrow=2)


# Id3
# ================================================================================
Vcam1.A <- plot_grid(wholeplot("Vcam1",seu),ncol=1)
Vcam1.B <- plot_grid(plotlist=spliplot("Vcam1",seu), nrow=2)


# Gli3
# ================================================================================
Gli3.A <- plot_grid(wholeplot("Gli3",seu),ncol=1)
Gli3.B <- plot_grid(plotlist=spliplot("Gli3",seu), nrow=2)

# Arl13b
# ================================================================================
Arl13b.A <- plot_grid(wholeplot("Arl13b",seu),ncol=1)
Arl13b.B <- plot_grid(plotlist=spliplot("Arl13b",seu), nrow=2)

# Ift88
# ================================================================================
Ift88.A <- plot_grid(wholeplot("Ift88",seu),ncol=1)
Ift88.B <- plot_grid(plotlist=spliplot("Ift88",seu), nrow=2)
# ================================================================================

#  printing Tek and Vcam1 and the others
pdf(paste0(resu, "wholecellsTekVcaGliArlIft.pdf"), width=14, height= 6)
plot_grid(Tek.A, Tek.B, ncol=2, rel_widths = c(2,5))
plot_grid(Vcam1.A, Vcam1.B, ncol=2, rel_widths = c(2,5))
plot_grid(Gli3.A, Gli3.B, ncol=2, rel_widths = c(2,5))
plot_grid(Arl13b.A, Arl13b.B, ncol=2, rel_widths = c(2,5))
plot_grid(Ift88.A, Ift88.B, ncol=2, rel_widths = c(2,5))
dev.off()


# ================================================================================

# Shh
# ================================================================================
Shh.A <- plot_grid(wholeplot("9530036O11Rik",seu),ncol=1)
Shh.B <- plot_grid(plotlist=spliplot("9530036O11Rik",seu), nrow=2)

#plot_grid(Shh.A, Shh.B, ncol=2, rel_widths = c(2,5))
# ================================================================================

# Ihh
# ================================================================================
Ihh.A <- plot_grid(wholeplot("Ihh",seu),ncol=1)
Ihh.B <- plot_grid(plotlist=spliplot("Ihh",seu), nrow=2)

# plot_grid(Ihh.A, Ihh.B, ncol=2, rel_widths = c(2,5))

# ================================================================================
# Dhh
# ================================================================================
Dhh.A <- plot_grid(wholeplot("Dhh",seu),ncol=1)
Dhh.B <- plot_grid(plotlist=spliplot("Dhh",seu), nrow=2)

# plot_grid(Dhh.A, Dhh.B, ncol=2, rel_widths = c(2,5))

pdf(paste0(resu,"wholecellsSIDhedgehog.pdf"), width=14, height= 6)
plot_grid(Shh.A, Shh.B, ncol=2, rel_widths = c(2,5))
plot_grid(Ihh.A, Ihh.B, ncol=2, rel_widths = c(2,5))
plot_grid(Dhh.A, Dhh.B, ncol=2, rel_widths = c(2,5))
dev.off()




