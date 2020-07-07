# taking doublets by William
# and contrast with clusters found
# --
# johaGL
library(tidyverse)
library(Seurat)

prloc="~/INMG_SingleCell/"
setwd(prloc)

load("rds/doubletsD0/plotdoublet.RData")  # plotdoublet
head(plotdoublet)
# load the seurat object and celltype vectors
integrated.musD0 <- readRDS("rds/integratedD0/integrated_seu_fitsne.rds")
veccelltype <- readRDS("rds/integratedD0/integr_celltype_Vector.rds")

integrated.musD0 <- RenameIdents(integrated.musD0, veccelltype )

integrated.musD0@meta.data$barcodes = rownames(integrated.musD0@meta.data)

# do dataframe to play around:
fool = data.frame(barcodes=integrated.musD0@meta.data$barcodes, 
                  orig.ident=integrated.musD0@meta.data$orig.ident)
# transform barcodes to homonymous respect to plotdoublet
fool$id=NA
fool$id[fool$orig.ident %in% "DellOrso" & str_detect(fool$barcodes,"-1_1")]=
  paste0("dorso1_wt1_",str_replace(fool$barcodes[fool$orig.ident %in% "DellOrso" & str_detect(fool$barcodes,"-1_1")],"_1",""))
fool$id[fool$orig.ident %in% "DellOrso" & str_detect(fool$barcodes,"-1_2")]=
  paste0("dorso_wt2_",str_replace(fool$barcodes[fool$orig.ident %in% "DellOrso" & str_detect(fool$barcodes,"-1_2")],"_2",""))
fool$id[fool$orig.ident %in% "Giordani" & str_detect(fool$barcodes, "wt1_i1_")]=
  paste0("gio_", fool$barcodes[fool$orig.ident %in% "Giordani" & str_detect(fool$barcodes,"wt1_i1_")])
fool$id[fool$orig.ident %in% "Giordani" & str_detect(fool$barcodes, "wt2_i2_")]=
  paste0("gio_", fool$barcodes[fool$orig.ident %in% "Giordani" & str_detect(fool$barcodes,"wt2_i2_")])
fool$id[fool$orig.ident %in% "DeMicheli"]= paste0("demichili_", fool$barcodes[fool$orig.ident %in% "DeMicheli"])

merged <- inner_join(fool, plotdoublet, by="id")
integrated.musD0@meta.data$doubletclass=NA
integrated.musD0@meta.data$doubletclass = merged$DFclass

thecols = c(rgb(0.9, 0.2, 0.3, 1),
  rgb(0.2, 0.5, 0.8, 0.5), 
   rgb(0.1, 0.7, 0.3, 0.1))

a <- DimPlot(integrated.musD0, group.by = "doubletclass", 
             cols = thecols)
b <- DimPlot(integrated.musD0, cols=
               colorRampPalette(brewer.pal(8,"Set2"))(length(levels(integrated.musD0@active.ident))),
             label=T, repel=T, size=3)
plot_grid(a,b)

# HOW WE EXPLORED SUBSTRINGS ALONG BARCODES VECTORS
# unique(substr(plotdoublet$id, 19,20))
# unique(substr(integrated.musD0@meta.data$barcodes, 19,20))
# unique(substr(plotdoublet$id, 1,3))
# unique(substr(integrated.musD0@meta.data$barcodes, 1,4))
# head(integrated.musD0@meta.data$barcodes[str_detect(integrated.musD0@meta.data$barcodes,"_2")==T], n=2)
# tail(plotdoublet[str_detect(plotdoublet$id, "gio")==T ,])
# plotdoublet$id[plotdoublet$groupe %in% c("dorso1","dorso2)")]

 