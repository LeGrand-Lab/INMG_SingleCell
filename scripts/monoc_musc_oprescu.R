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
library(inlmisc)

prloc = "~/INMG_SingleCell/"
datadir = "data/Oprescu/"
rdsdir= "rds/OprescuTimePoints/"
resu = "results/OprescuTimePoints/"
CELLTYPEcol = "celltype" # the name I give to this metadata column
REDUC = "umap"  #if desired use "tsne" ==> is FIt-SNE indeed !!

setwd(prloc)
print("this analysis runs on musc+sm object, markers-Clusters were verified by the boss")
seu <- readRDS(paste0(rdsdir,"muscpostSEU.rds"))

print("printing DimPlots : calculated clusters, then plus markers")
plotsNums = list(); clusby=c(CELLTYPEcol, "seurat_clusters")
# *
plotsNums = lapply(clusby, function(y){ 
  p = DimPlot(seu, reduction = REDUC, group.by = y, 
              label=T, repel=T, label.size = 3 ) })

markerstab = read.table(paste0(resu,"tablesMarkersSubPops/musc_checkMarkersAndClusters.txt"),
                        header = T,sep='\t')
markerstab$cluster = as.factor(markerstab$cluster)
tmpdf = data.frame(numclus= seu@meta.data$seurat_clusters)
tmpdf = left_join(tmpdf, markerstab,by=c("numclus"="cluster"))
tmpdf <- tmpdf %>% mutate(nb_mark = paste0(numclus," ",concatmarkers))
seu@meta.data$nb_mark = as.factor(tmpdf$nb_mark)
# *
plotplusmarkers = DimPlot(seu, reduction = REDUC, group.by = "concat", pt.size = 0.3,
  label=T, repel=T, label.size = 3 ) + theme(legend.text = element_text(size=8))

tmpdf <- tmpdf %>% mutate(nb_newtype = case_when(
  nb_mark == "0 Mb_Csrp3" ~ "myonuclei",
  nb_mark == "1 Amd1_Myh4" ~ "myonuclei",
  nb_mark == "2 Top2a_Hmgb2" ~ "MuSCrenew",
  nb_mark == "3 Crip1_Spp1" ~  "MuSCprol",
  nb_mark == "4 Meg3_Fos" ~ "Asc",
  nb_mark == "5 Myog_Cdkn1c" ~ "Myocytes.early",
  nb_mark ==  "6 Myl1_Mylpf" ~ "Myocytes.late",
  nb_mark == "7 mt-Nd2_Myh1" ~ "myonuclei",
  nb_mark ==   "8 Lyz2_Apoe" ~ "Imb",
  nb_mark == "9 Mpz_Pmp22" ~ "Mpz_Pmp22",  # 9: neuromuscular junction cells?
  TRUE ~ "Qsc"
) )
seu@meta.data$newtype = tmpdf$nb_newtype
# * 
plotNEWtypes = DimPlot(seu, reduction = REDUC, group.by = "newtype", pt.size = 0.3,
   label=T, repel=T, label.size = 3 )+ theme(legend.text = element_text(size=8))

# TODO : fix, it did not print
pdf(paste0(resu,"cartMusc_subclustersPrep.pdf"),width=12)
plot_grid( plotNEWtypes,plotplusmarkers,
           plotsNums[[2]], plotsNums[[1]],
  nrow= 2 ) + plot_annotation(title="MuSC and SC clustering, steps (inversed order), random colors")
dev.off()

print("PROBLEMATIC cluster 9, exclude to do monocle(it runs reductdim and sizef again)")
head(seu@active.ident)
seu.prepmono <- subset(seu, idents="9",invert=T)

#================================================================================
# do monocle:
# ================================================================================

print("doing monocle")
source("~/INMG_SingleCell/scripts/functions_stock.R", local=T)
cds <- domono(seu.prepmono, "SCT") # As RNA gives error :starting vector near the null space
saveRDS(cds,paste0(rdsdir,"musc_postMONOcds.rds"))

# ================================================================================  
#              do visuals
# ================================================================================
# general
# ================================================================================

# colors for type
daysorder = c("0.5 DPI", "2 DPI", "3.5 DPI", "5 DPI", "10 DPI", "21 DPI", "Noninjured")
dayscols = rev(viridis(length(daysorder), alpha=0.6))
names(dayscols) = daysorder## !!

# colors for "newtype" : 
MANUALCOLORS=c("orange2","lightblue","gold2","deepskyblue4","peachpuff2","aquamarine3","cadetblue4","violetred")

cds <- readRDS(paste0(rdsdir,"musc_postMONOcds.rds"))
# order cells
plot_cell_trajectory(cds, color_by="orig.ident")
cds <- orderCells(cds) # interacting with plot to define.

## !!
cartseu.plot <- DimPlot(seu.prepmono, group.by = "newtype", cols = MANUALCOLORS, 
                    label=T, repel=T,pt.size = 0.4, label.size = 4)
dpiseu.plot <- DimPlot(seu.prepmono, group.by = "orig.ident", cols = dayscols, 
                    repel=T,pt.size = 0.4, label.size = 4)


trajTypeplt <- plot_cell_trajectory(cds,color_by="newtype",cell_size = 1, alpha = 0.8) +
  scale_color_manual(values=MANUALCOLORS)

pseudpiplt <- plot_cell_trajectory(cds,color_by="orig.ident",cell_size = 1, alpha = 0.8) +
  scale_color_manual(values=dayscols)
pseuplot <- plot_cell_trajectory(cds,color_by="Pseudotime",cell_size = 0.8, alpha = 0.8)

# print !!
pdf(paste0(resu,"MUSC_monocleFirstStep.pdf"), width=13)
plot_grid(plot_grid(cartseu.plot,  trajTypeplt + theme(legend.position = "right"),
                    NULL ,ncol=3, rel_widths = c(2,2,1)),  
          plot_grid(dpiseu.plot, pseudpiplt, pseuplot, ncol=3), nrow=2
          ) + plot_annotation("Monocle : MuSC+SM subpopulations")
dev.off()

# ================================================================================

# markers on trajectories
# ================================================================================

domarkerplot <- function(cds, mymarker){
  plot <- plot_cell_trajectory(cds, cell_size = 0.8,
                               markers=mymarker, backbone_color = "gray",
                               use_color_gradient = T, alpha=0.6
  ) +   scale_color_gradient(low = "gray", high = "red") + labs(title=mymarker)+
    theme(legend.position = "none") 
  return(plot)
}


grem1sign = c("Grem1", "Nog", "Bmp2", "Bmp4", "Bmpr1a",
          "Bmpr1b", "Bmpr2", "Id1", "Id2", "Id3","Bmp7")

listgremplots = lapply(grem1sign, function(x) domarkerplot(cds, x)
)

# do a last plot to get legend only
tmp <- plot_cell_trajectory(cds, cell_size = 0.8,markers="Pax7",
                     use_color_gradient = T, alpha=0.6) +   scale_color_gradient(low = "gray", high = "red") 
legendalone <- cowplot::get_legend(tmp)
listgremplots[[length(grem1sign)+1]] <- legendalone

# known genes as well
knowngenes = c("Myod1","Six1", "Myf5", "Vcam1","Pax7","Myh3" )
listknown = lapply(knowngenes, function(x) domarkerplot(cds, x))
listknown[[length(listknown)+1]] <- legendalone

# printing all markers
pdf(paste0(resu,"MUSC_markersMonoc.pdf"), height = 12)
plot_grid(plotlist = listgremplots, nrow=4) +
 plot_annotation(title = "Gremlin signaling: Pseudotemporal trajectories (MuSC + myonuclei)")
plot_grid(plotlist = listknown, nrow=4) + 
  plot_annotation(title = "known genes, trajectories (MuSC + myonuclei)")
dev.off()

# ================================================================================

# ================================================================================
#                 do some other visuals
# ================================================================================
#  violin plots by marker by timepoint
# fix factors
seu.prepmono@meta.data$newtype <- factor(seu.prepmono@meta.data$newtype, 
                                         levels=c("Asc","Imb","MuSCprol","MuSCrenew", "Myocytes.early", 
                                                  "Myocytes.late","myonuclei","Qsc"))
seu.prepmono@meta.data$orig.ident <- factor(x=seu.prepmono@meta.data$orig.ident,
                                            levels=c( "0.5 DPI", "2 DPI", "3.5 DPI", "5 DPI", "10 DPI", "21 DPI", "Noninjured"))

listVlngrem = lapply(grem1sign, function(x){
  p <- VlnPlot(seu.prepmono, features=x, cols=MANUALCOLORS, group.by = "orig.ident",
               pt.size = 0 , split.by = "newtype")
  p + theme(legend.position = "none", axis.title.x=element_blank(), axis.title.y=element_text(size=8),
            axis.text=element_text(size=8))
})
tmp2 <- VlnPlot(seu.prepmono, features="Pax7",cols=MANUALCOLORS, group.by = "orig.ident",
                split.by = "newtype") 
legendalone2 <- cowplot::get_legend(tmp2)
listVlngrem[[length(listVlngrem)+1]] <- legendalone2
plot_grid(plotlist = listVlngrem, ncol=3) + plot_annotation("Signaling genes by celltype across time")


listtmp = VlnPlot(seu.prepmono, features=grem1sign, cols=MANUALCOLORS, group.by="newtype",
                  pt.size = 0)
listVlngremB = list()
for (i in 1:length(grem1sign)){
  if (i > length(grem1sign)-2){
    listVlngremB[[i]]  = listtmp[[i]] + theme(legend.position = "none", 
                                              axis.title.y=element_text(size=8),
                                              axis.text=element_text(size=8), axis.title.x=element_blank())
  }else{
    listVlngremB[[i]]  = listtmp[[i]] + theme(legend.position = "none", axis.title.x=element_blank(), 
                                              axis.title.y=element_text(size=8),
                                              axis.text.x=element_blank(), axis.text.y=element_text(size=8))
  }
}
tmp2 <- VlnPlot(seu.prepmono, features="Pax7",cols=MANUALCOLORS, split.by = "newtype") 
legendalone2 <- cowplot::get_legend(tmp2)
listVlngremB[[length(listVlngremB)+1]] <- legendalone2
plot_grid(plotlist = listVlngremB, nrow=4, rel_heights = c(1,1,1,1.5)) + 
  plot_annotation("Signaling genes by celltype (resume)")

pdf(paste0(resu,"MUSC_ViolinPlots1.pdf"),width=14, height = 17)
plot_grid(plotlist = listVlngrem, ncol=3) + plot_annotation("Signaling genes by celltype across time")
plot_grid(plotlist = listVlngremB, NULL,NULL,NULL, NULL, ncol=4) + plot_annotation("Signaling genes by celltype (resume)")
dev.off()
### **


## END
# ================================================================================

# ### logical ok but failed:
# timesets <- SplitObject(seu.prepmono,split.by="orig.ident")
# plotsgrid = list()
# for (t in 1:7){ j <- head(timesets[[t]]@meta.data$orig.ident, 1)
# print(j) }  # order obtained: Noninjured, 0.5 2 3 5 10 21
# 
# for (t in 1:7){ 
#   mytitle = head(timesets[[t]]@meta.data$orig.ident, 1)
#   k <- VlnPlot(seu.prepmono, features=grem1sign,  
#                pt.size = 0 ,  cols=MANUALCOLORS, group.by = "newtype",combine=F)
#   tmp = list()
#   for (i in 1:(length(k)-1)){
#     tmp[[i]] <- k[[i]] + theme(axis.text.x = element_blank(), axis.text.y =element_text(size=7), 
#                                axis.title.x= element_blank(),  axis.title.y= element_blank(),
#                                title = element_text(size=8),
#                                legend.position = "none") 
#   }
#   tmp[[length(k)]] <- k[[length(k)]] + theme(axis.text = element_text(size=7), 
#                                              axis.title.x= element_blank(), title = element_text(size=8),
#                                              legend.position = "none")
#   plotsgrid[[t]] <- plot_grid(plotlist=tmp, ncol=1) + plot_annotation(title=mytitle)
# }
# #  plotsgrid[[7]]
# X = c("Noninjured", "0.5 DPI", "2 DPI", "3.5 DPI", "5 DPI", "10 DPI", "21 DPI")
# 

# pdf(paste0(resu,"MUSC_ViolinPlots2.pdf"), width=14,height=14)
# plot_grid(plotlist=plotsgrid, labels=X, ncol=7)
# dev.off()

# pdf(paste0(resu,"MUSC_ViolinPlots2.pdf"), width=14,height=14)
# plot_grid(
#   plot_grid(NULL+labs(title="jop"),NULL,NULL, ncol=7),
#   plot_grid(
#     plotsgrid[[1]] + plot_annotation(title=X[1]),
#     plotsgrid[[2]] + plot_annotation(title=X[2]),
#     plotsgrid[[3]] + plot_annotation(title=X[3]),
#     plotsgrid[[4]] + plot_annotation(title=X[4]),
#     plotsgrid[[5]] + plot_annotation(title=X[5]),
#     plotsgrid[[6]] + plot_annotation(title=X[6]),
#     plotsgrid[[7]] + plot_annotation(title=X[7]),
#     ncol=7),
#   ncol=1, rel_heights = c(1,8)
# )
# dev.off()
# 



# ================================================================================
# NOT GOOD:
# # this supposes that colors have defined names
# 
# comparisontmp <- VlnPlot(seu.prepmono, features=grem1sign,  
#                          pt.size = 0 , split.by = "orig.ident", cols=dayscols, group.by = "newtype")
# comparisonli <- list()
# for (i in 1:length(grem1sign)){
#   if (i > length(grem1sign)-2){
#     comparisonli[[i]] = comparisontmp[[i]] + theme(legend.position = "none", axis.title.x=element_blank(),
#                                                    axis.text=element_text(size=7), title=element_text(size=8))
#   }else{
#     comparisonli[[i]] = comparisontmp[[i]] + theme(legend.position = "none", axis.title=element_blank(),
#                                                    axis.text.y=element_blank(), title=element_text(size=8),
#                                                    axis.text.x=element_text(size=7))
#     
#   }}
# 
# tmp2 <- VlnPlot(seu.prepmono, features="Pax7", split.by = "orig.ident", cols=dayscols)
# 
# legendalone2 <- cowplot::get_legend(tmp2)
# comparisonli[[length(comparisonli)+1]] <- legendalone2
# plot_grid(plotlist = comparisonli, ncol=3) + plot_annotation("For comparison purposes")