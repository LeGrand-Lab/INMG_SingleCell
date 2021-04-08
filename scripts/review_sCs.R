##
# johaGL 2021
#  Get seu object with only Sat cells (for Dr. Brun)
# (note: myonuclei excluded only for graphics)
##
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
sat.opr = "sat_seu_2021.rds"
CELLTYPEcol = "celltype" # the name I give to this metadata column
REDUC = "umap"  #if desired use "tsne" ==> is FIt-SNE indeed !!

setwd(prloc)
print(" markers-Clusters were verified by the boss")
seu <- readRDS(paste0(rdsdir,"muscpostSEU.rds"))

print("visualizing initial DimPlots : calculated clusters, plus markers")
plotsNums = list(); clusby=c(CELLTYPEcol, "seurat_clusters")
plotsNums = lapply(clusby, function(y){ 
  p = DimPlot(seu, reduction = REDUC, group.by = y, 
              label=T, repel=T, label.size = 3 ) })
plot_grid(plotlist=plotsNums)

markerstab = read.table(paste0(resu,"tablesMarkersSubPops/musc_checkMarkersAndClusters.txt"),
                        header = T,sep='\t')
markerstab$cluster = as.factor(markerstab$cluster)
tmpdf = data.frame(numclus= seu@meta.data$seurat_clusters)
tmpdf = left_join(tmpdf, markerstab,by=c("numclus"="cluster"))
tmpdf <- tmpdf %>% mutate(nb_mark = paste0(numclus," ",concatmarkers))
seu@meta.data$nb_mark = as.factor(tmpdf$nb_mark)
# *
plotplusmarkers = DimPlot(seu, reduction = REDUC, group.by = "nb_mark", pt.size = 0.3,
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
#pdf(paste0(resu,"cartMusc_subclustersPrep.pdf"),width=12)
plot_grid( plotNEWtypes,plotplusmarkers,
           plotsNums[[2]], plotsNums[[1]],
           nrow= 2 ) + plot_annotation(title="MuSC and SC clustering, steps (inversed order), random colors")
#dev.off()

print("PROBLEMATIC cluster 9, exclude it")
head(seu@active.ident)
seu.ok <- subset(seu, idents="9",invert=T)
seu.ok@meta.data$newtype <- factor(seu.ok@meta.data$newtype, 
                                         levels=c("Asc","Imb","MuSCprol","MuSCrenew", "Myocytes.early", 
                                                  "Myocytes.late","myonuclei","Qsc"))
seu.ok@meta.data$orig.ident <- factor(x=seu.ok@meta.data$orig.ident,
                                            levels=c( "0.5 DPI", "2 DPI", "3.5 DPI", "5 DPI", "10 DPI", "21 DPI", "Noninjured"))

MANUALCOLORS=c("orange2","lightblue","gold2","deepskyblue4","peachpuff2","aquamarine3","cadetblue4","violetred")

a <- DimPlot(seu.ok, reduction = REDUC, group.by = "newtype", pt.size = 0.4,
             label=T, repel=T, label.size = 3 ) + 
    theme(legend.text = element_text(size=8), axis.text=element_text(size=8),
          axis.title=element_text(size=8)) +
      scale_color_manual(values = MANUALCOLORS)  +
      labs(title="", subtitle="sub-populations")

b <- DimPlot(seu.ok, reduction = REDUC, group.by="orig.ident", pt.size = 0.4,
             label=T, repel=T, label.size = 3 )+ 
  theme(legend.text = element_text(size=8), axis.text=element_text(size=8),
        axis.title=element_text(size=8)) +
  scale_color_manual(values=rev(viridis_pal()(7))) + 
  labs(title="", subtitle="time-points")

abtitle <- ggdraw() + draw_label("Satellite cells and myonuclei")  
pdf(paste0(resu,"sCs_myo_forDrBrun.pdf"),width=10, height = 5)
plot_grid(abtitle,plot_grid(a,b),nrow=2, rel_heights = c(1,15))
dev.off()
#df4plot <- as.data.frame(seu.ok@reductions[["umap"]]@cell.embeddings)
#ggplot(df4plot) + geom_point(aes(UMAP_1, UMAP_2)) + theme_classic()

# =============================================================================
# save this seurat object as new one , FILTER OUT MYONUCLEI
# =============================================================================
seu.ok <- subset(seu.ok, idents=c(0,1,6,7), invert=T)

MANUALCOLORS2=c("orange2","lightblue","gold2","deepskyblue4","peachpuff2","violetred")
a2 <- DimPlot(seu.ok, reduction = REDUC, group.by = "newtype", pt.size = 0.4,
             label=T, repel=T, label.size = 3 ) + 
  theme(legend.text = element_text(size=8), axis.text=element_text(size=8),
        axis.title=element_text(size=8)) +
  scale_color_manual(values = MANUALCOLORS2)  +
  labs(title="", subtitle="sub-populations")

b2 <- DimPlot(seu.ok, reduction = REDUC, group.by="orig.ident", pt.size = 0.4,
             label=T, repel=T, label.size = 3 )+ 
  theme(legend.text = element_text(size=8), axis.text=element_text(size=8),
        axis.title=element_text(size=8)) +
  scale_color_manual(values=rev(viridis_pal()(7))) + 
  labs(title="", subtitle="time-points")

abtitle2 <- ggdraw() + draw_label("Satellite cells")  
pdf(paste0(resu,"SAT_forDrBrun_fitsne.pdf"),width=10, height = 5)
plot_grid(abtitle2,plot_grid(a2,b2),nrow=2, rel_heights = c(1,15))
dev.off()

saveRDS(seu.ok, file=paste0(rdsdir, sat.opr))





