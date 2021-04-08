##
# johaGL 2021
# Visualize only Sat cells (for Dr. Brun)
# (note: myonuclei included until new order)
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
library(ggthemes)
library(gridExtra)
prloc = "~/INMG_SingleCell/"
datadir = "data/Oprescu/"
rdsdir= "rds/OprescuTimePoints/"
resu = "results/OprescuTimePoints/"
sat.opr = "sat_seu_2021.rds"
CELLTYPEcol = "celltype" # the name I give to this metadata column
REDUC = "umap"  #if desired use "tsne" ==> is FIt-SNE indeed !!

setwd(prloc)
seu.ok <- readRDS(paste0(rdsdir, sat.opr))

markers <- c("Gli1", "Gli2", "Gli3", "Ptch1", "Smo", "Ihh", "Dhh")
shhgene <- "9530036O11Rik"

##########################
# DATA for FEATUREPLOTS FROM SCRATCH :D !!!
## **** new: custom function **
#########################
getcustom.dfFeat <- function(seu, REDUC, vec.markers){
  df4plot <- as.data.frame(seu@reductions[[REDUC]]@cell.embeddings)
  matcols <- colnames(seu@assays[["SCT"]])
  for (gene in vec.markers){
    vecexpr <- tryCatch({
      tmp <- as.vector(seu@assays[["SCT"]][gene, ])
      names(tmp) <- matcols
      tmp
    }, error=function(e){
      print(paste(gene," ==> looking into RNA assay (absent in SCT)"))
      tmp <- as.vector(seu@assays[["RNA"]][gene, ])
      names(tmp) <- matcols
      tmp
    })
    df4plot[[gene]] <- vecexpr[match(rownames(df4plot),names(vecexpr))]
  }
  return(df4plot)
}



# =============================================================================
# plot several markers into one single grid, SAT cells
# =============================================================================

colnames(df4plot)[colnames(df4plot)==shhgene] <- "Shh(9530036O11Rik)"
df4plot$bc <- rownames(df4plot)
dt_spread = reshape2::melt(data=df4plot,id=c("UMAP_1","UMAP_2", "bc"))
# cast???
# adjust colors to max value
maxval <- max(dt_spread$value)


# aestetics settings
mygray <- gray.colors(12, alpha=.2)[7] #moderately dark
somecols <- brewer_pal("seq", palette="YlGnBu")(9)[c(5,6,7,8,9)]
#transparent.colors = alpha(somecols, alpha = 0.9)
completecols <- c(mygray,somecols)
df4plot <- getcustom.dfFeat(seu.ok, REDUC, c(markers,shhgene))
# pdf
pdf(paste0(resu,"SAT_8markers_uniscale.pdf"),  width=10, height=10 )
ggplot(dt_spread, aes(UMAP_1,UMAP_2, colour= value)) + geom_point(size=1) +
  scale_colour_gradientn(colours = completecols, limits=c(0,maxval)) +
  theme_calc() +
  facet_wrap(~variable) + theme(strip.text.x = element_text(size = 15),
                                panel.grid.major = element_blank(), 
                                panel.grid.minor = element_blank())
dev.off() 

# =============================================================================
# plot markers in all the cells
# =============================================================================
filtered.seu <- readRDS(paste0(rdsdir, "filtered_opreFITSNEUMAP.rds"))
levels(filtered.seu@active.ident) #[1] "0"  "1"  "2"  "3"  "4"  "5"  ... "25" "26"
print("passing celltype names from saved rds vector (produced by scripts/impute_oprescuFULL.R) to seurat object")
new.cluster.ids <- readRDS(paste0(rdsdir,"opreFULL_transf_celltype_Vector"))
filtered.seu.imp <- RenameIdents(filtered.seu, new.cluster.ids)
filtered.seu.imp@meta.data[[CELLTYPEcol]] = filtered.seu.imp@active.ident
rm(filtered.seu)
#  to ease job, lets call "seu" the object 'filtered.seu.imp' and delete it
seu.big <- filtered.seu.imp; rm(filtered.seu.imp) 
#fix levels days injury:
seu.big@meta.data$orig.ident <- factor(x=seu.big@meta.data$orig.ident,
                levels=c( "0.5 DPI", "2 DPI", "3.5 DPI", "5 DPI", "10 DPI", "21 DPI", "Noninjured"))


df4plot.big <- getcustom.dfFeat(seu.big, REDUC, c(markers,shhgene))

colnames(df4plot.big)[colnames(df4plot.big)==shhgene] <- "Shh(9530036O11Rik)"

df4plot.big$bc <- rownames(df4plot.big)
HUGE_spread = reshape2::melt(data=df4plot.big,id=c("UMAP_1","UMAP_2", "bc"))

maxval.big <- max(HUGE_spread$value)

# aestetics settings MODIFIED mygray: 
mygray <- gray.colors(12, alpha=.2)[11]  # very LIGHT gray 
somecols <- brewer_pal("seq", palette="YlGnBu")(9)[c(5,6,7,8,9)]
completecols <- c(mygray,somecols)
# pdf
pdf(paste0(resu,"ALLpops_8markers_uniscale.pdf"),  width=10, height=10 )
ggplot(HUGE_spread, aes(UMAP_1,UMAP_2, colour= value)) + geom_point(size=0.2) +
  scale_colour_gradientn(colours = completecols, limits=c(0,maxval.big)) +
  theme_calc() +
  facet_wrap(~variable) + theme(strip.text.x = element_text(size = 15),
                                panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank())
dev.off() 



#################################"
# OLD VERSION:
#   ggl <- list()
# for (i in 1:length(markers)){
#   ggl[[i]] <- FeaturePlot(seu.ok, features=markers[i]) + 
#     theme(legend.position="none") +
#     theme_few(base_size=10) &
#     scale_colour_gradientn(colours = completecols) 
# }
# 
# ggl[[length(markers)+1]] <- FeaturePlot(seu.ok, features=shhgene) + 
#   labs(title="Shh (9530036O11Rik)") +
#   theme_few(base_size=10) & scale_colour_gradientn(colours = completecols)  
# 
# 
# pdf(paste0(resu,"SAT_8markerslocated.pdf"), width=15, height=14)
# plot_grid(ggl[[1]], ggl[[2]], ggl[[3]],
#           ggl[[4]], ggl[[5]], ggl[[6]],
#           ggl[[7]], ggl[[8]], nrow=3)
# dev.off()


# 
# FeaturePlot(seu.ok, features="Ptch1") &
#   scale_colour_gradientn(colours = rev(viridis_pal(alpha=.8, begin = 0.2,
#                                                    end=0.8,
#                                                    direction=-1,
#                                                    option="magma")(5)))
