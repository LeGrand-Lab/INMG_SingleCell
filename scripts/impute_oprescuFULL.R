# SPECIAL CELL TYPE IMPUTATION: **BY LABELING TRANSFER**
# FROM OPRESCU RAW MATRIX METADATA
# --
# JohaGL 
library(Seurat)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(cowplot)
library(vididis)

prloc = "~/INMG_SingleCell/"
datadir = "data/Oprescu/"
rdsdir= "rds/OprescuTimePoints/"
resu = "results/OprescuTimePoints/"
daysorder = c("0.5 DPI", "2 DPI", "3.5 DPI", "5 DPI", "10 DPI", "21 DPI", "Noninjured")
dayscols = viridis(length(daysorder), alpha=0.6)

setwd(prloc)

customLabelTransfer <- function(typecell.authors.df, seu){
  tmpdf = seu@meta.data %>% as.tibble() %>% select(barcode, seurat_clusters)
  colnames(tmpdf) = c("barcode", "calcSEUcluster")
  typecell.authors.df = typecell.authors.df %>% select(cellID,cluster)
  tmpdf <- left_join(tmpdf, typecell.authors.df, by=c("barcode"="cellID") )
  # finding the best matches between 'cluster' given by authors, and calcSEUcluster by seurat:
  bestmatch <- as_tibble(tmpdf) %>% group_by(calcSEUcluster) %>% 
    summarize(given_cellname=names(which.max(table(cluster))))
  new.cluster.ids <- arrange(bestmatch,as.integer(calcSEUcluster)) %>% pull(given_cellname)
  names(new.cluster.ids) <- seq( 0, (length(new.cluster.ids)-1) )
  return(new.cluster.ids)  # a string vector celltype, having numeric levels
}



typecell.authors.df <- read.table(paste0(datadir, "oprescu_ExtractedMetaData.txt"), sep="\t",  header=T, row.names = 1)
#                           timepoint   cluster       metacluster                  cellID labelDescription
# 21 DPI_TTTGGTTTCTACCAGA     21dpi Osr1_FAPs              FAPs 21 DPI_TTTGGTTTCTACCAGA               NA

seu <- readRDS(paste0(rdsdir, "filtered_opreFITSNEUMAP.rds"))

print("before transferring celltypes labels, checking CALCULATED seurat_clusters:")
musc.bc <-  WhichCells(seu2, idents=15)
myonu.bc <- WhichCells(seu2, idents=14)
faps.bc <- list()
faps.bc <- lapply(c(11,10,1,16,13,16,6,19,17), function(x) WhichCells(seu2, idents=x))

mu_myo <- doDimPlotHighlight(seu2, list(musc.bc,myonu.bc), c("royalblue","gold"), "tsne",
                                "14_Myonuclei and 15_MuSC"  )
faps_only <- doDimPlotHighlight(seu2, faps.bc, 
                                   viridis_pal(option="B")(length(faps.bc)), 
                                   "tsne",
                                   "FAPs distinct sub-clusters"  )

pdf(paste0(resu,"oprescu_calcClusters.pdf"), width=10)
plot_grid(mu_myo, faps_only)
dev.off()
  
print(all(seu@active.ident == seu@meta.data$seurat_clusters))
new.cluster.ids <- customLabelTransfer(typecell.authors.df, seu)
saveRDS(new.cluster.ids, paste0(rdsdir,"opreFULL_transf_celltype_Vector"))
seu <- RenameIdents(seu, new.cluster.ids)


#filtered.seu.markersU <- read.table(paste0(resu, "filtered_opreMARKERSboth.txt"), sep="\t",  header=T) 