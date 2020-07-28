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
library(inlmisc)
source("~/INMG_SingleCell/scripts/functions_stock.R", local=T)
prloc = "~/INMG_SingleCell/"
datadir = "data/Oprescu/"
rdsdir= "rds/OprescuTimePoints/"
resu = "results/OprescuTimePoints/"

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
#                         timepoint   cluster       metacluster     cellID            labelDescription
# 21 DPI_TTTGGTTTCTACCAGA     21dpi   Osr1_FAPs     FAPs       21 DPI_TTTGGTTTCTACCAGA               NA

filtered.seu <- readRDS(paste0(rdsdir, "filtered_opreFITSNEUMAP.rds"))

print("before transferring celltypes labels, checking CALCULATED seurat_clusters:")
musc.bc <-  WhichCells(filtered.seu, idents=15)
myonu.bc <- WhichCells(filtered.seu, idents=14)
faps.bc <- list()
faps.bc <- lapply(c(11,10,1,16,13,16,6,19,17), function(x) WhichCells(filtered.seu, idents=x))

mu_myo <- doDimPlotHighlight(filtered.seu, list(musc.bc,myonu.bc), c("cadetblue","orange2"), "tsne",
                                "15_MuSC and 14_Myonuclei "  )
faps_only <- doDimPlotHighlight(filtered.seu, faps.bc, 
                                  GetColors(length(faps.bc)+2,scheme="iridescent", rev=T, alpha=0.7)[1:length(faps.bc)], 
                                   "tsne",
                                   "FAPs distinct sub-clusters"  )
print(paste0("CALCULATED seurat_clusters plot saved into : ", resu, "xpreview_opre_calcClusters.pdf"))
pdf(paste0(resu,"xpreview_opre_calcClusters.pdf"), width=10)
plot_grid(mu_myo, faps_only)
dev.off()
  
print(paste0("imputed celltypes (leveled vector) saved here : ", rdsdir,
             "opreFULL_transf_celltype_Vector.rds"))
print(all(filtered.seu@active.ident == filtered.seu@meta.data$seurat_clusters))
new.cluster.ids <- customLabelTransfer(typecell.authors.df, filtered.seu)
saveRDS(new.cluster.ids, paste0(rdsdir,"opreFULL_transf_celltype_Vector.rds"))

# END

