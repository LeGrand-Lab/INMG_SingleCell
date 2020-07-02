#!/usr/bin/env Rscript
# --
# JohaGL

library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
source(file="~/INMG_SingleCell/scripts/functions_stock.R",local=T)

prloc="~/INMG_SingleCell/"
setwd(prloc)

resu = "results/OprescuD0/"
dir.create(resu,recursive = T)
rdsdir = "rds/OprescuD0/"
dir.create(rdsdir,recursive = T)

opre.data <- read.table("data/OprescuD0/oprescu_Noninjured_raw.txt", sep="\t",
                        header=T, row.names = 1)

opre <- CreateSeuratObject(opre.data, project="Oprescu", min.cells=3, min.features=200)

opre <- NormFeatScalePCA(opre, 5000, 10) # Oprescu(https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4120117):"
#Cells with more than 15% mitochondrial reads and less than 200 UMIs were excluded from downstream analysis "
# Joha: I tested 5% filtering, too stringent: 16886 features across 2472 samples.

pdf(paste0(resu,"preprocessSeu.pdf"))
VlnPlot(opre, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(opre, feature1="nCount_RNA", feature2 = "nFeature_RNA")
DimHeatmap(opre, dims=1:18, cells = 500, balanced = TRUE) 
ElbowPlot(opre)
dev.off()

opre <- KNNplusFITSNE(opre, 15, 0.5)

opre.markers <- FindAllMarkers(opre, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

topn <- opre.markers %>% group_by(cluster) %>% top_n(n = 4, wt = avg_logFC)

# print results (plots)

pdf(paste0(resu,"FITSNE.pdf"))
DimPlot(opre, reduction="tsne", label=TRUE, pt.size=0.5, 
        cols=definecolors(opre@active.ident), shape.by="orig.ident") + ggtitle("Oprescu D0")
dev.off()

pdf(paste0(resu,"HEATMAP.pdf"),width=13)
DoHeatmap(opre, features = topn$gene, group.colors=definecolors(opre@active.ident)) + NoLegend()
dev.off()

# save .rds object
saveRDS(opre,file=paste0(rdsdir,"opre_seu_fitsne.rds"))

#write table of all markers
write.table(opre.markers, paste0(resu,"ALLMARKERS_OprescuD0.txt"))
print("finished")

