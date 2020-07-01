# Impute cell types to INTEGRATED D0 job
# input: 
#   txt file ==> all specific markers from experience
#   txt file => reference markers
#   rds file 
# output : .rds file containing corresponding celltypes (vector)
#        and figures
# --
# JohaGL

library(dplyr)
library(Seurat)
library(ggplot2)
library(sctransform)
library(cowplot)
library(RColorBrewer)
source(file="~/INMG_SingleCell/scripts/functions_stock.R",local=T)

prloc="~/INMG_SingleCell/"
setwd(prloc)

# colors to distinguish author origin:
ctrans <- c(rgb(0.81, 0.15, 0.49, alpha = 0.4), 
            rgb(0.15, 0.34, 0.95, alpha=0.4),
            rgb(0.4, 0.77, 0.50, alpha=0.5))


delimiter = "\t" # TAB delimited in this integrated MARKERS
#   Markers groups and celltypes as discussed with Dr Le Grand:
ref = "refmarkers/newRefmarkersToCells_v1.txt" # >>check
refDF <- read.table(ref,sep="\t",header=TRUE)
head(refDF,n=2)

# integrated: paths
resu="results/integratedD0/"
rdsdir = "rds/integratedD0/"
seufile = "integrated_seu_fitsne.rds"
markersinresults = "ALLMARKERS_Pos_integratedD0.txt"
outsuffix="integr_celltype_Vector.rds"
integrtypeseu <- doCustomImputeCelltype(refDF, resu, markersinresults, delimiter,
                                         rdsdir, seufile, outsuffix)

# PLOTS
plotbyauthor <- DimPlot(object = muscle.integrated, cols=ctrans, reduction = "tsne", 
                        group.by = "orig.ident")
plotbycluster <- DimPlot(integrtypeseu, label=T, repel=T, cols= definecolors(integrtypeseu@active.ident)) +
  ggtitle("Integrated D0")

pdf(paste0(resu,"FITSNE_named_byauthors.pdf"), width=14)
plot_grid(plotbyauthor, plotbycluster)
dev.off()

pdf(paste0(resu,"FITSNE_interestingFeatures.pdf"))
FeaturePlot(muscle.integrated, features=c("Smoc2","Pi16", "Sfrp5", "Hic1"))
FeaturePlot(muscle.integrated, features= c("Pax7","Cd34", "Notch2"))  
dev.off()
