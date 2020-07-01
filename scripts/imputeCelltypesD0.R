# Impute cell types to D0 independent datasets
# input: 
#   txt file ==> all specific markers from experience
#   txt file => reference markers
#   rds file 
# output : .rds file containing corresponding celltypes (vector)
#         and figures
# --
# JohaGL
library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
source(file="~/INMG_SingleCell/scripts/functions_stock.R",local=T)

prloc="~/INMG_SingleCell/"
setwd(prloc)

#          *-* Global inputs *-*
delimiter = " " # espaced delimited results/../ALLMARKERS_....txt
#   Markers groups and celltypes as discussed with Dr Le Grand:
ref = "refmarkers/newRefmarkersToCells_v1.txt" # >>check
refDF <- read.table(ref,sep="\t",header=TRUE)
head(refDF,n=2)
#clusternames   gene
#1 B_lymphocytes.CD4   Cd74
#2 B_lymphocytes.CD4  Cd79a
#         *-* end Global inputs *-*

outsuffix="celltype_Vector.rds"

# ================================================================================
# imputing cell types to Giordani D0:
# ================================================================================
# 
resu = "results/GiordaniD0/"  # >> check
rdsdir = "rds/GiordaniD0/" # >> check
seufile = "gio_seu_fitsne.rds"   # >> check
markersinresults = "ALLMARKERS_GiordaniD0.txt"   # >> check

giotypeseu <- doCustomImputeCelltype(refDF, resu, markersinresults, delimiter,
                                     rdsdir, seufile, outsuffix)
pdf(paste0(resu, "FITSNE_named.pdf"))
DimPlot(giotypeseu, label=T,repel=T, 
        cols=definecolors(giotypeseu@active.ident))+
  ggtitle("Giordani D0")
dev.off()

# ================================================================================

# ================================================================================
# imputing cell types to DeMicheli D0:
# ================================================================================

resu2 = "results/DeMicheliD0/"  # >> check
rdsdir2 = "rds/DeMicheliD0/" # >> check
seufile2 = "dmizero_seu_fitsne.rds"   # >> check
markersinresults2 = "ALLMARKERS_DeMicheliD0.txt"   # >> check
delimiter = " " # space delimited results/../ALLMARKERS_....txt

demichtypeseu <- doCustomImputeCelltype(refDF,resu2,markersinresults2, delimiter,
                                        rdsdir2,seufile2, outsuffix)
pdf(paste0(resu2, "FITSNE_named.pdf"))
DimPlot(demichtypeseu, label=T,repel=T, 
        cols=definecolors(demichtypeseu@active.ident))+
  ggtitle("De Micheli D0")
dev.off()

# ================================================================================

# ================================================================================
# imputing cell types to DellOrso D0:
# ================================================================================
resu3 = "results/DellOrsoD0/"  # >> check
rdsdir3 = "rds/DellOrsoD0/" # >> check
seufile3 = "dorso_seu_fitsne.rds"   # >> check
markersinresults3 = "ALLMARKERS_DellOrsoD0.txt"   # >> check
delimiter = " " # space delimited results/../ALLMARKERS_....txt

dorsotypeseu <- doCustomImputeCelltype(refDF, resu3, markersinresults3, delimiter,
                                       rdsdir3, seufile3, outsuffix)


pdf(paste0(resu3, "FITSNE_named.pdf"))
DimPlot(dorsotypeseu, label=T,repel=T, 
        cols=definecolors(dorsotypeseu@active.ident))+
  ggtitle("DellOrso D0")
dev.off()
# ================================================================================
# END
# ================================================================================

# === extended version used first time for giordani in place of the local Function:
# markersDF <- read.table(paste0(resu, markersinresults ), 
#                         sep=delimiter,
#                         header=TRUE)
# seu <-readRDS(paste0(rdsdir,seufile))
# 
# # use function from functions_stock.R :
# matchedtypes <- customTransferLabels(markersDF,refDF, 6) #using 6top+ markers
# print(tail(matchedtypes,n=3))
# #  15                        16                        17 
# # "B_lymphocytes.CD4" "Neutrophyls_macrophages"               "Myonuclei" 
# 
# # check cell types have been imputed already or not
# unique(seu@active.ident)
# # impute cell types using 'matchedtypes' vector
# seu <- RenameIdents(seu, matchedtypes)
# unique(seu@active.ident)  # we see only 13 levels because homonyms are merged
# 
# mycols <- colorRampPalette(brewer.pal(8,"Set2"))(length(levels(seu@active.ident)))
# 
# pdf(paste0(resu, "FITSNE_named.pdf"))
# DimPlot(seu, label=T,repel=T, cols=mycols)
# dev.off()
# 
# # save 'matchedtypes', a vector compatible with ..._seu_fitsne.rds
# saveRDS(matchedtypes, file=paste0(rdsdir,outsuffix))

