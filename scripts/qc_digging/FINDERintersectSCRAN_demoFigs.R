# comparing dFinder vs scran::doublet
# ON SPLITTED EXPERIMENTS
# --
# JohaGL
library(VennDiagram)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(cowplot)

prloc="~/INMG_SingleCell/"
setwd(prloc)


scranres <- read.table("qcdoubl/spli_SCRAN/TABLE_DOUBLETS_SCRAN_splitted.txt",
                       sep="\t", header=T)
finderres <- read.table("qcdoubl/spli_FINDER/TABLE_DOUBLFINDER_SPLITTED.txt", 
                        sep="\t", header=T)
print("it is normal if scran results has many more cells (no seurat filters applied yet)")
paste("scran= ", dim(scranres)[1],"  Finder= ",dim(finderres)[1])

# if problems , fix to obtain homologated prefixes to barcodes:
#ATTENTION: barcodes on demicheli have already D0_A_ , or _D0_B_ ,  or D0_Cv3_

scranres <- scranres %>% mutate(newid = case_when(
  filerds == "dorsowt1" ~ 
    paste0(as.character(filerds),"_",as.character(barcode)),
  filerds == "dorsowt2" ~ 
    paste0(as.character(filerds),"_",as.character(barcode)),
  TRUE ~ as.character(barcode)
))
finderres <- finderres %>% mutate(newid = case_when(
  groupe == "dorso1_wt1" ~ paste0("dorsowt1","_",as.character(id)),
  groupe == "dorso_wt2" ~ paste0("dorsowt2","_",as.character(id)),
  TRUE ~ as.character(id)
))
print("checking that no duplicated rows exist ( distinct() ):")
dim(finderres %>% distinct()); dim(scranres %>% distinct())
comparedf <- left_join(finderres, scranres, by="newid")

print("equivalences in classification columns")
unique(comparedf$classific)
unique(comparedf$DFclass)
vennliA <- list(
  "scran_doub_lo" = comparedf$newid[comparedf$classific == "doublets low conf"],
  "Finder_doub_lo" = comparedf$newid[comparedf$DFclass == "Doublets - Low confidence"],
  "scran_doub_hi" = comparedf$newid[comparedf$classific == "doublets high conf"],
  "Finder_doub_hi" = comparedf$newid[comparedf$DFclass == "Doublets - High confidence"]
)

vennliB <- list(
  "scran_DOUBLET" = comparedf$newid[comparedf$classific == "doublets low conf" |
                            comparedf$classific == "doublets high conf"] ,
  "Finder_DOUBLET" = comparedf$newid[comparedf$DFclass == "Doublets - Low confidence" |
                            comparedf$DFclass ==  "Doublets - High confidence"] ,
  "scran_singlet" = comparedf$newid[comparedf$classific == "singlet"] ,
  "Finder_singlet" = comparedf$newid[comparedf$DFclass == "Singlet"]
)

a <- venn.diagram(vennliA, fill=0:3, alpha=0.3, filename=NULL, margin=0.1)
b <- venn.diagram(vennliB, fill=2:5, alpha=0.3, filename=NULL, margin=0.1)
grid.newpage()

print("getting barcodes in the intersection")
trueids_doubl = intersect(vennliB[["scran_DOUBLET"]], vennliB[["Finder_DOUBLET"]])

LIKELYDOUBL <- comparedf %>% filter(newid %in% trueids_doubl)

c <- ggplot(coso, aes(x=filerds, y=pANN)) + geom_boxplot() + 
  labs(title=paste0("DOUBLETS (intersection Finder and scran methods) n= ", dim(LIKELYDOUBL)[1]))
c2 <- ggplot(coso, aes(x=filerds, y=doublet_score)) + geom_boxplot()
print("doublets proportion in each sample:")
prop = table(LIKELYDOUBL$filerds) / ( table(comparedf$filerds) / 100 )
tprop = dim(LIKELYDOUBL)[1]/(dim(comparedf)[1]/100)
# D0_A_DeMich D0_B_DeMich D0_CvDeMich    dorsowt1    dorsowt2  GSM3520458  GSM3520459 
# 0.2702703   1.9438445   3.8647343   0.7522837   1.3066202   2.4089936   2.2005295 
d <- ggplot(as.data.frame(prop), aes(Var1,Freq)) + geom_col(fill="turquoise4") +  
  labs(y="proportion of calculated doublets by sample (%)", caption = paste("total %=",tprop)) + coord_flip()

umapinter1 <- ggplot(LIKELYDOUBL, aes(x=UMAP_1,y=UMAP_2, label=cluster)) + geom_point(aes(colour=factor(classific))) + 
  labs(title="scran::doubletCells") + theme(legend.position="none") + geom_text(check_overlap=T, size=3, alpha=0.7)
umapinter2 <- ggplot(LIKELYDOUBL, aes(UMAP_1,UMAP_2, label=cluster)) + geom_point(aes(colour=factor(DFclass))) +
  labs(title="DoubletFinder") + theme(legend.title=element_blank()) + geom_text(check_overlap=T, size=3, alpha=0.7)

clus1 <- ggplot(comparedf, aes(x=cluster,fill=classific)) + geom_bar(position="fill")+
  theme(legend.position="none")
clus2 <- ggplot(comparedf, aes(x=cluster,fill=DFclass)) + geom_bar(position="fill")


pdf("qcdoubl/FINDERintersectSCRAN.pdf", width=10)
plot_grid(a, b, ncol=2)
plot_grid(c,c2,d, ncol = 1, rel_widths = c(4,4,2))
plot_grid(umapinter1,umapinter2, clus1 , clus2, ncol=2, rel_widths = c(3,4))
dev.off()


