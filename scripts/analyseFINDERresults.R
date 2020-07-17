# comparing results DoubletFinder : integrated vs splitted analysis
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
rdsdir = "rds/doubletsD0spli_FINDER/"
outdir = "qcdoubl/spli_FINDER/"

integratedRDATAtablepath = "rds/doubletsD0/plotdoublet.Rdata"
load(integratedRDATAtablepath) # it is the RData  provided by Will: 'plotdoublet' object
head(plotdoublet)
is.data.frame(plotdoublet)
newsplires <- read.table(paste0(outdir, "TABLE_DOUBLFINDER_SPLITTED.txt"), sep="\t", header=T)
setdiff(colnames(newsplires),colnames(plotdoublet)) # its just TSNE added to newsplires, let 
print("checking groups")
unique(newsplires$groupe); unique(plotdoublet$groupe)

#prefixes in barcodes from 'plotdoublet' are:
#dorso1_wt1_; dorso_wt2_; gio_wt1_i1_; 	gio_wt2_i2_; demichili_D0_A_; demichili_D0_B_; demichili_D0_Cv3_
#ATTENTION: barcodes on demicheli have already D0_A_ , or _D0_B_ ,  or D0_Cv3_
newsplires <- newsplires %>% mutate(homoloid=case_when(
  groupe == "D0" ~ paste0("demichili_",id),
  TRUE ~ paste0(sample,"_",id)) ) 

length(intersect(newsplires$homoloid,plotdoublet$id))
#[1] 19143
dim(plotdoublet)
#[1] 19143    11

sel1 <- plotdoublet %>% select("id"=id, "DFclassINTE"=DFclass, "pANN_inte"=pANN)
sel2 <- newsplires %>% select("id"= homoloid, "DFclassSPLIT"=DFclass, "pANN_split"=pANN)

max(comparison$pANN_inte);max(comparison$pANN_split)
#[1] 0.001306051
#[1] 0
max(comparison$pANN_inte);max(comparison$pANN_split)
#[1] 0.6926426
#[1] 0.7219731

comparison = left_join(sel1, sel2, by="id")
dim(comparison)

h = 0.5 #known threshold high confidence
l = round(min(comparison$pANN_inte[comparison$DFclassINTE=="Doublets - Low confidence"]),2)
#[1] 0.3108402
a <- ggplot(comparison, aes(fill=DFclassINTE, x=pANN_inte)) + geom_histogram(alpha=0.35 , bins = 100) +
  ggtitle(paste0(
    "DoubletFinder on integrated object ( high conf , low conf  = ", h," , ", l,")" )) +
  theme(plot.title = element_text(size=9, face="bold"))

b <- ggplot(comparison, aes(fill=DFclassSPLIT, x=pANN_split)) + 
  geom_histogram(alpha=0.35 , bins = 100) +
  labs(title = paste0("DoubletFinder on SPLITTED object ( low conf varies across experiences* )" ),
       caption =  "* as expected, threshold is defined as quantile 92.5% of pANN, on each experience") +
  theme(plot.title = element_text(size=9, face="bold"), plot.caption = element_text(size=8))
griid <- plot_grid(a,b,nrow=2)
#ggdraw(add_sub(griid, "* as expected, threshold is defined dynamically by doubletFinder", lineheight = 0.05))


# Do Venn diagrams
vennlist1 <- list(
  "INTEdoublet_lo" = comparison$id[comparison$DFclassINTE=="Doublets - Low confidence"],
  "INTEdoublet_hi" = comparison$id[comparison$DFclassINTE=="Doublets - High confidence"],
  "SPLITdoublet_lo" = comparison$id[comparison$DFclassSPLIT=="Doublets - Low confidence"],
  "SPLITdoublet_hi" = comparison$id[comparison$DFclassSPLIT=="Doublets - High confidence"]
) # end list
myve1 <- venn.diagram(vennlist1, fill=0:3, alpha=0.3, filename=NULL, margin=0.1)
#grid.newpage()
#grid.draw(myve)
vennlist2 <- list(
  "INTEdoublet_lo" = comparison$id[comparison$DFclassINTE=="Doublets - Low confidence"],
 "INTEsinglet" = comparison$id[comparison$DFclassINTE=="Singlet"],
  "SPLITdoublet_lo" = comparison$id[comparison$DFclassSPLIT=="Doublets - Low confidence"],
 "SPLITsinglet" = comparison$id[comparison$DFclassSPLIT=="Singlet"]
) # end list
myve2 <- venn.diagram(vennlist2, fill=2:5, alpha=0.3, filename=NULL, margin=0.1)

vennlist3 <- list(
  "INTEdoublet_hi" = comparison$id[comparison$DFclassINTE=="Doublets - High confidence"],
  "INTEsinglet" = comparison$id[comparison$DFclassINTE=="Singlet"],
  "SPLITdoublet_hi" = comparison$id[comparison$DFclassSPLIT=="Doublets - High confidence"],
  "SPLITsinglet" = comparison$id[comparison$DFclassSPLIT=="Singlet"]
  ) # end list
myve3 <- venn.diagram(vennlist3, fill=4:7, alpha=0.3, filename=NULL, margin=0.1)

vennlist4 <- list(
  "INTEdoublet" = comparison$id[comparison$DFclassINTE=="Doublets - High confidence" |
                                  comparison$DFclassINTE=="Doublets - Low confidence"],
  "SPLITdoublet" = comparison$id[comparison$DFclassSPLIT=="Doublets - High confidence" |
                                   comparison$DFclassSPLIT=="Doublets - Low confidence"]
)
myve4 <- venn.diagram(vennlist4, fill=5:6, alpha=0.3, filename=NULL, margin=0.1)
vennsall <- plot_grid(myve1,myve2, myve3, myve4, nrow=2)

##
print("filtering cells discordantly classified by DoubletFinder")
discordant <- comparison %>% filter(DFclassINTE != DFclassSPLIT)
print("bringing the cluster numbers to this dataframe")
discordant <- left_join(discordant, plotdoublet, by="id")
mygeompoints <- 1
inpl <- ggplot(discordant, aes(x=UMAP_1,y=UMAP_2, label=cluster)) + geom_point(aes(colour=factor(DFclassINTE))) + 
  labs(title="Integrated Object") + theme(legend.position="none") + geom_text(check_overlap=T, size=3, alpha=0.7)

sppl <- ggplot(discordant, aes(UMAP_1,UMAP_2, label=cluster)) + geom_point(aes(colour=factor(DFclassSPLIT))) +
  labs(title="Separated experiments") + theme(legend.title=element_blank()) + geom_text(check_overlap=T, size=3, alpha=0.7)
inplsppl <- plot_grid(inpl,sppl,ncol=2, rel_widths = c(2,3.5))
title <- ggdraw() + draw_label(paste0("Discordantly classified cells, n=",
                                      as.character(dim(discordant)[1])), fontface='bold')
mygeompoins = plot_grid(title, inplsppl, ncol=1, rel_heights=c(0.1, 1))

##
### print plots
pdf(paste0(outdir,"plotsComparisonsFINDER.pdf"), width=10)
griid
vennsall
mygeompoints
dev.off()



# previous verson (NO labels)

# print("filtering cells discordantly classified by DoubletFinder")
# discordant <- comparison %>% filter(DFclassINTE != DFclassSPLIT)
# print("bringing the cluster numbers to this dataframe")
# discordant <- left_join(discordant, plotdoublet, by="id")
# mygeompoints <- 1
# inpl <- ggplot(discordant, aes(UMAP_1,UMAP_2)) + geom_point(aes(colour=factor(DFclassINTE))) + 
#   labs(title="Integrated Object") + theme(legend.position="none")
# 
# sppl <- ggplot(discordant, aes(UMAP_1,UMAP_2)) + geom_point(aes(colour=factor(DFclassSPLIT))) +
#   labs(title="Separated experiments") + theme(legend.title=element_blank())
# inplsppl <- plot_grid(inpl,sppl,ncol=2, rel_widths = c(2,3.5))
# title <- ggdraw() + draw_label(paste0("Discordantly classified cells, n=",
#                                       as.character(dim(discordant)[1])), fontface='bold')
# mygeompoints = plot_grid(title, inplsppl, ncol=1, rel_heights=c(0.1, 1))

