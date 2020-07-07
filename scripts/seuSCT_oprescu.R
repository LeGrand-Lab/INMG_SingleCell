#
#
# --
# Joha GL

library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(sctransform)
library(RColorBrewer)
library(viridis)
library(patchwork)
library(cowplot)

prloc = "~/INMG_SingleCell/"
prefix = "oprescu"
daysorder = c("0.5 DPI", "2 DPI", "3.5 DPI", "5 DPI", "10 DPI", "21 DPI", "Noninjured")
dayscols = viridis(length(daysorder), alpha=0.6)
resu = "results/Oprescu_" 

setwd(prloc)
dir.create(resu,recursive = T)