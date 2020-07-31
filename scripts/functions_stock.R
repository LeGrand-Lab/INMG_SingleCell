# Frequently used functions
# 
# --
# Joha 2020

doSCT_PCA_UMAP <- function(seu,numdim){
  seu <- SCTransform(seu, vars.to.regress = "percent.mt")
  seu <- RunPCA(seu)
  seu <- RunUMAP(seu, dims = 1:numdim)
  return(seu)
}

NormFeatScalePCA <- function(seu, nFeatRNA, percentmit){
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern= "^mt-")
  seu <- subset(seu, subset = nFeature_RNA > 200 & 
                  nFeature_RNA < nFeatRNA & percent.mt < percentmit)
  seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor=10000)
  seu <- FindVariableFeatures(seu, selection.method="vst", nfeatures= 2000)
  all.genes <- rownames(seu)
  seu <- ScaleData(seu, features = all.genes)
  seu <- RunPCA(seu,features=VariableFeatures(object = seu))  
  return(seu)
}

KNNplusFITSNE <- function(seu, nbDIM, resolu){
  #note: more resolution more groups,range 0.4-1.12
  seu <- FindNeighbors(seu, dims=1:nbDIM)
  seu <- FindClusters(seu, resolution=resolu)
  # Run FIt-SNE
  source("~/INMG_SingleCell/programsSC/FIt-SNE/fast_tsne.R", chdir=T)
  nbCELLS.div100 <- round(dim(seu)[2]/100) #nombre de cellules divisé par 100
  seu <- RunTSNE(object = seu,
                         perplexity=nbCELLS.div100, 
                         reduction="pca",
                         dims=1:nbDIM,
                         tsne.method = "FIt-SNE",
                         nthreads=4,
                         reduction.key="FItSNE_",
                         fast_tsne_path="~/INMG_SingleCell/programsSC/FIt-SNE/bin/fast_tsne",
                         max_iter=1000)
   return(seu)
}

definecolors <- function(celltype){ 
  # celltype is: seu@active.idents OR seu@meta.data$celltype
  print("for information: argument vector MUST HAVE LEVELS")
  c <- colorRampPalette(brewer.pal(8,"Set2"))(length(levels(celltype)))
  return(c)
}

colorsTransparent <- function(vectorlength){ 
  c <- colorRampPalette(c(rgb(1,1,1,0.5),rgb(1,0,0,0.5)))(vectorlength)
  return(c)
}

customTransferLabels <- function(markersDF, refDF, NBtop)  {
  # this function takes two dataframe (markers, Areference) and an integer
  # (i.e. how many top positive markers from markers tibble to use);
  # return vector string having celltypes, levels corresponding to cluster numbers:
  topmarkersDF <- markersDF %>% group_by(cluster) %>% top_n(n=NBtop,wt=avg_logFC)
  joined <- left_join(topmarkersDF, refDF, by = "gene") %>%
    select("cluster","gene","clusternames","avg_logFC") %>% distinct()
  bestmatch <- joined %>% group_by(cluster) %>%
    summarize(clusternames=names(which.max(table(clusternames))))
  vec = bestmatch$clusternames
  names(vec) = bestmatch$cluster
  return(vec)
}

doCustomImputeCelltype <- function(refDF, resu, markersinresults, delimiter,
                                   rdsdir, seufile, outsuffix){
  # Function for 'impute..' scripts, SAVES celltype_Vector.rds :
  seu <-readRDS(paste0(rdsdir,seufile))
  markersDF <- read.table(paste0(resu, markersinresults ), 
                          sep=delimiter,
                          header=TRUE)
  if(length(colnames(markersDF)) < 2) {
    print("ERROR: less than two columns in markersDF, delimiter is WRONG")
    return(1)
  }else {
    matchedtypes <- customTransferLabels(markersDF,refDF, 6) #using 6top+ markers
    seu <- RenameIdents(seu, matchedtypes)
    # save 'matchedtypes', a vector compatible with ..._seu_fitsne.rds
    print("saving celltypes (vector) specific to this seurat object into rds/")
    saveRDS(matchedtypes, file=paste0(rdsdir,outsuffix))
    return(seu) # returns seurat with clusters names
  }
}

NEWCustomImputeCelltype <- function(refDF, resu, markersinresults, delimiter,
                                   rdsdir, seufile, outsuffix, NBtop){
  # Function for 'impute..' scripts, SAVES celltype_Vector.rds :
  seu <-readRDS(paste0(rdsdir,seufile))

  markersDF <- read.table(paste0(resu, markersinresults ), 
                          sep=delimiter,
                          header=TRUE)
  if(length(colnames(markersDF)) < 2) {
    print("ERROR: less than two columns in markersDF, delimiter is WRONG")
    return(1)
  }else {
    matchedtypes <- customTransferLabels(markersDF,refDF, 6) #using 6top+ markers
    seu <- RenameIdents(seu, matchedtypes)
    # save 'matchedtypes', a vector compatible with ..._seu_fitsne.rds
    print("saving celltypes (vector) specific to this seurat object into rds/")
    saveRDS(matchedtypes, file=paste0(rdsdir,outsuffix))
    top.mk <- markersDF %>% group_by(cluster) %>% top_n(n=NBtop,wt=avg_logFC)
    mylist <- list(seu,top.mk)
    names(mylist) <- c("seu", "top.mk")
    return(list(seu, top.mk)) # returns seurat AND top markers
  }
}


getMatrixFromtxt <- function(filepath){
  m = read.table(filepath, sep="\t", header=T, row.names = 1)
  return(as.matrix(m))
}
getMatrixFrom10X <- function(dirpath){
  m = Read10X(data.dir=dirpath)
  return(as.matrix(m))
}
getMatrixFromCsv <- function(filepath){
  m = read.csv( filepath, sep=",", header=TRUE, row.names=1)
  return(as.matrix(m))
}

doSplitDeMicheli <- function(datadirD, filenm){
  dmizerol <- read.table(paste0(datadirD, filenm), sep="\t",
                         header=T, row.names=1)
  dim(dmizerol)
 print("NOTE: experiments on barcodes (colnames) prefixes")
  myexps = unique(substring(colnames(dmizerol), first=1, last=5)) #"D0_A_" "D0_B_" "D0_Cv"
  listofexps = sapply(myexps, function(x) {
    select(dmizerol, starts_with(x))
  })
  listofexps <- setNames(listofexps,myexps)
  nbcolssum = 0
  for (i in 1:length(listofexps)){ nbcolssum = nbcolssum + dim(listofexps[[i]])[2] }
  if( (dim(dmizerol)[2]) == nbcolssum){
    print(paste0("split operation sucessful, saving into ", datadirD))
    lapply( names(listofexps), function(k) { 
      write.table( listofexps[[k]],
                   paste0(datadirD, k,"DeMich.txt") , sep="\t", col.names = T);return(0)
    })
  }else{ print("errors in dimensions of splitted tables, check, nothing to save")
    return(1)}
}

knee_plot <- function(bc_rank) {
  knee_plt <- tibble(rank = bc_rank[["rank"]],
                     total = bc_rank[["total"]]) %>%
    distinct() %>%
    dplyr::filter(total > 0)
  annot <- tibble(inflection = metadata(bc_rank)[["inflection"]],
                  rank_cutoff = max(bc_rank$rank[bc_rank$total > metadata(bc_rank)[["inflection"]]]))
  p <- ggplot(knee_plt, aes(total, rank)) +
    geom_line() +
    geom_hline(aes(yintercept = rank_cutoff), data = annot, linetype = 2) +
    geom_vline(aes(xintercept = inflection), data = annot, linetype = 2) +
    scale_x_log10() +
    scale_y_log10() +
    annotation_logticks() +
    labs(y = "Rank", x = "Total UMIs")
  return(p)
}

runPrep_scran <- function(sce){
  rowData(sce) <- DataFrame(genes_names = rownames(sce))
  rowData(sce)$expressed <- scater::nexprs(sce,byrow=TRUE)>0
  per_cell <- perCellQCMetrics(sce[rowData(sce)$expressed,])
  colData(sce) <- cbind(colData(sce),per_cell)
  head(rowData(sce))
  sce <- scater::addPerFeatureQC(sce)
  sce <- computeSumFactors(sce) # by scran: scaling normalization (implements deconvolution)
  sce <- logNormCounts(sce)
  return(sce)
}
rundoublets_scran <- function(sce){
  dbl_dens <- doubletCells(sce) # *scran* doubletCells function
  sce$doublet_score <- 0
  sce$doublet_score <- log10(dbl_dens + 1)
  sce <- runTSNE(sce, perplexity=50, PCA=T, num_threads=4)
  return(sce)
}

doDimPlotHighlight <- function(seu, cellgrlist, mycolors, reduc, title){
  return(DimPlot(seu, label=T, repel=T, cells.highlight = cellgrlist ,
                 cols.highlight = mycolors, label.size= 3,
                 sizes.highlight = 0.2, pt.size=0.1, cols="darkgray",  
                reduction=reduc) + ggtitle(title) + 
           theme(legend.position="none"))
}

domono <- function(seu, experiment){
  print("two arguments to use this function: seurat_obj and experiment= RNA or SCT")
  gene_annotation <- data.frame(rownames(seu@assays[[experiment]]@data))
  rownames(gene_annotation) <- rownames(seu@assays[[experiment]]@data)
  colnames(gene_annotation) <- c("gene_short_name")
  gene_annotation = new("AnnotatedDataFrame",gene_annotation)
  cell_metadata = new("AnnotatedDataFrame",seu@meta.data)
  expression_matrix <- seu@assays[[experiment]]@data
  if (experiment == "RNA"){
    cds <- newCellDataSet(expression_matrix,
                          phenoData = cell_metadata,
                          featureData = gene_annotation,
                          expressionFamily=VGAM::negbinomial.size()) 
    # use option expressionFamily=VGAM::negbinomial.size() , when raw counts
  }else if (experiment == "SCT"){
    cds <- newCellDataSet(expression_matrix,
                          phenoData = cell_metadata,
                          featureData = gene_annotation) 
    #default is VGAM::vglmff(), use when input is normalized-scaled expression matrix
  }
  DelayedArray:::set_verbose_block_processing(TRUE)
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds, cores=4)
  cds <- preprocessCDS(cds, method="PCA", num_dim = 20,
                       norm_method = "log",
                       verbose = T, cores=4) #introduced in v3alpha, !yields ERROR if residualModelFormulaStr = "~Size_Factor + percent.mt",
  cds <- reduceDimension(cds, num_dim=20, reduction_method= 'tSNE', 
                         residualModelFormulaStr = "~Size_Factor + percent.mt",
                         verbose = T,cores=4)
  cds <- partitionCells(cds)
  cds <- learnGraph(cds, use_pca=TRUE, 
                    rge_method = "DDRTree", verbose=TRUE)
  return(cds)
}
