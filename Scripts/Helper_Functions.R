#ggplot constants
tagSize <- 6

#Helper Functions

### Opposite of %in%
`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

### For sanity checking if things are the same
SameElements <- function(a, b) return(identical(sort(a), sort(b))) #https://stackoverflow.com/questions/27912800/check-whether-two-vectors-contain-the-same-unordered-elements-in-r


### Modified from the DESeq2 plotPCA function
plotPCAcustom = function(object, intgroup="condition", ntop=500, returnData=T)
{
  # calculate the variance for each gene
  rv <- rowVars(assay(object))
  
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))
  
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
  
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=":"))
  } else {
    colData(object)[[intgroup]]
  }
  
  # assemble the data for the plot
  d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], group=group, intgroup.df, name=colnames(object))
  attr(d, "percentVar") <- percentVar[1:2]
  return(d)
}

### For determining TPM cutoffs of groups
groupMedian <- function(df,colName, groupName, sampleTable){
  return(rowMedians(df[,colnames(df) %in% sampleTable[sampleTable[,colName] == groupName,]$sample]))
}


### For returning DE results from DEseq2

de <- function(dds, cont, pval, logFC = 0){
  temp <- data.frame(results(dds, contrast =cont))
  temp <- merge(temp[temp$padj < pval,],ensemblGeneListH,by.x=0,by.y="ensembl_gene_id")
  temp <- temp[temp$Row.names %in% row.names(highTPM),]
  temp <- temp[abs(temp$log2FoldChange) >= logFC,]
  return(temp[complete.cases(temp) ==T,])
}


## Functions for Processing RcisTarget data from DE lists
TFpaste <- function(TFlist, returnMotifType = T){
  finalReturn <- c()
  if(TFlist == ""){
    finalReturn <- ""
  }
  else{
    testSplit <- unlist(strsplit(unlist(TFlist), "[.]"))
    testSplit2 <- list()
    for(i in 1:length(testSplit)){
      testSplit2[i] <- strsplit(testSplit[i], " [(]|[)]")
    }
    if(returnMotifType == T){
      for(i in 1:length(testSplit2)){
        finalReturn <- c(finalReturn, paste(unlist(strsplit(testSplit2[[i]][1], "; ")), strsplit(testSplit2[[i]][2], "; ")))
      }
      finalReturn <- paste(finalReturn, collapse = "; ")
    }
    else{
      finalReturn <- paste(unlist(testSplit2)[c(seq(1,length(unlist(testSplit2)),2))], collapse = "; ")
    }
  }
  return(finalReturn)
}




TFidentify <- function(geneList, geneFilter){
  geneList <- as.character(geneList)
  motifs_AUC <- calcAUC(geneList, motifRankings, nCores=1)
  motifEnrichmentTable <- addMotifAnnotation(motifs_AUC, nesThreshold=3,
                                             motifAnnot=motifAnnotations)
  motifEnrichmentTable_wGenes <- addSignificantGenes(motifEnrichmentTable,
                                                     rankings=motifRankings, 
                                                     geneSets=geneList)
  motifTest <- motifEnrichmentTable_wGenes
  motifTest$TF_highConf_genes <- ""
  motifTest$TF_lowConf_genes <- ""
  for(i in 1:nrow(motifTest)){
    motifTest$TF_highConf_genes[i] <- TFpaste(motifTest$TF_highConf[i], F)
    motifTest$TF_lowConf_genes[i] <- TFpaste(motifTest$TF_lowConf[i], F)
  }
  motifTest <- motifTest[motifTest$TF_highConf != "" | motifTest$TF_lowConf != "",]
  for(i in 1:nrow(motifTest)){
    motifTest$genes[i] <- paste(c(motifTest$TF_highConf_genes[i], motifTest$TF_lowConf_genes[i]), collapse = "; ")
  }
  motifTest$genes <- gsub("^\\s+|\\s+$", "", motifTest$genes)
  motifTest$genes <- gsub("^\\;+|\\;+$", "", motifTest$genes)
  entrezList <- strsplit(motifTest$genes, split = "; ")
  motifEntrez <- data.frame()
  for(i in 1:length(entrezList)){
    motifEntrez <- rbind(motifEntrez,data.frame(Matrix = motifTest[i,2],Gene = entrezList[[i]], NES = motifTest[i,3], nEnrGenes = motifTest[i,7], enrichedGenes = motifTest[i,9]))
  }
  #motifEntrez <- merge(motifEntrez, geneFilter, by.x = 2, by.y =7)
  motifEntrez <- motifEntrez[motifEntrez$Gene %in% geneFilter,] 
  return(motifEntrez)
}


# Modified for sc Lists
TFidentify2 <- function(geneList, geneFilter = NULL, activity, window = "Both"){
  if(window %in% c("Both", "500")){
  motifRankings <- motifRankings500
  geneList <- as.character(geneList)
  motifs_AUC <- calcAUC(geneList, motifRankings, nCores=1)
  motifEnrichmentTable <- addMotifAnnotation(motifs_AUC, nesThreshold=3,
                                             motifAnnot=motifAnnotations)
  motifEnrichmentTable_wGenes <- addSignificantGenes(motifEnrichmentTable,
                                                     rankings=motifRankings, 
                                                     geneSets=geneList)
  motifTest <- motifEnrichmentTable_wGenes
  motifTest$TF_highConf_genes <- ""
  motifTest$TF_lowConf_genes <- ""
  for(i in 1:nrow(motifTest)){
    motifTest$TF_highConf_genes[i] <- TFpaste(motifTest$TF_highConf[i], F)
    motifTest$TF_lowConf_genes[i] <- TFpaste(motifTest$TF_lowConf[i], F)
  }
  motifTest <- motifTest[motifTest$TF_highConf != "" | motifTest$TF_lowConf != "",]
  for(i in 1:nrow(motifTest)){
    motifTest$genes[i] <- paste(c(motifTest$TF_highConf_genes[i], motifTest$TF_lowConf_genes[i]), collapse = "; ")
  }
  motifTest$genes <- gsub("^\\s+|\\s+$", "", motifTest$genes)
  motifTest$genes <- gsub("^\\;+|\\;+$", "", motifTest$genes)
  entrezList <- strsplit(motifTest$genes, split = "; ")
  motifEntrez <- data.frame()
  for(i in 1:length(entrezList)){
    motifEntrez <- rbind(motifEntrez,data.frame(Matrix = motifTest[i,2],Gene = entrezList[[i]], NES = motifTest[i,3], nEnrGenes = motifTest[i,7], enrichedGenes = motifTest[i,9]))
  }
  if(!is.null(geneFilter)){
    motifEntrez <- motifEntrez[motifEntrez$Gene %in% geneFilter,] 
  }
  motifEntrez$window <- "500bp up/100bp down"
  temp <- motifEntrez
  }
  if(window %in% c("10K", "Both")){
  motifRankings <- motifRankings10K
  geneList <- as.character(geneList)
  motifs_AUC <- calcAUC(geneList, motifRankings, nCores=1)
  motifEnrichmentTable <- addMotifAnnotation(motifs_AUC, nesThreshold=3,
                                             motifAnnot=motifAnnotations)
  motifEnrichmentTable_wGenes <- addSignificantGenes(motifEnrichmentTable,
                                                     rankings=motifRankings, 
                                                     geneSets=geneList)
  motifTest <- motifEnrichmentTable_wGenes
  motifTest$TF_highConf_genes <- ""
  motifTest$TF_lowConf_genes <- ""
  for(i in 1:nrow(motifTest)){
    motifTest$TF_highConf_genes[i] <- TFpaste(motifTest$TF_highConf[i], F)
    motifTest$TF_lowConf_genes[i] <- TFpaste(motifTest$TF_lowConf[i], F)
  }
  motifTest <- motifTest[motifTest$TF_highConf != "" | motifTest$TF_lowConf != "",]
  for(i in 1:nrow(motifTest)){
    motifTest$genes[i] <- paste(c(motifTest$TF_highConf_genes[i], motifTest$TF_lowConf_genes[i]), collapse = "; ")
  }
  motifTest$genes <- gsub("^\\s+|\\s+$", "", motifTest$genes)
  motifTest$genes <- gsub("^\\;+|\\;+$", "", motifTest$genes)
  entrezList <- strsplit(motifTest$genes, split = "; ")
  motifEntrez <- data.frame()
  for(i in 1:length(entrezList)){
    motifEntrez <- rbind(motifEntrez,data.frame(Matrix = motifTest[i,2],Gene = entrezList[[i]], NES = motifTest[i,3], nEnrGenes = motifTest[i,7], enrichedGenes = motifTest[i,9]))
  }
  if(!is.null(geneFilter)){
    motifEntrez <- motifEntrez[motifEntrez$Gene %in% geneFilter,] 
  }
  motifEntrez$window <- "10Kb up/10kb down"
  }
  if(window == "Both"){
    motifEntrez <- rbind(temp, motifEntrez)
  }
  motifEntrez$activity <- activity
  return(motifEntrez)
}

#### Plotting Functions

DimPlotCustom <- function(seurat, group.by = "orig.ident", pt.size = 2, plotLegend = T){
  embeddings <- as.data.frame(seurat@reductions$umap@cell.embeddings)
  pTemp <- ggplot(data=embeddings, aes(x=UMAP_1, y = UMAP_2)) + geom_point(aes(size = 0.5))
  yRange <- ggplot_build(pTemp)$layout$panel_scales_y[[1]]$range$range
  xRange <- ggplot_build(pTemp)$layout$panel_scales_x[[1]]$range$range
  embeddings$group <- seurat@meta.data[,group.by]
  p <- ggplot(data=embeddings, aes(x=UMAP_1, y=UMAP_2)) + geom_point(aes(fill= group), size = pt.size, colour = "black", stroke = .1, shape = 21) + theme_classic() +
    xlim(xRange) + ylim(yRange) + ylab("UMAP 2") + xlab("UMAP 1") + theme(plot.tag = element_text(size = 12), legend.position = "bottom", legend.direction = "horizontal", text = element_text(size = labelFont), legend.text = element_text(size = 6), plot.margin = unit(c(0,0,0,0), "cm"), plot.title = element_text(hjust = 0.5, size = titleFont))
  if(plotLegend == F){
    p <- p + theme(legend.position = "none")
  }
  return(p)
}

FeaturePlotCustom <- function(seurat, genes, plot = T, tag = element_blank(), plotLegend = T, cellNum = T, pt.size = .25, labelFont = 6, titleFont = 8, sharedScale = "All", nrow = NULL, ncol = NULL, split.by = NULL, sideBySide = T, color0 = "grey60", colPalette = c("dodgerblue2", "gold", "red2"), assay = "RNA"){
  embeddings <- as.data.frame(seurat@reductions$umap@cell.embeddings)
  if(assay == "RNA"){
    expr <- seurat@assays$RNA@data[genes,, drop = F]
  } else if(assay == "SCENIC"){
    expr <- seurat@assays$SCENIC@data[genes,, drop = F]
  }
  scatter_col = c("grey60",colorRampPalette(c("dodgerblue2", "gold", "red2"))(max(expr)*100))
  if(is.null(split.by)){
    p <- lapply(genes, function(x) {ggplot2::ggplot(data=embeddings[match(names(sort(expr[x,], decreasing = F)), row.names(embeddings)),], aes(x=UMAP_1, y=UMAP_2)) + ggplot2::geom_point(aes(color= sort(expr[x,], decreasing = F)), size = pt.size) +
        ggplot2::labs(col="Expression", title = x) + ggplot2::theme_classic() +
        ggplot2::ylab(element_blank()) + ggplot2::xlab(element_blank()) + ggplot2::theme(plot.tag = element_text(size = 12), text = element_text(size = labelFont), legend.text = element_text(size = 6), plot.margin = unit(c(0,0,0,0), "cm"), plot.title = element_text(hjust = 0.5, size = titleFont)) +
        if(sharedScale %in% c("All", "Gene")){
          ggplot2::scale_colour_gradientn(colors = c(color0, colorRampPalette(colPalette)(100)), limits = c(0, max(expr)))
        } else if(sharedScale == "None"){
          ggplot2::scale_colour_gradientn(colors = c(color0, colorRampPalette(colPalette)(100)), limits = c(0, max(expr[x,])))
        }})
  } else {
    xLimits <- c(min(embeddings$UMAP_1), max(embeddings$UMAP_1))
    yLimits <- c(min(embeddings$UMAP_2), max(embeddings$UMAP_2))
    splitDF <- seurat@meta.data[drop = F,,split.by]
    splits <- unique(splitDF[,1])
    q <- lapply(splits, function(y) {
      lapply(genes, function(x) {
        ggplot2::ggplot(data=embeddings[match(names(sort(expr[x,colnames(expr) %in% row.names(splitDF[splitDF[,split.by] %in% y,, drop = F])], decreasing = F)), row.names(embeddings)),], aes(x=UMAP_1, y=UMAP_2)) + ggplot2::geom_point(aes(color= sort(expr[x,colnames(expr) %in% row.names(splitDF[splitDF[,split.by] %in% y,, drop = F])], decreasing = F)), size = pt.size) +
          ggplot2::labs(col="Expression", title = paste0(x, " - ", y)) + ggplot2::theme_classic() +
          ggplot2::ylab(element_blank()) + ggplot2::xlab(element_blank()) + ggplot2::theme(plot.tag = element_text(size = 12), text = element_text(size = labelFont), legend.text = element_text(size = 6), plot.margin = unit(c(0,0,0,0), "cm"), plot.title = element_text(hjust = 0.5, size = titleFont)) +
          ggplot2::xlim(xLimits) + ggplot2::ylim(yLimits) +
          if(sharedScale == "All"){
            ggplot2::scale_colour_gradientn(colors = c(color0, colorRampPalette(colPalette)(100)), limits = c(0, max(expr)))
          } else if(sharedScale == "Gene"){
            ggplot2::scale_colour_gradientn(colors = c(color0, colorRampPalette(colPalette)(100)), limits = c(0, max(expr[x,])))
          } else if(sharedScale == "None"){
            ggplot2::scale_colour_gradientn(colors = c(color0, colorRampPalette(colPalette)(100)))
          }
      })})
    if(sideBySide == T & !is.null(split.by)){
      p <- list()
      loopCounter <- 1
      for(i in 1:length(genes)){
        for(j in 1:length(splits)) {
          p[[loopCounter]] <-q[[j]][i]
          loopCounter <- loopCounter + 1
        }
      }
    }
    p <- unlist(p, recursive = F)
  }
  if(plot == T){
    if(is.null(ncol) & is.null(nrow) & !is.null(split.by)){
      ncol = length(splits)
    }
    if(plotLegend == F){
      patchwork::wrap_plots(p, nrow = nrow, ncol = ncol) & theme(legend.position = "none")
    } else {
      patchwork::wrap_plots(p, nrow = nrow, ncol = ncol)
    }
  } else {
    return(p)
  }
}



