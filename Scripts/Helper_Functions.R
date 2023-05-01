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
