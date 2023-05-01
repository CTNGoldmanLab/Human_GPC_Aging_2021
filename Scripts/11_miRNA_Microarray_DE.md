Differential Expression Analysis of Adult vs Fetal hGPC miRNAs
================
John Mariani
03/06/23

``` r
source("Scripts/Helper_Functions.R")
library(limma)
library(pd.mirna.3.0)
```

    ## Loading required package: Biostrings

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following object is masked from 'package:limma':
    ## 
    ##     plotMA

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
    ##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
    ##     Position, rank, rbind, Reduce, rownames, sapply, setdiff, sort,
    ##     table, tapply, union, unique, unsplit, which.max, which.min

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## Loading required package: XVector

    ## Loading required package: GenomeInfoDb

    ## 
    ## Attaching package: 'Biostrings'

    ## The following object is masked from 'package:base':
    ## 
    ##     strsplit

    ## Loading required package: RSQLite

    ## Loading required package: oligoClasses

    ## Welcome to oligoClasses version 1.60.0

    ## Loading required package: oligo

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## ================================================================================

    ## Welcome to oligo version 1.62.2

    ## ================================================================================

    ## 
    ## Attaching package: 'oligo'

    ## The following object is masked from 'package:limma':
    ## 
    ##     backgroundCorrect

    ## Loading required package: DBI

``` r
library(miRBaseConverter)
library(miRNAtap)
```

    ## Loading required package: AnnotationDbi

    ## 
    ## Attaching package: 'miRNAtap'

    ## The following object is masked from 'package:Biostrings':
    ## 
    ##     translate

``` r
library(tidyr)
```

    ## 
    ## Attaching package: 'tidyr'

    ## The following object is masked from 'package:S4Vectors':
    ## 
    ##     expand

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following object is masked from 'package:miRNAtap':
    ## 
    ##     select

    ## The following object is masked from 'package:AnnotationDbi':
    ## 
    ##     select

    ## The following object is masked from 'package:oligo':
    ## 
    ##     summarize

    ## The following object is masked from 'package:Biobase':
    ## 
    ##     combine

    ## The following objects are masked from 'package:Biostrings':
    ## 
    ##     collapse, intersect, setdiff, setequal, union

    ## The following object is masked from 'package:GenomeInfoDb':
    ## 
    ##     intersect

    ## The following object is masked from 'package:XVector':
    ## 
    ##     slice

    ## The following objects are masked from 'package:IRanges':
    ## 
    ##     collapse, desc, intersect, setdiff, slice, union

    ## The following objects are masked from 'package:S4Vectors':
    ## 
    ##     first, intersect, rename, setdiff, setequal, union

    ## The following objects are masked from 'package:BiocGenerics':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(plyr)
```

    ## ------------------------------------------------------------------------------

    ## You have loaded plyr after dplyr - this is likely to cause problems.
    ## If you need functions from both plyr and dplyr, please load plyr first, then dplyr:
    ## library(plyr); library(dplyr)

    ## ------------------------------------------------------------------------------

    ## 
    ## Attaching package: 'plyr'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     arrange, count, desc, failwith, id, mutate, rename, summarise,
    ##     summarize

    ## The following object is masked from 'package:oligo':
    ## 
    ##     summarize

    ## The following object is masked from 'package:XVector':
    ## 
    ##     compact

    ## The following object is masked from 'package:IRanges':
    ## 
    ##     desc

    ## The following object is masked from 'package:S4Vectors':
    ## 
    ##     rename

``` r
library(trqwe)
```

    ## 
    ## Attaching package: 'trqwe'

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     chop

    ## The following object is masked from 'package:oligo':
    ## 
    ##     se

``` r
library(xlsx)
library(ggplot2)
library(ggfortify)
library(ggrepel)
library(patchwork)
```

## Load

``` r
de_intersect <- read.delim("output/de_Adult_vs_Fetal_Intersect.txt")
supTable3a <- read.xlsx("Extended Data Tables/Extended Data Table 3 - Adult vs Fetal hGPC Bulk RNA-seq.xlsx", sheetName = "Adult vs Fetal hGPC DE")
adultRepressorNetwork <- read.table("output/adultRepressorNetwork.txt", header = T)
adultActivatorNetwork <- read.table("output/adultActivatorNetwork.txt", header = T)
fetalRepressorNetwork <- read.table("output/fetalRepressorNetwork.txt", header = T)
fetalActivatorNetwork <-  read.table("output/fetalActivatorNetwork.txt", header = T)

activators <- c("MYC", "NFIB", "STAT3", "HMGA2", "TEAD2")
```

## Load data and conduct differential expression with limma

``` r
#Function used throughout

dataFolder <- "data_for_import/miRNA/"

annotation <- read.csv("data_for_import/miRNA-3_0-st-v1.annotations.201405132.csv")


#Read in and run RMA on CEL files
fns <- list.celfiles(path = dataFolder)
fns
```

    ## [1] "20.3_wk_CD140a+_(miRNA-3_0).CEL" "20.5_wk_CD140a+_(miRNA-3_0).CEL"
    ## [3] "21_wk_CD140a+_(miRNA-3_0).CEL"   "21.1_wk_CD140a+_(miRNA-3_0).CEL"
    ## [5] "SID_372_(miRNA-3_0).CEL"         "SID348_(miRNA-3_0).CEL"         
    ## [7] "SID390_(miRNA-3_0).CEL"          "SID398_(miRNA-3_0).CEL"

``` r
Data <- read.celfiles(filenames=paste0(dataFolder,fns))
```

    ## Platform design info loaded.

    ## Reading in : data_for_import/miRNA/20.3_wk_CD140a+_(miRNA-3_0).CEL
    ## Reading in : data_for_import/miRNA/20.5_wk_CD140a+_(miRNA-3_0).CEL
    ## Reading in : data_for_import/miRNA/21_wk_CD140a+_(miRNA-3_0).CEL
    ## Reading in : data_for_import/miRNA/21.1_wk_CD140a+_(miRNA-3_0).CEL
    ## Reading in : data_for_import/miRNA/SID_372_(miRNA-3_0).CEL
    ## Reading in : data_for_import/miRNA/SID348_(miRNA-3_0).CEL
    ## Reading in : data_for_import/miRNA/SID390_(miRNA-3_0).CEL
    ## Reading in : data_for_import/miRNA/SID398_(miRNA-3_0).CEL

``` r
eset <- oligo::rma(Data)
```

    ## Background correcting
    ## Normalizing
    ## Calculating Expression

``` r
edata <- data.frame(exprs(eset))

humanMirAnnotation <- annotation[annotation$Species.Scientific.Name == "Homo sapiens" & annotation$Sequence.Type == "miRNA",]


# Filter on list of human mir probes
edata_mir <- edata[row.names(edata) %in% humanMirAnnotation$Probe.Set.ID,]

mirData <- data.frame(sample = names(edata_mir), group = c(rep("Fetal",4), rep("Adult",4)))


design <- model.matrix(~0+group, mirData)
colnames(design) <- make.names(colnames(design))

fitV <- lmFit(edata_mir, design)
fitV2 <- eBayes(fitV)
colnames(design) 
```

    ## [1] "groupAdult" "groupFetal"

``` r
cont.matrix <- makeContrasts(comparison = groupAdult - groupFetal,
                             levels=design)

fitV2 <- contrasts.fit(fitV2, cont.matrix)
fitV2 <- eBayes(fitV2)

mirFetalAdultAll <- topTable(fitV2, coef = 1,number = 1000000000, p.value = 0.01)


humanMirAnnotationDE <- humanMirAnnotation[humanMirAnnotation$Probe.Set.ID %in% row.names(mirFetalAdultAll),]
humanMirAnnotationDE[,12:13] <- miRNA_AccessionToName(humanMirAnnotationDE$Accession,targetVersion = "v22")

#Fails to find hsa-miR-3656 appropriately so is filled in manually
humanMirAnnotationDE[is.na(humanMirAnnotationDE$TargetName),]$TargetName <- "hsa-miR-3656"
humanMirAnnotationDE$mirName <- gsub(humanMirAnnotationDE$TargetName, pattern = "hsa-", replacement = "")

mirFetalAdultAll <- merge(mirFetalAdultAll,humanMirAnnotationDE, by.x=0,by.y=1)

### Make Sup Table 5a
supTable5a <- mirFetalAdultAll[,c(8:10,2,6)]
names(supTable5a) <- c("Probe_Set_Name", "Accession", "Transcript_ID_Array_Design", "Log2FC_Adult_vs_Fetal_GPC", "Adj_P_Val")
supTable5a <- supTable5a[order(supTable5a$Adj_P_Val, decreasing = F),]
write.xlsx(supTable5a, file = "Extended Data Tables/Extended Data Table 5 - Adult vs Fetal miRNA Microarray.xlsx", sheetName = "DE Adult vs Fetal GPC miRNAs", row.names = F)
```

## Find miRNA targets using miRNAtap

``` r
humanMirsDE <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("rank_product", "rank_final", "mir", "entrez"))
for(i in 1:length(mirFetalAdultAll$mirName)){
  tryCatch({
    #print(mirFetalAdultAll$mirName[i])
    temp <- as.data.frame(getPredictedTargets(as.character(mirFetalAdultAll$mirName[i]), species = 'hsa',method = 'geom', min_src = 2))
    temp$mir <- as.character(mirFetalAdultAll$mirName[i])
    temp$entrez <- row.names(temp)
    temp <- temp[,c("rank_product", "rank_final", "mir", "entrez")]
    humanMirsDE <- rbind(humanMirsDE, temp)
  }, error=function(e){})
}
```

    ## Loading required package: miRNAtap.db

``` r
table(mirFetalAdultAll$logFC > 0)
```

    ## 
    ## FALSE  TRUE 
    ##    33    23

``` r
mirPCA <- autoplot(prcomp(t(edata_mir)), data = mirData, size = 1) + geom_point(shape = 21, size = 5, aes(fill =  mirData$group)) + theme_bw() + ylim(-.8,.9) + theme(legend.position = "bottom", legend.direction = "horizontal") + scale_fill_manual(values = c("#C40000", "#00008B"))

mirPCA
```

![](10_miRNA_Microarray_DE_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
#write.table(mirFetalAdult, "output/de_mirFetalAdultAll.txt", sep ="\t", row.names = F, quote = F)
```

## Determine proper direction miRs

``` r
filename="data_for_import/ensembl_miR_list.csv"
if(file.exists(filename)){
  ensemblGeneListMir <- read.csv(filename)} else{
    marth <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = 'http://jan2019.archive.ensembl.org/', ensemblRedirect = T)
    ensemblGeneListMir <- getBM(attributes = c("entrezgene", "external_gene_name", "gene_biotype", "description"), filters = "entrezgene",unique(row.names(humanMirsDE)), mart = marth)
    write.csv(ensemblGeneListMir, filename, row.names = F)
  } 

nrow(humanMirsDE[unique(humanMirsDE$entrez),])
```

    ## [1] 10919

``` r
humanMirsLabeled <- merge(humanMirsDE, ensemblGeneListMir, by.x = "entrez", by.y = "entrezgene")




#### 
mirFetalAdultLabeled <- merge(mirFetalAdultAll, humanMirsDE, by.x = 20, by.y = 3)
mirFetalAdultLabeled <- merge(mirFetalAdultLabeled, ensemblGeneListMir, by.x = "entrez", by.y = "entrezgene")

adult_vs_fetal_mir <- de_intersect
adult_vs_fetal_mir <- adult_vs_fetal_mir[,c(1,3,7,8,9,10)]
names(adult_vs_fetal_mir) <- c("ensembl_id", "a2b5_FC", "a2b5_padj", "external_gene_name", "gene_biotype", "description")

adult_vs_fetal_mir <- merge(adult_vs_fetal_mir, mirFetalAdultLabeled, by.x = "external_gene_name", by.y = "external_gene_name")
adult_vs_fetal_mir$X <- NULL
adult_vs_fetal_mir <- adult_vs_fetal_mir[!duplicated(adult_vs_fetal_mir),]


mirDownGeneUpFinal<- adult_vs_fetal_mir[adult_vs_fetal_mir$a2b5_FC > 0 & adult_vs_fetal_mir$logFC < 0,]
mirUpGeneDownFinal<- adult_vs_fetal_mir[adult_vs_fetal_mir$a2b5_FC< 0 & adult_vs_fetal_mir$logFC > 0,]

mirProperDirection <- rbind(mirDownGeneUpFinal,mirUpGeneDownFinal)


### Make Sup Table 5b
supTable5b <- mirProperDirection
supTable5b <- supTable5b[,c(8,10,1,3)]
names(supTable5b) <- c("miR_Name", "miR_Log2FC_Adult_vs_Fetal_GPC", "External_Gene_Name", "Target_Gene_Log2FC_Adult_vs_Fetal_A2B5_GPC")
write.xlsx(supTable5b, file = "Extended Data Tables/Extended Data Table 5 - Adult vs Fetal miRNA microarray.xlsx", sheetName = "Predicted miRNA targets", row.names = F, append = T)
```

## Make miRNA Dot Plot

``` r
mirDotPlot <- mirProperDirection %>%
  group_by(mirName) %>%
  dplyr::summarise(logFCmean = mean(a2b5_FC,na.rm = TRUE), mirLogFC = mean(logFC), n = n())

mirDotPlotFetal <- mirDotPlot[mirDotPlot$mirLogFC < 0,]
mirDotPlotAdult <- mirDotPlot[mirDotPlot$mirLogFC > 0,]

mean(mirDotPlotFetal$n)
```

    ## [1] 36.34375

``` r
mean(mirDotPlotAdult$n)
```

    ## [1] 46.3913

``` r
sd(mirDotPlotFetal$n)
```

    ## [1] 24.53255

``` r
sd(mirDotPlotAdult$n)
```

    ## [1] 37.81623

``` r
signif(length(unique(mirDownGeneUpFinal$external_gene_name)) / length(de_intersect[de_intersect$log2FoldChange > 0,]$external_gene_name)*100,3)
```

    ## [1] 48.8

``` r
signif(length(unique(mirUpGeneDownFinal$external_gene_name)) / length(de_intersect[de_intersect$log2FoldChange < 0,]$external_gene_name)*100, 3)
```

    ## [1] 39.9

``` r
length(unique(mirUpGeneDownFinal$external_gene_name))
```

    ## [1] 663

``` r
miRinteresting <- c("miR-9-3p", "miR-9-5p", "miR-193b-3p", "miR-338-5p", "miR-24-3p", "miR-193a-5p", "miR-31-5p", "miR-584-5p", "miR-330-3p", "miR-409-3p", "miR-379-5p", "miR-432-5p", "miR-219a-2-3p")

mirDotPlot$label <- ""
mirDotPlot[mirDotPlot$mirName %in% miRinteresting,]$label <- mirDotPlot[mirDotPlot$mirName %in% miRinteresting,]$mirName



mirGGdotPlot <- ggplot(mirDotPlot, aes(y = logFCmean, x = mirLogFC, label = label)) + geom_point(shape = 21, aes(fill = mirLogFC, stroke = 1, size = n),  alpha = 0.5) +  scale_size(range = c(0.5, 12)) + theme_bw() + geom_text_repel() +
  scale_fill_gradient2(midpoint = 0, low = "#00008B", mid = "lightgrey", high = "#C40000", space = "Lab", guide = guide_colourbar(direction = "horizontal", title = "Adult vs. Fetal miRNA log2FC", title.position = "top")) + theme(legend.position = "bottom") + xlab("miRNA Log2FC") + ylab("Average Predicted Target Log2FC")


mirGGdotPlot
```

    ## Warning: ggrepel: 1 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](10_miRNA_Microarray_DE_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->
\### miR HM

``` r
mirHM <- edata_mir %>% set_rownames(row.names(.))  %>%
  as_tibble(rownames = "row") %>%
  pivot_longer(-row, names_to = "Sample", values_to = "Intensity")


#interestingAFmiRtargets <- mirProperDirection[mirProperDirection$external_gene_name %in% afHMgenes,]

mirHM <- mirHM[mirHM$row %in% mirFetalAdultAll$Row.names, ]

mirHM$group <- mapvalues(mirHM$Sample, mirData$sample, as.character(mirData$group))

mirBar <- mirFetalAdultAll
mirBar$logFC <- as.numeric(mirBar$logFC)
mirBar$Row.names <- as.factor(mirBar$Row.names)
mirBar$direction <- "Fetal"
mirBar[mirBar$logFC > 0,]$direction <- "Adult"


#Grab SE from limma
mirSE <- fitV2$stdev.unscaled * fitV2$sigma
colnames(mirSE) <- "StandardError"
mirBar <- merge(mirBar, mirSE, by.x = "Row.names", by.y = 0)
mirBar <- mirBar[order(mirBar$logFC),]

mirBar$shortName <- gsub("hsa-", "", mirBar$Transcript.ID.Array.Design.)
mirBar$shortName <- factor(mirBar$shortName, levels = mirBar$shortName)
mirHM$row <- factor(mirHM$row, levels = mirBar$Row.names)




m1 <- ggplot(mirHM, aes(row, Sample)) + theme_bw() + geom_tile(aes(fill = mirHM$Intensity)) + theme(legend.position = "bottom", legend.direction = "horizontal", axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())  + scale_fill_gradientn(colours = c("#009900","#fffcbd","#ff2020")) + facet_grid(vars(group),scales = "free", space = "free", switch = "y") + scale_y_discrete(expand = c(0,0))

m2 <- ggplot(mirBar, aes(x = shortName, y = logFC, fill = direction)) + 
   geom_errorbar(aes(ymin=logFC-StandardError, ymax=logFC+StandardError), width=.2,
                 position=position_dodge(.9)) +
  geom_col(colour = "black") + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.x = element_blank(), legend.position = "none") + scale_fill_manual(values = c("#C40000", "#00008B"))


m2
```

![](10_miRNA_Microarray_DE_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
(((mirPCA /mirGGdotPlot) + plot_layout(heights = c(1,1.75))) | (m2 / m1)) + plot_layout(widths = c(1,2))
```

    ## Warning: Use of `mirHM$Intensity` is discouraged.
    ## ℹ Use `Intensity` instead.

    ## Warning: ggrepel: 12 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](10_miRNA_Microarray_DE_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->

``` r
#ggsave("panels/mirSignaling.pdf", units = "in", dpi = 300, width = 16, height = 8, device = NULL)
```

## Make miRNA Network

``` r
miRsCurated <- read.csv("data_for_import/miRNAtap_Curated_Network.csv")

mirAdultEdges <- mirProperDirection[mirProperDirection$mirName %in% miRsCurated$Adult_miR_Network & mirProperDirection$external_gene_name %in% miRsCurated$Adult_miR_Network,]

mirAdultEdges <- mirAdultEdges[,c(8,1)]

mirFetalEdges <- mirProperDirection[mirProperDirection$mirName %in% miRsCurated$Fetal_miR_Network & mirProperDirection$external_gene_name %in% miRsCurated$Fetal_miR_Network,]

mirFetalEdges <- mirFetalEdges[,c(8,1)]

mirFetalEdges <- mirFetalEdges[!duplicated(mirFetalEdges[, c("mirName", "external_gene_name")]),]

mirAdultEdges <- mirAdultEdges[!duplicated(mirAdultEdges[, c("mirName", "external_gene_name")]),]

mirAdultNodes <- data.frame(node = unique(c(mirAdultEdges$mirName, mirAdultEdges$external_gene_name)))
mirFetalNodes <- data.frame(node = unique(c(mirFetalEdges$mirName, mirFetalEdges$external_gene_name)))

mirAdultNodes$type <- ifelse(mirAdultNodes$node %in% mirFetalAdultLabeled$mirName, "miRNA", "Gene Target")
mirFetalNodes$type <- ifelse(mirFetalNodes$node %in% mirFetalAdultLabeled$mirName, "miRNA", "Gene Target")

mirAdultNodes$enriched <- ifelse(mirAdultNodes$node %in% mirFetalAdultLabeled$mirName, "Adult", "Fetal")
mirFetalNodes$enriched <- ifelse(mirFetalNodes$node %in% mirFetalAdultLabeled$mirName, "Fetal", "Adult")

write.csv(mirAdultEdges, "output/mirFetalEdges.csv", quote = F, row.names = F)
write.csv(mirFetalEdges, "output/mirAdultEdges.csv", quote = F, row.names = F)
write.csv(mirAdultNodes, "output/mirAdultNodes.csv", quote = F, row.names = F)
write.csv(mirFetalNodes, "output/mirFetalNodes.csv", quote = F, row.names = F)
```

### TrasnmiR

``` r
TransmiR <- read.delim("data_for_import/transmiR.tsv", sep = "\t", header = F)
humanMirAnnotationDE$strippedName <- gsub(x = humanMirAnnotationDE$TargetName, pattern = "-5p", replacement = "")
humanMirAnnotationDE$strippedName <- gsub(x = humanMirAnnotationDE$strippedName, pattern = "-3p", replacement = "")
humanMirAnnotationDE$strippedNameLower <- tolower(humanMirAnnotationDE$strippedName)


TransmiRDE <- TransmiR[TransmiR$V2 %in% tolower(humanMirAnnotationDE$strippedName),]
leftout <- humanMirAnnotationDE[tolower(humanMirAnnotationDE$strippedName) %not in% TransmiR$V2,]


### hsa-miR-3656 is not in the db.
### hsa-miR-378d can be in two different loci

TransmiRDE <- TransmiRDE[TransmiRDE$V1 %in% de_intersect$external_gene_name,]
unique(TransmiRDE$V1)
```

    ##  [1] "AHR"    "ARNT"   "BCL11A" "BCL6"   "BMPR1A" "BRD3"   "CCND2"  "CTBP2" 
    ##  [9] "E2F6"   "EHMT2"  "ELF1"   "ELL2"   "EOMES"  "EPAS1"  "ETV1"   "ETV4"  
    ## [17] "EZH2"   "FGFR1"  "FOXK1"  "FOXM1"  "FOXO1"  "FOXP1"  "GLI2"   "HDAC2" 
    ## [25] "HES1"   "HEY1"   "IL1B"   "INTS11" "KDM1A"  "KDM4A"  "KDM4C"  "MAX"   
    ## [33] "MAZ"    "MBD3"   "MEF2A"  "MYC"    "NFE2L2" "NIPBL"  "NME2"   "NR2F1" 
    ## [41] "NR2F2"  "NR3C1"  "OTX2"   "PBX1"   "PCGF2"  "PRDM5"  "RUNX1"  "RUNX2" 
    ## [49] "SMAD1"  "SNAI2"  "SOX2"   "SOX4"   "SOX9"   "STAT1"  "STAT3"  "TCF3"  
    ## [57] "TEAD1"  "TEAD2"  "TEAD4"  "TFAP2C" "THAP11" "TP53"   "TRIM24" "TRIM28"
    ## [65] "YAP1"   "ZEB1"

``` r
TransmirDEcollapsed <- TransmiRDE[,1:2]
TransmirDEcollapsed <- TransmirDEcollapsed[!duplicated(TransmirDEcollapsed),]
TransmirDEcollapsed <- merge(TransmirDEcollapsed, humanMirAnnotationDE, by.x = "V2", by.y = "strippedNameLower")
TransmirDEcollapsed <- merge(TransmirDEcollapsed, mirFetalAdultAll, by.x = "mirName", by.y = "mirName")
TransmirDEcollapsed <- merge(TransmirDEcollapsed, de_intersect, by.x = "V1", by.y = "external_gene_name")
TransmirDEcollapsed <- TransmirDEcollapsed[,c("mirName", "V1", "logFC", "log2FoldChange")]
TransmirDEcollapsed <- TransmirDEcollapsed[order(TransmirDEcollapsed$logFC, decreasing = F),]
TransmirDEcollapsed$mirName <- factor(TransmirDEcollapsed$mirName, levels = unique(TransmirDEcollapsed$mirName))
TransmirDEcollapsed$dotColor <- ifelse(TransmirDEcollapsed$log2FoldChange < 0, "Fetal", "Adult")
TransmirLabelColor <- ifelse(TransmirDEcollapsed[!duplicated(TransmirDEcollapsed$mirName),]$logFC < 0,  "#00008B","#C40000")

### Sup Table 5c
supTable5c <- TransmirDEcollapsed[,c(1:3)]
supTable5c <- merge(supTable5c, supTable3a, by.x = 2, by.y = "External_Gene_Name")
supTable5c <- supTable5c[,c(2,1,3,5,7)]
names(supTable5c)[1:3] <- c("mirName", "External_Gene_Name", "mir_Log2FC")
write.xlsx(supTable5c, file = "Extended Data Tables/Extended Data Table 5 - Adult vs Fetal miRNA microarray.xlsx", sheetName = "TransmiR Predictions", row.names = F, append = T)
```

## Plot TF-miRNA pairs

``` r
ggplot(TransmirDEcollapsed, aes(x=mirName, y=log2FoldChange)) + geom_violin(trim=TRUE) +  geom_jitter(aes(color = dotColor), shape=16, position=position_jitter(0.2)) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color = TransmirLabelColor)) + scale_color_manual(values = c("#C40000", "#00008B")) + xlab("miR Name") + ylab("Log2FC of Upstream TF")
```

    ## Warning: Vectorized input to `element_text()` is not officially supported.
    ## ℹ Results may be unexpected or may change in future versions of ggplot2.

    ## Warning: Groups with fewer than two data points have been dropped.
    ## Groups with fewer than two data points have been dropped.
    ## Groups with fewer than two data points have been dropped.
    ## Groups with fewer than two data points have been dropped.
    ## Groups with fewer than two data points have been dropped.

![](10_miRNA_Microarray_DE_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

## Make Transmir Network

``` r
humanMirAnnotationDEup <- humanMirAnnotationDE[humanMirAnnotationDE$Transcript.ID.Array.Design. %in% mirFetalAdultAll[mirFetalAdultAll$logFC > 0,]$Transcript.ID.Array.Design.,]

humanMirAnnotationDEdown <- humanMirAnnotationDE[humanMirAnnotationDE$Transcript.ID.Array.Design. %in% mirFetalAdultAll[mirFetalAdultAll$logFC < 0,]$Transcript.ID.Array.Design.,]

mirAdultRepressors <- TransmirDEcollapsed[TransmirDEcollapsed$mirName %in% humanMirAnnotationDEdown$mirName & TransmirDEcollapsed$V1 %in% adultRepressorNetwork$Source,]
mirAdultRepressors$type <- "Adult Repressor"
mirFetalActivators <- TransmirDEcollapsed[TransmirDEcollapsed$mirName %in% humanMirAnnotationDEdown$mirName  & TransmirDEcollapsed$V1 %in% fetalActivatorNetwork$Source,]
mirFetalActivators$type <- "Fetal Activator"
mirFetalRepressors <- TransmirDEcollapsed[TransmirDEcollapsed$mirName %in% humanMirAnnotationDEup$mirName  & TransmirDEcollapsed$V1 %in% fetalRepressorNetwork$Source,]
mirFetalRepressors$type <- "Fetal Repressor"
mirAdultActivators <- TransmirDEcollapsed[TransmirDEcollapsed$mirName %in% humanMirAnnotationDEup$mirName  & TransmirDEcollapsed$V1 %in% adultActivatorNetwork$Source,]
mirAdultActivators$type <- "Adult Activator"



mirTFedgesFetalMiRs <- rbind(mirAdultRepressors, mirFetalActivators)
#write.table(mirTFedgesFetalMiRs, "output/mirTFedgesFetalMiRs.txt", sep = "\t", quote = F, row.names = F)

mirTFedgesAdultMiRs <- rbind(mirFetalRepressors, mirAdultActivators)
#write.table(mirTFedgesAdultMiRs, "output/mirTFedgesAdultMiRs.txt", sep = "\t", quote = F, row.names = F)
```

``` r
sessionInfo()
```

    ## R version 4.2.3 (2023-03-15)
    ## Platform: aarch64-apple-darwin20 (64-bit)
    ## Running under: macOS Ventura 13.2.1
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] miRNAtap.db_0.99.10     patchwork_1.1.2         ggrepel_0.9.3          
    ##  [4] ggfortify_0.4.16        ggplot2_3.4.2           xlsx_0.6.5             
    ##  [7] trqwe_0.1               plyr_1.8.8              dplyr_1.1.1            
    ## [10] tidyr_1.3.0             miRNAtap_1.32.0         AnnotationDbi_1.60.2   
    ## [13] miRBaseConverter_1.22.0 pd.mirna.3.0_3.12.0     DBI_1.1.3              
    ## [16] oligo_1.62.2            Biobase_2.58.0          oligoClasses_1.60.0    
    ## [19] RSQLite_2.3.1           Biostrings_2.66.0       GenomeInfoDb_1.34.9    
    ## [22] XVector_0.38.0          IRanges_2.32.0          S4Vectors_0.36.2       
    ## [25] BiocGenerics_0.44.0     limma_3.54.2           
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] bitops_1.0-7                matrixStats_0.63.0         
    ##  [3] bit64_4.0.5                 httr_1.4.5                 
    ##  [5] rprojroot_2.0.3             tools_4.2.3                
    ##  [7] utf8_1.2.3                  R6_2.5.1                   
    ##  [9] affyio_1.68.0               colorspace_2.1-0           
    ## [11] withr_2.5.0                 tidyselect_1.2.0           
    ## [13] gridExtra_2.3               bit_4.0.5                  
    ## [15] compiler_4.2.3              preprocessCore_1.60.2      
    ## [17] chron_2.3-60                cli_3.6.1                  
    ## [19] DelayedArray_0.24.0         labeling_0.4.2             
    ## [21] scales_1.2.1                stringr_1.5.0              
    ## [23] digest_0.6.31               rmarkdown_2.21             
    ## [25] pkgconfig_2.0.3             htmltools_0.5.5            
    ## [27] MatrixGenerics_1.10.0       highr_0.10                 
    ## [29] fastmap_1.1.1               rlang_1.1.0                
    ## [31] rstudioapi_0.14             farver_2.1.1               
    ## [33] generics_0.1.3              RCurl_1.98-1.12            
    ## [35] magrittr_2.0.3              GenomeInfoDbData_1.2.9     
    ## [37] Matrix_1.5-4                Rcpp_1.0.10                
    ## [39] munsell_0.5.0               fansi_1.0.4                
    ## [41] proto_1.0.0                 lifecycle_1.0.3            
    ## [43] sqldf_0.4-11                stringi_1.7.12             
    ## [45] yaml_2.3.7                  SummarizedExperiment_1.28.0
    ## [47] zlibbioc_1.44.0             grid_4.2.3                 
    ## [49] affxparser_1.70.0           blob_1.2.4                 
    ## [51] crayon_1.5.2                lattice_0.21-8             
    ## [53] splines_4.2.3               xlsxjars_0.6.1             
    ## [55] KEGGREST_1.38.0             knitr_1.42                 
    ## [57] pillar_1.9.0                GenomicRanges_1.50.2       
    ## [59] codetools_0.2-19            glue_1.6.2                 
    ## [61] evaluate_0.20               BiocManager_1.30.20        
    ## [63] png_0.1-8                   vctrs_0.6.1                
    ## [65] foreach_1.5.2               gtable_0.3.3               
    ## [67] purrr_1.0.1                 gsubfn_0.7                 
    ## [69] cachem_1.0.7                xfun_0.38                  
    ## [71] ff_4.0.9                    tibble_3.2.1               
    ## [73] rJava_1.0-6                 iterators_1.0.14           
    ## [75] memoise_2.0.1
