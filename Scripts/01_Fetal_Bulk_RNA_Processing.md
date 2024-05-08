Processing of Fetal and Adult GPC Bulk RNA-seq
================
John Mariani
3/6/2023

## Load in Libraries

``` r
library(readr)
library(tximport)
library(biomaRt)
library(DESeq2)
source(file = "Scripts/Helper_Functions.R")
```

## Read in RSEM gene output

``` r
temp = list.files(path = "data_for_import/genes", pattern="genes.results")
temp <- temp[c(1:3,15:34)]

names(temp) <- substr(temp,1,nchar(temp)-19)


txi.rsem <- tximport(paste0("data_for_import/genes/",temp), type = "rsem")
```

    ## It looks like you are importing RSEM genes.results files, setting txIn=FALSE

    ## reading in files with read_tsv

    ## 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23

``` r
for(i in 1:3){
  colnames(txi.rsem[[i]]) <- names(temp)
}
```

## Grab Ensembl 95 Gene annotations from biomaRt unless you’ve already done so

``` r
filename="data_for_import/ensemblGeneList.csv"
if(file.exists(filename)){
  ensemblGeneListH <- read.csv(filename)} else{
    marth <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = 'http://jan2019.archive.ensembl.org/', ensemblRedirect = T)
    ensemblGeneListH <- getBM(attributes = c("ensembl_gene_id","external_gene_name", "gene_biotype", "description"), filters = "ensembl_gene_id",values = row.names(txi.rsem$counts), mart = marth)
    write.csv(ensemblGeneListH, filename, row.names = F)
  }
```

## Read in sample information

``` r
sampleTableFull <- read.csv("data_for_import/sampleTableFull.csv")
```

## Preprocessing

``` r
txi.rsem$length[txi.rsem$length == 0] <- 1


#Annotate the abundance dataframe
TPM <- merge(txi.rsem$abundance, ensemblGeneListH,by.x=0,by.y="ensembl_gene_id")
write.table(TPM, "output/TPM.txt", quote = F, row.names = F, sep = "\t")



lowTPMfull <- data.frame(row.names = ensemblGeneListH$ensembl_gene_id)

for(i in unique(sampleTableFull$line)){
  lowTPMfull[,i] <- groupMedian(txi.rsem$abundance, "line", i, sampleTableFull)
}
lowTPMfull$external_gene_name <- ensemblGeneListH$external_gene_name

tpmCutoff <- 1
highTPM<- lowTPMfull[apply(lowTPMfull[,1:ncol(lowTPMfull)-1], 1, max)>tpmCutoff, ]
```

## Save data for other scripts

``` r
saveRDS(txi.rsem, file = "RDS/txi.rsem.rds")
saveRDS(highTPM, file = "RDS/highTPM.rds")
```

## Session Info

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
    ##  [1] DESeq2_1.38.3               SummarizedExperiment_1.28.0
    ##  [3] Biobase_2.58.0              MatrixGenerics_1.10.0      
    ##  [5] matrixStats_0.63.0          GenomicRanges_1.50.2       
    ##  [7] GenomeInfoDb_1.34.9         IRanges_2.32.0             
    ##  [9] S4Vectors_0.36.2            BiocGenerics_0.44.0        
    ## [11] biomaRt_2.54.1              tximport_1.26.1            
    ## [13] readr_2.1.4                
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] httr_1.4.5             vroom_1.6.1            bit64_4.0.5           
    ##  [4] BiocFileCache_2.6.1    blob_1.2.4             GenomeInfoDbData_1.2.9
    ##  [7] yaml_2.3.7             progress_1.2.2         pillar_1.9.0          
    ## [10] RSQLite_2.3.1          lattice_0.21-8         glue_1.6.2            
    ## [13] digest_0.6.31          RColorBrewer_1.1-3     XVector_0.38.0        
    ## [16] colorspace_2.1-0       htmltools_0.5.5        Matrix_1.5-4          
    ## [19] XML_3.99-0.14          pkgconfig_2.0.3        zlibbioc_1.44.0       
    ## [22] xtable_1.8-4           scales_1.2.1           tzdb_0.3.0            
    ## [25] BiocParallel_1.32.6    tibble_3.2.1           annotate_1.76.0       
    ## [28] KEGGREST_1.38.0        generics_0.1.3         ggplot2_3.4.2         
    ## [31] cachem_1.0.7           cli_3.6.1              magrittr_2.0.3        
    ## [34] crayon_1.5.2           memoise_2.0.1          evaluate_0.20         
    ## [37] fansi_1.0.4            xml2_1.3.3             tools_4.2.3           
    ## [40] prettyunits_1.1.1      hms_1.1.3              lifecycle_1.0.3       
    ## [43] stringr_1.5.0          locfit_1.5-9.7         munsell_0.5.0         
    ## [46] DelayedArray_0.24.0    AnnotationDbi_1.60.2   Biostrings_2.66.0     
    ## [49] compiler_4.2.3         rlang_1.1.0            grid_4.2.3            
    ## [52] RCurl_1.98-1.12        rstudioapi_0.14        rappdirs_0.3.3        
    ## [55] bitops_1.0-7           rmarkdown_2.21         gtable_0.3.3          
    ## [58] codetools_0.2-19       DBI_1.1.3              curl_5.0.0            
    ## [61] R6_2.5.1               knitr_1.42             dplyr_1.1.1           
    ## [64] fastmap_1.1.1          bit_4.0.5              utf8_1.2.3            
    ## [67] filelock_1.0.2         rprojroot_2.0.3        stringi_1.7.12        
    ## [70] parallel_4.2.3         Rcpp_1.0.10            vctrs_0.6.1           
    ## [73] geneplotter_1.76.0     png_0.1-8              dbplyr_2.3.2          
    ## [76] tidyselect_1.2.0       xfun_0.38
