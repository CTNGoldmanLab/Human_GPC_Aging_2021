Processing of In Vitro Aging Overexperssion scRNA-seq
================
John Mariani
1/24/2023

``` r
library(Seurat)
library(ggplot2)
library(plyr)
library(patchwork)

`%not in%` <- function(x, table) is.na(match(x, table, nomatch = NA_integer_))


options(future.globals.maxSize = 16000 * 1024^2)
```

# Load in Count Matrices and create seurat objects

``` r
sampleList <- list.files("Matrices/Overexpression/")
sampleList
```

    ## [1] "000113_Sample_E2F6_Rep1" "000114_Sample_GFP_Rep1" 
    ## [3] "000115_Sample_ZNF_Rep1"

``` r
# Read in filtered matrices 
raw <- sapply(sampleList, function(x) {print(x) ; Read10X(paste0("Matrices/Overexpression/",x,"/star_hs/filtered"))})
```

    ## [1] "000113_Sample_E2F6_Rep1"
    ## [1] "000114_Sample_GFP_Rep1"
    ## [1] "000115_Sample_ZNF_Rep1"

``` r
sets <- length(raw)

#  Convert to list of seurat objects and filter for quality
ObjectsH <- sapply(c(1:sets), function(x) CreateSeuratObject(raw[[x]], project = sampleList[x]))
```

## Make Quality Violin Plots

``` r
for (i in 1:sets) {
  ObjectsH[[i]] <- PercentageFeatureSet(ObjectsH[[i]], pattern = "^MT-", col.name = "percent.mt")
}

saveRDS(ObjectsH, "RDS/ObjectsH_Pre.rds")

for (i in 1:sets) {
  ObjectsH[[i]] <- subset(x = ObjectsH[[i]], subset = nFeature_RNA > 500 & percent.mt < 15)
}

VlnPlot(ObjectsH[[3]], c("nFeature_RNA", "nCount_RNA", "percent.mt"))
```

![](Aging_TF_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
saveRDS(ObjectsH, "RDS/ObjectsH_Post.rds")
```

## Merge datasets and compute Cell Cycle Scores

``` r
merged <- merge(ObjectsH[[1]], y = ObjectsH[2:length(ObjectsH)])
```

    ## Warning in CheckDuplicateCellNames(object.list = objects): Some cell names are
    ## duplicated across objects provided. Renaming to enforce unique cell names.

``` r
merged <- NormalizeData(merged)

s.genes = cc.genes$s.genes
s.genes[s.genes == "MLF1IP"] <- "CENPU"
g2m.genes = cc.genes$g2m.genes
g2m.genes[g2m.genes == "FAM64A"] <- "PIMREG"
g2m.genes[g2m.genes == "HN1"] <- "JPT1"

merged <- CellCycleScoring(merged, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
merged$CC.difference <- merged$S.Score - merged$G2M.Score
```

## Integrate Data and Save

``` r
integrationList <- SplitObject(merged, split.by = "orig.ident")


for(i in 1:length(integrationList)){
  integrationList[[i]] <- SCTransform(integrationList[[i]], verbose = T, vars.to.regress = c("percent.mt", "nCount_RNA", "CC.difference"), conserve.memory = F)
}
```

    ## Calculating cell attributes from input UMI matrix: log_umi

    ## Variance stabilizing transformation of count matrix of size 14768 by 4517

    ## Model formula is y ~ log_umi

    ## Get Negative Binomial regression parameters per gene

    ## Using 2000 genes, 4517 cells

    ##   |                                                                              |                                                                      |   0%  |                                                                              |==================                                                    |  25%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================================| 100%

    ## Found 57 outliers - those will be ignored in fitting/regularization step

    ## Second step: Get residuals using fitted parameters for 14768 genes

    ##   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   3%  |                                                                              |=====                                                                 |   7%  |                                                                              |=======                                                               |  10%  |                                                                              |=========                                                             |  13%  |                                                                              |============                                                          |  17%  |                                                                              |==============                                                        |  20%  |                                                                              |================                                                      |  23%  |                                                                              |===================                                                   |  27%  |                                                                              |=====================                                                 |  30%  |                                                                              |=======================                                               |  33%  |                                                                              |==========================                                            |  37%  |                                                                              |============================                                          |  40%  |                                                                              |==============================                                        |  43%  |                                                                              |=================================                                     |  47%  |                                                                              |===================================                                   |  50%  |                                                                              |=====================================                                 |  53%  |                                                                              |========================================                              |  57%  |                                                                              |==========================================                            |  60%  |                                                                              |============================================                          |  63%  |                                                                              |===============================================                       |  67%  |                                                                              |=================================================                     |  70%  |                                                                              |===================================================                   |  73%  |                                                                              |======================================================                |  77%  |                                                                              |========================================================              |  80%  |                                                                              |==========================================================            |  83%  |                                                                              |=============================================================         |  87%  |                                                                              |===============================================================       |  90%  |                                                                              |=================================================================     |  93%  |                                                                              |====================================================================  |  97%  |                                                                              |======================================================================| 100%

    ## Computing corrected count matrix for 14768 genes

    ##   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   3%  |                                                                              |=====                                                                 |   7%  |                                                                              |=======                                                               |  10%  |                                                                              |=========                                                             |  13%  |                                                                              |============                                                          |  17%  |                                                                              |==============                                                        |  20%  |                                                                              |================                                                      |  23%  |                                                                              |===================                                                   |  27%  |                                                                              |=====================                                                 |  30%  |                                                                              |=======================                                               |  33%  |                                                                              |==========================                                            |  37%  |                                                                              |============================                                          |  40%  |                                                                              |==============================                                        |  43%  |                                                                              |=================================                                     |  47%  |                                                                              |===================================                                   |  50%  |                                                                              |=====================================                                 |  53%  |                                                                              |========================================                              |  57%  |                                                                              |==========================================                            |  60%  |                                                                              |============================================                          |  63%  |                                                                              |===============================================                       |  67%  |                                                                              |=================================================                     |  70%  |                                                                              |===================================================                   |  73%  |                                                                              |======================================================                |  77%  |                                                                              |========================================================              |  80%  |                                                                              |==========================================================            |  83%  |                                                                              |=============================================================         |  87%  |                                                                              |===============================================================       |  90%  |                                                                              |=================================================================     |  93%  |                                                                              |====================================================================  |  97%  |                                                                              |======================================================================| 100%

    ## Calculating gene attributes

    ## Wall clock passed: Time difference of 1.405293 mins

    ## Determine variable features

    ## Place corrected count matrix in counts slot

    ## Regressing out percent.mt, nCount_RNA, CC.difference

    ## Centering data matrix

    ## Set default assay to SCT

    ## Calculating cell attributes from input UMI matrix: log_umi

    ## Variance stabilizing transformation of count matrix of size 15088 by 4526

    ## Model formula is y ~ log_umi

    ## Get Negative Binomial regression parameters per gene

    ## Using 2000 genes, 4526 cells

    ##   |                                                                              |                                                                      |   0%  |                                                                              |==================                                                    |  25%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================================| 100%

    ## There are 1 estimated thetas smaller than 1e-07 - will be set to 1e-07

    ## Found 81 outliers - those will be ignored in fitting/regularization step

    ## Second step: Get residuals using fitted parameters for 15088 genes

    ##   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   3%  |                                                                              |=====                                                                 |   6%  |                                                                              |=======                                                               |  10%  |                                                                              |=========                                                             |  13%  |                                                                              |===========                                                           |  16%  |                                                                              |==============                                                        |  19%  |                                                                              |================                                                      |  23%  |                                                                              |==================                                                    |  26%  |                                                                              |====================                                                  |  29%  |                                                                              |=======================                                               |  32%  |                                                                              |=========================                                             |  35%  |                                                                              |===========================                                           |  39%  |                                                                              |=============================                                         |  42%  |                                                                              |================================                                      |  45%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |======================================                                |  55%  |                                                                              |=========================================                             |  58%  |                                                                              |===========================================                           |  61%  |                                                                              |=============================================                         |  65%  |                                                                              |===============================================                       |  68%  |                                                                              |==================================================                    |  71%  |                                                                              |====================================================                  |  74%  |                                                                              |======================================================                |  77%  |                                                                              |========================================================              |  81%  |                                                                              |===========================================================           |  84%  |                                                                              |=============================================================         |  87%  |                                                                              |===============================================================       |  90%  |                                                                              |=================================================================     |  94%  |                                                                              |====================================================================  |  97%  |                                                                              |======================================================================| 100%

    ## Computing corrected count matrix for 15088 genes

    ##   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   3%  |                                                                              |=====                                                                 |   6%  |                                                                              |=======                                                               |  10%  |                                                                              |=========                                                             |  13%  |                                                                              |===========                                                           |  16%  |                                                                              |==============                                                        |  19%  |                                                                              |================                                                      |  23%  |                                                                              |==================                                                    |  26%  |                                                                              |====================                                                  |  29%  |                                                                              |=======================                                               |  32%  |                                                                              |=========================                                             |  35%  |                                                                              |===========================                                           |  39%  |                                                                              |=============================                                         |  42%  |                                                                              |================================                                      |  45%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |======================================                                |  55%  |                                                                              |=========================================                             |  58%  |                                                                              |===========================================                           |  61%  |                                                                              |=============================================                         |  65%  |                                                                              |===============================================                       |  68%  |                                                                              |==================================================                    |  71%  |                                                                              |====================================================                  |  74%  |                                                                              |======================================================                |  77%  |                                                                              |========================================================              |  81%  |                                                                              |===========================================================           |  84%  |                                                                              |=============================================================         |  87%  |                                                                              |===============================================================       |  90%  |                                                                              |=================================================================     |  94%  |                                                                              |====================================================================  |  97%  |                                                                              |======================================================================| 100%

    ## Calculating gene attributes

    ## Wall clock passed: Time difference of 1.384985 mins

    ## Determine variable features

    ## Place corrected count matrix in counts slot

    ## Regressing out percent.mt, nCount_RNA, CC.difference

    ## Centering data matrix

    ## Set default assay to SCT

    ## Calculating cell attributes from input UMI matrix: log_umi

    ## Variance stabilizing transformation of count matrix of size 15365 by 5584

    ## Model formula is y ~ log_umi

    ## Get Negative Binomial regression parameters per gene

    ## Using 2000 genes, 5000 cells

    ##   |                                                                              |                                                                      |   0%  |                                                                              |==================                                                    |  25%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================================| 100%

    ## Found 84 outliers - those will be ignored in fitting/regularization step

    ## Second step: Get residuals using fitted parameters for 15365 genes

    ##   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   3%  |                                                                              |=====                                                                 |   6%  |                                                                              |=======                                                               |  10%  |                                                                              |=========                                                             |  13%  |                                                                              |===========                                                           |  16%  |                                                                              |==============                                                        |  19%  |                                                                              |================                                                      |  23%  |                                                                              |==================                                                    |  26%  |                                                                              |====================                                                  |  29%  |                                                                              |=======================                                               |  32%  |                                                                              |=========================                                             |  35%  |                                                                              |===========================                                           |  39%  |                                                                              |=============================                                         |  42%  |                                                                              |================================                                      |  45%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |======================================                                |  55%  |                                                                              |=========================================                             |  58%  |                                                                              |===========================================                           |  61%  |                                                                              |=============================================                         |  65%  |                                                                              |===============================================                       |  68%  |                                                                              |==================================================                    |  71%  |                                                                              |====================================================                  |  74%  |                                                                              |======================================================                |  77%  |                                                                              |========================================================              |  81%  |                                                                              |===========================================================           |  84%  |                                                                              |=============================================================         |  87%  |                                                                              |===============================================================       |  90%  |                                                                              |=================================================================     |  94%  |                                                                              |====================================================================  |  97%  |                                                                              |======================================================================| 100%

    ## Computing corrected count matrix for 15365 genes

    ##   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   3%  |                                                                              |=====                                                                 |   6%  |                                                                              |=======                                                               |  10%  |                                                                              |=========                                                             |  13%  |                                                                              |===========                                                           |  16%  |                                                                              |==============                                                        |  19%  |                                                                              |================                                                      |  23%  |                                                                              |==================                                                    |  26%  |                                                                              |====================                                                  |  29%  |                                                                              |=======================                                               |  32%  |                                                                              |=========================                                             |  35%  |                                                                              |===========================                                           |  39%  |                                                                              |=============================                                         |  42%  |                                                                              |================================                                      |  45%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |======================================                                |  55%  |                                                                              |=========================================                             |  58%  |                                                                              |===========================================                           |  61%  |                                                                              |=============================================                         |  65%  |                                                                              |===============================================                       |  68%  |                                                                              |==================================================                    |  71%  |                                                                              |====================================================                  |  74%  |                                                                              |======================================================                |  77%  |                                                                              |========================================================              |  81%  |                                                                              |===========================================================           |  84%  |                                                                              |=============================================================         |  87%  |                                                                              |===============================================================       |  90%  |                                                                              |=================================================================     |  94%  |                                                                              |====================================================================  |  97%  |                                                                              |======================================================================| 100%

    ## Calculating gene attributes

    ## Wall clock passed: Time difference of 1.600314 mins

    ## Determine variable features

    ## Place corrected count matrix in counts slot

    ## Regressing out percent.mt, nCount_RNA, CC.difference

    ## Centering data matrix

    ## Set default assay to SCT

``` r
features <- SelectIntegrationFeatures(object.list = integrationList, nfeatures = 3000)
features <- features[features %not in% c("eGFP", "ZNF274", "E2F6")]

integrationList <- PrepSCTIntegration(object.list = integrationList, anchor.features = features, verbose = T)
integrated <- FindIntegrationAnchors(object.list = integrationList, normalization.method = "SCT", anchor.features = features, verbose = T)
```

    ## Finding all pairwise anchors

    ## Running CCA

    ## Merging objects

    ## Finding neighborhoods

    ## Finding anchors

    ##  Found 11022 anchors

    ## Filtering anchors

    ##  Retained 8910 anchors

    ## Running CCA

    ## Merging objects

    ## Finding neighborhoods

    ## Finding anchors

    ##  Found 12245 anchors

    ## Filtering anchors

    ##  Retained 9598 anchors

    ## Running CCA

    ## Merging objects

    ## Finding neighborhoods

    ## Finding anchors

    ##  Found 12847 anchors

    ## Filtering anchors

    ##  Retained 10513 anchors

``` r
integrated <- IntegrateData(anchorset = integrated, normalization.method = "SCT", verbose = T)
```

    ## Merging dataset 2 into 3

    ## Extracting anchors for merged samples

    ## Finding integration vectors

    ## Finding integration vector weights

    ## Integrating data

    ## Merging dataset 1 into 3 2

    ## Extracting anchors for merged samples

    ## Finding integration vectors

    ## Finding integration vector weights

    ## Integrating data

``` r
DefaultAssay(integrated) <- "integrated"
integrated <- RunPCA(integrated, verbose = FALSE)
integrated <- RunUMAP(integrated, dims = 1:50, verbose = FALSE, umap.method = "umap-learn")
integrated <- FindNeighbors(integrated, dims = 1:50, verbose = FALSE)
integrated <- FindClusters(integrated, verbose = FALSE, resolution = .25)

DimPlot(integrated)
```

![](Aging_TF_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
integrated$ogClusters <- Idents(integrated)

saveRDS(integrated, "RDS/integrated_OE.rds")
```

## Session Info

``` r
sessionInfo()
```

    ## R version 4.1.1 (2021-08-10)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Red Hat Enterprise Linux Server 7.9 (Maipo)
    ## 
    ## Matrix products: default
    ## BLAS:   /gpfs/fs1/sfw2/r/4.1.1/b1/lib64/R/lib/libRblas.so
    ## LAPACK: /gpfs/fs1/sfw2/r/4.1.1/b1/lib64/R/lib/libRlapack.so
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] patchwork_1.1.1    plyr_1.8.6         ggplot2_3.3.5      sp_1.4-6          
    ## [5] SeuratObject_4.1.0 Seurat_4.1.1      
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] Rtsne_0.16            colorspace_2.0-2      deldir_1.0-6         
    ##   [4] ellipsis_0.3.2        ggridges_0.5.3        rprojroot_2.0.2      
    ##   [7] rgdal_1.5-28          rstudioapi_0.13       spatstat.data_3.0-0  
    ##  [10] farver_2.1.0          leiden_0.3.9          listenv_0.8.0        
    ##  [13] ggrepel_0.9.1         fansi_1.0.2           codetools_0.2-18     
    ##  [16] splines_4.1.1         knitr_1.37            polyclip_1.10-0      
    ##  [19] jsonlite_1.7.3        ica_1.0-2             cluster_2.1.2        
    ##  [22] png_0.1-7             rgeos_0.5-9           uwot_0.1.11          
    ##  [25] shiny_1.7.1           sctransform_0.3.3     spatstat.sparse_3.0-0
    ##  [28] compiler_4.1.1        httr_1.4.2            assertthat_0.2.1     
    ##  [31] Matrix_1.4-0          fastmap_1.1.0         lazyeval_0.2.2       
    ##  [34] cli_3.3.0             later_1.3.0           htmltools_0.5.2      
    ##  [37] tools_4.1.1           igraph_1.3.0          gtable_0.3.0         
    ##  [40] glue_1.6.2            RANN_2.6.1            reshape2_1.4.4       
    ##  [43] dplyr_1.0.8           rappdirs_0.3.3        Rcpp_1.0.8           
    ##  [46] scattermore_0.8       vctrs_0.3.8           nlme_3.1-155         
    ##  [49] progressr_0.10.0      lmtest_0.9-39         spatstat.random_3.0-1
    ##  [52] xfun_0.29             stringr_1.4.0         globals_0.14.0       
    ##  [55] mime_0.12             miniUI_0.1.1.1        lifecycle_1.0.1      
    ##  [58] irlba_2.3.5           goftest_1.2-3         future_1.24.0        
    ##  [61] MASS_7.3-55           zoo_1.8-11            scales_1.1.1         
    ##  [64] spatstat.core_2.4-2   promises_1.2.0.1      spatstat.utils_3.0-1 
    ##  [67] parallel_4.1.1        RColorBrewer_1.1-2    yaml_2.3.4           
    ##  [70] reticulate_1.24       pbapply_1.5-0         gridExtra_2.3        
    ##  [73] rpart_4.1.16          stringi_1.7.6         highr_0.9            
    ##  [76] rlang_1.0.1           pkgconfig_2.0.3       matrixStats_0.61.0   
    ##  [79] evaluate_0.15         lattice_0.20-45       ROCR_1.0-11          
    ##  [82] purrr_0.3.4           tensor_1.5            labeling_0.4.2       
    ##  [85] htmlwidgets_1.5.4     cowplot_1.1.1         tidyselect_1.1.1     
    ##  [88] here_1.0.1            parallelly_1.30.0     RcppAnnoy_0.0.19     
    ##  [91] magrittr_2.0.2        R6_2.5.1              generics_0.1.2       
    ##  [94] DBI_1.1.2             withr_2.4.3           mgcv_1.8-38          
    ##  [97] pillar_1.7.0          fitdistrplus_1.1-6    survival_3.2-13      
    ## [100] abind_1.4-5           tibble_3.1.6          future.apply_1.8.1   
    ## [103] crayon_1.5.1          KernSmooth_2.23-20    utf8_1.2.2           
    ## [106] spatstat.geom_3.0-3   plotly_4.10.0         rmarkdown_2.11       
    ## [109] grid_4.1.1            data.table_1.14.2     digest_0.6.29        
    ## [112] xtable_1.8-4          tidyr_1.2.0           httpuv_1.6.5         
    ## [115] munsell_0.5.0         viridisLite_0.4.0
