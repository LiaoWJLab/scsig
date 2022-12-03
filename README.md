
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mysc

### 1.Introduction

mysc is an R package to perform single cell RNAseq data analysis. Main
advantages: 1. standard single cell pipeline: standard\_sc() 2. doublet
detection: doublet\_detect() 3. easy dimension plot: dong\_dimplot() 4.
find marker genes of clusters: dong\_find\_markers() 5. training cell
annotation model for reference: training\_scRef() 6. identifying tumor
cells:

### 2.Installation

The package is not yet on CRAN. You can install from Github:

### 3.Usage

Main documentation is on the `mysc` function in the package:

``` r
library('mysc')
library('IOBR')
library('mydb')
```

#### 3.1 Example

``` r
#' for-single-cell-data-derived-from-10X-genomics
data_path<-"H:/03-NSCLC/13-NSCLC-scRNA-Immunotherapy/3-Data/1-sc-data/1-eset/TIA2/filtered_feature_bc_matrix"

help(standard_sc)
#> starting httpd help server ... done
sce<-standard_sc(eset               = NULL, 
                 file_type          = "10X", 
                 data_path          = data_path,  
                 project            = "TIA2", 
                 nPCs               = 30,
                 res                = 1,  
                 verbose            = FALSE, 
                 index              = 1,   #default resolution
                 plot               = TRUE,  
                 minFeature         = 20, #filtering option
                 minCount           = 10,
                 percent.mt         = 20,  # set to TRUE if user want to identify marker genes of each clusters
                 findmarkers        = FALSE, 
                 already_normalized = FALSE)
#> Warning: Feature names cannot have underscores ('_'), replacing with dashes
#> ('-')
#> An object of class Seurat 
#> 60623 features across 2059 samples within 1 assay 
#> Active assay: RNA (60623 features, 0 variable features)
#>  >>>  Data will be deposite in H:/18-Github/mysc/1-TIA2/
#> [1] ">>> Step-1: Quality control"
#> Mitochondrial genes have been removed
#> Ribosomal genes have been removed
#>                    orig.ident nCount_RNA nFeature_RNA percent.mt percent.ribo
#> AAACCTGAGATCACGG-1       TIA2        920          441          0            0
#> AAACCTGGTAGGAGTC-1       TIA2       1164          817          0            0
#> AAACCTGGTGCAGGTA-1       TIA2        599          106          0            0
#> AAACGGGAGGGCTCTC-1       TIA2       3593          453          0            0
#> AAACGGGAGTTCCACA-1       TIA2        889          499          0            0
#> AAACGGGCATCTCGCT-1       TIA2        893          505          0            0
#> Warning in cor(x = data[, 1], y = data[, 2]): 标准差为零
#> Warning: CombinePlots is being deprecated. Plots should now be combined using
#> the patchwork system.
#> Warning in SingleExIPlot(type = type, data = data[, x, drop = FALSE], idents =
#> idents, : All cells have the same value of percent.ribo.
#> Warning in SingleExIPlot(type = type, data = data[, x, drop = FALSE], idents =
#> idents, : All cells have the same value of percent.mt.
#> Warning in SingleExIPlot(type = type, data = data[, x, drop = FALSE], idents =
#> idents, : All cells have the same value of percent.ribo.
#> >>> For cell subset: Default parameters are : minFeature = 2000, minCount = 1000, percent.mt = 20
#>  >>>  After filtering cells with low features, low count and high expression of mitochondrial genes
#> An object of class Seurat 
#> 60623 features across 2059 samples within 1 assay 
#> Active assay: RNA (60623 features, 0 variable features)
#> [1] ">>> Step-2: Data normalization and dimension reduction"
#> [1] ">>> Step-3: Find clusters"
#> Computing nearest neighbor graph
#> Computing SNN
#> Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
#> 
#> Number of nodes: 2059
#> Number of edges: 76073
#> 
#> Running Louvain algorithm...
#> Maximum modularity in 10 random starts: 0.9454
#> Number of communities: 7
#> Elapsed time: 0 seconds
#>   0   1   2   3   4   5   6 
#> 626 362 318 229 187 177 160 
#>   0   1   2   3   4   5   6   7   8   9  10  11  12  13 
#> 385 237 229 176 170 169 146 143 114 104  59  47  45  35 
#>   0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15 
#> 344 243 176 167 146 146 142 121 114 107 104  62  59  47  46  35
#> Warning: The following arguments are not used: do.fast
#> Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
#> To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
#> This message will be shown once per session
```

<img src="man/figuresunnamed-chunk-4-1.png" width="100%" />

    #> [1] "'#7FC97F', '#BEAED4', '#FDC086', '#386CB0', '#F0027F', '#BF5B17', '#1B9E77', '#7570B3', '#66A61E', '#E6AB02', '#E64B35FF', '#4DBBD5FF', '#F39B7FFF', '#91D1C2FF', '#DC0000FF', '#374E55FF', '#DF8F44FF', '#B24745FF', '#3B4992FF', '#631879FF', '#A20056FF', '#0073C2FF', '#EFC000FF', '#CD534CFF', '#7AA6DCFF', '#003C67FF', '#8F7700FF', '#3B3B3BFF', '#A73030FF', '#377EB8', '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#66C2A5', '#FC8D62', '#8DA0CB', '#E78AC3', '#A6D854', '#FFD92F', '#E5C494', '#8DD3C7', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#FCCDE5', '#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C', '#FB9A99', '#E31A1C', '#FDBF6F', '#FF7F00', '#224444', '#e68a00', '#33adff', '#a6a6a6', '#439373', '#67001F', '#1B9E77', '#FB9A99', '#92C5DE', '#33A02C', '#92C5DE', '#4393C3', '#B2DF8A', '#CAB2D6', '#56B4E9', '#BC3C29FF'"

<img src="man/figuresunnamed-chunk-4-2.png" width="100%" /><img src="man/figuresunnamed-chunk-4-3.png" width="100%" /><img src="man/figuresunnamed-chunk-4-4.png" width="100%" />

``` r
#' for-single-cell-data-with-raw-count-as-input--------------------------------
data("sc_tnbc")
sce<-standard_sc(eset               = sc_tnbc, 
                 file_type          = "10X", 
                 data_path          = NULL,
                 project            = "TNBC", 
                 nPCs               = 30,
                 res                = 1,  
                 verbose            = FALSE, 
                 index              = 1,
                 plot               = TRUE,  
                 minFeature         = 20, 
                 minCount           = 10,
                 percent.mt         = 20,  
                 findmarkers        = FALSE, # set to TRUE if user want to identify marker genes of each clusters
                 already_normalized = FALSE)
#> An object of class Seurat 
#> 33694 features across 1097 samples within 1 assay 
#> Active assay: RNA (33694 features, 0 variable features)
#>  >>>  Data will be deposite in H:/18-Github/mysc/1-TNBC/
#> [1] ">>> Step-1: Quality control"
#> Mitochondrial genes have been removed
#> Ribosomal genes have been removed
#>                  orig.ident nCount_RNA nFeature_RNA percent.mt percent.ribo
#> AAACCTGCACCTTGTC       TNBC      13976         3590          0            0
#> AAACGGGAGTCCTCCT       TNBC       8732         2629          0            0
#> AAACGGGTCCAGAGGA       TNBC      17138         4166          0            0
#> AAAGATGCAGTTTACG       TNBC       9519         1703          0            0
#> AAAGCAACAGGAATGC       TNBC      10285         2778          0            0
#> AAAGCAATCGGAATCT       TNBC       8436         2822          0            0
#> Warning in cor(x = data[, 1], y = data[, 2]): 标准差为零
#> Warning: CombinePlots is being deprecated. Plots should now be combined using
#> the patchwork system.
#> Warning in SingleExIPlot(type = type, data = data[, x, drop = FALSE], idents =
#> idents, : All cells have the same value of percent.ribo.
#> Warning in SingleExIPlot(type = type, data = data[, x, drop = FALSE], idents =
#> idents, : All cells have the same value of percent.mt.
#> Warning in SingleExIPlot(type = type, data = data[, x, drop = FALSE], idents =
#> idents, : All cells have the same value of percent.ribo.
#> >>> For cell subset: Default parameters are : minFeature = 2000, minCount = 1000, percent.mt = 20
#>  >>>  After filtering cells with low features, low count and high expression of mitochondrial genes
#> An object of class Seurat 
#> 33694 features across 1097 samples within 1 assay 
#> Active assay: RNA (33694 features, 0 variable features)
#> [1] ">>> Step-2: Data normalization and dimension reduction"
#> [1] ">>> Step-3: Find clusters"
#> Computing nearest neighbor graph
#> Computing SNN
#> Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
#> 
#> Number of nodes: 1097
#> Number of edges: 37333
#> 
#> Running Louvain algorithm...
#> Maximum modularity in 10 random starts: 0.9356
#> Number of communities: 6
#> Elapsed time: 0 seconds
#>   0   1   2   3   4   5 
#> 437 303 226  56  54  21 
#>   0   1   2   3   4   5   6   7   8 
#> 292 259 180 145  71  54  46  29  21 
#>   0   1   2   3   4   5   6   7   8   9  10 
#> 231 164 141 139 116 110  64  54  29  28  21
#> Warning: The following arguments are not used: do.fast
```

<img src="man/figuresunnamed-chunk-4-5.png" width="100%" />

    #> [1] "'#7FC97F', '#BEAED4', '#FDC086', '#386CB0', '#F0027F', '#BF5B17', '#1B9E77', '#7570B3', '#66A61E', '#E6AB02', '#E64B35FF', '#4DBBD5FF', '#F39B7FFF', '#91D1C2FF', '#DC0000FF', '#374E55FF', '#DF8F44FF', '#B24745FF', '#3B4992FF', '#631879FF', '#A20056FF', '#0073C2FF', '#EFC000FF', '#CD534CFF', '#7AA6DCFF', '#003C67FF', '#8F7700FF', '#3B3B3BFF', '#A73030FF', '#377EB8', '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#66C2A5', '#FC8D62', '#8DA0CB', '#E78AC3', '#A6D854', '#FFD92F', '#E5C494', '#8DD3C7', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#FCCDE5', '#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C', '#FB9A99', '#E31A1C', '#FDBF6F', '#FF7F00', '#224444', '#e68a00', '#33adff', '#a6a6a6', '#439373', '#67001F', '#1B9E77', '#FB9A99', '#92C5DE', '#33A02C', '#92C5DE', '#4393C3', '#B2DF8A', '#CAB2D6', '#56B4E9', '#BC3C29FF'"

<img src="man/figuresunnamed-chunk-4-6.png" width="100%" /><img src="man/figuresunnamed-chunk-4-7.png" width="100%" /><img src="man/figuresunnamed-chunk-4-8.png" width="100%" />

``` r
head(sce@meta.data)
#>                  orig.ident nCount_RNA nFeature_RNA percent.mt percent.ribo
#> AAACCTGCACCTTGTC       TNBC      13976         3590          0            0
#> AAACGGGAGTCCTCCT       TNBC       8732         2629          0            0
#> AAACGGGTCCAGAGGA       TNBC      17138         4166          0            0
#> AAAGATGCAGTTTACG       TNBC       9519         1703          0            0
#> AAAGCAACAGGAATGC       TNBC      10285         2778          0            0
#> AAAGCAATCGGAATCT       TNBC       8436         2822          0            0
#>                  RNA_snn_res.0.2 seurat_clusters RNA_snn_res.0.8 RNA_snn_res.1
#> AAACCTGCACCTTGTC               3               9               3             9
#> AAACGGGAGTCCTCCT               2               5               2             5
#> AAACGGGTCCAGAGGA               0               0               0             0
#> AAAGATGCAGTTTACG               1               1               1             1
#> AAAGCAACAGGAATGC               0               0               0             0
#> AAAGCAATCGGAATCT               1               3               4             3
```

\#\#find marker genes of each clusters

``` r
dong_find_markers(sce        = sce, 
                  group      = "RNA_snn_res.1", 
                  top        = 10,
                  verbose    = FALSE,
                  fig.type   = "pdf",
                  pt.size    = 0.5, 
                  cols       = "normal", 
                  seed       = 54321, 
                  show_col   = F, 
                  max_cols   = 36,
                  show_genes = 10,
                  hwidth     = 13, 
                  hheight    = 8.2,
                  show_plot  = T, 
                  path       = NULL )
#> >>> Idents of Seurat object is: RNA_snn_res.1
#> 
#>   0   1   2   3   4   5   6   7   8   9  10 
#> 231 164 141 139 116 110  64  54  29  28  21
#> Warning in DoHeatmap(sce, features = top10$gene, angle = 60, size = 3.5, : The
#> following features were omitted as they were not found in the scale.data slot
#> for the RNA assay: SLC29A1, CEBPD
#> Scale for 'fill' is already present. Adding another scale for 'fill', which
#> will replace the existing scale.
#> >>> Processing cluster: 0
#>               p_val avg_log2FC pct.1 pct.2    p_val_adj
#> IFRD1  4.536179e-83   1.591789 0.978 0.721 1.528420e-78
#> SQSTM1 8.382420e-81   1.329052 1.000 0.940 2.824373e-76
#> CEBPD  8.756984e-81   1.752058 0.991 0.915 2.950578e-76
#> SBSPON 4.691434e-76   1.840927 0.896 0.328 1.580732e-71
#> CLU    7.005458e-74   2.136658 0.987 0.777 2.360419e-69
#> MDFI   1.507978e-71   1.511643 0.987 0.658 5.080981e-67
#> >>> Processing cluster: 1
#>               p_val avg_log2FC pct.1 pct.2    p_val_adj
#> PCSK1N 1.107386e-81  2.9368162 1.000 0.708 3.731228e-77
#> CALML5 4.530106e-78  3.1175341 1.000 0.887 1.526374e-73
#> CFB    7.360834e-75  1.2726108 0.835 0.192 2.480159e-70
#> SEZ6L2 1.834793e-72  0.7422809 0.671 0.100 6.182153e-68
#> S100A6 4.723056e-72  2.1794245 1.000 0.997 1.591387e-67
#> CRLF1  4.244460e-65  0.5915279 0.573 0.071 1.430128e-60
#> >>> Processing cluster: 10
#>                  p_val avg_log2FC pct.1 pct.2     p_val_adj
#> ATP6V0D2 7.748189e-119  1.2349968 0.810 0.010 2.610675e-114
#> RUFY4    4.422410e-117  0.9240462 0.714 0.007 1.490087e-112
#> SIGLEC15 2.321180e-106  2.1729594 0.952 0.022 7.820984e-102
#> STRIP2   4.588293e-103  0.7451051 0.571 0.004  1.545979e-98
#> DPP4     3.257959e-101  2.0573517 0.810 0.015  1.097737e-96
#> SLC9B2    2.241746e-89  2.5858379 1.000 0.035  7.553339e-85
#> >>> Processing cluster: 2
#>                 p_val avg_log2FC pct.1 pct.2    p_val_adj
#> AZGP1    7.216641e-55  2.1730921 0.993 0.682 2.431575e-50
#> TMEM106C 1.991554e-54  1.5285037 1.000 0.704 6.710342e-50
#> MARCKSL1 3.420926e-51  1.3516377 1.000 0.868 1.152647e-46
#> TCEB2    7.210931e-51  0.9801733 1.000 0.987 2.429651e-46
#> CLPSL1   2.557722e-50  1.4948290 0.858 0.323 8.617987e-46
#> IDH1     4.199267e-48  1.5013079 0.993 0.666 1.414901e-43
#> >>> Processing cluster: 3
#>               p_val avg_log2FC pct.1 pct.2    p_val_adj
#> CLDN6  2.867708e-91  1.7841916 0.878 0.151 9.662457e-87
#> SMIM22 2.010032e-82  1.4009411 0.885 0.178 6.772602e-78
#> PCAT19 1.307326e-77  0.6016770 0.683 0.076 4.404904e-73
#> MAGIX  2.590624e-71  1.4688294 0.971 0.370 8.728848e-67
#> KIF1A  2.735859e-67  0.3446962 0.590 0.062 9.218203e-63
#> RAB6B  5.700943e-66  0.6366021 0.655 0.096 1.920876e-61
#> >>> Processing cluster: 4
#>                p_val avg_log2FC pct.1 pct.2    p_val_adj
#> ALOX5AP 3.101915e-96   2.557231 0.888 0.130 1.045159e-91
#> CYBB    1.150280e-87   1.734240 0.922 0.138 3.875753e-83
#> LCP1    1.994185e-87   1.606348 0.905 0.151 6.719207e-83
#> TGFBI   3.171163e-85   2.540471 0.948 0.202 1.068492e-80
#> SPI1    7.847299e-82   1.922512 0.940 0.178 2.644069e-77
#> CD53    2.424690e-80   1.398227 0.914 0.151 8.169751e-76
#> >>> Processing cluster: 5
#>                p_val avg_log2FC pct.1 pct.2     p_val_adj
#> MS4A4A 5.903920e-111   1.819220 0.855 0.076 1.989267e-106
#> MSR1   2.024414e-107   2.035037 0.909 0.109 6.821060e-103
#> MS4A7  3.836149e-106   2.641641 0.982 0.151 1.292552e-101
#> CD163  1.221016e-102   1.738753 0.891 0.095  4.114092e-98
#> MS4A6A 4.380433e-102   2.756329 0.945 0.136  1.475943e-97
#> FCGR3A  8.570964e-99   1.974184 0.836 0.087  2.887900e-94
#> >>> Processing cluster: 6
#>                 p_val avg_log2FC pct.1 pct.2    p_val_adj
#> SLC26A7  7.480742e-43  1.0404211 0.672 0.098 2.520561e-38
#> TOX      1.066455e-41  0.9500629 0.953 0.241 3.593313e-37
#> CNGA1    6.929179e-35  0.4857899 0.641 0.108 2.334717e-30
#> PLA2G4A  7.754693e-34  1.1210204 0.906 0.275 2.612866e-29
#> FAM198A  2.598654e-32  0.4799967 0.812 0.187 8.755904e-28
#> C10orf99 4.434466e-32  0.9888397 0.250 0.012 1.494149e-27
#> >>> Processing cluster: 7
#>                p_val avg_log2FC pct.1 pct.2     p_val_adj
#> CFH    1.549368e-162   1.824041 0.852 0.012 5.220439e-158
#> NTM    1.774126e-160   1.510111 0.741 0.004 5.977740e-156
#> PDGFRB 6.010230e-159   2.649824 0.907 0.018 2.025087e-154
#> FBN1   1.553514e-150   2.063519 0.852 0.016 5.234409e-146
#> ECM2   5.919022e-149   1.210730 0.722 0.006 1.994355e-144
#> CHN1   7.750026e-148   1.813266 0.852 0.017 2.611294e-143
#> >>> Processing cluster: 8
#>              p_val avg_log2FC pct.1 pct.2    p_val_adj
#> RIBC2 3.870992e-71  0.3593246 0.655 0.021 1.304292e-66
#> MCM10 3.446968e-64  0.8768738 0.897 0.057 1.161421e-59
#> DTL   9.232544e-59  0.7859208 0.828 0.052 3.110813e-54
#> EXO1  3.112930e-58  0.5442979 0.759 0.041 1.048871e-53
#> CDC45 8.124437e-54  0.6153790 0.793 0.054 2.737448e-49
#> CDC6  2.542558e-48  0.4410393 0.724 0.047 8.566895e-44
#> >>> Processing cluster: 9
#>               p_val avg_log2FC pct.1 pct.2    p_val_adj
#> NEK2   1.803522e-61  1.9820619 1.000 0.084 6.076787e-57
#> PBK    2.064247e-60  1.5950544 1.000 0.086 6.955274e-56
#> CENPA  4.710395e-58  1.1745086 0.964 0.077 1.587120e-53
#> CDC25C 3.077076e-57  0.4231027 0.607 0.023 1.036790e-52
#> NMU    1.576128e-54  1.6298499 0.786 0.051 5.310606e-50
#> KIF20A 6.266232e-54  0.9329978 0.786 0.051 2.111344e-49
```

<img src="man/figuresunnamed-chunk-6-1.png" width="100%" />

## Doublet detection by DoubletFinder or scds packages

``` r
#find-doublet------------------------------------------------------------------
sce_df<-doublet_detect(sce                  = sce,
                         already_normalized = T,
                         eset               = NULL,
                         file_type          = "10X",
                         data_path          = NULL,
                         project            = "TNBC",
                         index              = 1,
                         check_data         = TRUE,
                         filter_data        = FALSE,
                         minFeature         = 2000,
                         minCount           = 1000,
                         percent.mt         =20,
                         method             = "doubletfinder")
#> An object of class Seurat 
#> 33694 features across 1097 samples within 1 assay 
#> Active assay: RNA (33694 features, 2000 variable features)
#>  3 dimensional reductions calculated: pca, tsne, umap
#> [1] ">>> Quality control"
#> [1] "Creating artificial doublets for pN = 5%"
#> [1] "Creating Seurat object..."
#> [1] "Normalizing Seurat object..."
#> [1] "Finding variable genes..."
#> [1] "Scaling data..."
#> [1] "Running PCA..."
#> [1] "Calculating PC distance matrix..."
#> [1] "Defining neighborhoods..."
#> [1] "Computing pANN across all pK..."
#> [1] "pK = 0.01..."
#> [1] "pK = 0.02..."
#> [1] "pK = 0.03..."
#> [1] "pK = 0.04..."
#> [1] "pK = 0.05..."
#> [1] "pK = 0.06..."
#> [1] "pK = 0.07..."
#> [1] "pK = 0.08..."
#> [1] "pK = 0.09..."
#> [1] "pK = 0.1..."
#> [1] "pK = 0.11..."
#> [1] "pK = 0.12..."
#> [1] "pK = 0.13..."
#> [1] "pK = 0.14..."
#> [1] "pK = 0.15..."
#> [1] "pK = 0.16..."
#> [1] "pK = 0.17..."
#> [1] "pK = 0.18..."
#> [1] "pK = 0.19..."
#> [1] "pK = 0.2..."
#> [1] "pK = 0.21..."
#> [1] "pK = 0.22..."
#> [1] "pK = 0.23..."
#> [1] "pK = 0.24..."
#> [1] "pK = 0.25..."
#> [1] "pK = 0.26..."
#> [1] "pK = 0.27..."
#> [1] "pK = 0.28..."
#> [1] "pK = 0.29..."
#> [1] "pK = 0.3..."
#> [1] "Creating artificial doublets for pN = 10%"
#> [1] "Creating Seurat object..."
#> [1] "Normalizing Seurat object..."
#> [1] "Finding variable genes..."
#> [1] "Scaling data..."
#> [1] "Running PCA..."
#> [1] "Calculating PC distance matrix..."
#> [1] "Defining neighborhoods..."
#> [1] "Computing pANN across all pK..."
#> [1] "pK = 0.01..."
#> [1] "pK = 0.02..."
#> [1] "pK = 0.03..."
#> [1] "pK = 0.04..."
#> [1] "pK = 0.05..."
#> [1] "pK = 0.06..."
#> [1] "pK = 0.07..."
#> [1] "pK = 0.08..."
#> [1] "pK = 0.09..."
#> [1] "pK = 0.1..."
#> [1] "pK = 0.11..."
#> [1] "pK = 0.12..."
#> [1] "pK = 0.13..."
#> [1] "pK = 0.14..."
#> [1] "pK = 0.15..."
#> [1] "pK = 0.16..."
#> [1] "pK = 0.17..."
#> [1] "pK = 0.18..."
#> [1] "pK = 0.19..."
#> [1] "pK = 0.2..."
#> [1] "pK = 0.21..."
#> [1] "pK = 0.22..."
#> [1] "pK = 0.23..."
#> [1] "pK = 0.24..."
#> [1] "pK = 0.25..."
#> [1] "pK = 0.26..."
#> [1] "pK = 0.27..."
#> [1] "pK = 0.28..."
#> [1] "pK = 0.29..."
#> [1] "pK = 0.3..."
#> [1] "Creating artificial doublets for pN = 15%"
#> [1] "Creating Seurat object..."
#> [1] "Normalizing Seurat object..."
#> [1] "Finding variable genes..."
#> [1] "Scaling data..."
#> [1] "Running PCA..."
#> [1] "Calculating PC distance matrix..."
#> [1] "Defining neighborhoods..."
#> [1] "Computing pANN across all pK..."
#> [1] "pK = 0.01..."
#> [1] "pK = 0.02..."
#> [1] "pK = 0.03..."
#> [1] "pK = 0.04..."
#> [1] "pK = 0.05..."
#> [1] "pK = 0.06..."
#> [1] "pK = 0.07..."
#> [1] "pK = 0.08..."
#> [1] "pK = 0.09..."
#> [1] "pK = 0.1..."
#> [1] "pK = 0.11..."
#> [1] "pK = 0.12..."
#> [1] "pK = 0.13..."
#> [1] "pK = 0.14..."
#> [1] "pK = 0.15..."
#> [1] "pK = 0.16..."
#> [1] "pK = 0.17..."
#> [1] "pK = 0.18..."
#> [1] "pK = 0.19..."
#> [1] "pK = 0.2..."
#> [1] "pK = 0.21..."
#> [1] "pK = 0.22..."
#> [1] "pK = 0.23..."
#> [1] "pK = 0.24..."
#> [1] "pK = 0.25..."
#> [1] "pK = 0.26..."
#> [1] "pK = 0.27..."
#> [1] "pK = 0.28..."
#> [1] "pK = 0.29..."
#> [1] "pK = 0.3..."
#> [1] "Creating artificial doublets for pN = 20%"
#> [1] "Creating Seurat object..."
#> [1] "Normalizing Seurat object..."
#> [1] "Finding variable genes..."
#> [1] "Scaling data..."
#> [1] "Running PCA..."
#> [1] "Calculating PC distance matrix..."
#> [1] "Defining neighborhoods..."
#> [1] "Computing pANN across all pK..."
#> [1] "pK = 0.01..."
#> [1] "pK = 0.02..."
#> [1] "pK = 0.03..."
#> [1] "pK = 0.04..."
#> [1] "pK = 0.05..."
#> [1] "pK = 0.06..."
#> [1] "pK = 0.07..."
#> [1] "pK = 0.08..."
#> [1] "pK = 0.09..."
#> [1] "pK = 0.1..."
#> [1] "pK = 0.11..."
#> [1] "pK = 0.12..."
#> [1] "pK = 0.13..."
#> [1] "pK = 0.14..."
#> [1] "pK = 0.15..."
#> [1] "pK = 0.16..."
#> [1] "pK = 0.17..."
#> [1] "pK = 0.18..."
#> [1] "pK = 0.19..."
#> [1] "pK = 0.2..."
#> [1] "pK = 0.21..."
#> [1] "pK = 0.22..."
#> [1] "pK = 0.23..."
#> [1] "pK = 0.24..."
#> [1] "pK = 0.25..."
#> [1] "pK = 0.26..."
#> [1] "pK = 0.27..."
#> [1] "pK = 0.28..."
#> [1] "pK = 0.29..."
#> [1] "pK = 0.3..."
#> [1] "Creating artificial doublets for pN = 25%"
#> [1] "Creating Seurat object..."
#> [1] "Normalizing Seurat object..."
#> [1] "Finding variable genes..."
#> [1] "Scaling data..."
#> [1] "Running PCA..."
#> [1] "Calculating PC distance matrix..."
#> [1] "Defining neighborhoods..."
#> [1] "Computing pANN across all pK..."
#> [1] "pK = 0.01..."
#> [1] "pK = 0.02..."
#> [1] "pK = 0.03..."
#> [1] "pK = 0.04..."
#> [1] "pK = 0.05..."
#> [1] "pK = 0.06..."
#> [1] "pK = 0.07..."
#> [1] "pK = 0.08..."
#> [1] "pK = 0.09..."
#> [1] "pK = 0.1..."
#> [1] "pK = 0.11..."
#> [1] "pK = 0.12..."
#> [1] "pK = 0.13..."
#> [1] "pK = 0.14..."
#> [1] "pK = 0.15..."
#> [1] "pK = 0.16..."
#> [1] "pK = 0.17..."
#> [1] "pK = 0.18..."
#> [1] "pK = 0.19..."
#> [1] "pK = 0.2..."
#> [1] "pK = 0.21..."
#> [1] "pK = 0.22..."
#> [1] "pK = 0.23..."
#> [1] "pK = 0.24..."
#> [1] "pK = 0.25..."
#> [1] "pK = 0.26..."
#> [1] "pK = 0.27..."
#> [1] "pK = 0.28..."
#> [1] "pK = 0.29..."
#> [1] "pK = 0.3..."
#> [1] "Creating artificial doublets for pN = 30%"
#> [1] "Creating Seurat object..."
#> [1] "Normalizing Seurat object..."
#> [1] "Finding variable genes..."
#> [1] "Scaling data..."
#> [1] "Running PCA..."
#> [1] "Calculating PC distance matrix..."
#> [1] "Defining neighborhoods..."
#> [1] "Computing pANN across all pK..."
#> [1] "pK = 0.01..."
#> [1] "pK = 0.02..."
#> [1] "pK = 0.03..."
#> [1] "pK = 0.04..."
#> [1] "pK = 0.05..."
#> [1] "pK = 0.06..."
#> [1] "pK = 0.07..."
#> [1] "pK = 0.08..."
#> [1] "pK = 0.09..."
#> [1] "pK = 0.1..."
#> [1] "pK = 0.11..."
#> [1] "pK = 0.12..."
#> [1] "pK = 0.13..."
#> [1] "pK = 0.14..."
#> [1] "pK = 0.15..."
#> [1] "pK = 0.16..."
#> [1] "pK = 0.17..."
#> [1] "pK = 0.18..."
#> [1] "pK = 0.19..."
#> [1] "pK = 0.2..."
#> [1] "pK = 0.21..."
#> [1] "pK = 0.22..."
#> [1] "pK = 0.23..."
#> [1] "pK = 0.24..."
#> [1] "pK = 0.25..."
#> [1] "pK = 0.26..."
#> [1] "pK = 0.27..."
#> [1] "pK = 0.28..."
#> [1] "pK = 0.29..."
#> [1] "pK = 0.3..."
```

<img src="man/figuresunnamed-chunk-7-1.png" width="100%" />

    #> NULL
    #> [1] "Creating 366 artificial doublets..."
    #> [1] "Creating Seurat object..."
    #> [1] "Normalizing Seurat object..."
    #> [1] "Finding variable genes..."
    #> [1] "Scaling data..."
    #> [1] "Running PCA..."
    #> [1] "Calculating PC distance matrix..."
    #> [1] "Computing pANN..."
    #> [1] "Classifying doublets.."
    #> Doublet Singlet 
    #>      61    1036

``` r
head(sce_df@meta.data)
#>                  orig.ident nCount_RNA nFeature_RNA percent.mt percent.ribo
#> AAACCTGCACCTTGTC       TNBC      13976         3590          0            0
#> AAACGGGAGTCCTCCT       TNBC       8732         2629          0            0
#> AAACGGGTCCAGAGGA       TNBC      17138         4166          0            0
#> AAAGATGCAGTTTACG       TNBC       9519         1703          0            0
#> AAAGCAACAGGAATGC       TNBC      10285         2778          0            0
#> AAAGCAATCGGAATCT       TNBC       8436         2822          0            0
#>                  RNA_snn_res.0.2 seurat_clusters RNA_snn_res.0.8 RNA_snn_res.1
#> AAACCTGCACCTTGTC               3               9               3             9
#> AAACGGGAGTCCTCCT               2               5               2             5
#> AAACGGGTCCAGAGGA               0               0               0             0
#> AAAGATGCAGTTTACG               1               1               1             1
#> AAAGCAACAGGAATGC               0               0               0             0
#> AAAGCAATCGGAATCT               1               3               4             3
#>                  pANN_0.25_0.01_82 DF.classifications_0.25_0.01_82
#> AAACCTGCACCTTGTC        0.33333333                         Doublet
#> AAACGGGAGTCCTCCT        0.06666667                         Singlet
#> AAACGGGTCCAGAGGA        0.13333333                         Singlet
#> AAAGATGCAGTTTACG        0.13333333                         Singlet
#> AAAGCAACAGGAATGC        0.06666667                         Singlet
#> AAAGCAATCGGAATCT        0.20000000                         Singlet
#>                  DoubletFinder_final
#> AAACCTGCACCTTGTC             Singlet
#> AAACGGGAGTCCTCCT             Singlet
#> AAACGGGTCCAGAGGA             Singlet
#> AAAGATGCAGTTTACG             Singlet
#> AAAGCAACAGGAATGC             Singlet
#> AAAGCAATCGGAATCT             Singlet
```

``` r
p1<-dong_dimplot(sce = sce_df,
                 reduction = "umap",
                 groups = "DoubletFinder_final",
                 split.by = NULL,
                 label=T,
                 label.size = 5,
                 pt.size = 0.5,
                 cols = "normal", seed = 54321, show_col = F, max_cols = 36,
                 width = 8, height = 8, w_index = 7, w_add = 2,
                 max_category = 14,
                 show_plot = T,
                 path = "./",
                 index = paste0("1-DoubletFinder"),
                 legend.position = "right",
                 legend.direction = "vertical",
                 legend.size = 10)
#> Warning in dir.create(file_store): '.'已存在
#> >>> Processing group:: DoubletFinder_final
#> >>>> theme: classic. Options: light, bw, classic and classic2
#> >>>>>> Figure name is:: 1-DoubletFinder-1-umap-DoubletFinder_final.pdf
```

<img src="man/figuresunnamed-chunk-9-1.png" width="100%" />

``` r
# cnv-estimation---------------------------------------------------------------
copykat_res<-copykat_plus(eset        = sc_tnbc, 
                          file_type   = "10X", 
                          data_path   = NULL,
                          index       = 1, 
                          project     = "TNBC", 
                          id_type     = "symbol",
                          check_data  = TRUE,
                          filter_data = FALSE,
                          minFeature  = 2000,
                          minCount    = 1000,
                          percent.mt  = 20)
#> An object of class Seurat 
#> 33694 features across 1097 samples within 1 assay 
#> Active assay: RNA (33694 features, 0 variable features)
#>  >>>  Data will be deposite in H:/18-Github/mysc/1-TNBC/
#> [1] ">>> Quality control"
#> Mitochondrial genes have been removed
#> Ribosomal genes have been removed
#>                  orig.ident nCount_RNA nFeature_RNA percent.mt percent.ribo
#> AAACCTGCACCTTGTC       TNBC      13976         3590          0            0
#> AAACGGGAGTCCTCCT       TNBC       8732         2629          0            0
#> AAACGGGTCCAGAGGA       TNBC      17138         4166          0            0
#> AAAGATGCAGTTTACG       TNBC       9519         1703          0            0
#> AAAGCAACAGGAATGC       TNBC      10285         2778          0            0
#> AAAGCAATCGGAATCT       TNBC       8436         2822          0            0
#> Warning in cor(x = data[, 1], y = data[, 2]): 标准差为零
#> Warning: CombinePlots is being deprecated. Plots should now be combined using
#> the patchwork system.
#> Warning in SingleExIPlot(type = type, data = data[, x, drop = FALSE], idents =
#> idents, : All cells have the same value of percent.ribo.
#> Warning in SingleExIPlot(type = type, data = data[, x, drop = FALSE], idents =
#> idents, : All cells have the same value of percent.mt.
```

<img src="man/figuresunnamed-chunk-10-1.png" width="100%" /><img src="man/figuresunnamed-chunk-10-2.png" width="100%" />

    #>               AAACCTGCACCTTGTC AAACGGGAGTCCTCCT AAACGGGTCCAGAGGA
    #> RP11-34P13.3                 0                0                0
    #> FAM138A                      0                0                0
    #> OR4F5                        0                0                0
    #> RP11-34P13.7                 0                0                0
    #> RP11-34P13.8                 0                0                0
    #> RP11-34P13.14                0                0                0
    #> RP11-34P13.9                 0                0                0
    #> FO538757.3                   0                0                0
    #> FO538757.2                   1                0                0
    #> AP006222.2                   0                0                0
    #> [1] "running copykat v1.0.4"
    #> [1] "step1: read and filter data ..."
    #> [1] "33694 genes, 1097 cells in raw data"
    #> [1] "12230 genes past LOW.DR filtering"
    #> [1] "step 2: annotations gene coordinates ..."
    #> [1] "start annotation ..."
    #> [1] "step 3: smoothing data with dlm ..."
    #> [1] "step 4: measuring baselines ..."
    #> number of iterations= 347 
    #> number of iterations= 730 
    #> number of iterations= 246 
    #> number of iterations= 770 
    #> number of iterations= 790 
    #> number of iterations= 500 
    #> [1] "step 5: segmentation..."
    #> [1] "step 6: convert to genomic bins..."
    #> [1] "step 7: adjust baseline ..."
    #> [1] "step 8: final prediction ..."
    #> [1] "step 9: saving results..."
    #> [1] "step 10: ploting heatmap ..."

<img src="man/figuresunnamed-chunk-10-3.png" width="100%" />

    #> Time difference of 14.28615 mins

``` r
sce<-readRDS("H:/01-1-GC-Project/11-GC-CAF/2-CAFs-project/4-scRNAseq/0-data/0-Integrated-seurat-sce-of-18-samples.rds")

head(sce@meta.data)
#>                      orig.ident nCount_RNA nFeature_RNA percent.mt percent.ribo
#> AAACCTGAGACTGTAA-1_1    5846_t1       2795          939          0            0
#> AAACCTGCACATGGGA-1_1    5846_t1       1918          730          0            0
#> AAACCTGCACTCTGTC-1_1    5846_t1        655          369          0            0
#> AAACCTGGTAAATGAC-1_1    5846_t1        670          368          0            0
#> AAACCTGGTAGTGAAT-1_1    5846_t1       5936         1600          0            0
#> AAACCTGTCCGTACAA-1_1    5846_t1       1066          557          0            0
#>                           id id_patient tissue_type Age Sex grade2
#> AAACCTGAGACTGTAA-1_1 5846_t1      p5846       tumor  72   M     MP
#> AAACCTGCACATGGGA-1_1 5846_t1      p5846       tumor  72   M     MP
#> AAACCTGCACTCTGTC-1_1 5846_t1      p5846       tumor  72   M     MP
#> AAACCTGGTAAATGAC-1_1 5846_t1      p5846       tumor  72   M     MP
#> AAACCTGGTAGTGAAT-1_1 5846_t1      p5846       tumor  72   M     MP
#> AAACCTGTCCGTACAA-1_1 5846_t1      p5846       tumor  72   M     MP
#>                                                    Grade MSI_subtype
#> AAACCTGAGACTGTAA-1_1 Moderately to poorly differentiated         MSS
#> AAACCTGCACATGGGA-1_1 Moderately to poorly differentiated         MSS
#> AAACCTGCACTCTGTC-1_1 Moderately to poorly differentiated         MSS
#> AAACCTGGTAAATGAC-1_1 Moderately to poorly differentiated         MSS
#> AAACCTGGTAGTGAAT-1_1 Moderately to poorly differentiated         MSS
#> AAACCTGTCCGTACAA-1_1 Moderately to poorly differentiated         MSS
#>                       NT_response Neoadjuvant_treatment
#> AAACCTGAGACTGTAA-1_1 Na<U+00EF>ve                  None
#> AAACCTGCACATGGGA-1_1 Na<U+00EF>ve                  None
#> AAACCTGCACTCTGTC-1_1 Na<U+00EF>ve                  None
#> AAACCTGGTAAATGAC-1_1 Na<U+00EF>ve                  None
#> AAACCTGGTAGTGAAT-1_1 Na<U+00EF>ve                  None
#> AAACCTGTCCGTACAA-1_1 Na<U+00EF>ve                  None
#>                      Pathologic_T_AJCC_7th_edition
#> AAACCTGAGACTGTAA-1_1                           pT3
#> AAACCTGCACATGGGA-1_1                           pT3
#> AAACCTGCACTCTGTC-1_1                           pT3
#> AAACCTGGTAAATGAC-1_1                           pT3
#> AAACCTGGTAGTGAAT-1_1                           pT3
#> AAACCTGTCCGTACAA-1_1                           pT3
#>                      Pathologic_N_AJCC_7th_edition
#> AAACCTGAGACTGTAA-1_1                           pN1
#> AAACCTGCACATGGGA-1_1                           pN1
#> AAACCTGCACTCTGTC-1_1                           pN1
#> AAACCTGGTAAATGAC-1_1                           pN1
#> AAACCTGGTAGTGAAT-1_1                           pN1
#> AAACCTGTCCGTACAA-1_1                           pN1
#>                      Reported_recurrence_or_death_in_12_months_follow_up
#> AAACCTGAGACTGTAA-1_1                                                  No
#> AAACCTGCACATGGGA-1_1                                                  No
#> AAACCTGCACTCTGTC-1_1                                                  No
#> AAACCTGGTAAATGAC-1_1                                                  No
#> AAACCTGGTAGTGAAT-1_1                                                  No
#> AAACCTGTCCGTACAA-1_1                                                  No
#>                      Anatomical_region Other_features_of_interest
#> AAACCTGAGACTGTAA-1_1       GEJ+stomach               Not_reported
#> AAACCTGCACATGGGA-1_1       GEJ+stomach               Not_reported
#> AAACCTGCACTCTGTC-1_1       GEJ+stomach               Not_reported
#> AAACCTGGTAAATGAC-1_1       GEJ+stomach               Not_reported
#> AAACCTGGTAGTGAAT-1_1       GEJ+stomach               Not_reported
#> AAACCTGTCCGTACAA-1_1       GEJ+stomach               Not_reported
#>                           Diagnosis Surgical_procedure
#> AAACCTGAGACTGTAA-1_1 Adenocarcinoma  Total gastrectomy
#> AAACCTGCACATGGGA-1_1 Adenocarcinoma  Total gastrectomy
#> AAACCTGCACTCTGTC-1_1 Adenocarcinoma  Total gastrectomy
#> AAACCTGGTAAATGAC-1_1 Adenocarcinoma  Total gastrectomy
#> AAACCTGGTAGTGAAT-1_1 Adenocarcinoma  Total gastrectomy
#> AAACCTGTCCGTACAA-1_1 Adenocarcinoma  Total gastrectomy
#>                                                                      Histological_type
#> AAACCTGAGACTGTAA-1_1 Invasive adenocarcarcinoma with mucinous and signet ring features
#> AAACCTGCACATGGGA-1_1 Invasive adenocarcarcinoma with mucinous and signet ring features
#> AAACCTGCACTCTGTC-1_1 Invasive adenocarcarcinoma with mucinous and signet ring features
#> AAACCTGGTAAATGAC-1_1 Invasive adenocarcarcinoma with mucinous and signet ring features
#> AAACCTGGTAGTGAAT-1_1 Invasive adenocarcarcinoma with mucinous and signet ring features
#> AAACCTGTCCGTACAA-1_1 Invasive adenocarcarcinoma with mucinous and signet ring features
#>                            cell_barcode condition final_celltype
#> AAACCTGAGACTGTAA-1_1 AAACCTGAGACTGTAA-2     tumor              B
#> AAACCTGCACATGGGA-1_1 AAACCTGCACATGGGA-2     tumor           mast
#> AAACCTGCACTCTGTC-1_1 AAACCTGCACTCTGTC-2     tumor     epithelial
#> AAACCTGGTAAATGAC-1_1 AAACCTGGTAAATGAC-2     tumor     epithelial
#> AAACCTGGTAGTGAAT-1_1 AAACCTGGTAGTGAAT-2     tumor     macrophage
#> AAACCTGTCCGTACAA-1_1 AAACCTGTCCGTACAA-2     tumor     epithelial
#>                      cluster_celltype index integrated_snn_res.1
#> AAACCTGAGACTGTAA-1_1              B_8  -1_1                   11
#> AAACCTGCACATGGGA-1_1          31_mast  -1_1                   21
#> AAACCTGCACTCTGTC-1_1     epithelial_5  -1_1                   17
#> AAACCTGGTAAATGAC-1_1     epithelial_8  -1_1                    7
#> AAACCTGGTAGTGAAT-1_1    macrophage_c5  -1_1                    8
#> AAACCTGTCCGTACAA-1_1     epithelial_3  -1_1                    4
#>                      seurat_clusters
#> AAACCTGAGACTGTAA-1_1              11
#> AAACCTGCACATGGGA-1_1              21
#> AAACCTGCACTCTGTC-1_1              17
#> AAACCTGGTAAATGAC-1_1               7
#> AAACCTGGTAGTGAAT-1_1               8
#> AAACCTGTCCGTACAA-1_1               4
sce_pred<-training_sc_anno(sce          = sce,
                           group        = "final_celltype", 
                           remove_cell  = NULL,
                           subset_cell  = TRUE,
                           propotion    =  0.2,
                           dims         = 45,
                           model        = "mda",
                           cores        = 1,
                           verbose      = FALSE, 
                           fig.type     = "pdf",
                           pt.size      = 1, 
                           cols         = "normal",
                           seed         = 54321,
                           show_col     = F, 
                           max_cols     = 36,
                           show_plot    = T,
                           path         = NULL )
#> >>> Idents of Seurat object is: final_celltype
#> 
#>           B         CD4         CD8          DC endothelial  epithelial 
#>         900        1241        4560          55        1447        8217 
#> fibroblasts  macrophage        mast          NK   pericytes      plasma 
#>        2430        1940         329        1251         544        2076 
#>        Treg 
#>         951
#> >>> Count of each cell type:
#>           Var1 Freq
#> 1            B  900
#> 2          CD4 1241
#> 3          CD8 4560
#> 4           DC   55
#> 5  endothelial 1447
#> 6   epithelial 8217
#> 7  fibroblasts 2430
#> 8   macrophage 1940
#> 9         mast  329
#> 10          NK 1251
#> 11   pericytes  544
#> 12      plasma 2076
#> 13        Treg  951
#> >>> Cell types included in the study:
#>  [1] "B"           "CD4"         "CD8"         "DC"          "endothelial"
#>  [6] "epithelial"  "fibroblasts" "macrophage"  "mast"        "NK"         
#> [11] "pericytes"   "plasma"      "Treg"       
#> 
#>           B         CD4         CD8          DC endothelial  epithelial 
#>         180         248         912          11         289        1643 
#> fibroblasts  macrophage        mast          NK   pericytes      plasma 
#>         486         388          66         250         109         415 
#>        Treg 
#>         190 
#>           B         CD4         CD8          DC endothelial  epithelial 
#>         180         248         912          11         289        1643 
#> fibroblasts  macrophage        mast          NK   pericytes      plasma 
#>         486         388          66         250         109         415 
#>        Treg 
#>         190
#> >>> Subseting cells in each subset randomly: 20%...
#> >>> Final training data:
#> 
#>           B         CD4         CD8          DC endothelial  epithelial 
#>         180         248         912          11         289        1643 
#> fibroblasts  macrophage        mast          NK   pericytes      plasma 
#>         486         388          66         250         109         415 
#>        Treg 
#>         190
#> Centering and scaling data matrix
#> 10:12:37 UMAP embedding parameters a = 0.9922 b = 1.112
#> Found more than one class "dist" in cache; using the first, from namespace 'BiocGenerics'
#> Also defined by 'spam'
#> 10:12:37 Read 5187 rows and found 45 numeric columns
#> 10:12:37 Using Annoy for neighbor search, n_neighbors = 30
#> Found more than one class "dist" in cache; using the first, from namespace 'BiocGenerics'
#> Also defined by 'spam'
#> 10:12:37 Building Annoy index with metric = cosine, n_trees = 50
#> 0%   10   20   30   40   50   60   70   80   90   100%
#> [----|----|----|----|----|----|----|----|----|----|
#> **************************************************|
#> 10:12:37 Writing NN index file to temp file C:\Users\70890\AppData\Local\Temp\RtmpYJGbf1\file49149411ffc
#> 10:12:37 Searching Annoy index using 1 thread, search_k = 3000
#> 10:12:39 Annoy recall = 100%
#> 10:12:40 Commencing smooth kNN distance calibration using 1 thread
#> 10:12:42 Initializing from normalized Laplacian + noise
#> 10:12:42 Commencing optimization for 500 epochs, with 217390 positive edges
#> 10:12:56 Optimization finished
#> >>> Processing group:: final_celltype
#> >>>> theme: classic. Options: light, bw, classic and classic2
#> >>>>>> Figure name is:: 1-Reference-data--1-umap-final_celltype.pdf
```

<img src="man/figuresunnamed-chunk-11-1.png" width="100%" />

    #> o  Extracting feature space for each cell type...
    #> DONE!
    #> o  Training models for each cell type...
    #> 载入需要的程辑包：lattice
    #> 
    #> 载入程辑包：'caret'
    #> 
    #> The following object is masked from 'package:survival':
    #> 
    #>     cluster
    #> 
    #> The following object is masked from 'package:purrr':
    #> 
    #>     lift

<img src="man/figuresunnamed-chunk-11-2.png" width="100%" />

    #> DONE!

<img src="man/figuresunnamed-chunk-11-3.png" width="100%" />

``` r
sce <-scPred:: scPredict(sce, sce_pred,
                        recompute_alignment = FALSE,
                        threshold = 0.65)
#> o  Matching reference with new dataset...
#>   - 2000 features present in reference loadings
#>   - 1319 features shared between reference and new dataset
#>   - 65.95% of features in the reference are present in new dataset
#> o  Aligning new data to reference...
#> Harmony 1/20
#> Harmony 2/20
#> Harmony 3/20
#> Harmony 4/20
#> Harmony 5/20
#> Harmony 6/20
#> Harmony 7/20
#> Harmony 8/20
#> Harmony 9/20
#> Harmony converged after 9 iterations
#> o  Classifying cells...
#> DONE!
table(sce$scpred_prediction)
#> 
#>           B         CD4         CD8          DC endothelial  epithelial 
#>        1055         685        5375          21        1314        8066 
#> fibroblasts  macrophage        mast          NK   pericytes      plasma 
#>        2374        1726         327        1115         500        2029 
#>        Treg  unassigned 
#>         850         504
p1<-dong_dimplot(sce = sce,
                 reduction = "umap",
                 groups = "scpred_prediction",
                 split.by = NULL,
                 label=T,
                 label.size = 5,
                 pt.size = 0.5,
                 cols = "normal", seed = 54321, show_col = F, max_cols = 36,
                 width = 8, height = 8, w_index = 7, w_add = 2,
                 max_category = 14,
                 show_plot = F,
                 path = "./",
                 index = paste0("1-predict-data"),
                 legend.position = "right",
                 legend.direction = "vertical",
                 legend.size = 10)
#> Warning in dir.create(file_store): '.'已存在
#> >>> Processing group:: scpred_prediction
#> >>>> theme: classic. Options: light, bw, classic and classic2
#> >>>>>> Figure name is:: 1-predict-data-1-umap-scpred_prediction.pdf
p2<-dong_dimplot(sce = sce,
                 reduction = "umap",
                 groups = "final_celltype",
                 split.by = NULL,
                 label=T,
                 label.size = 5,
                 pt.size = 0.5,
                 cols = "normal", seed = 54321, show_col = F, max_cols = 36,
                 width = 8, height = 8, w_index = 7, w_add = 2,
                 max_category = 14,
                 show_plot = F,
                 path = "./",
                 index = paste0("2-predict-data"),
                 legend.position = "right",
                 legend.direction = "vertical",
                 legend.size = 10)
#> Warning in dir.create(file_store): '.'已存在
#> >>> Processing group:: final_celltype
#> >>>> theme: classic. Options: light, bw, classic and classic2
#> >>>>>> Figure name is:: 2-predict-data-1-umap-final_celltype.pdf
```

``` r
p1+p2
```

<img src="man/figuresunnamed-chunk-13-1.png" width="100%" />

### Citation

If you use mysc in published research, please cite:

### Contact

E-mail any questions to <dongqiangzeng0808@gmail.com>
