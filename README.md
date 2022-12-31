
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scsig

### 1.Introduction

scsig is an R package to perform single cell RNAseq data analysis. Main
advantages: 1. standard single cell pipeline: standard_sc() 2. doublet
detection: doublet_detect() 3. easy dimension plot: dong_dimplot() 4.
find marker genes of clusters: dong_find_markers() 5. training cell
annotation model for reference: training_scRef() 6. identifying tumor
cells: copykat_plus

### 2.Installation

The package is not yet on CRAN. You can install from Github:

### 3.Usage

Main documentation is on the `scsig` function in the package:

``` r
library('scsig')
library('IOBR')
library('mydb')
library(Seurat)
```

#### 3.1 Example

``` r
#' for-single-cell-data-derived-from-10X-genomics
data_path<-"E:/03-NSCLC/13-NSCLC-scRNA-Immunotherapy/3-Data/1-sc-data/1-eset/TIA2/filtered_feature_bc_matrix"

help(standard_sc)
#> 打开httpd帮助服务器… 好了
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
                 cutoff_percent_mt  = 25,
                 minCount           = 10,
                 findmarkers        = FALSE, 
                 already_normalized = FALSE,
                 save_path          = "./test")
#> Warning: Feature names cannot have underscores ('_'), replacing with dashes
#> ('-')
#> An object of class Seurat 
#> 60623 features across 2059 samples within 1 assay 
#> Active assay: RNA (60623 features, 0 variable features)
#> >>>>=== Palette option for random: 1: palette1; 2: palette2; 3: palette3;  4: palette4
#> [1] "'#5f75ae', '#64a841', '#e5486e', '#de8e06', '#b5aa0f', '#7ba39d', '#b15928', '#6a3d9a', '#cab2d6', '#374E55FF', '#00A1D5FF', '#6A6599FF', '#80796BFF', '#e31a1c', '#fb9a99', '#1f78b4', '#a6cee3', '#008280FF', '#3C5488FF', '#8F7700FF', '#666666', '#A20056FF', '#fdbf6f', '#E78AC3', '#b2df8a', '#386CB0', '#CD534CFF', '#008B45FF', '#7AA6DCFF', '#00A087FF', '#A73030FF', '#631879FF', '#003C67FF'"
#>  >>>-----  Data will be deposite in E:/18-Github/scsig/./test/
#> [1] ">>>----- Step-1: Quality control"
#> [1] "Mitochondrial genes: "
#>  [1] "MT-TF"   "MT-RNR1" "MT-TV"   "MT-RNR2" "MT-TL1"  "MT-ND1"  "MT-TI"  
#>  [8] "MT-TQ"   "MT-TM"   "MT-ND2"  "MT-TW"   "MT-TA"   "MT-TN"   "MT-TC"  
#> [15] "MT-TY"   "MT-CO1"  "MT-TS1"  "MT-TD"   "MT-CO2"  "MT-TK"   "MT-ATP8"
#> [22] "MT-ATP6" "MT-CO3"  "MT-TG"   "MT-ND3"  "MT-TR"   "MT-ND4L" "MT-ND4" 
#> [29] "MT-TH"   "MT-TS2"  "MT-TL2"  "MT-ND5"  "MT-ND6"  "MT-TE"   "MT-CYB" 
#> [36] "MT-TT"   "MT-TP"
#> Ribosomal genes have been removed
#> species must be one of `human` or `mouse`.
#> Warning: CombinePlots is being deprecated. Plots should now be combined using
#> the patchwork system.
#> >>>------ For cell subset: Default parameters are : minFeature = 2000, minCount = 1000, cutoff_percent_mt = 25
#> >>>------ After filtering cells with low features, low count and high expression of mitochondrial genes
#> An object of class Seurat 
#> 60623 features across 1649 samples within 1 assay 
#> Active assay: RNA (60623 features, 0 variable features)
#> [1] ">>>------ Step-2: Data normalization and dimension reduction"
#> >>>------ Cell cycle scoring...
```

<img src="man/figuresunnamed-chunk-4-1.png" width="100%" />

    #>                    orig.ident nCount_RNA nFeature_RNA percent.mt percent.rp
    #> AAACCTGAGATCACGG-1       TIA2        920          441   3.913043  40.326087
    #> AAACCTGGTAGGAGTC-1       TIA2       1164          817  10.395189   1.030928
    #> AAACCTGGTGCAGGTA-1       TIA2        599          106  14.357262   2.504174
    #> AAACGGGAGGGCTCTC-1       TIA2       3593          453   2.894517   7.180629
    #> AAACGGGAGTTCCACA-1       TIA2        889          499  21.034871  13.160855
    #> AAACGGGCATCTCGCT-1       TIA2        893          505  21.836506  16.909295
    #>                         S.Score    G2M.Score Phase
    #> AAACCTGAGATCACGG-1 -0.037260496 -0.033503464    G1
    #> AAACCTGGTAGGAGTC-1 -0.068185303 -0.083362112    G1
    #> AAACCTGGTGCAGGTA-1 -0.008520917 -0.001408458    G1
    #> AAACGGGAGGGCTCTC-1 -0.018726516 -0.018105992    G1
    #> AAACGGGAGTTCCACA-1 -0.031509755  0.011596915   G2M
    #> AAACGGGCATCTCGCT-1  0.032450378  0.056456534   G2M
    #> ------------------------------------------------------
    #> [1] ">>>------ Step-3: Find clusters"
    #> Computing nearest neighbor graph
    #> Computing SNN
    #> Resolution is 0.4
    #> >>>-- Findclusters when resolution = RNA_snn_res.0.4 ...
    #>   0   1   2   3   4   5   6   7 
    #> 497 274 238 175 147 134 114  70
    #> Resolution is 0.8
    #> >>>-- Findclusters when resolution = RNA_snn_res.0.8 ...
    #>   0   1   2   3   4   5   6   7   8   9  10 
    #> 296 208 200 175 147 144 133 131 114  71  30
    #> Resolution is 1.2
    #> >>>-- Findclusters when resolution = RNA_snn_res.1.2 ...
    #>   0   1   2   3   4   5   6   7   8   9  10  11  12  13 
    #> 285 210 174 145 144 135 114 105  71  69  63  62  42  30
    #> Resolution is 1.6
    #> >>>-- Findclusters when resolution = RNA_snn_res.1.6 ...
    #>   0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15 
    #> 230 146 144 144 133 120 117 114 105  71  69  63  62  59  42  30
    #> Resolution is 2
    #> >>>-- Findclusters when resolution = RNA_snn_res.2 ...
    #>   0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16 
    #> 227 146 144 121 120 117 114 105 104  71  69  62  62  58  57  42  30
    #> Scale for 'edge_colour' is already present. Adding another scale for
    #> 'edge_colour', which will replace the existing scale.
    #> >>>-- Findclusters when resolution = RNA_snn_res.1 ...
    #>   0   1   2   3   4   5   6   7   8   9  10  11  12 
    #> 288 208 175 147 147 144 133 114  71  69  62  61  30
    #> >>>-- These results will be used to draw dimension plot....
    #> ------------------------------------------------------
    #> Warning: The following arguments are not used: do.fast
    #> Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
    #> To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
    #> This message will be shown once per session

<img src="man/figuresunnamed-chunk-4-2.png" width="100%" /><img src="man/figuresunnamed-chunk-4-3.png" width="100%" /><img src="man/figuresunnamed-chunk-4-4.png" width="100%" /><img src="man/figuresunnamed-chunk-4-5.png" width="100%" />

    #> ------------------------------------------------------
    #> >>> Processing method:: pca

<img src="man/figuresunnamed-chunk-4-6.png" width="100%" />

    #> >>> Processing group:: Phase
    #> >>>>Options for `theme`: light, bw, classic and classic2
    #> >>>>Options for 'legend.position' : none, left, right, bottom, top
    #> >>>>Options for 'legend.direction' : horizontal, vertical

<img src="man/figuresunnamed-chunk-4-7.png" width="100%" />

    #> ------------------------------------------------------
    #> >>> Processing method:: tsne

<img src="man/figuresunnamed-chunk-4-8.png" width="100%" />

    #> >>> Processing group:: Phase
    #> >>>>Options for `theme`: light, bw, classic and classic2
    #> >>>>Options for 'legend.position' : none, left, right, bottom, top
    #> >>>>Options for 'legend.direction' : horizontal, vertical

<img src="man/figuresunnamed-chunk-4-9.png" width="100%" />

    #> ------------------------------------------------------
    #> >>> Processing method:: umap

<img src="man/figuresunnamed-chunk-4-10.png" width="100%" />

    #> >>> Processing group:: Phase
    #> >>>>Options for `theme`: light, bw, classic and classic2
    #> >>>>Options for 'legend.position' : none, left, right, bottom, top
    #> >>>>Options for 'legend.direction' : horizontal, vertical

<img src="man/figuresunnamed-chunk-4-11.png" width="100%" /><img src="man/figuresunnamed-chunk-4-12.png" width="100%" />

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
                 findmarkers        = FALSE, # set to TRUE if user want to identify marker genes of each clusters
                 already_normalized = FALSE,
                 save_path          = "./test")
#> An object of class Seurat 
#> 33694 features across 1097 samples within 1 assay 
#> Active assay: RNA (33694 features, 0 variable features)
#> >>>>=== Palette option for random: 1: palette1; 2: palette2; 3: palette3;  4: palette4
#> [1] "'#5f75ae', '#64a841', '#e5486e', '#de8e06', '#b5aa0f', '#7ba39d', '#b15928', '#6a3d9a', '#cab2d6', '#374E55FF', '#00A1D5FF', '#6A6599FF', '#80796BFF', '#e31a1c', '#fb9a99', '#1f78b4', '#a6cee3', '#008280FF', '#3C5488FF', '#8F7700FF', '#666666', '#A20056FF', '#fdbf6f', '#E78AC3', '#b2df8a', '#386CB0', '#CD534CFF', '#008B45FF', '#7AA6DCFF', '#00A087FF', '#A73030FF', '#631879FF', '#003C67FF'"
#>  >>>-----  Data will be deposite in E:/18-Github/scsig/./test/
#> [1] ">>>----- Step-1: Quality control"
#> [1] "Mitochondrial genes: "
#>  [1] "MT-ND1"  "MT-ND2"  "MT-CO1"  "MT-CO2"  "MT-ATP8" "MT-ATP6" "MT-CO3" 
#>  [8] "MT-ND3"  "MT-ND4L" "MT-ND4"  "MT-ND5"  "MT-ND6"  "MT-CYB"
#> Ribosomal genes have been removed
#> species must be one of `human` or `mouse`.
#> Warning: CombinePlots is being deprecated. Plots should now be combined using
#> the patchwork system.
#> >>>------ For cell subset: Default parameters are : minFeature = 2000, minCount = 1000, cutoff_percent_mt = 25
#> >>>------ After filtering cells with low features, low count and high expression of mitochondrial genes
#> An object of class Seurat 
#> 33694 features across 1097 samples within 1 assay 
#> Active assay: RNA (33694 features, 0 variable features)
#> [1] ">>>------ Step-2: Data normalization and dimension reduction"
#> >>>------ Cell cycle scoring...
```

<img src="man/figuresunnamed-chunk-4-13.png" width="100%" />

    #>                  orig.ident nCount_RNA nFeature_RNA percent.mt percent.rp
    #> AAACCTGCACCTTGTC       TNBC      13976         3590  16.463938  10.997424
    #> AAACGGGAGTCCTCCT       TNBC       8732         2629   5.817682   9.757215
    #> AAACGGGTCCAGAGGA       TNBC      17138         4166  10.228731   9.738593
    #> AAAGATGCAGTTTACG       TNBC       9519         1703   5.200126  16.482824
    #> AAAGCAACAGGAATGC       TNBC      10285         2778  18.395722  13.349538
    #> AAAGCAATCGGAATCT       TNBC       8436         2822  16.797060   5.263158
    #>                      S.Score   G2M.Score Phase
    #> AAACCTGCACCTTGTC  0.07723506  0.42077041   G2M
    #> AAACGGGAGTCCTCCT -0.07738382 -0.13026069    G1
    #> AAACGGGTCCAGAGGA -0.03538864 -0.09533497    G1
    #> AAAGATGCAGTTTACG -0.04374442 -0.09006494    G1
    #> AAAGCAACAGGAATGC -0.07575323 -0.07794991    G1
    #> AAAGCAATCGGAATCT  0.12829573  0.35263693   G2M
    #> ------------------------------------------------------
    #> [1] ">>>------ Step-3: Find clusters"
    #> Computing nearest neighbor graph
    #> Computing SNN
    #> Resolution is 0.4
    #> >>>-- Findclusters when resolution = RNA_snn_res.0.4 ...
    #>   0   1   2   3   4   5 
    #> 434 303 230  55  54  21
    #> Resolution is 0.8
    #> >>>-- Findclusters when resolution = RNA_snn_res.0.8 ...
    #>   0   1   2   3   4   5   6   7   8   9 
    #> 281 180 168 153 136  54  49  28  27  21
    #> Resolution is 1.2
    #> >>>-- Findclusters when resolution = RNA_snn_res.1.2 ...
    #>   0   1   2   3   4   5   6   7   8   9  10  11 
    #> 214 167 149 136 117  91  71  54  28  27  22  21
    #> Resolution is 1.6
    #> >>>-- Findclusters when resolution = RNA_snn_res.1.6 ...
    #>   0   1   2   3   4   5   6   7   8   9  10  11  12  13 
    #> 169 148 143 116  91  78  65  65  54  49  43  28  27  21
    #> Resolution is 2
    #> >>>-- Findclusters when resolution = RNA_snn_res.2 ...
    #>   0   1   2   3   4   5   6   7   8   9  10  11  12  13 
    #> 166 148 142 113  93  78  68  66  54  49  44  28  27  21
    #> Scale for 'edge_colour' is already present. Adding another scale for
    #> 'edge_colour', which will replace the existing scale.
    #> >>>-- Findclusters when resolution = RNA_snn_res.1 ...
    #>   0   1   2   3   4   5   6   7   8   9  10 
    #> 219 180 168 148 136  67  54  49  28  27  21
    #> >>>-- These results will be used to draw dimension plot....
    #> ------------------------------------------------------
    #> Warning: The following arguments are not used: do.fast

<img src="man/figuresunnamed-chunk-4-14.png" width="100%" /><img src="man/figuresunnamed-chunk-4-15.png" width="100%" /><img src="man/figuresunnamed-chunk-4-16.png" width="100%" /><img src="man/figuresunnamed-chunk-4-17.png" width="100%" />

    #> ------------------------------------------------------
    #> >>> Processing method:: pca

<img src="man/figuresunnamed-chunk-4-18.png" width="100%" />

    #> >>> Processing group:: Phase
    #> >>>>Options for `theme`: light, bw, classic and classic2
    #> >>>>Options for 'legend.position' : none, left, right, bottom, top
    #> >>>>Options for 'legend.direction' : horizontal, vertical

<img src="man/figuresunnamed-chunk-4-19.png" width="100%" />

    #> ------------------------------------------------------
    #> >>> Processing method:: tsne

<img src="man/figuresunnamed-chunk-4-20.png" width="100%" />

    #> >>> Processing group:: Phase
    #> >>>>Options for `theme`: light, bw, classic and classic2
    #> >>>>Options for 'legend.position' : none, left, right, bottom, top
    #> >>>>Options for 'legend.direction' : horizontal, vertical

<img src="man/figuresunnamed-chunk-4-21.png" width="100%" />

    #> ------------------------------------------------------
    #> >>> Processing method:: umap

<img src="man/figuresunnamed-chunk-4-22.png" width="100%" />

    #> >>> Processing group:: Phase
    #> >>>>Options for `theme`: light, bw, classic and classic2
    #> >>>>Options for 'legend.position' : none, left, right, bottom, top
    #> >>>>Options for 'legend.direction' : horizontal, vertical

<img src="man/figuresunnamed-chunk-4-23.png" width="100%" /><img src="man/figuresunnamed-chunk-4-24.png" width="100%" />

``` r
head(sce@meta.data)
#>                  orig.ident nCount_RNA nFeature_RNA percent.mt percent.rp
#> AAACCTGCACCTTGTC       TNBC      13976         3590  16.463938  10.997424
#> AAACGGGAGTCCTCCT       TNBC       8732         2629   5.817682   9.757215
#> AAACGGGTCCAGAGGA       TNBC      17138         4166  10.228731   9.738593
#> AAAGATGCAGTTTACG       TNBC       9519         1703   5.200126  16.482824
#> AAAGCAACAGGAATGC       TNBC      10285         2778  18.395722  13.349538
#> AAAGCAATCGGAATCT       TNBC       8436         2822  16.797060   5.263158
#>                      S.Score   G2M.Score Phase RNA_snn_res.0.4 seurat_clusters
#> AAACCTGCACCTTGTC  0.07723506  0.42077041   G2M               3               8
#> AAACGGGAGTCCTCCT -0.07738382 -0.13026069    G1               2               1
#> AAACGGGTCCAGAGGA -0.03538864 -0.09533497    G1               0               0
#> AAAGATGCAGTTTACG -0.04374442 -0.09006494    G1               1               2
#> AAAGCAACAGGAATGC -0.07575323 -0.07794991    G1               0               5
#> AAAGCAATCGGAATCT  0.12829573  0.35263693   G2M               1               4
#>                  RNA_snn_res.0.8 RNA_snn_res.1.2 RNA_snn_res.1.6 RNA_snn_res.2
#> AAACCTGCACCTTGTC               7               8              11            11
#> AAACGGGAGTCCTCCT               1               4               3             3
#> AAACGGGTCCAGAGGA               0               0               6             7
#> AAAGATGCAGTTTACG               2               1               0             0
#> AAAGCAACAGGAATGC               0               6               5             5
#> AAAGCAATCGGAATCT               4               3              10            10
#>                  RNA_snn_res.1
#> AAACCTGCACCTTGTC             8
#> AAACGGGAGTCCTCCT             1
#> AAACGGGTCCAGAGGA             0
#> AAAGATGCAGTTTACG             2
#> AAAGCAACAGGAATGC             5
#> AAAGCAATCGGAATCT             4
```

\##find marker genes of each clusters

``` r
dong_find_markers(sce        = sce, 
                  group      = "RNA_snn_res.1", 
                  assay      = "RNA", 
                  verbose    = FALSE,
                  fig.type   = "pdf",
                  pt.size    = 0.5, 
                  cols       = "normal", 
                  seed       = 54321, 
                  show_col   = F, 
                  show_genes = 10,
                  hwidth     = 13, 
                  hheight    = 10,
                  show_plot  = T, 
                  path       = "./test" )
#> >>>---Assay used to estimation:
#> [1] "RNA"
#> >>> Idents of Seurat object is: RNA_snn_res.1
#> 
#>   0   1  10   2   3   4   5   6   7   8   9 
#> 219 180  21 168 148 136  67  54  49  28  27
#> >>>>=== Palette option for random: 1: palette1; 2: palette2; 3: palette3;  4: palette4
#>               p_val avg_log2FC pct.1 pct.2    p_val_adj cluster   gene
#> IFRD1  3.772950e-83   1.559992 0.963 0.315 1.131885e-79       0  IFRD1
#> CEBPD  8.688609e-81   1.428352 0.968 0.387 2.606583e-77       0  CEBPD
#> SAT1   8.133918e-77   1.386217 0.959 0.412 2.440175e-73       0   SAT1
#> CLU    1.793408e-72   1.417723 0.945 0.354 5.380225e-69       0    CLU
#> MMP7   8.665662e-72   1.514401 0.886 0.246 2.599699e-68       0   MMP7
#> SBSPON 1.308976e-70   1.400558 0.836 0.190 3.926929e-67       0 SBSPON
#> Scale for 'fill' is already present. Adding another scale for 'fill', which
#> will replace the existing scale.
#> >>> Idents of Seurat object is: RNA_snn_res.1
#> 
#>   0   1  10   2   3   4   5   6   7   8   9 
#> 219 180  21 168 148 136  67  54  49  28  27
#> >>>--- Results of DEs..
#>                         p_val  avg_log2FC pct.1 pct.2     p_val_adj cluster
#> IFRD1            3.772950e-83 1.559992098 0.963 0.315  1.131885e-79       0
#> CEBPD            8.688609e-81 1.428352286 0.968 0.387  2.606583e-77       0
#> SAT1             8.133918e-77 1.386217176 0.959 0.412  2.440175e-73       0
#> CLU              1.793408e-72 1.417722957 0.945 0.354  5.380225e-69       0
#> MMP7             8.665662e-72 1.514400588 0.886 0.246  2.599699e-68       0
#> SBSPON           1.308976e-70 1.400558488 0.836 0.190  3.926929e-67       0
#> RGS16            1.420488e-70 1.425266700 0.785 0.211  4.261463e-67       0
#> MDFI             7.946212e-69 1.362323639 0.936 0.342  2.383864e-65       0
#> MAOB             3.088222e-66 1.420113326 0.744 0.182  9.264667e-63       0
#> SOD3             1.715063e-64 1.367957423 0.918 0.295  5.145189e-61       0
#> CNTNAP3B         2.567220e-63 1.360174583 0.699 0.175  7.701659e-60       0
#> C1orf186         2.944844e-62 1.293212803 0.872 0.296  8.834533e-59       0
#> PHLDA1           9.566478e-62 1.394727848 0.845 0.304  2.869943e-58       0
#> MGP              1.347189e-61 1.208254549 0.982 0.449  4.041568e-58       0
#> GLS              1.788027e-61 1.296018831 0.977 0.339  5.364080e-58       0
#> HIST2H2BE        8.672002e-61 1.325082323 0.849 0.322  2.601600e-57       0
#> FBXO32           5.026748e-60 1.287911408 0.913 0.362  1.508024e-56       0
#> SLC43A3          1.719911e-59 1.262193654 0.918 0.346  5.159734e-56       0
#> PGM2L1           7.057559e-59 1.351603178 0.817 0.272  2.117268e-55       0
#> CHI3L1           1.199778e-57 1.218972893 0.927 0.393  3.599334e-54       0
#> TSPAN12          4.173623e-57 1.265473224 0.795 0.239  1.252087e-53       0
#> HSPA1B           6.596440e-57 1.232143990 0.936 0.354  1.978932e-53       0
#> TOB1             2.993929e-54 1.232668215 0.863 0.366  8.981788e-51       0
#> CHI3L2           6.360197e-53 1.299222841 0.808 0.293  1.908059e-49       0
#> TMEM61           7.203171e-53 1.249486098 0.790 0.268  2.160951e-49       0
#> BAMBI            1.223380e-52 1.199265554 0.877 0.364  3.670140e-49       0
#> SLC29A1          3.582348e-52 1.135602024 0.909 0.421  1.074704e-48       0
#> CLDN4            6.073088e-52 1.107878156 0.932 0.484  1.821926e-48       0
#> SFRP1            4.168156e-51 1.143015118 0.886 0.371  1.250447e-47       0
#> PLA2G4A          4.418644e-51 1.111917264 0.694 0.200  1.325593e-47       0
#> S100P            4.720082e-51 1.162690742 0.877 0.368  1.416025e-47       0
#> TIFA             1.417369e-50 1.278832613 0.804 0.331  4.252108e-47       0
#> NNMT             1.501401e-50 1.074234602 0.868 0.320  4.504204e-47       0
#> HSPA1A           4.863034e-50 1.137473006 0.922 0.403  1.458910e-46       0
#> CDC25B           1.280576e-48 1.196672057 0.813 0.375  3.841728e-45       0
#> TUBB2B           2.086562e-47 1.066091797 0.749 0.244  6.259686e-44       0
#> SFN              3.955321e-46 1.078065714 0.795 0.277  1.186596e-42       0
#> NCOA7            6.038263e-46 1.175883259 0.749 0.284  1.811479e-42       0
#> S100A1           6.516236e-46 1.127712078 0.950 0.407  1.954871e-42       0
#> EFNA1            1.513951e-45 1.061804494 0.890 0.424  4.541852e-42       0
#> BGN              1.853730e-45 1.026281677 0.881 0.396  5.561190e-42       0
#> CRISPLD1         5.608405e-45 1.080939798 0.858 0.355  1.682521e-41       0
#> OVOS2            5.868989e-45 1.204459302 0.758 0.274  1.760697e-41       0
#> TIMM10           7.497512e-45 1.038819378 0.932 0.384  2.249253e-41       0
#> FAM3C            9.375787e-44 1.095794553 0.836 0.344  2.812736e-40       0
#> GCSH             1.267788e-43 1.153690504 0.790 0.362  3.803364e-40       0
#> QPCT             4.697679e-42 1.122736862 0.731 0.274  1.409304e-38       0
#> SOX4             1.185534e-41 0.969062544 0.904 0.499  3.556603e-38       0
#> SNAI2            2.033597e-41 1.130112821 0.703 0.253  6.100791e-38       0
#> NRP2             5.270891e-41 1.214809000 0.740 0.341  1.581267e-37       0
#> RHOBTB3          2.810528e-40 1.073763375 0.776 0.378  8.431585e-37       0
#> MLLT11           1.144886e-39 1.066585939 0.703 0.239  3.434659e-36       0
#> TUBB2A           1.738548e-39 0.992486184 0.826 0.335  5.215644e-36       0
#> HIST1H2AE        8.942485e-39 0.920020976 0.817 0.312  2.682745e-35       0
#> PCOLCE2          2.566431e-38 0.903985264 0.726 0.240  7.699292e-35       0
#> EFNB2            2.378771e-37 1.046761235 0.416 0.088  7.136313e-34       0
#> EPHX2            2.746253e-37 1.070176892 0.584 0.191  8.238758e-34       0
#> FZD5             2.093218e-36 1.067760114 0.534 0.166  6.279653e-33       0
#> GPM6B            2.380867e-36 0.994626600 0.763 0.367  7.142602e-33       0
#> MGST1            5.004018e-36 0.978881650 0.881 0.366  1.501205e-32       0
#> CDKN1A           6.806004e-36 1.078042702 0.644 0.238  2.041801e-32       0
#> FAM89A           8.520912e-36 1.067307268 0.667 0.264  2.556274e-32       0
#> SDC4             1.446519e-35 1.003292281 0.767 0.400  4.339558e-32       0
#> PDZK1IP1         2.301219e-35 0.988541418 0.772 0.394  6.903657e-32       0
#> TUBA1A           7.181325e-35 0.940532977 0.822 0.407  2.154397e-31       0
#> HSPB1            1.600178e-34 0.957110811 0.813 0.424  4.800535e-31       0
#> PIM1             2.583987e-34 0.985475790 0.758 0.379  7.751960e-31       0
#> ARID5B           1.104897e-33 0.991856726 0.712 0.343  3.314690e-30       0
#> GMPR             1.177409e-33 0.962835142 0.580 0.200  3.532228e-30       0
#> PAM              1.835960e-33 0.954192478 0.767 0.400  5.507880e-30       0
#> IDH1             6.283752e-33 0.874931556 0.776 0.326  1.885126e-29       0
#> SLC6A15          8.096510e-33 0.924975927 0.361 0.068  2.428953e-29       0
#> MATN3            9.739107e-33 0.863559885 0.653 0.241  2.921732e-29       0
#> TNFSF13B         1.051880e-32 0.978770256 0.776 0.367  3.155640e-29       0
#> PRELP            1.251866e-32 0.975739147 0.530 0.165  3.755599e-29       0
#> DLX5             1.764271e-32 1.011783456 0.489 0.151  5.292812e-29       0
#> TMX2             2.397885e-32 0.849699656 0.826 0.390  7.193656e-29       0
#> FXYD3            2.633664e-32 0.886908804 0.845 0.416  7.900993e-29       0
#> TM4SF1           3.716359e-32 0.872348884 0.900 0.522  1.114908e-28       0
#> IGFBP7           6.201609e-32 0.691152277 0.813 0.322  1.860483e-28       0
#> PTGS2            1.099658e-30 0.925716898 0.461 0.134  3.298975e-27       0
#> DEFB1            1.420200e-30 0.876100036 0.470 0.138  4.260599e-27       0
#> SSRP1            1.452904e-30 0.800282551 0.913 0.397  4.358711e-27       0
#> ACOX2            1.930621e-30 0.846627274 0.333 0.062  5.791862e-27       0
#> ISG20            3.510338e-30 0.919357558 0.717 0.318  1.053102e-26       0
#> SLC25A37         5.842457e-30 0.863620115 0.758 0.362  1.752737e-26       0
#> SLPI             7.162164e-30 0.939338531 0.630 0.269  2.148649e-26       0
#> TNFRSF21         8.965120e-30 0.910699229 0.735 0.351  2.689536e-26       0
#> HILPDA           1.873251e-29 0.942726147 0.721 0.327  5.619752e-26       0
#> SOHLH1           1.937274e-29 0.809480180 0.648 0.246  5.811823e-26       0
#> IMPA2            2.731604e-29 0.883508742 0.767 0.448  8.194813e-26       0
#> HIST1H1C         4.461323e-29 0.842744672 0.744 0.378  1.338397e-25       0
#> KLF10            6.347973e-29 0.904509895 0.712 0.363  1.904392e-25       0
#> SNX22            1.896300e-28 0.829402567 0.635 0.243  5.688901e-25       0
#> MIA              1.974482e-28 0.761919587 0.731 0.322  5.923445e-25       0
#> GAS6             2.989243e-28 0.821275790 0.731 0.372  8.967728e-25       0
#> C10orf10         3.734138e-28 0.785726860 0.740 0.312  1.120241e-24       0
#> PROM1            1.000624e-27 0.864370564 0.740 0.437  3.001871e-24       0
#> TPD52L1          2.242688e-27 0.842424654 0.767 0.404  6.728064e-24       0
#> PLOD2            7.146965e-27 0.822363442 0.749 0.416  2.144090e-23       0
#> METTL7A          1.207430e-26 0.835258753 0.667 0.297  3.622289e-23       0
#> S100B            1.498262e-26 0.788584491 0.735 0.345  4.494786e-23       0
#> MYLK             2.024009e-26 0.806629798 0.708 0.379  6.072027e-23       0
#> LGALS3           4.014314e-26 0.779474788 0.808 0.464  1.204294e-22       0
#> RPL39L           4.682508e-26 0.832744898 0.731 0.408  1.404752e-22       0
#> CITED4           1.017776e-25 0.823593876 0.744 0.424  3.053329e-22       0
#> CKS2             1.254453e-25 0.582022452 0.790 0.364  3.763358e-22       0
#> SERPINE2         3.515351e-25 0.891530885 0.676 0.337  1.054605e-21       0
#> EPHX1            3.808113e-25 0.757342203 0.772 0.396  1.142434e-21       0
#> DBI              4.468752e-25 0.745287324 0.808 0.435  1.340626e-21       0
#> DSC3             4.564227e-25 0.797410042 0.607 0.260  1.369268e-21       0
#> C2orf82          6.739529e-25 0.778626982 0.808 0.433  2.021859e-21       0
#> NDRG2            1.195627e-24 0.758328632 0.721 0.385  3.586881e-21       0
#> MYO1B            1.754724e-24 0.772115683 0.717 0.341  5.264173e-21       0
#> CTHRC1           2.461201e-24 0.691244174 0.795 0.397  7.383602e-21       0
#> ZBTB10           3.373414e-23 0.829630014 0.726 0.371  1.012024e-19       0
#> HIST1H2BG        4.043896e-23 0.734385107 0.721 0.335  1.213169e-19       0
#> MMP2             4.324921e-23 0.699297035 0.521 0.186  1.297476e-19       0
#> LCN2             1.428433e-22 0.717334173 0.644 0.325  4.285299e-19       0
#> ROPN1B           2.695571e-22 0.711369745 0.708 0.382  8.086714e-19       0
#> GATA3            6.291415e-22 0.767276065 0.457 0.172  1.887424e-18       0
#> TINCR            1.132228e-21 0.857732768 0.566 0.270  3.396684e-18       0
#> TMEM158          1.364976e-21 0.863188766 0.621 0.292  4.094929e-18       0
#> AC005152.3       2.416992e-21 0.738301110 0.708 0.383  7.250976e-18       0
#> B3GNT7           3.153810e-21 0.737058996 0.621 0.278  9.461431e-18       0
#> SOCS3            4.209903e-21 0.740858207 0.676 0.385  1.262971e-17       0
#> BTG2             6.519972e-21 0.658895909 0.717 0.385  1.955992e-17       0
#> PPP1R1B          1.053197e-20 0.714575295 0.648 0.329  3.159592e-17       0
#> NDRG1            1.301653e-20 0.598150662 0.735 0.379  3.904958e-17       0
#> NAB1             1.618205e-20 0.742098219 0.530 0.221  4.854615e-17       0
#> CSRP1            2.007994e-20 0.685897344 0.717 0.397  6.023981e-17       0
#> KLK1             2.350812e-20 0.745330160 0.338 0.100  7.052435e-17       0
#> MAP1B            2.365277e-20 0.691207498 0.694 0.361  7.095832e-17       0
#> IFITM1           1.353873e-19 0.565005149 0.603 0.260  4.061619e-16       0
#> STAT1            1.640716e-19 0.657863531 0.703 0.380  4.922148e-16       0
#> GPRC5A           2.045899e-19 0.707710435 0.416 0.156  6.137697e-16       0
#> MYL9             2.425016e-19 0.661024261 0.781 0.427  7.275048e-16       0
#> IDI1             3.051989e-19 0.702571807 0.676 0.412  9.155967e-16       0
#> SCIN             6.033598e-19 0.660296146 0.388 0.133  1.810079e-15       0
#> MAFF             1.693120e-18 0.772641568 0.603 0.315  5.079359e-15       0
#> ELF3             7.957180e-18 0.638687562 0.749 0.441  2.387154e-14       0
#> C1R              1.050360e-17 0.493126047 0.648 0.301  3.151079e-14       0
#> AGT              1.222292e-17 0.482464211 0.589 0.313  3.666877e-14       0
#> AQP3             4.449435e-17 0.725041985 0.635 0.317  1.334830e-13       0
#> A2M              4.878270e-17 0.613972690 0.685 0.441  1.463481e-13       0
#> SERPINA5         4.979894e-17 0.604724828 0.598 0.279  1.493968e-13       0
#> SAA1             5.914639e-17 0.606629070 0.557 0.268  1.774392e-13       0
#> PMP22            7.253757e-17 0.650106914 0.795 0.478  2.176127e-13       0
#> SCRG1            2.492430e-16 0.601521487 0.607 0.335  7.477289e-13       0
#> DDIT3            3.403920e-16 0.698131892 0.598 0.310  1.021176e-12       0
#> KCNMB1           1.020658e-15 0.684657509 0.575 0.297  3.061973e-12       0
#> TOX              1.137847e-15 0.484568155 0.507 0.212  3.413540e-12       0
#> TMEM79           2.737024e-15 0.678122540 0.621 0.343  8.211072e-12       0
#> C1QL4            3.728415e-15 0.620746481 0.502 0.238  1.118525e-11       0
#> SDC2             7.720423e-15 0.550405099 0.699 0.423  2.316127e-11       0
#> ST3GAL1          8.442908e-15 0.595530543 0.626 0.341  2.532872e-11       0
#> PROX1            1.085949e-14 0.636524244 0.397 0.180  3.257847e-11       0
#> MAL2             1.146628e-14 0.545748324 0.721 0.444  3.439884e-11       0
#> TBC1D1           1.390697e-14 0.529787661 0.712 0.400  4.172092e-11       0
#> CKS1B            2.001815e-14 0.559035426 0.813 0.442  6.005444e-11       0
#> HES1             2.064653e-14 0.572454374 0.658 0.407  6.193958e-11       0
#> SOX8             2.190841e-14 0.543970884 0.644 0.346  6.572523e-11       0
#> RAD21            2.381904e-14 0.506573553 0.785 0.479  7.145711e-11       0
#> C1S              3.383181e-14 0.429813009 0.603 0.325  1.014954e-10       0
#> CDKN2A           3.445524e-14 0.640485607 0.562 0.300  1.033657e-10       0
#> VWA5A            5.288639e-14 0.671609808 0.457 0.229  1.586592e-10       0
#> SERTAD4          7.560385e-14 0.535571077 0.689 0.426  2.268116e-10       0
#> COL11A1          8.539204e-14 0.387604356 0.630 0.327  2.561761e-10       0
#> CTGF             1.096822e-13 0.490315354 0.584 0.327  3.290466e-10       0
#> CYP26B1          1.236106e-13 0.523176160 0.306 0.110  3.708317e-10       0
#> DTNB             1.496486e-13 0.579854168 0.589 0.312  4.489459e-10       0
#> SLC9A3R2         1.975579e-13 0.499422119 0.671 0.399  5.926737e-10       0
#> HIBCH            2.334968e-13 0.491984635 0.676 0.369  7.004904e-10       0
#> MARCKSL1         2.546472e-13 0.606934965 0.877 0.507  7.639416e-10       0
#> MAPK8IP2         2.703596e-13 0.645790964 0.384 0.183  8.110789e-10       0
#> CACYBP           2.849370e-13 0.527419538 0.813 0.469  8.548111e-10       0
#> NEGR1            4.736331e-13 0.512867489 0.388 0.172  1.420899e-09       0
#> HMGB1            5.214878e-13 0.412134609 0.749 0.438  1.564463e-09       0
#> SH3BGR           6.769144e-13 0.479647088 0.584 0.289  2.030743e-09       0
#> SDC1             7.911680e-13 0.530965254 0.758 0.477  2.373504e-09       0
#> FMO2             1.187497e-12 0.332697644 0.571 0.284  3.562491e-09       0
#> FXYD6            1.463838e-12 0.516398900 0.680 0.452  4.391515e-09       0
#> ANXA1            1.540196e-12 0.539571223 0.685 0.465  4.620588e-09       0
#> ACAN             1.767271e-12 0.498530662 0.256 0.087  5.301812e-09       0
#> C8orf46          1.831447e-12 0.563527588 0.502 0.264  5.494341e-09       0
#> ABCA5            2.512167e-12 0.606068496 0.411 0.208  7.536500e-09       0
#> MFAP2            3.731790e-12 0.489787977 0.658 0.428  1.119537e-08       0
#> CPED1            4.802931e-12 0.501883309 0.479 0.238  1.440879e-08       0
#> ROPN1            5.562422e-12 0.487013074 0.548 0.280  1.668727e-08       0
#> BOC              8.238101e-12 0.561520681 0.416 0.212  2.471430e-08       0
#> CCDC80           1.168518e-11 0.343683584 0.438 0.206  3.505553e-08       0
#> NEAT1            1.577208e-11 0.445954840 0.817 0.475  4.731624e-08       0
#> ENPP5            1.579939e-11 0.549167164 0.548 0.303  4.739817e-08       0
#> ARL6IP1          1.737999e-11 0.358694072 0.712 0.432  5.213997e-08       0
#> TACSTD2          2.427797e-11 0.531817004 0.644 0.437  7.283391e-08       0
#> GLRX             2.960483e-11 0.521589165 0.658 0.418  8.881448e-08       0
#> C1orf56          3.367150e-11 0.516753532 0.639 0.397  1.010145e-07       0
#> CP               3.608626e-11 0.445446872 0.616 0.355  1.082588e-07       0
#> TEKT3            4.348062e-11 0.506815710 0.735 0.462  1.304419e-07       0
#> COL9A3           4.438220e-11 0.404435470 0.635 0.408  1.331466e-07       0
#> SERPINB5         5.151870e-11 0.428809469 0.644 0.384  1.545561e-07       0
#> NOTUM            5.544840e-11 0.536933329 0.251 0.095  1.663452e-07       0
#> IFI6             6.012654e-11 0.435448897 0.621 0.386  1.803796e-07       0
#> CRYAB            6.338370e-11 0.546941677 0.763 0.524  1.901511e-07       0
#> PEG10            6.496811e-11 0.470280192 0.676 0.443  1.949043e-07       0
#> KRT15            1.030321e-10 0.538033152 0.470 0.267  3.090963e-07       0
#> PLPP3            1.065981e-10 0.406173323 0.301 0.129  3.197942e-07       0
#> GOLM1            1.211898e-10 0.456573763 0.671 0.443  3.635693e-07       0
#> TAGLN            1.220007e-10 0.467516047 0.635 0.401  3.660022e-07       0
#> FKBP4            1.489570e-10 0.439847928 0.703 0.453  4.468711e-07       0
#> TMEM106C         1.697788e-10 0.358673266 0.635 0.374  5.093365e-07       0
#> FAM84A           2.193344e-10 0.584261447 0.283 0.124  6.580031e-07       0
#> HRCT1            2.358669e-10 0.485901071 0.612 0.346  7.076007e-07       0
#> KRT7             2.400096e-10 0.597847517 0.840 0.567  7.200287e-07       0
#> PLOD3            2.908375e-10 0.481314304 0.589 0.392  8.725124e-07       0
#> RP11-357H14.17   3.151018e-10 0.483307653 0.361 0.175  9.453053e-07       0
#> CD24             3.682743e-10 0.627607882 0.918 0.592  1.104823e-06       0
#> SMTN             3.949216e-10 0.471144959 0.630 0.415  1.184765e-06       0
#> ERRFI1           4.237294e-10 0.511178208 0.603 0.388  1.271188e-06       0
#> SLC26A2          6.962083e-10 0.386122426 0.571 0.317  2.088625e-06       0
#> NET1             9.269293e-10 0.466676590 0.635 0.445  2.780788e-06       0
#> ITGB4            1.658648e-09 0.440610547 0.534 0.312  4.975944e-06       0
#> SLBP             1.685679e-09 0.378349608 0.680 0.459  5.057036e-06       0
#> GSDMC            1.980814e-09 0.515933278 0.393 0.216  5.942443e-06       0
#> KLF2             2.100521e-09 0.440623192 0.557 0.319  6.301562e-06       0
#> VTCN1            2.347033e-09 0.422485945 0.648 0.450  7.041098e-06       0
#> SERTAD4-AS1      2.792587e-09 0.443748339 0.603 0.408  8.377761e-06       0
#> LNX1             2.853983e-09 0.402816743 0.548 0.314  8.561948e-06       0
#> CLDN3            2.957870e-09 0.456107125 0.644 0.506  8.873610e-06       0
#> ACTL6A           4.762641e-09 0.407165598 0.658 0.452  1.428792e-05       0
#> DNM3OS           4.955710e-09 0.365827241 0.416 0.213  1.486713e-05       0
#> SOX9             6.628722e-09 0.472551407 0.621 0.379  1.988617e-05       0
#> CYP1B1           6.873691e-09 0.344862494 0.324 0.154  2.062107e-05       0
#> PPFIA1           9.112426e-09 0.449192431 0.658 0.484  2.733728e-05       0
#> IFIH1            9.226436e-09 0.494865604 0.393 0.213  2.767931e-05       0
#> LEMD1            1.017417e-08 0.427930458 0.635 0.468  3.052250e-05       0
#> CARHSP1          1.167290e-08 0.357641753 0.671 0.456  3.501871e-05       0
#> CLCN4            1.209899e-08 0.361963993 0.539 0.319  3.629698e-05       0
#> IER5             1.218124e-08 0.483532225 0.571 0.394  3.654371e-05       0
#> LINC01436        1.249303e-08 0.474518741 0.466 0.286  3.747909e-05       0
#> CTNND2           1.858791e-08 0.417081344 0.406 0.224  5.576372e-05       0
#> SELM             2.001064e-08 0.422750896 0.717 0.497  6.003191e-05       0
#> HSPH1            2.474822e-08 0.409447454 0.658 0.461  7.424466e-05       0
#> PTRF             2.632567e-08 0.253849801 0.525 0.292  7.897700e-05       0
#> SUN3             3.904669e-08 0.397350304 0.361 0.194  1.171401e-04       0
#> SQLE             4.330356e-08 0.471283081 0.562 0.419  1.299107e-04       0
#> KLF5             4.897222e-08 0.462418067 0.397 0.237  1.469167e-04       0
#> JHDM1D-AS1       5.390441e-08 0.444366837 0.461 0.289  1.617132e-04       0
#> AZGP1            5.727733e-08 0.241325586 0.548 0.338  1.718320e-04       0
#> ATF3             9.917755e-08 0.427069051 0.589 0.423  2.975327e-04       0
#> ACTA2            1.019649e-07 0.356680658 0.539 0.318  3.058946e-04       0
#> SPON2            1.345380e-07 0.312805403 0.511 0.306  4.036139e-04       0
#> TSPAN2           2.047881e-07 0.399446570 0.306 0.157  6.143644e-04       0
#> LAT2             2.057230e-07 0.364341806 0.484 0.281  6.171691e-04       0
#> NFKBIE           2.359130e-07 0.431179388 0.575 0.454  7.077390e-04       0
#> VANGL1           2.404358e-07 0.383840457 0.653 0.508  7.213075e-04       0
#> AP1M2            2.489947e-07 0.361116434 0.635 0.428  7.469842e-04       0
#> KRT18            2.774144e-07 0.407882739 0.662 0.527  8.322431e-04       0
#> TMSB15A          3.624471e-07 0.417973637 0.311 0.172  1.087341e-03       0
#> TRIB1            4.927133e-07 0.428470759 0.575 0.372  1.478140e-03       0
#> MBNL1-AS1        5.152052e-07 0.434863883 0.447 0.282  1.545616e-03       0
#> TTC39A           6.773256e-07 0.439606340 0.384 0.245  2.031977e-03       0
#> IFT172           7.739230e-07 0.409443524 0.447 0.287  2.321769e-03       0
#> CYP39A1          1.063189e-06 0.374581817 0.511 0.331  3.189566e-03       0
#> RP11-25K19.1     1.126529e-06 0.378016530 0.329 0.188  3.379588e-03       0
#> NCCRP1           1.200569e-06 0.383992141 0.384 0.230  3.601707e-03       0
#> HEY2             1.380682e-06 0.444606547 0.338 0.206  4.142047e-03       0
#> MALL             1.474776e-06 0.343663023 0.315 0.174  4.424329e-03       0
#> CCT5             1.491629e-06 0.326950298 0.708 0.489  4.474886e-03       0
#> HIST1H4E         1.535757e-06 0.368701446 0.297 0.166  4.607272e-03       0
#> ALDH1B1          1.572012e-06 0.433610805 0.498 0.319  4.716036e-03       0
#> NR2F2            1.613419e-06 0.319993098 0.575 0.438  4.840256e-03       0
#> ICAM1            1.690204e-06 0.307243791 0.534 0.334  5.070613e-03       0
#> DNM3             1.813633e-06 0.369663708 0.402 0.247  5.440898e-03       0
#> KLHDC3           1.912671e-06 0.417719901 0.731 0.539  5.738013e-03       0
#> RAB3IP           2.021385e-06 0.498278189 0.420 0.286  6.064155e-03       0
#> CIART            2.027292e-06 0.398249563 0.406 0.263  6.081875e-03       0
#> NKD2             2.111743e-06 0.471291465 0.283 0.159  6.335229e-03       0
#> ITGA10           2.351376e-06 0.411589491 0.333 0.199  7.054127e-03       0
#> ANXA2R           2.532931e-06 0.425585702 0.374 0.236  7.598792e-03       0
#> MYO6             3.098967e-06 0.329348808 0.594 0.452  9.296902e-03       0
#> PRSS33           3.483589e-06 0.283773272 0.502 0.378  1.045077e-02       0
#> KCTD1            3.776360e-06 0.377974276 0.507 0.336  1.132908e-02       0
#> RDH10            4.295487e-06 0.310956550 0.557 0.374  1.288646e-02       0
#> PRPS2            4.611984e-06 0.313725744 0.557 0.362  1.383595e-02       0
#> LAMB3            5.297453e-06 0.499302553 0.356 0.227  1.589236e-02       0
#> ELF5             5.561959e-06 0.352303561 0.548 0.405  1.668588e-02       0
#> FBXO2            6.589082e-06 0.370923276 0.511 0.334  1.976724e-02       0
#> NFKBIA           7.823244e-06 0.187071227 0.575 0.379  2.346973e-02       0
#> S100A6           8.082354e-06 0.258290781 0.616 0.420  2.424706e-02       0
#> TPM2             8.339087e-06 0.194819074 0.594 0.399  2.501726e-02       0
#> PRSS8            8.966309e-06 0.346172068 0.626 0.497  2.689893e-02       0
#> TFAP2A           9.593431e-06 0.315391786 0.607 0.400  2.878029e-02       0
#> HACD1            1.143103e-05 0.257174912 0.589 0.372  3.429308e-02       0
#> NFATC1           1.210873e-05 0.370947034 0.521 0.370  3.632618e-02       0
#> TPM1             1.482324e-05 0.344984996 0.671 0.535  4.446971e-02       0
#> FHL1             1.604050e-05 0.324888532 0.311 0.179  4.812149e-02       0
#> PFKP             1.609715e-05 0.343311447 0.580 0.401  4.829146e-02       0
#> MT2A             1.782448e-05 0.289516039 0.543 0.379  5.347344e-02       0
#> SORBS2           2.480045e-05 0.289171353 0.616 0.466  7.440134e-02       0
#> FGF13            2.571557e-05 0.320110164 0.502 0.346  7.714672e-02       0
#> CALD1            2.790859e-05 0.160842769 0.539 0.396  8.372578e-02       0
#> GGCT             2.826431e-05 0.371607987 0.575 0.467  8.479292e-02       0
#> SYT8             3.100998e-05 0.290187765 0.493 0.350  9.302994e-02       0
#> KRT8             3.109326e-05 0.345655777 0.676 0.526  9.327979e-02       0
#> RP11-19E11.1     3.588312e-05 0.292132853 0.539 0.399  1.076494e-01       0
#> PEG3             4.158112e-05 0.274075447 0.562 0.413  1.247434e-01       0
#> TUBB4B           4.285396e-05 0.291211666 0.648 0.484  1.285619e-01       0
#> MSRB3            4.417612e-05 0.230245695 0.365 0.226  1.325283e-01       0
#> CDC42EP1         4.544444e-05 0.323222057 0.603 0.424  1.363333e-01       0
#> S100A10          4.644296e-05 0.349160241 0.676 0.523  1.393289e-01       0
#> FBLN2            4.669809e-05 0.300985658 0.676 0.494  1.400943e-01       0
#> POSTN            5.214653e-05 0.076535150 0.406 0.208  1.564396e-01       0
#> SLC40A1          5.332646e-05 0.114845257 0.256 0.137  1.599794e-01       0
#> FAM46B           5.470363e-05 0.414750464 0.379 0.269  1.641109e-01       0
#> WDR34            6.277590e-05 0.281283964 0.607 0.479  1.883277e-01       0
#> BACE2            7.415540e-05 0.260595434 0.616 0.437  2.224662e-01       0
#> LMCD1            8.030932e-05 0.199649560 0.333 0.213  2.409280e-01       0
#> NPR3             9.267960e-05 0.392114317 0.324 0.220  2.780388e-01       0
#> IL32             1.019077e-04 0.270455247 0.575 0.451  3.057230e-01       0
#> PBX1             1.109515e-04 0.255574054 0.589 0.443  3.328544e-01       0
#> TUBB             1.393849e-04 0.302853344 0.621 0.508  4.181547e-01       0
#> PRSS22           1.475726e-04 0.345785234 0.507 0.364  4.427177e-01       0
#> SLC39A14         1.890736e-04 0.341754659 0.251 0.157  5.672209e-01       0
#> CA8              2.281186e-04 0.150898535 0.397 0.247  6.843559e-01       0
#> CTNNB1           2.404909e-04 0.248876690 0.589 0.462  7.214726e-01       0
#> OPRK1            2.641240e-04 0.377292393 0.279 0.188  7.923719e-01       0
#> IL17B            3.040906e-04 0.346828878 0.269 0.171  9.122719e-01       0
#> CSF3R            3.070025e-04 0.280108656 0.534 0.384  9.210075e-01       0
#> LAPTM4B          3.136325e-04 0.294690359 0.689 0.522  9.408975e-01       0
#> SOD2             3.200010e-04 0.143202761 0.543 0.400  9.600030e-01       0
#> PDGFRA           4.218266e-04 0.184621014 0.466 0.312  1.000000e+00       0
#> HEBP2            4.354306e-04 0.313597022 0.689 0.525  1.000000e+00       0
#> PYCR1            4.433429e-04 0.210324226 0.557 0.440  1.000000e+00       0
#> TSPAN5           4.605033e-04 0.261630831 0.388 0.263  1.000000e+00       0
#> MYBL1            5.139284e-04 0.297110154 0.347 0.248  1.000000e+00       0
#> NQO1             5.210380e-04 0.258419415 0.562 0.415  1.000000e+00       0
#> SPINT1           5.304187e-04 0.230345269 0.639 0.479  1.000000e+00       0
#> C6orf132         5.327927e-04 0.224473610 0.575 0.433  1.000000e+00       0
#> RTN4RL2          6.714819e-04 0.214717072 0.274 0.177  1.000000e+00       0
#> PLSCR1           7.926947e-04 0.249788815 0.493 0.427  1.000000e+00       0
#> CDK4             1.030787e-03 0.219358084 0.594 0.510  1.000000e+00       0
#> STMN1            1.062519e-03 0.253766275 0.662 0.502  1.000000e+00       0
#> HSP90AB1         1.081385e-03 0.555402625 0.941 0.606  1.000000e+00       0
#> MESP1            1.137044e-03 0.224367480 0.557 0.431  1.000000e+00       0
#> PTS              1.378415e-03 0.215252247 0.562 0.467  1.000000e+00       0
#> LIMCH1           1.493769e-03 0.225586192 0.521 0.378  1.000000e+00       0
#> DNTTIP1          1.533648e-03 0.206368832 0.607 0.478  1.000000e+00       0
#> PIR              1.581471e-03 0.267959819 0.429 0.328  1.000000e+00       0
#> SYNGR1           1.635452e-03 0.234999168 0.484 0.343  1.000000e+00       0
#> APP              1.743825e-03 0.230516573 0.635 0.473  1.000000e+00       0
#> ABL2             1.821361e-03 0.288592508 0.457 0.350  1.000000e+00       0
#> TM4SF1-AS1       1.926068e-03 0.172750425 0.429 0.311  1.000000e+00       0
#> TRAF3IP3         1.952107e-03 0.134173483 0.416 0.278  1.000000e+00       0
#> GPNMB            1.967344e-03 0.155910074 0.571 0.423  1.000000e+00       0
#> MFSD6            2.081845e-03 0.307107001 0.356 0.268  1.000000e+00       0
#> CAPN2            2.251340e-03 0.230527234 0.607 0.507  1.000000e+00       0
#> OVOL1            2.621174e-03 0.330656947 0.301 0.222  1.000000e+00       0
#> HMGA1            2.742425e-03 0.238170955 0.598 0.505  1.000000e+00       0
#> MB               3.316974e-03 0.261422150 0.352 0.262  1.000000e+00       0
#> S100A9           3.681765e-03 0.132736213 0.420 0.349  1.000000e+00       0
#> EPS8             3.693355e-03 0.163036377 0.297 0.207  1.000000e+00       0
#> VASN             3.802273e-03 0.201779557 0.557 0.445  1.000000e+00       0
#> NUDT4            3.970280e-03 0.168025950 0.502 0.368  1.000000e+00       0
#> NSG1             4.156893e-03 0.229428580 0.256 0.175  1.000000e+00       0
#> COL11A2          4.717472e-03 0.122871313 0.311 0.210  1.000000e+00       0
#> UBE2L6           4.763329e-03 0.244258939 0.479 0.344  1.000000e+00       0
#> PDLIM4           4.936921e-03 0.242334816 0.333 0.247  1.000000e+00       0
#> LAMB1            5.260351e-03 0.133483003 0.251 0.172  1.000000e+00       0
#> ABCC3            5.813125e-03 0.151068067 0.283 0.190  1.000000e+00       0
#> FOXC1            6.136866e-03 0.208225417 0.438 0.330  1.000000e+00       0
#> KIAA0040         6.397137e-03 0.254096668 0.352 0.273  1.000000e+00       0
#> INSIG1           6.507185e-03 0.214375113 0.475 0.371  1.000000e+00       0
#> BARD1            6.726283e-03 0.218670921 0.406 0.312  1.000000e+00       0
#> KCNQ1OT1         8.085349e-03 0.154232979 0.575 0.476  1.000000e+00       0
#> SNHG25           8.555355e-03 0.159653175 0.479 0.358  1.000000e+00       0
#> SELENBP1         8.736057e-03 0.250041180 0.420 0.335  1.000000e+00       0
#> PTPRF            8.983512e-03 0.136344824 0.553 0.444  1.000000e+00       0
#> GRB14            9.000440e-03 0.234343395 0.292 0.222  1.000000e+00       0
#> ERBB3            9.943854e-03 0.151642321 0.571 0.423  1.000000e+00       0
#> CD163           2.265804e-151 2.097133318 0.839 0.040 6.797411e-148       1
#> MS4A7           1.269907e-139 2.179377122 0.911 0.076 3.809722e-136       1
#> CYBB            2.875015e-139 2.079585000 0.906 0.076 8.625044e-136       1
#> FPR3            5.081822e-136 2.072075033 0.856 0.068 1.524547e-132       1
#> SLCO2B1         1.028761e-132 2.016369316 0.756 0.043 3.086282e-129       1
#> MSR1            1.424992e-132 2.080712488 0.811 0.063 4.274976e-129       1
#> SPI1            7.851604e-129 2.022431852 0.950 0.091 2.355481e-125       1
#> SRGN            2.533615e-128 2.237337600 1.000 0.097 7.600844e-125       1
#> C1QC            1.393465e-127 2.264643243 0.967 0.069 4.180395e-124       1
#> PLTP            2.444559e-125 2.050292721 0.883 0.091 7.333678e-122       1
#> PILRA           1.041725e-124 2.042424075 0.794 0.067 3.125175e-121       1
#> C3AR1           5.538652e-124 1.861659649 0.750 0.046 1.661595e-120       1
#> LYZ             3.773477e-122 1.996056813 0.944 0.095 1.132043e-118       1
#> LY96            1.076064e-121 2.040979831 0.944 0.130 3.228191e-118       1
#> CD53            6.583455e-121 1.932236375 0.878 0.094 1.975036e-117       1
#> LST1            1.080199e-120 1.929482713 0.922 0.095 3.240596e-117       1
#> FCGR1A          1.339522e-120 1.851754012 0.711 0.039 4.018566e-117       1
#> LILRB4          4.189542e-120 1.924232741 0.767 0.061 1.256863e-116       1
#> CD4             7.381740e-120 1.957507889 0.844 0.092 2.214522e-116       1
#> PYCARD          1.368481e-119 1.957686500 0.939 0.115 4.105444e-116       1
#> C5AR1           2.132272e-119 1.938728248 0.711 0.046 6.396817e-116       1
#> SAMSN1          4.129424e-119 1.889216143 0.761 0.059 1.238827e-115       1
#> CD68            1.367330e-118 2.199663762 0.989 0.111 4.101989e-115       1
#> HNMT            2.360067e-118 2.025060326 0.956 0.148 7.080202e-115       1
#> PLEK            3.149716e-118 1.880890118 0.822 0.079 9.449147e-115       1
#> MS4A6A          4.236927e-118 1.879844362 0.839 0.073 1.271078e-114       1
#> FCGR2A          5.618979e-118 2.084456788 0.861 0.100 1.685694e-114       1
#> LAPTM5          2.074347e-117 2.174537180 0.994 0.104 6.223040e-114       1
#> MARCH1          5.702576e-117 1.882434234 0.756 0.060 1.710773e-113       1
#> RGS1            1.654160e-116 1.994754069 0.917 0.081 4.962481e-113       1
#> LINC01094       1.242329e-115 1.876926993 0.672 0.037 3.726987e-112       1
#> TGFBI           1.438261e-115 1.929164722 0.928 0.118 4.314784e-112       1
#> AIF1            1.928002e-115 2.152238576 0.972 0.089 5.784007e-112       1
#> HLA-DQA1        1.156042e-114 1.942396503 0.911 0.089 3.468127e-111       1
#> MPEG1           1.695196e-114 1.851431356 0.717 0.052 5.085589e-111       1
#> ADAP2           6.108722e-114 1.880783109 0.778 0.071 1.832616e-110       1
#> HLA-DQA2        2.851108e-113 1.863985437 0.889 0.095 8.553323e-110       1
#> NPL             9.585038e-113 2.032107472 0.789 0.092 2.875511e-109       1
#> FYB             1.166580e-112 2.048994510 0.933 0.162 3.499739e-109       1
#> GMFG            4.994904e-112 1.918178914 0.944 0.126 1.498471e-108       1
#> FCER1G          2.053295e-110 2.222969431 0.994 0.113 6.159886e-107       1
#> SLC7A7          3.247417e-110 1.859194435 0.739 0.064 9.742251e-107       1
#> C1QA            6.803750e-110 2.244403224 0.944 0.073 2.041125e-106       1
#> TYROBP          2.395954e-109 2.232765642 0.989 0.096 7.187863e-106       1
#> FCGR3A          1.159189e-108 1.818030208 0.700 0.052 3.477566e-105       1
#> ACP5            1.804193e-108 1.557887551 0.939 0.107 5.412579e-105       1
#> HLA-DRB5        2.375502e-108 1.887471277 0.917 0.119 7.126507e-105       1
#> TREM2           1.560363e-107 1.771776367 0.611 0.028 4.681090e-104       1
#> KCTD12          5.724560e-107 1.897112950 0.811 0.099 1.717368e-103       1
#> PTPRC           3.729105e-106 1.810602564 0.856 0.107 1.118731e-102       1
#> ARHGAP18        3.881405e-106 1.915788452 0.861 0.118 1.164422e-102       1
#> HLA-DMB         4.420929e-106 1.873543830 0.856 0.099 1.326279e-102       1
#> IL18            2.786976e-105 1.927452790 0.789 0.097 8.360927e-102       1
#> GPR34           4.080015e-105 1.804962550 0.683 0.055 1.224004e-101       1
#> C1QB            6.428258e-105 2.194336689 0.928 0.074 1.928477e-101       1
#> DMXL2           2.680100e-104 1.838027005 0.678 0.057 8.040300e-101       1
#> HCLS1           1.117097e-102 1.817489099 0.800 0.096  3.351291e-99       1
#> KYNU            5.245320e-102 1.696003449 0.722 0.070  1.573596e-98       1
#> MS4A4A          9.796332e-102 1.750316956 0.667 0.052  2.938899e-98       1
#> CTSS            9.893530e-102 2.091655825 0.983 0.182  2.968059e-98       1
#> STAB1           5.240174e-101 1.814040575 0.733 0.081  1.572052e-97       1
#> CD86            9.226827e-101 1.700679331 0.694 0.061  2.768048e-97       1
#> FERMT3          1.593519e-100 1.777443599 0.789 0.099  4.780558e-97       1
#> CD14            3.158103e-100 2.146146590 0.967 0.110  9.474308e-97       1
#> C1orf162         8.511910e-99 1.772865067 0.700 0.070  2.553573e-95       1
#> HLA-DPA1         1.016416e-98 1.985953357 0.944 0.108  3.049249e-95       1
#> LY86             1.107582e-98 1.705271649 0.689 0.064  3.322747e-95       1
#> ITGB2            1.324244e-98 1.907956455 0.956 0.134  3.972732e-95       1
#> APOC1            2.068199e-98 2.323871091 0.972 0.070  6.204596e-95       1
#> RP11-1143G9.4    7.390642e-98 1.623942240 0.756 0.085  2.217193e-94       1
#> GIMAP4           2.296191e-97 1.756301091 0.694 0.069  6.888572e-94       1
#> HLA-DMA          4.566355e-97 1.949758957 0.950 0.156  1.369907e-93       1
#> CD84             7.243627e-97 1.620395562 0.744 0.086  2.173088e-93       1
#> CCL3             9.534689e-97 1.917088904 0.872 0.080  2.860407e-93       1
#> PLA2G7           1.631096e-96 1.665642660 0.539 0.022  4.893288e-93       1
#> ZEB2             4.918246e-96 1.701104903 0.800 0.110  1.475474e-92       1
#> MPP1             6.433474e-96 1.829926784 0.722 0.092  1.930042e-92       1
#> FTL              1.241253e-95 2.227685143 0.994 0.190  3.723759e-92       1
#> IGSF6            1.302723e-95 1.666581444 0.656 0.056  3.908169e-92       1
#> SAMHD1           3.664622e-95 1.708953330 0.828 0.113  1.099387e-91       1
#> IL7R             5.705275e-95 1.560077180 0.728 0.077  1.711583e-91       1
#> LAIR1            1.428356e-94 1.789540809 0.717 0.085  4.285069e-91       1
#> SLAMF8           2.301295e-94 1.704947940 0.589 0.038  6.903884e-91       1
#> CXCL3            1.534054e-93 1.715899336 0.644 0.048  4.602161e-90       1
#> IFI30            4.170329e-93 1.610763369 0.722 0.084  1.251099e-89       1
#> HCK              8.396885e-93 1.684355308 0.706 0.079  2.519066e-89       1
#> CTSB             1.438331e-92 2.168242091 0.983 0.206  4.314993e-89       1
#> CARD16           3.716011e-92 1.758392683 0.794 0.121  1.114803e-88       1
#> BMP2K            1.359352e-91 1.841549657 0.783 0.120  4.078056e-88       1
#> HLA-DRB1         1.779033e-91 1.922686895 0.950 0.140  5.337098e-88       1
#> HLA-DQB1         2.747419e-91 1.838762827 0.939 0.150  8.242257e-88       1
#> CYBA             3.352394e-91 1.953003056 1.000 0.287  1.005718e-87       1
#> APOE             4.101584e-91 2.288455782 0.978 0.094  1.230475e-87       1
#> RNASE1           6.932652e-91 2.168898719 0.928 0.125  2.079796e-87       1
#> CXCL2            8.161956e-91 1.705006648 0.639 0.060  2.448587e-87       1
#> LCP1             9.729568e-91 1.603137288 0.789 0.100  2.918870e-87       1
#> EVI2B            1.028590e-90 1.539924968 0.717 0.081  3.085769e-87       1
#> CSF1R            2.561024e-90 1.603350481 0.794 0.108  7.683073e-87       1
#> DAB2             9.862693e-90 1.789173175 0.839 0.160  2.958808e-86       1
#> CLEC2B           1.031577e-89 1.537268427 0.872 0.141  3.094732e-86       1
#> CLEC4E           1.692231e-89 1.568028239 0.528 0.026  5.076693e-86       1
#> CD48             3.082843e-89 1.631650311 0.689 0.082  9.248530e-86       1
#> CD74             2.555967e-88 2.092976595 1.000 0.119  7.667900e-85       1
#> MAFB             5.813067e-87 1.900309883 0.900 0.197  1.743920e-83       1
#> HLA-DRA          1.132897e-86 2.033273813 0.961 0.121  3.398690e-83       1
#> IL4I1            5.201584e-86 1.712231943 0.678 0.080  1.560475e-82       1
#> CD37             8.805182e-86 1.672019197 0.694 0.091  2.641555e-82       1
#> HSD17B11         1.299637e-85 1.607592750 0.783 0.124  3.898912e-82       1
#> BST2             8.503279e-85 1.472819243 0.878 0.148  2.550984e-81       1
#> HLA-DPB1         7.899891e-84 1.966057344 0.928 0.105  2.369967e-80       1
#> CAPG             8.022122e-84 1.905369919 0.972 0.282  2.406637e-80       1
#> CTSC             1.042665e-83 1.925582418 0.911 0.210  3.127996e-80       1
#> IRF8             1.082676e-83 1.616712269 0.517 0.032  3.248027e-80       1
#> AKR1B1           2.554827e-83 1.875231555 0.878 0.198  7.664481e-80       1
#> HAVCR2           3.098176e-83 1.577681521 0.589 0.053  9.294528e-80       1
#> IQGAP2           2.499839e-82 1.575667167 0.567 0.047  7.499517e-79       1
#> ADAMDEC1         6.845129e-81 1.569639402 0.567 0.050  2.053539e-77       1
#> GPX1             7.171263e-81 1.635805523 0.989 0.397  2.151379e-77       1
#> SLA              7.787544e-81 1.527606797 0.606 0.059  2.336263e-77       1
#> ALOX5AP          8.707546e-81 1.477772291 0.717 0.086  2.612264e-77       1
#> NCKAP1L          1.731812e-80 1.561901626 0.644 0.075  5.195436e-77       1
#> PTAFR            1.842006e-80 1.571740974 0.583 0.055  5.526019e-77       1
#> GPR65            6.105056e-80 1.588063812 0.583 0.058  1.831517e-76       1
#> FCGRT            6.800766e-80 1.876573810 0.961 0.255  2.040230e-76       1
#> EMP3             1.069772e-79 1.692376228 0.922 0.174  3.209315e-76       1
#> HPSE             1.530571e-79 1.501262699 0.444 0.015  4.591714e-76       1
#> UBB              1.635276e-79 1.709372565 0.972 0.164  4.905828e-76       1
#> HMOX1            6.538270e-79 1.856498783 0.833 0.213  1.961481e-75       1
#> MYO1F            7.946244e-79 1.618261907 0.783 0.146  2.383873e-75       1
#> RGS2             2.278563e-78 1.717508398 0.850 0.165  6.835689e-75       1
#> SLC8A1           8.085132e-78 1.522315503 0.550 0.047  2.425539e-74       1
#> LRRC25           9.352198e-78 1.557594144 0.539 0.045  2.805659e-74       1
#> PLAUR            1.036398e-77 1.767719343 0.933 0.219  3.109195e-74       1
#> EPB41L3          2.335929e-77 1.687371949 0.656 0.098  7.007788e-74       1
#> NCF2             4.029619e-77 1.507488540 0.533 0.044  1.208886e-73       1
#> XIST             7.816709e-77 1.444525136 0.833 0.149  2.345013e-73       1
#> LGALS9           1.025905e-76 1.686473239 0.744 0.136  3.077716e-73       1
#> MARCKS           1.191269e-76 1.842383633 0.944 0.252  3.573808e-73       1
#> SPP1             1.550822e-76 1.844647552 0.833 0.076  4.652465e-73       1
#> C1orf54          2.161538e-76 1.692005207 0.767 0.142  6.484613e-73       1
#> C15orf48         5.125087e-76 1.712277287 0.839 0.159  1.537526e-72       1
#> RNASE6           1.428786e-75 1.350992083 0.650 0.077  4.286358e-72       1
#> CXCR4            1.183914e-74 1.605061152 0.706 0.116  3.551742e-71       1
#> B2M              4.639472e-74 1.578276233 1.000 0.337  1.391841e-70       1
#> FCGR2B           6.411367e-74 1.373286891 0.583 0.061  1.923410e-70       1
#> NCF4             5.032241e-73 1.489965089 0.556 0.057  1.509672e-69       1
#> SNX10            9.108182e-73 1.538937969 0.750 0.141  2.732455e-69       1
#> AP1S2            3.711459e-72 1.785244217 0.917 0.269  1.113438e-68       1
#> EBI3             4.539844e-72 1.382090336 0.467 0.029  1.361953e-68       1
#> PLXND1           8.856293e-72 1.705057755 0.744 0.166  2.656888e-68       1
#> VSIG4            2.012593e-71 1.432637929 0.489 0.037  6.037778e-68       1
#> IGSF21           2.066291e-71 1.459936495 0.517 0.046  6.198872e-68       1
#> LSP1             1.528980e-70 1.543043109 0.817 0.146  4.586939e-67       1
#> CXCL8            3.821431e-70 1.652748287 0.628 0.077  1.146429e-66       1
#> CCR1             1.639432e-69 1.353066966 0.550 0.057  4.918297e-66       1
#> ITGAX            2.409315e-69 1.459047973 0.483 0.039  7.227944e-66       1
#> ARRB2            5.653747e-69 1.664298612 0.750 0.190  1.696124e-65       1
#> CFD              8.581476e-69 1.486724356 0.533 0.059  2.574443e-65       1
#> CTSL             1.323059e-68 1.849718142 0.894 0.289  3.969177e-65       1
#> RNASET2          4.190783e-68 1.693181797 0.917 0.286  1.257235e-64       1
#> DUSP1            4.721566e-68 1.679854712 0.939 0.274  1.416470e-64       1
#> CD83             1.201427e-67 1.441889496 0.689 0.129  3.604281e-64       1
#> LPAR6            2.787039e-67 1.409770561 0.561 0.068  8.361117e-64       1
#> LILRB1           3.120066e-67 1.425362338 0.478 0.040  9.360199e-64       1
#> TFEC             3.615417e-67 1.438852850 0.522 0.055  1.084625e-63       1
#> NR1H3            2.657296e-66 1.670245724 0.589 0.097  7.971887e-63       1
#> KCNMA1           3.319015e-66 1.497069131 0.500 0.053  9.957045e-63       1
#> OTOA             4.175512e-66 1.508543526 0.489 0.049  1.252654e-62       1
#> HLA-A            5.348931e-66 1.503611083 0.972 0.324  1.604679e-62       1
#> ATF5             9.933773e-66 1.646362152 0.700 0.157  2.980132e-62       1
#> TCF4             4.735485e-65 1.334286823 0.789 0.173  1.420645e-61       1
#> SLC2A3           5.430113e-65 1.386374061 0.717 0.135  1.629034e-61       1
#> CSTA             1.100804e-64 1.367225758 0.533 0.063  3.302411e-61       1
#> GATM             2.849653e-64 1.336041995 0.450 0.035  8.548959e-61       1
#> ARHGDIB          3.571722e-64 1.462964074 0.944 0.363  1.071517e-60       1
#> EVI2A            3.892757e-64 1.308911034 0.594 0.085  1.167827e-60       1
#> NCF1             4.823957e-64 1.577588510 0.639 0.118  1.447187e-60       1
#> EFHD2            7.792747e-64 1.649099838 0.883 0.238  2.337824e-60       1
#> GPR183           8.606674e-64 1.389386762 0.767 0.177  2.582002e-60       1
#> CCL4             8.915866e-64 1.693479287 0.733 0.107  2.674760e-60       1
#> IFI27L2          9.998560e-64 1.463342247 0.772 0.183  2.999568e-60       1
#> SLC16A10         2.168568e-63 1.409087929 0.483 0.049  6.505705e-60       1
#> LIMS1            2.476198e-63 1.671692543 0.906 0.260  7.428594e-60       1
#> SDS              1.566459e-61 1.483710937 0.550 0.084  4.699378e-58       1
#> FOLR2            5.076802e-61 1.291833764 0.428 0.034  1.523041e-57       1
#> CYTH4            5.435247e-61 1.401616899 0.522 0.069  1.630574e-57       1
#> LPXN             9.561082e-61 1.431083338 0.567 0.092  2.868325e-57       1
#> LILRB2           1.802283e-60 1.345054567 0.461 0.045  5.406850e-57       1
#> HTRA1            2.743729e-60 1.177946807 0.744 0.132  8.231186e-57       1
#> PLAU             3.710378e-60 1.208054204 0.661 0.117  1.113113e-56       1
#> MNDA             3.834340e-60 1.178971761 0.544 0.072  1.150302e-56       1
#> CMKLR1           7.702205e-60 1.349531923 0.367 0.019  2.310662e-56       1
#> ADAM8            9.467736e-60 1.337376549 0.528 0.071  2.840321e-56       1
#> GPNMB1           1.084968e-59 1.506291299 0.917 0.361  3.254905e-56       1
#> CST3             5.277789e-59 1.257251811 0.983 0.337  1.583337e-55       1
#> SOD21            5.971496e-59 1.467679359 0.922 0.332  1.791449e-55       1
#> ARL4C            6.050428e-59 1.572044556 0.833 0.258  1.815128e-55       1
#> APOC2            6.187777e-59 1.323273518 0.339 0.013  1.856333e-55       1
#> ID2              6.838419e-59 1.538664137 0.889 0.275  2.051526e-55       1
#> TGFB1            5.795105e-58 1.256861954 0.667 0.126  1.738531e-54       1
#> GPR84            6.019331e-58 1.259924815 0.406 0.031  1.805799e-54       1
#> PLXNC1           1.437997e-57 1.289199695 0.444 0.044  4.313992e-54       1
#> GPX3             3.138318e-56 1.443509867 0.633 0.134  9.414954e-53       1
#> HLA-DOA          3.229975e-56 1.349751861 0.411 0.036  9.689924e-53       1
#> TNFAIP8          1.216892e-55 1.397962891 0.639 0.144  3.650677e-52       1
#> OSM              1.288481e-55 1.255272633 0.411 0.036  3.865442e-52       1
#> ABCA1            1.505423e-55 1.605639057 0.683 0.179  4.516269e-52       1
#> EID1             2.991430e-55 0.985808565 0.800 0.164  8.974289e-52       1
#> IL10RA           4.309929e-55 1.386364246 0.533 0.088  1.292979e-51       1
#> SLC16A3          1.184818e-54 1.400420557 0.911 0.333  3.554455e-51       1
#> SGK1             1.538991e-54 1.442766171 0.811 0.241  4.616974e-51       1
#> PKIB             2.750180e-54 1.254616007 0.628 0.120  8.250539e-51       1
#> NCEH1            3.503204e-54 1.307983570 0.356 0.023  1.050961e-50       1
#> ALOX5            4.065259e-54 1.205053879 0.506 0.071  1.219578e-50       1
#> CCL3L3           4.646573e-54 1.256375203 0.494 0.067  1.393972e-50       1
#> SH2B3            6.185132e-54 1.325815340 0.461 0.058  1.855540e-50       1
#> SERPINF1         9.776414e-54 0.822294823 0.756 0.141  2.932924e-50       1
#> CCL5             1.055245e-53 1.501412724 0.650 0.149  3.165736e-50       1
#> MT1H             1.353146e-53 1.387156787 0.433 0.045  4.059438e-50       1
#> NAIP             2.397388e-53 1.353380565 0.483 0.069  7.192163e-50       1
#> ALDH2            4.259863e-53 1.428582973 0.739 0.212  1.277959e-49       1
#> RAB31            5.439366e-53 1.363331810 0.833 0.245  1.631810e-49       1
#> VMO1             1.072057e-52 1.353797689 0.417 0.046  3.216170e-49       1
#> RASSF4           1.655464e-52 1.467325662 0.794 0.258  4.966391e-49       1
#> OLR1             3.892376e-52 1.209902819 0.372 0.029  1.167713e-48       1
#> PARVG            8.123142e-52 1.256932049 0.428 0.048  2.436943e-48       1
#> GK               1.272228e-51 1.340535890 0.544 0.106  3.816683e-48       1
#> TRPV2            1.708451e-51 1.257305821 0.456 0.059  5.125354e-48       1
#> FMNL1            3.881642e-51 1.397291064 0.583 0.125  1.164493e-47       1
#> CPVL             4.074483e-51 1.367650959 0.622 0.148  1.222345e-47       1
#> EPSTI1           4.798729e-51 1.220962782 0.478 0.069  1.439619e-47       1
#> CLEC5A           1.074621e-50 1.213662463 0.372 0.032  3.223863e-47       1
#> HCST             1.279240e-50 1.436936193 0.600 0.136  3.837721e-47       1
#> ADAM28           1.751210e-50 1.173206101 0.467 0.063  5.253630e-47       1
#> CCL18            1.762125e-50 1.211533350 0.361 0.031  5.286374e-47       1
#> FAM198B          1.827314e-50 1.229985650 0.506 0.082  5.481942e-47       1
#> C2               3.392590e-50 1.488304751 0.639 0.161  1.017777e-46       1
#> SLC11A1          5.037061e-50 1.393858146 0.433 0.059  1.511118e-46       1
#> BHLHE41          5.776332e-50 1.215000921 0.400 0.041  1.732900e-46       1
#> MEF2C            1.584100e-49 1.176023547 0.667 0.162  4.752301e-46       1
#> BCAT1            1.937596e-49 1.297230555 0.456 0.065  5.812788e-46       1
#> TBXAS1           6.828323e-49 1.355220521 0.539 0.111  2.048497e-45       1
#> IFI16            3.485606e-48 1.137134642 0.744 0.195  1.045682e-44       1
#> STX11            4.757787e-48 1.199471732 0.328 0.023  1.427336e-44       1
#> SLC39A8          6.968645e-48 1.456973481 0.606 0.149  2.090593e-44       1
#> DOCK10           1.180862e-47 1.200138910 0.417 0.051  3.542585e-44       1
#> ABI3             2.048995e-47 1.142782198 0.500 0.084  6.146984e-44       1
#> PRDM1            3.649663e-47 1.156356915 0.517 0.094  1.094899e-43       1
#> DSE              4.691678e-47 1.344932933 0.672 0.184  1.407503e-43       1
#> AOAH             6.711505e-47 1.139214901 0.367 0.035  2.013452e-43       1
#> CYTIP            8.943846e-47 1.093497274 0.467 0.071  2.683154e-43       1
#> BTK              2.168037e-46 1.103615931 0.406 0.047  6.504110e-43       1
#> PTGS1            2.345855e-46 1.216842141 0.389 0.045  7.037566e-43       1
#> CACNA2D4         2.530708e-46 1.236178881 0.406 0.051  7.592125e-43       1
#> ADGRE2           2.621162e-46 1.241129611 0.394 0.048  7.863485e-43       1
#> NRP1             3.718134e-46 1.320507595 0.639 0.173  1.115440e-42       1
#> PECAM1           3.934705e-46 1.110459447 0.350 0.032  1.180412e-42       1
#> CLEC4A           3.980889e-46 1.172858083 0.478 0.077  1.194267e-42       1
#> LGMN             7.857728e-46 1.661809738 0.800 0.301  2.357319e-42       1
#> AQP9             9.234031e-46 1.194414459 0.278 0.013  2.770209e-42       1
#> DOK2             2.705986e-45 1.110131743 0.411 0.052  8.117957e-42       1
#> FAM26F           4.865350e-45 1.330629490 0.578 0.141  1.459605e-41       1
#> PTPRE            7.711065e-45 1.344332989 0.683 0.206  2.313320e-41       1
#> CASP1            1.207992e-44 1.241488792 0.572 0.141  3.623976e-41       1
#> IL10             1.332120e-44 1.079103296 0.278 0.014  3.996360e-41       1
#> MCOLN2           1.575378e-44 1.109952056 0.289 0.016  4.726134e-41       1
#> SEPP1            2.882367e-44 1.575393126 0.717 0.197  8.647100e-41       1
#> SLC1A3           3.157346e-44 1.307810802 0.522 0.118  9.472037e-41       1
#> MMP9             4.049272e-44 1.148620335 0.706 0.142  1.214782e-40       1
#> LXN              5.553579e-44 1.052444327 0.428 0.060  1.666074e-40       1
#> GNA15            7.155805e-44 1.145611979 0.494 0.096  2.146742e-40       1
#> BCL2A1           9.473988e-44 1.097548593 0.350 0.036  2.842196e-40       1
#> PDE4B            1.014945e-43 1.260883656 0.422 0.067  3.044836e-40       1
#> PARP14           4.148096e-43 1.191972027 0.583 0.148  1.244429e-39       1
#> CXorf21          5.045468e-43 1.002892024 0.389 0.048  1.513641e-39       1
#> CD93             7.572912e-43 1.042027827 0.350 0.036  2.271874e-39       1
#> SERPINA1         1.586751e-42 1.122096292 0.806 0.243  4.760252e-39       1
#> CD40             3.535827e-42 1.111426801 0.394 0.053  1.060748e-38       1
#> NR4A2            5.332113e-42 1.299903395 0.828 0.314  1.599634e-38       1
#> DLEU7            7.617882e-42 1.132495553 0.250 0.011  2.285364e-38       1
#> IL1B             1.260690e-41 1.132745661 0.417 0.058  3.782071e-38       1
#> CCL4L2           2.120454e-41 1.272371737 0.517 0.106  6.361361e-38       1
#> TLR4             2.173954e-41 1.148134644 0.344 0.038  6.521861e-38       1
#> CCL2             6.877955e-41 1.154886861 0.456 0.085  2.063387e-37       1
#> MAF              3.489703e-40 1.267201320 0.539 0.145  1.046911e-36       1
#> GAPLINC          9.011542e-40 1.115938975 0.339 0.039  2.703463e-36       1
#> ZNF331           2.493891e-39 1.424046245 0.544 0.153  7.481673e-36       1
#> VIM              3.074412e-39 0.999115287 0.883 0.475  9.223237e-36       1
#> AC092484.1       4.848439e-39 1.045268202 0.322 0.034  1.454532e-35       1
#> MT1M             5.459519e-39 1.123681435 0.289 0.025  1.637856e-35       1
#> CD300A           8.382411e-39 1.059271410 0.356 0.046  2.514723e-35       1
#> SIGLEC10         2.513883e-38 1.099333928 0.256 0.016  7.541650e-35       1
#> LIPA             2.616747e-38 1.447428163 0.683 0.267  7.850240e-35       1
#> MVP              3.372912e-38 1.164783906 0.678 0.228  1.011873e-34       1
#> TNFAIP3          4.099982e-38 1.251938849 0.578 0.168  1.229995e-34       1
#> NFKBIA1          1.911622e-37 1.223423703 0.800 0.344  5.734865e-34       1
#> OAS1             2.006945e-37 1.070066544 0.456 0.095  6.020835e-34       1
#> ZFP36            3.437829e-37 1.169316421 0.828 0.335  1.031349e-33       1
#> SLC2A5           3.925847e-37 1.051094514 0.278 0.024  1.177754e-33       1
#> TNFRSF4          1.931583e-36 1.039874689 0.311 0.036  5.794750e-33       1
#> SERPINB9         2.134913e-36 1.028013824 0.333 0.044  6.404738e-33       1
#> CXCL16           2.170862e-36 1.231467893 0.800 0.369  6.512587e-33       1
#> CD300E           2.628833e-36 1.033629163 0.289 0.028  7.886500e-33       1
#> CASP4            3.594342e-36 1.161703241 0.656 0.227  1.078303e-32       1
#> IL2RG            4.831829e-36 0.997712534 0.383 0.063  1.449549e-32       1
#> G0S2             6.645321e-36 1.322483827 0.494 0.121  1.993596e-32       1
#> HSD17B4          1.434547e-35 1.049143515 0.406 0.075  4.303641e-32       1
#> WAS              1.828600e-35 1.056362881 0.439 0.093  5.485799e-32       1
#> VASH1            1.869632e-35 1.001589370 0.367 0.057  5.608897e-32       1
#> PTGER4           2.303282e-35 1.154478501 0.411 0.084  6.909846e-32       1
#> IFI35            3.988623e-35 1.033076195 0.456 0.101  1.196587e-31       1
#> OLFML2B          4.128578e-35 0.720286735 0.500 0.107  1.238573e-31       1
#> APOBEC3C         7.777393e-35 1.044800905 0.372 0.063  2.333218e-31       1
#> FGL2             5.519149e-34 0.908771971 0.344 0.051  1.655745e-30       1
#> GIMAP7           1.123625e-33 1.048306589 0.339 0.052  3.370875e-30       1
#> HS3ST1           1.996741e-33 1.183673205 0.472 0.131  5.990223e-30       1
#> CKLF             2.008676e-33 0.929936688 0.811 0.361  6.026027e-30       1
#> FBP1             2.551073e-33 1.020564193 0.289 0.034  7.653218e-30       1
#> GALM             5.119010e-33 1.110523033 0.467 0.121  1.535703e-29       1
#> HLA-B            1.006338e-32 0.917147963 0.883 0.425  3.019014e-29       1
#> CLECL1           1.879344e-31 0.936382418 0.322 0.049  5.638033e-28       1
#> OSCAR            1.745583e-30 0.877497578 0.294 0.039  5.236748e-27       1
#> TREM1            7.604767e-30 0.969623963 0.256 0.029  2.281430e-26       1
#> SGPP1            2.549987e-29 0.937598044 0.267 0.034  7.649962e-26       1
#> FN1              4.088818e-29 0.874649012 0.617 0.162  1.226645e-25       1
#> BASP1            4.297837e-29 1.150780215 0.633 0.241  1.289351e-25       1
#> SPARC            5.081264e-29 0.397610618 0.639 0.154  1.524379e-25       1
#> TNFAIP2          8.463095e-29 1.173092926 0.572 0.230  2.538929e-25       1
#> CXCL1            2.203043e-28 0.836008934 0.289 0.044  6.609129e-25       1
#> PSMB9            4.264564e-28 0.915729730 0.550 0.180  1.279369e-24       1
#> RGS10            4.508234e-28 0.794909653 0.778 0.395  1.352470e-24       1
#> POU2F2           4.751078e-28 0.989036006 0.311 0.055  1.425323e-24       1
#> GBP1             6.436771e-28 0.861862093 0.333 0.062  1.931031e-24       1
#> SELPLG           1.889541e-27 0.978313862 0.322 0.061  5.668624e-24       1
#> ETV5             5.667972e-27 1.047885477 0.456 0.142  1.700392e-23       1
#> THBD             9.068410e-27 1.005947268 0.294 0.052  2.720523e-23       1
#> SAMD9L           1.352665e-26 0.873060332 0.317 0.059  4.057995e-23       1
#> C10orf54         1.797707e-26 0.919486543 0.417 0.115  5.393122e-23       1
#> HBEGF            2.285055e-26 0.946891586 0.328 0.067  6.855166e-23       1
#> UCP2             2.510207e-26 1.033823663 0.656 0.302  7.530621e-23       1
#> ABCC31           5.611428e-26 0.962558933 0.478 0.156  1.683428e-22       1
#> CD36             9.648723e-26 0.782540636 0.272 0.043  2.894617e-22       1
#> IER3             9.649722e-26 1.012815547 0.739 0.334  2.894917e-22       1
#> FGR              1.086683e-24 0.911611078 0.256 0.040  3.260048e-21       1
#> COQ2             1.713931e-24 1.031245739 0.422 0.136  5.141793e-21       1
#> PPM1N            2.016113e-24 0.767810364 0.267 0.044  6.048340e-21       1
#> IRF7             2.311687e-24 0.859219791 0.283 0.051  6.935061e-21       1
#> IL1RN            2.339238e-24 0.721342601 0.267 0.043  7.017714e-21       1
#> MT1G             2.608896e-24 1.203509948 0.528 0.180  7.826689e-21       1
#> SQRDL            2.681633e-24 0.936194456 0.728 0.364  8.044900e-21       1
#> PTPN6            3.075172e-24 1.087353852 0.606 0.299  9.225517e-21       1
#> PDPN             3.337422e-24 0.814026545 0.256 0.039  1.001227e-20       1
#> GBP2             7.383460e-24 0.890255964 0.511 0.182  2.215038e-20       1
#> CELF2            8.413001e-24 0.976007783 0.583 0.250  2.523900e-20       1
#> PSTPIP1          9.015476e-24 1.005117752 0.378 0.108  2.704643e-20       1
#> ERO1A            1.815777e-23 0.896668259 0.483 0.171  5.447332e-20       1
#> CD33             1.010934e-22 0.789521495 0.278 0.052  3.032802e-19       1
#> FABP5            1.275429e-22 0.845200386 0.717 0.385  3.826288e-19       1
#> PSMB8            4.935517e-22 0.853794124 0.633 0.282  1.480655e-18       1
#> MYO1G            1.044937e-21 0.840344035 0.250 0.045  3.134812e-18       1
#> CORO1A           1.114691e-21 0.884106386 0.628 0.260  3.344073e-18       1
#> PTPN22           1.904667e-21 0.665138196 0.261 0.048  5.714000e-18       1
#> BID              5.397609e-21 1.002352629 0.650 0.349  1.619283e-17       1
#> CSF2RA           1.058817e-20 0.847730799 0.400 0.138  3.176452e-17       1
#> MITF             1.436053e-20 0.906701203 0.444 0.171  4.308159e-17       1
#> SAMD9            3.522932e-20 0.794165059 0.328 0.086  1.056880e-16       1
#> FUCA1            4.820784e-19 1.134601518 0.544 0.282  1.446235e-15       1
#> PLEK2            1.084911e-18 0.741155600 0.278 0.068  3.254732e-15       1
#> RAC2             1.728678e-18 0.601685831 0.339 0.098  5.186033e-15       1
#> MGLL             1.817303e-17 0.973244142 0.456 0.214  5.451910e-14       1
#> LGALS3BP         4.758861e-17 0.547466845 0.467 0.174  1.427658e-13       1
#> IFI44            1.761449e-16 0.759315456 0.300 0.089  5.284346e-13       1
#> TRIM22           2.664297e-16 0.706714223 0.256 0.064  7.992891e-13       1
#> COLEC12          3.707896e-16 0.789677328 0.300 0.094  1.112369e-12       1
#> LINC00936        5.394112e-16 0.933559662 0.400 0.174  1.618233e-12       1
#> S100A91          1.145859e-15 0.595749908 0.672 0.302  3.437577e-12       1
#> ICAM11           4.548277e-15 0.769373581 0.606 0.328  1.364483e-11       1
#> TSC22D3          4.839046e-15 0.731878057 0.450 0.204  1.451714e-11       1
#> MT1F             5.159059e-15 1.112242298 0.550 0.284  1.547718e-11       1
#> SEPT6            6.779257e-15 0.641779390 0.361 0.134  2.033777e-11       1
#> ANPEP            7.143709e-15 0.720010186 0.522 0.249  2.143113e-11       1
#> LIMD2            9.573858e-15 0.793772559 0.439 0.202  2.872157e-11       1
#> FOXP1            2.955251e-14 0.698646294 0.456 0.201  8.865754e-11       1
#> CSTB             3.551574e-14 0.568760678 0.694 0.449  1.065472e-10       1
#> CYP27A1          3.985552e-14 0.909438000 0.328 0.130  1.195666e-10       1
#> RHOB             4.458858e-14 0.929174215 0.622 0.368  1.337657e-10       1
#> ECM1             8.692523e-14 0.484898101 0.289 0.089  2.607757e-10       1
#> PLSCR11          8.722458e-14 0.656090566 0.717 0.386  2.616738e-10       1
#> VAMP5            1.045500e-13 0.814586202 0.444 0.209  3.136499e-10       1
#> CD109            1.109722e-13 0.618938819 0.278 0.088  3.329165e-10       1
#> XAF1             1.136217e-13 0.721872190 0.261 0.082  3.408650e-10       1
#> FABP3            1.360537e-13 0.898035081 0.272 0.095  4.081612e-10       1
#> NFKBID           1.444444e-13 0.758166350 0.322 0.125  4.333333e-10       1
#> FRMD4A           5.137972e-13 0.787888904 0.561 0.335  1.541392e-09       1
#> BTG21            5.269624e-13 0.617345086 0.694 0.403  1.580887e-09       1
#> SLC25A5          1.007705e-12 0.546648922 0.694 0.457  3.023116e-09       1
#> DUSP6            1.638828e-12 0.780658867 0.561 0.335  4.916484e-09       1
#> CRIP1            2.738934e-12 0.524837199 0.333 0.130  8.216802e-09       1
#> FILIP1L          2.842756e-12 0.345780033 0.294 0.097  8.528269e-09       1
#> TCIRG1           8.829882e-12 0.618439389 0.567 0.338  2.648965e-08       1
#> LAT21            9.166238e-12 0.665222160 0.517 0.284  2.749871e-08       1
#> ID3              6.783991e-11 0.372295569 0.339 0.131  2.035197e-07       1
#> SCD              1.096781e-10 0.631958987 0.594 0.373  3.290344e-07       1
#> C3               1.825169e-10 0.626532707 0.422 0.218  5.475508e-07       1
#> FAM105A          3.189892e-10 0.602951599 0.467 0.274  9.569675e-07       1
#> NKG7             4.331094e-10 0.754384760 0.306 0.145  1.299328e-06       1
#> TNFSF10          1.016516e-09 0.504726315 0.267 0.105  3.049549e-06       1
#> CTSK             1.133383e-09 0.250392075 0.367 0.142  3.400149e-06       1
#> HLA-F            2.000362e-09 0.635141737 0.517 0.324  6.001087e-06       1
#> CXCL12           2.173376e-09 0.564289364 0.383 0.195  6.520127e-06       1
#> CITED2           3.225228e-09 0.579079439 0.589 0.370  9.675684e-06       1
#> TAP1             5.414258e-09 0.570675478 0.439 0.264  1.624277e-05       1
#> DUSP2            7.365111e-09 0.506109653 0.317 0.156  2.209533e-05       1
#> MX1              1.440752e-08 0.627425917 0.272 0.126  4.322255e-05       1
#> IFI61            2.518949e-08 0.446830530 0.633 0.394  7.556848e-05       1
#> SORL1            2.893820e-08 0.612671902 0.317 0.170  8.681459e-05       1
#> CD82             3.132037e-08 0.636013784 0.583 0.401  9.396112e-05       1
#> GNPTAB           4.366392e-08 0.548525122 0.478 0.323  1.309917e-04       1
#> RARRES1          6.314925e-08 0.434212654 0.456 0.241  1.894477e-04       1
#> PARVB            1.821161e-07 0.593785885 0.511 0.365  5.463484e-04       1
#> ACSL1            2.085486e-07 0.637124641 0.378 0.252  6.256457e-04       1
#> SAT11            4.263621e-07 0.363181796 0.733 0.480  1.279086e-03       1
#> RUNX3            5.242095e-07 0.461985091 0.250 0.121  1.572628e-03       1
#> S100A4           7.827932e-07 0.323252539 0.700 0.433  2.348380e-03       1
#> UBE2L61          9.082202e-07 0.533122996 0.522 0.341  2.724660e-03       1
#> COTL1            1.142051e-06 0.355181638 0.650 0.473  3.426152e-03       1
#> SLC25A19         1.197471e-06 0.641125017 0.256 0.135  3.592414e-03       1
#> FAM162A          2.210285e-06 0.387642667 0.633 0.432  6.630855e-03       1
#> ISG15            2.304187e-06 0.513464986 0.572 0.338  6.912561e-03       1
#> IFITM2           9.850243e-06 0.237788585 0.639 0.347  2.955073e-02       1
#> GCHFR            1.754352e-05 0.679524044 0.328 0.226  5.263056e-02       1
#> EPB41L2          2.869658e-05 0.423655922 0.372 0.250  8.608974e-02       1
#> DAPP1            3.757341e-05 0.477639164 0.311 0.206  1.127202e-01       1
#> MT1X             4.361332e-05 0.697057376 0.528 0.375  1.308400e-01       1
#> CAMK1            6.807683e-05 0.522577004 0.256 0.157  2.042305e-01       1
#> NRP21            9.478243e-05 0.178530453 0.633 0.378  2.843473e-01       1
#> MTSS1            9.846397e-05 0.472293048 0.389 0.291  2.953919e-01       1
#> SCCPDH           1.669454e-04 0.436164937 0.428 0.324  5.008363e-01       1
#> PPIF             1.993955e-04 0.437717157 0.483 0.396  5.981865e-01       1
#> CCDC85B          2.671860e-04 0.277061103 0.556 0.397  8.015581e-01       1
#> IFIH11           3.584384e-04 0.380480821 0.333 0.232  1.000000e+00       1
#> HLA-C            3.760018e-04 0.247261373 0.667 0.502  1.000000e+00       1
#> ARPC4            4.239630e-04 0.182294169 0.600 0.441  1.000000e+00       1
#> HIST1H4C         4.547754e-04 0.035208022 0.528 0.309  1.000000e+00       1
#> LMO4             4.916025e-04 0.363176363 0.428 0.325  1.000000e+00       1
#> ANKRD37          8.214114e-04 0.336994496 0.361 0.267  1.000000e+00       1
#> CD52             8.316715e-04 0.491552595 0.322 0.237  1.000000e+00       1
#> DNAJC15          9.980769e-04 0.278493071 0.467 0.340  1.000000e+00       1
#> TUBA1C           1.223188e-03 0.235356762 0.567 0.455  1.000000e+00       1
#> REEP4            1.259595e-03 0.389404930 0.372 0.288  1.000000e+00       1
#> C1orf21          1.408177e-03 0.433246519 0.256 0.180  1.000000e+00       1
#> IGKC             1.833586e-03 0.049363656 0.211 0.357  1.000000e+00       1
#> ARID5A           1.902528e-03 0.366564309 0.383 0.304  1.000000e+00       1
#> MT2A1            2.104308e-03 0.327846571 0.461 0.402  1.000000e+00       1
#> LGALS31          3.554714e-03 0.240569410 0.628 0.514  1.000000e+00       1
#> UPP1             3.573307e-03 0.344070923 0.383 0.313  1.000000e+00       1
#> MTHFD2           3.664804e-03 0.288064944 0.511 0.403  1.000000e+00       1
#> SOCS31           4.359613e-03 0.226859739 0.561 0.420  1.000000e+00       1
#> HBB              6.520677e-03 0.111571385 0.511 0.395  1.000000e+00       1
#> IFITM3           9.271228e-03 0.077144767 0.561 0.427  1.000000e+00       1
#> ST3GAL11         9.548284e-03 0.286390625 0.467 0.384  1.000000e+00       1
#> PCSK1N           1.290257e-82 1.954099485 0.982 0.304  3.870770e-79       2
#> CALML5           6.381724e-81 1.959466191 1.000 0.214  1.914517e-77       2
#> CFB              4.616001e-79 1.600838910 0.833 0.168  1.384800e-75       2
#> SEZ6L2           2.384760e-74 1.551716911 0.673 0.097  7.154279e-71       2
#> S100A61          3.319170e-74 1.628334506 0.970 0.367  9.957511e-71       2
#> APOD             4.931076e-65 1.819955699 0.851 0.159  1.479323e-61       2
#> S100A16          2.757717e-64 1.411321257 0.833 0.205  8.273152e-61       2
#> CRLF1            5.674772e-64 1.303683647 0.565 0.070  1.702432e-60       2
#> RBM3             1.057009e-63 1.438084382 0.976 0.352  3.171026e-60       2
#> PLP2             1.396215e-63 1.560352360 1.000 0.263  4.188645e-60       2
#> UBA1             6.080662e-60 1.521379804 0.929 0.319  1.824199e-56       2
#> BRK1             9.953265e-60 1.537480351 0.911 0.311  2.985980e-56       2
#> IFITM31          1.254088e-57 1.428224302 0.988 0.352  3.762265e-54       2
#> CDK16            2.341490e-55 1.389920652 0.958 0.354  7.024471e-52       2
#> KLK6             3.247577e-55 1.398216116 0.613 0.121  9.742731e-52       2
#> NDUFB11          3.506351e-55 1.325371167 0.970 0.352  1.051905e-51       2
#> MAGIX            4.954805e-53 1.220962270 0.869 0.270  1.486442e-49       2
#> KRT19            4.306958e-52 1.261608210 0.940 0.452  1.292087e-48       2
#> KRT17            4.418184e-52 1.420097283 0.798 0.237  1.325455e-48       2
#> APP1             2.154490e-51 1.230723762 0.923 0.429  6.463469e-48       2
#> PKP1             1.705213e-50 1.312866366 0.512 0.083  5.115638e-47       2
#> PPDPF            1.750953e-50 1.118553761 0.940 0.512  5.252858e-47       2
#> FMO21            1.213678e-48 1.454332863 0.774 0.263  3.641035e-45       2
#> S100A14          1.465223e-48 1.322536690 0.869 0.352  4.395670e-45       2
#> PDXK             3.326717e-48 1.226498097 0.929 0.414  9.980152e-45       2
#> ATP9A            3.701625e-48 1.282823068 0.786 0.273  1.110488e-44       2
#> KLK8             3.506263e-47 1.140662038 0.381 0.038  1.051879e-43       2
#> WDR13            1.078581e-46 1.246038128 0.905 0.355  3.235744e-43       2
#> CDH1             1.125649e-46 1.278627663 0.869 0.374  3.376947e-43       2
#> GABRP            1.206760e-46 1.324847265 0.845 0.395  3.620281e-43       2
#> PDIA6            1.654827e-46 1.182089388 0.899 0.452  4.964481e-43       2
#> EPPK1            6.496747e-46 1.150638476 0.631 0.146  1.949024e-42       2
#> CD241            2.576060e-45 1.011782034 0.988 0.597  7.728180e-42       2
#> BPIFB1           2.945983e-45 1.002455870 0.315 0.022  8.837949e-42       2
#> LINC01615        3.254380e-45 1.157055449 0.464 0.071  9.763139e-42       2
#> SUN1             2.313162e-44 1.251371627 0.881 0.384  6.939487e-41       2
#> SEC13            3.251161e-44 1.307655621 0.821 0.368  9.753483e-41       2
#> RPUSD3           4.789211e-44 1.259706719 0.845 0.350  1.436763e-40       2
#> UXT              7.578818e-44 1.132396014 0.923 0.379  2.273645e-40       2
#> DPM1             9.390941e-44 1.148389233 0.917 0.423  2.817282e-40       2
#> PQBP1            1.102863e-43 1.124301504 0.917 0.358  3.308590e-40       2
#> CUL7             3.252069e-43 1.198237312 0.768 0.273  9.756207e-40       2
#> PPP1R14C         4.355206e-43 1.149241007 0.702 0.215  1.306562e-39       2
#> GRXCR1           8.670015e-43 1.082298378 0.393 0.048  2.601004e-39       2
#> NUPR2            1.201985e-42 1.211209839 0.798 0.304  3.605955e-39       2
#> THUMPD3-AS1      1.576957e-42 1.288387976 0.833 0.381  4.730872e-39       2
#> SETD5            1.706951e-42 1.264466416 0.815 0.344  5.120852e-39       2
#> HEBP21           7.451852e-42 1.073630258 0.946 0.488  2.235556e-38       2
#> MEIS3            4.009822e-41 0.892400673 0.357 0.038  1.202947e-37       2
#> NPDC1            4.195373e-41 1.076980300 0.768 0.264  1.258612e-37       2
#> TIMP1            6.631051e-41 1.031282948 0.833 0.294  1.989315e-37       2
#> SMIM22           1.531437e-40 0.881805870 0.673 0.175  4.594312e-37       2
#> SOX18            1.690733e-40 0.941901422 0.530 0.110  5.072200e-37       2
#> FAM92A1          4.533433e-40 1.095802724 0.702 0.224  1.360030e-36       2
#> TACSTD21         4.217921e-39 1.112519651 0.887 0.405  1.265376e-35       2
#> METRN            4.939557e-39 1.153349972 0.887 0.400  1.481867e-35       2
#> TIMM17B          5.618021e-39 1.079089312 0.899 0.388  1.685406e-35       2
#> FNDC4            1.211102e-38 0.856070085 0.381 0.050  3.633307e-35       2
#> PDIA3            3.022925e-38 1.056226628 0.899 0.441  9.068774e-35       2
#> HDAC6            3.892145e-38 1.059129252 0.774 0.283  1.167644e-34       2
#> KRT81            4.052300e-38 1.087333900 0.875 0.466  1.215690e-34       2
#> NEAT11           4.209190e-38 0.927059674 0.881 0.482  1.262757e-34       2
#> DNAJC151         5.373054e-38 1.008987549 0.786 0.284  1.611916e-34       2
#> S100A41          5.771181e-38 1.170053731 0.815 0.416  1.731354e-34       2
#> KCNG1            8.694397e-38 1.184783707 0.643 0.220  2.608319e-34       2
#> FADS3            1.509182e-37 1.013429675 0.690 0.224  4.527546e-34       2
#> FSCN1            2.194676e-37 1.031278825 0.732 0.251  6.584027e-34       2
#> C6orf15          3.483412e-37 0.953755299 0.375 0.054  1.045024e-33       2
#> KLK7             4.818316e-37 0.966328686 0.315 0.033  1.445495e-33       2
#> ST6GAL1          9.057411e-37 1.010396121 0.649 0.189  2.717223e-33       2
#> PPP1R14B         1.045967e-36 1.088709243 0.851 0.410  3.137902e-33       2
#> IRF2BP2          1.622758e-36 1.090286840 0.863 0.423  4.868275e-33       2
#> PTPRF1           2.018675e-36 1.104373926 0.839 0.398  6.056024e-33       2
#> KANK4            2.458604e-36 0.979816943 0.333 0.040  7.375813e-33       2
#> LRP5             3.458088e-36 1.016705976 0.708 0.238  1.037427e-32       2
#> FBLN21           4.068051e-36 1.084730551 0.845 0.474  1.220415e-32       2
#> SPINT11          4.818744e-36 1.067575660 0.845 0.451  1.445623e-32       2
#> TADA3            7.130558e-36 1.174175066 0.833 0.399  2.139167e-32       2
#> MAD2L1BP         1.603359e-35 1.134865859 0.810 0.378  4.810078e-32       2
#> KIAA0513         1.665430e-35 1.035443866 0.518 0.122  4.996290e-32       2
#> CT83             1.692267e-35 0.929815195 0.506 0.111  5.076800e-32       2
#> PRICKLE1         3.321669e-35 1.059575337 0.679 0.240  9.965007e-32       2
#> APMAP            3.672821e-35 1.100194274 0.845 0.450  1.101846e-31       2
#> C6orf1321        5.016151e-35 1.075307598 0.845 0.392  1.504845e-31       2
#> TMEM132A         5.057663e-35 1.111353777 0.815 0.399  1.517299e-31       2
#> C1GALT1          8.784307e-35 1.193183246 0.750 0.354  2.635292e-31       2
#> ECEL1            1.699650e-34 0.946757178 0.482 0.105  5.098950e-31       2
#> PTGFRN           2.583597e-34 1.013744894 0.839 0.366  7.750792e-31       2
#> KRT6B            2.677381e-34 0.892803424 0.363 0.055  8.032143e-31       2
#> LINC00176        2.706504e-34 0.902377556 0.339 0.045  8.119513e-31       2
#> GPC1             4.521088e-34 1.076547378 0.732 0.310  1.356327e-30       2
#> TPT1-AS1         4.660613e-34 1.251777958 0.750 0.341  1.398184e-30       2
#> SLC35A2          5.309784e-34 1.036105630 0.845 0.376  1.592935e-30       2
#> C2orf821         7.449495e-34 1.026552818 0.845 0.447  2.234848e-30       2
#> KRBOX4           1.076938e-33 0.953439708 0.625 0.182  3.230814e-30       2
#> S100A2           1.970506e-33 1.129270045 0.613 0.213  5.911517e-30       2
#> IL17RC           2.591822e-33 1.112775319 0.726 0.304  7.775466e-30       2
#> PXDN             3.254946e-33 1.065498255 0.762 0.329  9.764839e-30       2
#> CHPF             3.684985e-33 1.057834747 0.815 0.382  1.105496e-29       2
#> TUNAR            4.789449e-33 0.908174286 0.327 0.044  1.436835e-29       2
#> DNAAF1           6.142114e-33 0.931988187 0.482 0.111  1.842634e-29       2
#> TRIB2            9.590956e-33 0.927196206 0.613 0.185  2.877287e-29       2
#> IL17RE           1.183240e-32 1.146091511 0.667 0.256  3.549720e-29       2
#> TRIM58           1.806883e-32 0.863077175 0.363 0.057  5.420649e-29       2
#> CACNA1A          2.740700e-32 1.070069481 0.655 0.248  8.222101e-29       2
#> PHACTR1          4.073163e-32 1.070248810 0.625 0.229  1.221949e-28       2
#> DSP              5.565082e-32 1.022944657 0.857 0.446  1.669525e-28       2
#> SHISA9           6.324755e-32 0.983415337 0.518 0.139  1.897427e-28       2
#> GPKOW            6.978272e-32 0.967833379 0.821 0.369  2.093481e-28       2
#> NAALADL2         7.293643e-32 1.053188008 0.613 0.210  2.188093e-28       2
#> VHL              8.901653e-32 1.033387560 0.702 0.288  2.670496e-28       2
#> FTSJ1            2.944611e-31 0.944291391 0.780 0.319  8.833832e-28       2
#> HOMER2           8.555605e-31 0.963090783 0.565 0.169  2.566682e-27       2
#> WDR45            8.847294e-31 0.922825575 0.875 0.414  2.654188e-27       2
#> PCDHB9           1.248196e-30 0.938864282 0.583 0.174  3.744587e-27       2
#> NOA1             1.559912e-30 1.028215062 0.738 0.340  4.679735e-27       2
#> SECTM1           1.678612e-30 0.878149437 0.571 0.172  5.035837e-27       2
#> EMC3             2.500377e-30 1.057186044 0.774 0.402  7.501131e-27       2
#> TRDC             3.009658e-30 0.902733736 0.304 0.042  9.028973e-27       2
#> PVRL4            4.351760e-30 0.960139424 0.768 0.335  1.305528e-26       2
#> TTLL3            4.996311e-30 0.907460524 0.500 0.132  1.498893e-26       2
#> PRRT3            6.716118e-30 0.849416226 0.548 0.159  2.014835e-26       2
#> IRF2BPL          6.864461e-30 1.010480947 0.810 0.407  2.059338e-26       2
#> OTUD5            8.171290e-30 0.997909667 0.804 0.405  2.451387e-26       2
#> CALB2            1.109177e-29 0.883834197 0.679 0.254  3.327532e-26       2
#> JAGN1            1.494715e-29 1.083765326 0.792 0.418  4.484144e-26       2
#> CHMP2B           1.621589e-29 1.036302863 0.780 0.428  4.864767e-26       2
#> TARBP1           1.818932e-29 0.877303458 0.679 0.252  5.456796e-26       2
#> LAMC2            5.438918e-29 0.793020792 0.393 0.080  1.631675e-25       2
#> SLC2A11          6.158521e-29 0.883890682 0.613 0.206  1.847556e-25       2
#> CTSF             1.576884e-28 1.016433213 0.798 0.409  4.730652e-25       2
#> KRT23            2.735937e-28 0.957384351 0.780 0.460  8.207811e-25       2
#> ITGB8            2.786782e-28 0.959691002 0.667 0.283  8.360345e-25       2
#> FUT2             3.390671e-28 0.909229064 0.292 0.043  1.017201e-24       2
#> FLNA             4.255114e-28 0.994675528 0.792 0.412  1.276534e-24       2
#> PTPN14           5.575642e-28 0.967643833 0.696 0.302  1.672692e-24       2
#> CEBPB            7.620687e-28 0.909727937 0.833 0.426  2.286206e-24       2
#> MISP             1.544382e-27 0.844301248 0.363 0.074  4.633147e-24       2
#> CAMK1D           2.146457e-27 0.852524803 0.565 0.175  6.439371e-24       2
#> FGF131           2.884762e-27 0.893459840 0.708 0.318  8.654285e-24       2
#> TPBG             3.031187e-27 0.935895131 0.708 0.323  9.093561e-24       2
#> MYO5B            4.582618e-27 0.941223210 0.631 0.271  1.374786e-23       2
#> LDOC1            5.301416e-27 0.959065530 0.744 0.363  1.590425e-23       2
#> SDR16C5          8.431991e-27 0.882633011 0.440 0.112  2.529597e-23       2
#> PRRX2            8.890554e-27 0.656661638 0.518 0.141  2.667166e-23       2
#> SEC11C           1.352817e-26 0.925830464 0.815 0.422  4.058450e-23       2
#> FKBP9            1.507420e-26 0.986122212 0.726 0.375  4.522259e-23       2
#> C1S1             1.643573e-26 0.727483948 0.726 0.318  4.930719e-23       2
#> NPTXR            1.648158e-26 0.917701434 0.554 0.185  4.944473e-23       2
#> LINC00960        1.906199e-26 0.828155159 0.530 0.164  5.718596e-23       2
#> ARAF             1.934139e-26 0.902081021 0.738 0.347  5.802418e-23       2
#> AGT1             1.983801e-26 0.901918952 0.696 0.309  5.951404e-23       2
#> TRPM8            2.310089e-26 0.634645624 0.536 0.154  6.930268e-23       2
#> VEGFA            3.342623e-26 1.105026784 0.685 0.315  1.002787e-22       2
#> KCNK12           3.849082e-26 0.841240872 0.274 0.040  1.154725e-22       2
#> MUC5B            4.197519e-26 0.904477477 0.554 0.188  1.259256e-22       2
#> ZNF674-AS1       6.053449e-26 0.706259905 0.452 0.118  1.816035e-22       2
#> NT5E             6.651837e-26 0.862687208 0.524 0.165  1.995551e-22       2
#> LAD1             7.531075e-26 0.918111020 0.762 0.375  2.259323e-22       2
#> ARMC1            9.407764e-26 0.743095863 0.661 0.236  2.822329e-22       2
#> SIPA1L2          9.525530e-26 0.731325616 0.554 0.174  2.857659e-22       2
#> SDC11            2.951759e-25 0.905729577 0.827 0.480  8.855277e-22       2
#> ACTR3B           3.276704e-25 0.900080796 0.613 0.255  9.830113e-22       2
#> LINC00518        5.874046e-25 0.723212397 0.369 0.081  1.762214e-21       2
#> AC022007.5       6.032957e-25 0.876519138 0.536 0.175  1.809887e-21       2
#> MYO61            7.221759e-25 0.932111478 0.780 0.426  2.166528e-21       2
#> CLDN6            1.432320e-24 0.686783652 0.554 0.175  4.296961e-21       2
#> WFDC2            1.752250e-24 0.701625569 0.387 0.093  5.256751e-21       2
#> SAA11            1.840752e-24 0.866612574 0.661 0.265  5.522256e-21       2
#> ARPC41           1.936766e-24 0.896846026 0.774 0.411  5.810299e-21       2
#> SEMA6A           2.741722e-24 0.792782359 0.333 0.069  8.225166e-21       2
#> GPR27            2.914567e-24 0.726968948 0.315 0.059  8.743702e-21       2
#> MFI2             2.942647e-24 0.886833627 0.643 0.285  8.827940e-21       2
#> CRYAB1           4.462558e-24 0.805023401 0.857 0.520  1.338767e-20       2
#> KRT16            5.875472e-24 0.716756240 0.411 0.109  1.762641e-20       2
#> LSR              6.318116e-24 0.829175256 0.833 0.457  1.895435e-20       2
#> CLCN41           7.544701e-24 0.934612233 0.655 0.310  2.263410e-20       2
#> PORCN            1.225349e-23 0.677767212 0.429 0.111  3.676048e-20       2
#> THUMPD3          1.463985e-23 0.904581455 0.696 0.338  4.391956e-20       2
#> PRAF2            3.328354e-23 0.837224550 0.685 0.335  9.985061e-20       2
#> MTSS1L           3.497437e-23 0.903918464 0.696 0.355  1.049231e-19       2
#> EGLN1            4.023329e-23 0.737494222 0.619 0.260  1.206999e-19       2
#> IGSF3            5.585077e-23 0.885791519 0.750 0.398  1.675523e-19       2
#> SYNCRIP          5.586783e-23 0.791298072 0.804 0.463  1.676035e-19       2
#> NRBP2            7.281676e-23 0.808408767 0.655 0.288  2.184503e-19       2
#> IRX4             7.562675e-23 0.720589313 0.339 0.075  2.268803e-19       2
#> CHST7            1.230283e-22 0.743963282 0.387 0.100  3.690848e-19       2
#> MAL21            1.444681e-22 0.896671958 0.750 0.454  4.334044e-19       2
#> MLLT4            2.216508e-22 0.860238413 0.708 0.365  6.649525e-19       2
#> IFIT2            2.257379e-22 0.724871427 0.423 0.123  6.772138e-19       2
#> FBXO9            2.411620e-22 0.778908413 0.714 0.362  7.234861e-19       2
#> DDR1             3.039934e-22 0.819545953 0.756 0.383  9.119802e-19       2
#> SLC15A1          3.134490e-22 0.712693601 0.250 0.040  9.403470e-19       2
#> NDUFAF6          4.388669e-22 0.862297473 0.750 0.405  1.316601e-18       2
#> THBS1            4.591737e-22 0.790472678 0.655 0.298  1.377521e-18       2
#> FLNB             7.038140e-22 0.823576376 0.720 0.378  2.111442e-18       2
#> CLDN7            7.487579e-22 0.811864533 0.845 0.498  2.246274e-18       2
#> SLC52A1          1.261857e-21 0.892516239 0.440 0.144  3.785572e-18       2
#> PEG31            2.841414e-21 0.820310411 0.720 0.393  8.524243e-18       2
#> SLC2A12          5.496400e-21 0.749790384 0.375 0.101  1.648920e-17       2
#> KCNQ1OT11        5.982261e-21 0.836233989 0.798 0.441  1.794678e-17       2
#> NMB              6.350994e-21 0.814550743 0.571 0.243  1.905298e-17       2
#> IL321            7.207894e-21 0.804500882 0.780 0.421  2.162368e-17       2
#> RASD1            8.292657e-21 0.609145093 0.470 0.155  2.487797e-17       2
#> PADI2            8.880541e-21 0.717570270 0.673 0.311  2.664162e-17       2
#> FZD8             1.045364e-20 0.701975245 0.298 0.065  3.136092e-17       2
#> KLF51            1.397161e-20 0.726176804 0.560 0.216  4.191484e-17       2
#> MMP15            1.522667e-20 0.820561871 0.696 0.360  4.568002e-17       2
#> MAP2             1.556484e-20 0.793592439 0.488 0.172  4.669451e-17       2
#> FBXO21           2.080135e-20 0.614937195 0.673 0.314  6.240405e-17       2
#> EFNA11           2.862531e-20 0.748221967 0.857 0.455  8.587593e-17       2
#> BCAS4            3.461976e-20 0.741898431 0.601 0.259  1.038593e-16       2
#> FXYD31           5.392578e-20 0.765385300 0.768 0.453  1.617773e-16       2
#> CD55             6.120524e-20 0.764047064 0.768 0.421  1.836157e-16       2
#> KRT71            6.178207e-20 0.760222321 0.869 0.577  1.853462e-16       2
#> RRP36            6.705121e-20 0.780524463 0.768 0.425  2.011536e-16       2
#> ALDOC            7.136370e-20 0.704378446 0.256 0.048  2.140911e-16       2
#> RGN              9.222575e-20 0.634485493 0.304 0.069  2.766773e-16       2
#> ANPEP1           1.156272e-19 0.434620134 0.589 0.240  3.468815e-16       2
#> TBC1D25          1.568702e-19 0.675726984 0.637 0.262  4.706107e-16       2
#> KRT86            1.866856e-19 1.043512932 0.500 0.228  5.600567e-16       2
#> ATP1B1           1.923282e-19 0.734185924 0.756 0.413  5.769845e-16       2
#> LEMD11           2.798607e-19 0.739087683 0.821 0.443  8.395821e-16       2
#> OGG1             2.924069e-19 0.702612156 0.560 0.221  8.772207e-16       2
#> MORN3            4.410334e-19 0.617424947 0.405 0.122  1.323100e-15       2
#> MAFK             4.949777e-19 0.703775544 0.452 0.156  1.484933e-15       2
#> GBA              5.059464e-19 0.778321985 0.762 0.463  1.517839e-15       2
#> EGLN3            6.452137e-19 0.727652066 0.339 0.093  1.935641e-15       2
#> MANEAL           1.379253e-18 0.600753717 0.274 0.058  4.137759e-15       2
#> MUC1             1.453265e-18 0.793446534 0.607 0.330  4.359794e-15       2
#> TMEM63A          1.588830e-18 0.723330235 0.702 0.368  4.766491e-15       2
#> SLC2A1           1.613048e-18 0.891582445 0.500 0.215  4.839145e-15       2
#> TNNC2            1.668299e-18 0.683103931 0.250 0.051  5.004896e-15       2
#> STARD10          1.789015e-18 0.707393658 0.708 0.381  5.367046e-15       2
#> KRT14            2.005467e-18 0.604762544 0.262 0.055  6.016401e-15       2
#> C1R1             2.196536e-18 0.480978850 0.619 0.325  6.589608e-15       2
#> TBC1D7           2.635535e-18 0.783680275 0.643 0.323  7.906606e-15       2
#> PVRL1            3.303197e-18 0.538776286 0.506 0.186  9.909592e-15       2
#> GAN              5.507192e-18 0.706489487 0.571 0.255  1.652158e-14       2
#> COL6A2           6.204744e-18 0.497679906 0.726 0.385  1.861423e-14       2
#> TLDC1            6.898662e-18 0.713073479 0.685 0.366  2.069599e-14       2
#> ZNF654           6.913678e-18 0.627133760 0.482 0.175  2.074103e-14       2
#> SULF2            7.019345e-18 0.754563169 0.619 0.288  2.105803e-14       2
#> AHNAK2           7.148124e-18 0.709884133 0.500 0.196  2.144437e-14       2
#> SAA2             7.583696e-18 0.673124859 0.369 0.113  2.275109e-14       2
#> GATA6            7.970449e-18 0.637925184 0.464 0.169  2.391135e-14       2
#> TFAP2A1          1.013744e-17 0.760141716 0.690 0.396  3.041231e-14       2
#> TJP1             1.238242e-17 0.585711361 0.595 0.259  3.714726e-14       2
#> SLC26A21         1.640032e-17 0.836330718 0.583 0.328  4.920096e-14       2
#> FOXC11           1.845177e-17 0.732240238 0.607 0.306  5.535532e-14       2
#> ILF2             2.137684e-17 0.623693348 0.714 0.399  6.413053e-14       2
#> COL9A31          2.728332e-17 0.762038941 0.696 0.409  8.184996e-14       2
#> KLK10            3.306814e-17 0.739381121 0.250 0.056  9.920441e-14       2
#> PSORS1C1         4.091681e-17 0.614475654 0.423 0.145  1.227504e-13       2
#> AQP31            4.170262e-17 0.593127999 0.649 0.332  1.251079e-13       2
#> PRRT3-AS1        4.474331e-17 0.625546456 0.310 0.083  1.342299e-13       2
#> GS1-124K5.4      5.715806e-17 0.635304094 0.429 0.151  1.714742e-13       2
#> CGGBP1           5.876478e-17 0.756340081 0.673 0.411  1.762943e-13       2
#> RBP1             7.146320e-17 0.671082659 0.613 0.307  2.143896e-13       2
#> GALNT2           8.422414e-17 0.707493835 0.649 0.341  2.526724e-13       2
#> RGS20            1.256143e-16 0.578871340 0.315 0.086  3.768430e-13       2
#> MYOF             1.694140e-16 0.662569217 0.560 0.249  5.082420e-13       2
#> ITGB41           1.828854e-16 0.687928045 0.601 0.312  5.486563e-13       2
#> CCDC167          3.062651e-16 0.678022594 0.726 0.467  9.187952e-13       2
#> EMP1             3.115115e-16 0.708176086 0.708 0.422  9.345346e-13       2
#> ELK1             3.607793e-16 0.548136431 0.577 0.242  1.082338e-12       2
#> ANXA3            3.820917e-16 0.551099265 0.292 0.076  1.146275e-12       2
#> HSPB2            5.790328e-16 0.702966558 0.607 0.316  1.737098e-12       2
#> TMEM208          7.491519e-16 0.643026409 0.750 0.462  2.247456e-12       2
#> DSEL             8.052078e-16 0.614068684 0.613 0.302  2.415623e-12       2
#> SORBS21          8.647094e-16 0.681135152 0.708 0.457  2.594128e-12       2
#> NTN1             8.947327e-16 0.634843028 0.512 0.219  2.684198e-12       2
#> GSDMC1           1.044094e-15 0.602899351 0.500 0.207  3.132281e-12       2
#> TMEM200A         1.338143e-15 0.603670721 0.268 0.068  4.014428e-12       2
#> CRABP2           1.370239e-15 0.497535298 0.423 0.175  4.110718e-12       2
#> OVOL11           1.374150e-15 0.536318602 0.482 0.194  4.122449e-12       2
#> HYOU1            1.530610e-15 0.558093979 0.619 0.318  4.591829e-12       2
#> H2AFY2           1.650963e-15 0.497274740 0.464 0.171  4.952889e-12       2
#> FBXL16           1.685661e-15 0.638101446 0.369 0.126  5.056982e-12       2
#> RNF8             1.719438e-15 0.629099415 0.655 0.360  5.158315e-12       2
#> SOX91            1.726328e-15 0.653038662 0.661 0.385  5.178983e-12       2
#> RTP4             1.766251e-15 0.594158158 0.530 0.228  5.298753e-12       2
#> CCDC64B          1.994943e-15 0.631383938 0.702 0.399  5.984829e-12       2
#> SLC6A14          2.251892e-15 0.493008365 0.280 0.073  6.755677e-12       2
#> LPIN1            2.871789e-15 0.689374522 0.571 0.294  8.615367e-12       2
#> PAX1             2.926296e-15 0.663984869 0.339 0.110  8.778887e-12       2
#> GATA6-AS1        3.120068e-15 0.484225298 0.274 0.071  9.360204e-12       2
#> ERRFI11          3.529717e-15 0.667230140 0.661 0.390  1.058915e-11       2
#> PGBD5            3.545799e-15 0.553037634 0.345 0.109  1.063740e-11       2
#> SERPINB51        4.339033e-15 0.731365870 0.667 0.394  1.301710e-11       2
#> BRPF1            4.990405e-15 0.439884752 0.351 0.109  1.497122e-11       2
#> MFSD3            9.139976e-15 0.589903176 0.661 0.357  2.741993e-11       2
#> ABHD11-AS1       9.736602e-15 0.587032802 0.256 0.066  2.920981e-11       2
#> HSPA5            9.991289e-15 0.616168504 0.738 0.471  2.997387e-11       2
#> RP11-161M6.2     1.005908e-14 0.594823216 0.506 0.219  3.017724e-11       2
#> RIC3             1.008536e-14 0.608997767 0.619 0.357  3.025608e-11       2
#> MTMR14           1.430911e-14 0.665756293 0.661 0.372  4.292733e-11       2
#> CLDN31           1.438354e-14 0.653874320 0.810 0.483  4.315063e-11       2
#> ZNRF2            1.504523e-14 0.506987729 0.429 0.157  4.513569e-11       2
#> ARFGEF3          1.682782e-14 0.567220483 0.387 0.143  5.048346e-11       2
#> DDIT4            2.253900e-14 0.646136228 0.655 0.404  6.761701e-11       2
#> LSM5             2.310727e-14 0.552029950 0.762 0.429  6.932182e-11       2
#> NBL1             2.512901e-14 0.407424868 0.643 0.360  7.538702e-11       2
#> EDN1             3.490553e-14 0.581238485 0.357 0.125  1.047166e-10       2
#> ACSL11           3.494342e-14 0.549789813 0.518 0.228  1.048303e-10       2
#> RHOV             3.556498e-14 0.545831023 0.435 0.172  1.066949e-10       2
#> ZNF217           3.574343e-14 0.654336944 0.595 0.319  1.072303e-10       2
#> ARHGAP29         5.870874e-14 0.570966177 0.673 0.396  1.761262e-10       2
#> FAM208B          7.489985e-14 0.625646334 0.667 0.366  2.246995e-10       2
#> NEBL             7.927205e-14 0.598767599 0.542 0.269  2.378161e-10       2
#> CLMN             8.513381e-14 0.641733651 0.679 0.404  2.554014e-10       2
#> RP11-783K16.5    9.151979e-14 0.566805685 0.387 0.144  2.745594e-10       2
#> SCGB3A1          1.018985e-13 0.612720974 0.571 0.309  3.056954e-10       2
#> LYPD3            1.108998e-13 0.460310729 0.310 0.096  3.326995e-10       2
#> LDHA             1.136360e-13 0.660788358 0.708 0.459  3.409081e-10       2
#> MMP11            1.208875e-13 0.088627115 0.476 0.217  3.626625e-10       2
#> KCNE5            1.436075e-13 0.370671870 0.417 0.153  4.308225e-10       2
#> LA16c-380H5.5    1.927444e-13 0.486637019 0.470 0.208  5.782333e-10       2
#> RP11-160O5.1     2.196026e-13 0.623302398 0.280 0.087  6.588079e-10       2
#> MCAM             2.579433e-13 0.459126332 0.589 0.329  7.738299e-10       2
#> EBP              2.758611e-13 0.535834443 0.690 0.414  8.275832e-10       2
#> IRX5             3.392699e-13 0.635033116 0.470 0.219  1.017810e-09       2
#> MALL1            4.585508e-13 0.636665396 0.399 0.167  1.375653e-09       2
#> GAPDH            4.694692e-13 0.621411728 0.673 0.482  1.408408e-09       2
#> NAT8L            4.860189e-13 0.571783342 0.452 0.194  1.458057e-09       2
#> LCN21            5.149018e-13 0.464061244 0.625 0.346  1.544705e-09       2
#> DOK5             6.036589e-13 0.364844806 0.333 0.111  1.810977e-09       2
#> PDZD2            6.159988e-13 0.662490934 0.292 0.096  1.847996e-09       2
#> CDH3             7.655471e-13 0.599642291 0.637 0.366  2.296641e-09       2
#> SLC7A5           7.841935e-13 0.545878020 0.542 0.286  2.352581e-09       2
#> ELF31            7.862278e-13 0.587685809 0.714 0.464  2.358683e-09       2
#> OAZ3             9.209788e-13 0.513525860 0.435 0.181  2.762936e-09       2
#> NDRG11           9.381699e-13 0.691857139 0.625 0.419  2.814510e-09       2
#> NFASC            1.045689e-12 0.447922370 0.256 0.072  3.137067e-09       2
#> TMPRSS13         1.144001e-12 0.552365240 0.518 0.245  3.432002e-09       2
#> EIF3B            1.255542e-12 0.598460622 0.726 0.443  3.766627e-09       2
#> MAP3K8           1.332182e-12 0.540441373 0.607 0.318  3.996547e-09       2
#> ERBB31           1.456158e-12 0.610865581 0.661 0.414  4.368473e-09       2
#> SRGAP3           1.589995e-12 0.501421532 0.327 0.114  4.769986e-09       2
#> GLRX1            1.683959e-12 0.641760364 0.667 0.429  5.051878e-09       2
#> KRTCAP3          1.802204e-12 0.571594706 0.708 0.443  5.406613e-09       2
#> STC2             1.987528e-12 0.531779876 0.429 0.182  5.962584e-09       2
#> VGLL1            3.786155e-12 0.509272700 0.589 0.335  1.135847e-08       2
#> SLC25A371        4.922759e-12 0.595515524 0.637 0.406  1.476828e-08       2
#> ADAMTS9          5.545248e-12 0.472072347 0.321 0.114  1.663574e-08       2
#> PPFIA11          6.012983e-12 0.555964793 0.744 0.478  1.803895e-08       2
#> ABLIM1           1.241583e-11 0.470570412 0.274 0.089  3.724750e-08       2
#> IRX3             1.310853e-11 0.578808583 0.655 0.417  3.932558e-08       2
#> CCT6A            1.419243e-11 0.541129881 0.714 0.485  4.257729e-08       2
#> SUSD2            1.431862e-11 0.386782619 0.262 0.082  4.295585e-08       2
#> ST14             1.523417e-11 0.584521171 0.685 0.463  4.570250e-08       2
#> GPRC5A1          1.640789e-11 0.539522449 0.393 0.174  4.922366e-08       2
#> SPON21           1.872858e-11 0.430066667 0.571 0.307  5.618573e-08       2
#> FSTL1            1.970103e-11 0.276679330 0.417 0.186  5.910310e-08       2
#> TM4SF1-AS11      2.064699e-11 0.471240445 0.560 0.294  6.194098e-08       2
#> EGFL7            2.157073e-11 0.457829837 0.304 0.112  6.471220e-08       2
#> BARX1            2.368326e-11 0.521102447 0.649 0.364  7.104977e-08       2
#> SLC9A7           2.369465e-11 0.290058808 0.405 0.156  7.108396e-08       2
#> RASGRP1          2.450807e-11 0.425112245 0.268 0.086  7.352420e-08       2
#> PRSS23           2.788707e-11 0.021047037 0.387 0.170  8.366120e-08       2
#> HLA-C1           2.853572e-11 0.556493560 0.714 0.495  8.560715e-08       2
#> ACTG2            2.910933e-11 0.555326344 0.595 0.358  8.732800e-08       2
#> CTSV             3.273121e-11 0.560179221 0.530 0.298  9.819363e-08       2
#> C3orf38          3.734090e-11 0.672266345 0.542 0.353  1.120227e-07       2
#> KLF13            4.914191e-11 0.456843707 0.494 0.243  1.474257e-07       2
#> LGALSL           8.714281e-11 0.460654544 0.429 0.193  2.614284e-07       2
#> PROX11           9.514281e-11 0.445571365 0.429 0.186  2.854284e-07       2
#> LINC00707        9.736312e-11 0.488510530 0.274 0.096  2.920894e-07       2
#> RAI14            1.001252e-10 0.523717437 0.560 0.334  3.003757e-07       2
#> PROM11           1.047032e-10 0.507023414 0.762 0.450  3.141096e-07       2
#> MTURN            1.076730e-10 0.519576783 0.494 0.252  3.230191e-07       2
#> CLUL1            1.277901e-10 0.511466726 0.339 0.137  3.833703e-07       2
#> YES1             1.332090e-10 0.534566140 0.583 0.328  3.996270e-07       2
#> SCD1             1.392575e-10 0.508875838 0.601 0.375  4.177725e-07       2
#> CENPV            1.491983e-10 0.350352202 0.351 0.138  4.475949e-07       2
#> KLHDC31          1.515733e-10 0.573821030 0.768 0.543  4.547200e-07       2
#> CRNDE            1.520891e-10 0.523898913 0.696 0.479  4.562672e-07       2
#> LOXL1            2.040105e-10 0.240871026 0.405 0.182  6.120316e-07       2
#> ALCAM            2.065219e-10 0.410129363 0.524 0.310  6.195657e-07       2
#> CAPN21           2.763543e-10 0.534870600 0.726 0.491  8.290628e-07       2
#> TMEM154          2.929746e-10 0.458186603 0.262 0.090  8.789238e-07       2
#> MSLN             2.932164e-10 0.511830804 0.387 0.175  8.796492e-07       2
#> S100A8           3.151502e-10 0.593452400 0.470 0.308  9.454505e-07       2
#> TCF7L1           3.210403e-10 0.445049848 0.524 0.264  9.631208e-07       2
#> FKBP10           3.308977e-10 0.519988035 0.649 0.420  9.926930e-07       2
#> CDC25B1          4.627394e-10 0.460469356 0.679 0.423  1.388218e-06       2
#> TMPRSS3          4.999665e-10 0.467331623 0.470 0.229  1.499900e-06       2
#> MIR4435-2HG      5.302647e-10 0.403240998 0.619 0.390  1.590794e-06       2
#> ELF51            5.569844e-10 0.465795061 0.637 0.397  1.670953e-06       2
#> HOPX             7.025927e-10 0.210993899 0.298 0.109  2.107778e-06       2
#> CELF4            7.472321e-10 0.489104763 0.512 0.318  2.241696e-06       2
#> HSPA6            7.936976e-10 0.666654853 0.399 0.196  2.381093e-06       2
#> SLC6A8           8.454919e-10 0.489738260 0.363 0.159  2.536476e-06       2
#> RP11-554I8.2     9.700762e-10 0.421556332 0.476 0.241  2.910229e-06       2
#> SOX41            1.041067e-09 0.542888345 0.798 0.540  3.123200e-06       2
#> SDF2L1           1.090383e-09 0.443045255 0.690 0.439  3.271150e-06       2
#> VTCN11           1.091411e-09 0.463793482 0.696 0.452  3.274234e-06       2
#> CXADR            1.190359e-09 0.481372466 0.577 0.335  3.571077e-06       2
#> MIR210HG         1.541703e-09 0.411254333 0.339 0.141  4.625108e-06       2
#> PCAT19           1.621065e-09 0.404929053 0.310 0.125  4.863196e-06       2
#> THBS2            1.628828e-09 0.313435573 0.458 0.253  4.886483e-06       2
#> PHYH             1.907537e-09 0.499877120 0.673 0.477  5.722612e-06       2
#> SLPI1            2.000450e-09 0.329590878 0.530 0.307  6.001351e-06       2
#> COL6A1           2.174851e-09 0.347987900 0.643 0.416  6.524554e-06       2
#> RAB25            2.261600e-09 0.410505910 0.589 0.366  6.784801e-06       2
#> DBNDD1           2.378552e-09 0.491928518 0.554 0.320  7.135657e-06       2
#> VASN1            2.480502e-09 0.473838158 0.673 0.431  7.441506e-06       2
#> SQRDL1           2.953121e-09 0.464273611 0.595 0.393  8.859363e-06       2
#> CASC15           3.691702e-09 0.389025301 0.250 0.091  1.107511e-05       2
#> UFD1L            3.753040e-09 0.491419755 0.661 0.469  1.125912e-05       2
#> PDP1             4.864258e-09 0.537825944 0.554 0.346  1.459277e-05       2
#> PDGFA            5.877262e-09 0.328808163 0.387 0.175  1.763178e-05       2
#> KLHL35           5.944626e-09 0.451896813 0.548 0.348  1.783388e-05       2
#> INHBB            6.133190e-09 0.512886827 0.310 0.133  1.839957e-05       2
#> PDLIM41          6.237994e-09 0.425279562 0.458 0.229  1.871398e-05       2
#> MESP2            7.922498e-09 0.400568517 0.482 0.256  2.376749e-05       2
#> IFIT1            8.457059e-09 0.212865065 0.292 0.114  2.537118e-05       2
#> NET11            8.745671e-09 0.450053315 0.679 0.448  2.623701e-05       2
#> RP11-400K9.4     9.071183e-09 0.402300947 0.250 0.095  2.721355e-05       2
#> WEE1             9.194811e-09 0.421422175 0.542 0.333  2.758443e-05       2
#> DAPP11           1.063845e-08 0.327107103 0.405 0.191  3.191535e-05       2
#> CLDN1            1.064138e-08 0.282917295 0.345 0.157  3.192415e-05       2
#> CAMK2N1          1.070094e-08 0.390236043 0.304 0.129  3.210283e-05       2
#> MMP71            1.463740e-08 0.369718561 0.577 0.337  4.391220e-05       2
#> RARRES11         1.556493e-08 0.347525306 0.435 0.248  4.669478e-05       2
#> MB1              1.564938e-08 0.411861065 0.464 0.247  4.694814e-05       2
#> MEX3A            1.634760e-08 0.495523281 0.548 0.324  4.904281e-05       2
#> S100A101         2.020097e-08 0.503596538 0.655 0.535  6.060292e-05       2
#> PODXL2           2.138686e-08 0.405997181 0.631 0.433  6.416059e-05       2
#> AC005152.31      2.186742e-08 0.418529022 0.619 0.417  6.560227e-05       2
#> PLXNA2           2.201501e-08 0.405854209 0.476 0.253  6.604503e-05       2
#> NEDD9            2.255078e-08 0.431530736 0.440 0.227  6.765234e-05       2
#> IFITM11          2.696218e-08 0.362772891 0.500 0.297  8.088655e-05       2
#> ART3             3.412548e-08 0.401350088 0.381 0.198  1.023764e-04       2
#> AARD             3.580262e-08 0.392355092 0.631 0.422  1.074078e-04       2
#> OCLN             4.564663e-08 0.394092487 0.470 0.272  1.369399e-04       2
#> CLDN41           4.676845e-08 0.495927699 0.792 0.534  1.403054e-04       2
#> HSP90AB11        4.866473e-08 0.604406035 0.917 0.629  1.459942e-04       2
#> GLYATL2          5.057347e-08 0.358634120 0.429 0.230  1.517204e-04       2
#> SYCE1L           5.890765e-08 0.377105914 0.446 0.235  1.767230e-04       2
#> TEKT31           6.274848e-08 0.447425457 0.690 0.485  1.882454e-04       2
#> TP53BP2          6.533055e-08 0.431973537 0.637 0.410  1.959917e-04       2
#> ORM2             7.410543e-08 0.105531113 0.250 0.095  2.223163e-04       2
#> TINCR1           7.886516e-08 0.334059307 0.518 0.295  2.365955e-04       2
#> CELF21           8.503262e-08 0.314931682 0.488 0.271  2.550979e-04       2
#> FASN             9.439755e-08 0.366325755 0.452 0.248  2.831927e-04       2
#> IFT1721          1.059666e-07 0.329688993 0.512 0.284  3.178998e-04       2
#> CSF3R1           1.158102e-07 0.453051499 0.560 0.388  3.474306e-04       2
#> PSCA             1.183010e-07 0.458931819 0.440 0.245  3.549030e-04       2
#> S100P1           1.391529e-07 0.387155286 0.673 0.433  4.174586e-04       2
#> SMTN1            1.579567e-07 0.413689685 0.613 0.429  4.738700e-04       2
#> SLC2A4RG         1.733568e-07 0.381130110 0.625 0.411  5.200703e-04       2
#> CAPG1            1.791779e-07 0.208955769 0.595 0.360  5.375337e-04       2
#> DNTTIP11         2.434020e-07 0.545342205 0.643 0.479  7.302061e-04       2
#> RANBP1           2.452794e-07 0.349629308 0.702 0.475  7.358383e-04       2
#> HLA-B1           3.793571e-07 0.426926109 0.690 0.466  1.138071e-03       2
#> SELM1            4.091436e-07 0.448536634 0.679 0.516  1.227431e-03       2
#> ISG201           4.615746e-07 0.287332396 0.548 0.370  1.384724e-03       2
#> NCCRP11          4.801674e-07 0.263481297 0.417 0.233  1.440502e-03       2
#> CSTB1            5.152298e-07 0.407581214 0.655 0.460  1.545690e-03       2
#> SMOC1            6.220081e-07 0.282687713 0.452 0.262  1.866024e-03       2
#> PSTPIP2          6.563779e-07 0.282092340 0.512 0.327  1.969134e-03       2
#> NPR31            7.011169e-07 0.362827488 0.393 0.213  2.103351e-03       2
#> ANKRD371         8.178824e-07 0.338548222 0.440 0.254  2.453647e-03       2
#> S100A92          8.930560e-07 0.462336097 0.488 0.340  2.679168e-03       2
#> LINC00152        9.208104e-07 0.338260844 0.571 0.405  2.762431e-03       2
#> MBD2             1.077887e-06 0.386378770 0.613 0.450  3.233662e-03       2
#> CHAF1B           1.163458e-06 0.323735104 0.327 0.161  3.490375e-03       2
#> CHI3L21          1.362248e-06 0.253654362 0.536 0.370  4.086745e-03       2
#> CAMK11           1.470967e-06 0.116709298 0.321 0.146  4.412901e-03       2
#> PLOD21           1.636707e-06 0.378336793 0.607 0.460  4.910121e-03       2
#> MTL5             1.805373e-06 0.291800115 0.583 0.378  5.416120e-03       2
#> EHD1             1.884342e-06 0.256215292 0.458 0.260  5.653025e-03       2
#> RDH101           1.939070e-06 0.374005203 0.524 0.390  5.817211e-03       2
#> PAM1             2.055513e-06 0.379654776 0.625 0.446  6.166539e-03       2
#> HES11            2.209461e-06 0.394584343 0.601 0.431  6.628384e-03       2
#> NR2F21           2.357152e-06 0.334307219 0.637 0.435  7.071455e-03       2
#> LMO41            2.541524e-06 0.355924368 0.488 0.315  7.624572e-03       2
#> MLPH             2.855162e-06 0.216613083 0.280 0.127  8.565487e-03       2
#> PRSS81           3.341306e-06 0.393056165 0.690 0.492  1.002392e-02       2
#> RCAN1            3.646703e-06 0.293784682 0.464 0.337  1.094011e-02       2
#> AIF1L            4.105392e-06 0.365671744 0.595 0.447  1.231618e-02       2
#> TTYH1            4.270043e-06 0.334589405 0.607 0.414  1.281013e-02       2
#> NANOS1           4.524569e-06 0.240765963 0.423 0.243  1.357371e-02       2
#> CTA-293F17.1     4.919646e-06 0.290113696 0.452 0.291  1.475894e-02       2
#> SDC41            5.068865e-06 0.321299287 0.643 0.442  1.520659e-02       2
#> SAC3D1           5.361921e-06 0.170471212 0.464 0.249  1.608576e-02       2
#> FADS2            6.196328e-06 0.340955200 0.440 0.255  1.858898e-02       2
#> KRT181           6.312425e-06 0.402513139 0.720 0.524  1.893728e-02       2
#> ELN              7.798080e-06 0.267168430 0.327 0.173  2.339424e-02       2
#> GCNT1            8.798727e-06 0.346771087 0.470 0.333  2.639618e-02       2
#> BACE21           9.598060e-06 0.360322118 0.607 0.449  2.879418e-02       2
#> THEM6            9.852393e-06 0.299325823 0.458 0.280  2.955718e-02       2
#> TPM11            1.025797e-05 0.413253550 0.744 0.530  3.077391e-02       2
#> ZBTB101          1.108056e-05 0.304202178 0.571 0.419  3.324169e-02       2
#> PMAIP1           1.135196e-05 0.334510318 0.482 0.321  3.405587e-02       2
#> GRB10            1.154205e-05 0.266398608 0.440 0.249  3.462616e-02       2
#> FH               1.200660e-05 0.333636437 0.619 0.447  3.601979e-02       2
#> BNIP3            1.202343e-05 0.372585496 0.583 0.459  3.607030e-02       2
#> MAP2K3           1.277851e-05 0.327895074 0.512 0.376  3.833554e-02       2
#> SNX8             1.514674e-05 0.271061199 0.554 0.369  4.544023e-02       2
#> SCARB1           2.342756e-05 0.376489747 0.589 0.445  7.028267e-02       2
#> TMEM25           2.412298e-05 0.228978552 0.470 0.278  7.236895e-02       2
#> ZG16B            2.527418e-05 0.225923031 0.488 0.348  7.582253e-02       2
#> IRAK1            2.645021e-05 0.286628892 0.565 0.382  7.935063e-02       2
#> ADM              3.355457e-05 0.303557874 0.321 0.194  1.006637e-01       2
#> POLD2            3.771673e-05 0.260056499 0.613 0.423  1.131502e-01       2
#> RPS3             4.076405e-05 0.205630205 0.661 0.440  1.222921e-01       2
#> HMGB3            4.089486e-05 0.294793592 0.643 0.468  1.226846e-01       2
#> CKB              4.420673e-05 0.253298039 0.637 0.423  1.326202e-01       2
#> ALDH1B11         4.759543e-05 0.268550102 0.482 0.332  1.427863e-01       2
#> MESP11           4.826306e-05 0.250578191 0.595 0.431  1.447892e-01       2
#> HLA-A1           5.520301e-05 0.195178773 0.613 0.397  1.656090e-01       2
#> NFIB             5.680807e-05 0.274011552 0.595 0.437  1.704242e-01       2
#> BCL9L            5.794534e-05 0.248299751 0.375 0.217  1.738360e-01       2
#> HES4             6.616312e-05 0.171559102 0.536 0.384  1.984894e-01       2
#> MX11             6.632357e-05 0.026059434 0.268 0.129  1.989707e-01       2
#> COL2A1           7.391285e-05 0.402755273 0.423 0.304  2.217386e-01       2
#> TUBA4A           7.725824e-05 0.230329719 0.464 0.339  2.317747e-01       2
#> RNASET21         7.890726e-05 0.084368279 0.560 0.358  2.367218e-01       2
#> CREB5            9.840216e-05 0.275013666 0.363 0.219  2.952065e-01       2
#> KCTD11           1.045202e-04 0.189676154 0.488 0.349  3.135607e-01       2
#> SMYD2            1.086319e-04 0.299422379 0.595 0.429  3.258956e-01       2
#> GAL              1.128560e-04 0.083575464 0.381 0.266  3.385681e-01       2
#> CENPQ            1.157593e-04 0.235643816 0.506 0.362  3.472780e-01       2
#> KLRG2            1.284616e-04 0.232740958 0.333 0.191  3.853847e-01       2
#> LY6E             1.360836e-04 0.306397211 0.583 0.459  4.082507e-01       2
#> PHLDA3           1.425339e-04 0.305817617 0.476 0.333  4.276016e-01       2
#> KLF101           1.669031e-04 0.275852212 0.536 0.414  5.007093e-01       2
#> ZNF83            1.894518e-04 0.226030073 0.446 0.283  5.683554e-01       2
#> PRR15L           2.062831e-04 0.152416600 0.327 0.184  6.188494e-01       2
#> TNFAIP21         2.315652e-04 0.099757231 0.405 0.265  6.946957e-01       2
#> PTPRS            2.494706e-04 0.286626485 0.321 0.195  7.484119e-01       2
#> TMEM139          2.685819e-04 0.216865184 0.470 0.324  8.057457e-01       2
#> TJP3             2.695699e-04 0.202413956 0.304 0.174  8.087097e-01       2
#> RBBP7            2.785517e-04 0.284663075 0.631 0.503  8.356551e-01       2
#> MBP              2.855605e-04 0.238871511 0.554 0.411  8.566815e-01       2
#> SPEG             3.115378e-04 0.313236184 0.345 0.217  9.346134e-01       2
#> HLA-F1           3.331278e-04 0.155765068 0.470 0.335  9.993834e-01       2
#> TM4SF11          3.489305e-04 0.371085158 0.744 0.571  1.000000e+00       2
#> MGP1             4.457367e-04 0.370917434 0.714 0.526  1.000000e+00       2
#> PHKG1            4.758141e-04 0.245195933 0.351 0.223  1.000000e+00       2
#> SELENBP11        4.802161e-04 0.267979220 0.470 0.330  1.000000e+00       2
#> PPA1             5.730121e-04 0.248272612 0.625 0.479  1.000000e+00       2
#> TRIB11           5.807406e-04 0.208907809 0.542 0.390  1.000000e+00       2
#> ADAM15           5.831394e-04 0.258787890 0.619 0.492  1.000000e+00       2
#> TRIB3            5.945220e-04 0.205567817 0.417 0.300  1.000000e+00       2
#> CP1              6.107266e-04 0.195228306 0.494 0.392  1.000000e+00       2
#> RAB6B            6.269806e-04 0.125067668 0.268 0.149  1.000000e+00       2
#> CMSS1            6.292445e-04 0.120329340 0.446 0.315  1.000000e+00       2
#> PAK1IP1          6.955724e-04 0.231138767 0.536 0.398  1.000000e+00       2
#> MAP7D3           7.899351e-04 0.297044470 0.393 0.276  1.000000e+00       2
#> SPHK1            8.017782e-04 0.258674705 0.530 0.421  1.000000e+00       2
#> SCCPDH1          8.233886e-04 0.155513369 0.446 0.322  1.000000e+00       2
#> PITX1            9.076478e-04 0.183933428 0.571 0.395  1.000000e+00       2
#> MFAP21           1.025644e-03 0.217963666 0.554 0.460  1.000000e+00       2
#> TNFRSF12A        1.393033e-03 0.243604454 0.571 0.422  1.000000e+00       2
#> NFATC11          1.532045e-03 0.183564981 0.494 0.383  1.000000e+00       2
#> RRAGD            1.617360e-03 0.086506381 0.339 0.201  1.000000e+00       2
#> MAFF1            1.769889e-03 0.185943746 0.458 0.357  1.000000e+00       2
#> SOX81            1.874404e-03 0.152462283 0.512 0.386  1.000000e+00       2
#> ADHFE1           2.319700e-03 0.242015742 0.315 0.210  1.000000e+00       2
#> PDIA4            2.379784e-03 0.254954906 0.601 0.495  1.000000e+00       2
#> USP18            2.418933e-03 0.140788741 0.274 0.166  1.000000e+00       2
#> RAB11FIP1        2.441468e-03 0.145047463 0.446 0.307  1.000000e+00       2
#> KRT82            2.597019e-03 0.283276190 0.690 0.532  1.000000e+00       2
#> PTP4A3           2.616342e-03 0.191816521 0.589 0.466  1.000000e+00       2
#> PMP221           2.755021e-03 0.263199032 0.667 0.519  1.000000e+00       2
#> PRSS221          3.077799e-03 0.135153032 0.506 0.372  1.000000e+00       2
#> AQP5             3.240920e-03 0.076868575 0.554 0.417  1.000000e+00       2
#> QDPR             3.249699e-03 0.181547444 0.548 0.426  1.000000e+00       2
#> C1orf116         3.310115e-03 0.226448725 0.357 0.240  1.000000e+00       2
#> PLEKHB1          3.312182e-03 0.213011175 0.613 0.471  1.000000e+00       2
#> RAB30-AS1        3.470268e-03 0.146503676 0.530 0.421  1.000000e+00       2
#> NUCKS1           3.516876e-03 0.201558439 0.631 0.485  1.000000e+00       2
#> NNMT1            3.533425e-03 0.176856766 0.518 0.413  1.000000e+00       2
#> RGS101           4.360631e-03 0.125647558 0.595 0.433  1.000000e+00       2
#> B4GALT1          5.208879e-03 0.159410999 0.506 0.413  1.000000e+00       2
#> C10orf101        5.408858e-03 0.235185821 0.464 0.385  1.000000e+00       2
#> TTF2             5.697773e-03 0.079857960 0.423 0.324  1.000000e+00       2
#> COL4A2           6.094699e-03 0.017238771 0.423 0.334  1.000000e+00       2
#> RMI2             6.218437e-03 0.072945159 0.357 0.235  1.000000e+00       2
#> HES6             6.230367e-03 0.173993302 0.345 0.229  1.000000e+00       2
#> C1orf561         6.310782e-03 0.173134486 0.524 0.432  1.000000e+00       2
#> SERPINH1         6.481283e-03 0.248163513 0.637 0.524  1.000000e+00       2
#> RAB3IP1          6.641383e-03 0.091886953 0.411 0.295  1.000000e+00       2
#> TMEM791          6.718439e-03 0.097485089 0.482 0.383  1.000000e+00       2
#> MGLL1            6.846430e-03 0.053651981 0.357 0.235  1.000000e+00       2
#> PFKP1            7.219781e-03 0.241227754 0.494 0.426  1.000000e+00       2
#> PHLDA2           8.457503e-03 0.106883138 0.399 0.310  1.000000e+00       2
#> GRB141           8.934882e-03 0.243312268 0.310 0.223  1.000000e+00       2
#> TMEM106C1        9.634129e-58 1.550635288 0.959 0.342  2.890239e-54       3
#> MARCKSL11        1.079709e-54 1.231641057 0.973 0.519  3.239128e-51       3
#> AZGP11           1.356838e-52 1.530683259 0.939 0.293  4.070513e-49       3
#> IDH11            1.519149e-50 1.427358054 0.959 0.331  4.557448e-47       3
#> CLPSL1           2.063225e-50 1.544419024 0.797 0.263  6.189676e-47       3
#> IGFBP2           1.193092e-46 1.323517193 0.878 0.342  3.579275e-43       3
#> FRMD3            2.417946e-46 1.331882481 0.662 0.158  7.253838e-43       3
#> PDGFRA1          4.522484e-46 1.326326111 0.824 0.268  1.356745e-42       3
#> SDC21            1.114392e-45 1.327411145 0.926 0.408  3.343177e-42       3
#> PCOLCE21         1.400276e-44 1.210306758 0.777 0.269  4.200827e-41       3
#> AQP51            2.049444e-44 1.419441960 0.865 0.371  6.148333e-41       3
#> PRSS331          2.604867e-44 1.289945123 0.899 0.326  7.814600e-41       3
#> COL2A11          5.875869e-44 1.196525635 0.797 0.248  1.762761e-40       3
#> LEFTY2           6.896814e-43 1.143225646 0.615 0.135  2.069044e-39       3
#> SSRP11           1.325174e-42 1.172785841 0.966 0.428  3.975522e-39       3
#> HIBCH1           1.296761e-41 1.266758316 0.899 0.357  3.890284e-38       3
#> MGST11           1.940918e-41 1.251107270 0.980 0.389  5.822753e-38       3
#> LNX11            1.070256e-40 1.251346222 0.824 0.289  3.210769e-37       3
#> MIA1             1.349339e-40 1.182166087 0.865 0.332  4.048018e-37       3
#> PYCR11           3.822048e-40 1.237691884 0.878 0.398  1.146614e-36       3
#> SYT81            1.557414e-39 1.303856681 0.811 0.311  4.672242e-36       3
#> SCRG11           1.663023e-39 1.196274669 0.878 0.313  4.989068e-36       3
#> SLC9A3R21        2.833326e-39 1.199824656 0.905 0.383  8.499978e-36       3
#> DBI1             8.864107e-39 1.133265687 0.919 0.446  2.659232e-35       3
#> SOHLH11          1.546717e-38 1.151462315 0.777 0.256  4.640150e-35       3
#> S100B1           3.709223e-38 1.176465704 0.851 0.356  1.112767e-34       3
#> TNFRSF211        8.443501e-38 1.172515753 0.865 0.359  2.533050e-34       3
#> FKBP41           1.650638e-37 1.138403428 0.905 0.440  4.951914e-34       3
#> ROPN1B1          1.660683e-37 1.145683879 0.905 0.375  4.982048e-34       3
#> ARL6IP11         3.147403e-37 1.041934176 0.959 0.414  9.442208e-34       3
#> PEG101           3.159497e-37 1.165817322 0.878 0.429  9.478492e-34       3
#> HIST1H2AE1       6.183159e-37 1.177679474 0.878 0.340  1.854948e-33       3
#> TBC1D11          6.205628e-37 1.217804631 0.865 0.399  1.861688e-33       3
#> CYP39A11         2.481615e-36 1.133729134 0.797 0.300  7.444845e-33       3
#> EFHD1            5.565737e-36 1.191051067 0.838 0.369  1.669721e-32       3
#> SERTAD4-AS11     1.770931e-35 1.168599585 0.858 0.383  5.312792e-32       3
#> S100A11          1.829095e-35 1.118355052 0.980 0.443  5.487286e-32       3
#> ACTA21           2.292031e-35 1.118222447 0.784 0.296  6.876093e-32       3
#> CTNND21          3.206089e-35 1.081979289 0.669 0.197  9.618268e-32       3
#> HSP90AB12        5.688443e-35 0.929762274 0.973 0.626  1.706533e-31       3
#> RAC3             5.845250e-35 1.116517766 0.595 0.163  1.753575e-31       3
#> SERTAD41         1.424519e-34 1.111117467 0.899 0.413  4.273556e-31       3
#> COL11A11         1.444166e-34 0.939761196 0.818 0.320  4.332498e-31       3
#> MYL91            2.006192e-34 1.093186899 0.899 0.435  6.018575e-31       3
#> LAPTM4B1         2.867599e-34 1.079926914 0.912 0.499  8.602798e-31       3
#> PBX11            4.396400e-34 1.133959570 0.838 0.415  1.318920e-30       3
#> TPM21            9.899325e-34 0.944236953 0.912 0.364  2.969798e-30       3
#> HAPLN1           1.253971e-33 1.034571661 0.520 0.122  3.761914e-30       3
#> IDH1-AS1         3.953597e-33 1.092223990 0.426 0.081  1.186079e-29       3
#> SNHG251          4.444961e-33 1.161526617 0.777 0.320  1.333488e-29       3
#> TIMM101          1.726112e-32 1.031622009 0.953 0.421  5.178336e-29       3
#> EXTL1            1.989655e-32 1.013680070 0.615 0.178  5.968965e-29       3
#> FBXL22           2.522176e-32 1.129445397 0.534 0.144  7.566529e-29       3
#> TMX21            4.816205e-32 1.010193375 0.919 0.408  1.444861e-28       3
#> GDPD2            6.520669e-32 1.143500226 0.439 0.093  1.956201e-28       3
#> C1QL41           1.238908e-31 0.909985317 0.709 0.226  3.716724e-28       3
#> SYNGR11          1.252781e-31 1.052291450 0.791 0.306  3.758344e-28       3
#> GLS1             1.284130e-31 1.032266056 0.932 0.394  3.852391e-28       3
#> KRT83            1.690957e-31 1.011787420 0.919 0.499  5.072871e-28       3
#> CSRP11           1.324350e-30 1.055915614 0.858 0.399  3.973051e-27       3
#> ACTA1            1.726007e-30 0.976806131 0.446 0.098  5.178022e-27       3
#> CA6              3.488243e-30 1.003107821 0.399 0.077  1.046473e-26       3
#> PLEKHB11         3.847735e-30 1.062615214 0.845 0.438  1.154320e-26       3
#> SFRP11           6.887209e-30 1.020558518 0.824 0.419  2.066163e-26       3
#> CKS1B1           7.600852e-30 0.943697306 0.959 0.447  2.280255e-26       3
#> STMN11           8.220011e-30 0.950807583 0.899 0.477  2.466003e-26       3
#> FAM3C1           2.565241e-29 0.978184292 0.858 0.377  7.695724e-26       3
#> NT5DC2           2.673880e-29 1.064949866 0.797 0.424  8.021639e-26       3
#> WFDC1            4.481317e-29 1.036494631 0.365 0.066  1.344395e-25       3
#> ANGPT1           5.199004e-29 1.061291124 0.500 0.135  1.559701e-25       3
#> FABP7            6.064579e-29 0.944337216 0.264 0.030  1.819374e-25       3
#> DKK1             1.742273e-28 1.059225724 0.459 0.113  5.226819e-25       3
#> PPP1R12A         1.821704e-28 1.103303128 0.804 0.409  5.465111e-25       3
#> CA81             4.120856e-28 0.840731299 0.662 0.217  1.236257e-24       3
#> ITGA101          4.749947e-28 0.859355877 0.581 0.171  1.424984e-24       3
#> LIMCH11          7.308061e-28 1.036625865 0.764 0.351  2.192418e-24       3
#> PLOD31           7.399414e-28 1.075954118 0.757 0.380  2.219824e-24       3
#> HRCT11           9.515855e-28 0.992409488 0.777 0.340  2.854757e-24       3
#> CPED11           1.416450e-27 0.957303315 0.649 0.230  4.249350e-24       3
#> ROPN11           2.137284e-27 0.963488134 0.723 0.273  6.411852e-24       3
#> NFIB1            2.895773e-27 1.067503299 0.777 0.412  8.687319e-24       3
#> C16orf74         7.680675e-27 1.082966278 0.534 0.174  2.304202e-23       3
#> HIST1H2BG1       1.541187e-26 0.950098786 0.797 0.352  4.623560e-23       3
#> TAGLN1           2.047407e-26 0.955235507 0.811 0.391  6.142220e-23       3
#> RP11-357H14.171  2.984965e-26 0.838013488 0.547 0.160  8.954894e-23       3
#> SPDEF            3.040827e-26 0.912192374 0.284 0.041  9.122481e-23       3
#> PHGDH            6.198332e-26 0.894348816 0.770 0.354  1.859500e-22       3
#> HACD11           7.120027e-26 0.935973192 0.777 0.359  2.136008e-22       3
#> PDIA41           1.605472e-25 0.958115516 0.845 0.459  4.816416e-22       3
#> SNHG19           1.663089e-25 1.018317371 0.743 0.371  4.989267e-22       3
#> NEGR11           2.110977e-25 0.931698359 0.534 0.165  6.332931e-22       3
#> CRABP1           2.980245e-25 0.930000204 0.784 0.354  8.940736e-22       3
#> SH3BGR1          3.253879e-25 0.909917040 0.709 0.292  9.761636e-22       3
#> MFAP22           4.420004e-25 0.937760255 0.818 0.420  1.326001e-21       3
#> PPP1R1B1         5.696404e-25 0.924915698 0.750 0.337  1.708921e-21       3
#> MYOZ1            6.774488e-25 0.886304654 0.615 0.224  2.032346e-21       3
#> ENPP51           1.035033e-24 0.907436735 0.730 0.293  3.105100e-21       3
#> LDHB             1.182194e-24 0.885016229 0.824 0.450  3.546581e-21       3
#> DTNB1            1.994958e-24 0.940145503 0.723 0.312  5.984874e-21       3
#> TUBB2B1          4.167498e-24 0.807800522 0.716 0.287  1.250250e-20       3
#> ITGA6            6.896855e-24 0.943388073 0.649 0.246  2.069057e-20       3
#> IGFBP5           7.882221e-24 0.994844492 0.581 0.229  2.364666e-20       3
#> AIF1L1           1.168752e-23 0.903627358 0.791 0.419  3.506256e-20       3
#> DSC31            1.388538e-23 0.766017045 0.723 0.268  4.165613e-20       3
#> RGCC             2.231547e-23 0.898536038 0.736 0.344  6.694641e-20       3
#> CACYBP1          3.062446e-23 0.860061723 0.878 0.485  9.187337e-20       3
#> NDRG21           3.207538e-23 0.868452972 0.797 0.398  9.622615e-20       3
#> SERPINH11        3.924576e-23 0.846932557 0.885 0.488  1.177373e-19       3
#> AP1M21           5.460651e-23 0.890524846 0.818 0.415  1.638195e-19       3
#> MATN31           7.675552e-23 0.704017419 0.655 0.272  2.302666e-19       3
#> NME1             9.005795e-23 0.888340599 0.791 0.454  2.701738e-19       3
#> TTYH11           9.586965e-23 0.897206611 0.730 0.399  2.876090e-19       3
#> MT1E             1.072084e-22 0.787389023 0.750 0.342  3.216253e-19       3
#> MYO1B1           1.126041e-22 0.772378058 0.770 0.360  3.378122e-19       3
#> SUN31            1.984938e-22 0.828471188 0.541 0.178  5.954814e-19       3
#> DNM3OS1          2.444730e-22 0.876898291 0.574 0.203  7.334190e-19       3
#> PMP222           2.648783e-22 0.822670118 0.858 0.492  7.946349e-19       3
#> ANO1             3.216497e-22 0.780862824 0.426 0.116  9.649490e-19       3
#> MYLK1            4.430509e-22 0.827414419 0.791 0.391  1.329153e-18       3
#> HSPA1A1          4.627626e-22 0.853137866 0.899 0.446  1.388288e-18       3
#> CTGF1            7.507732e-22 0.790056830 0.696 0.329  2.252320e-18       3
#> VDR              8.150858e-22 0.802695086 0.432 0.121  2.445258e-18       3
#> TUBB2A1          9.870322e-22 0.836388542 0.797 0.376  2.961097e-18       3
#> CKB1             2.761010e-21 0.759898448 0.811 0.400  8.283029e-18       3
#> WIF1             3.813066e-21 0.775924871 0.365 0.093  1.143920e-17       3
#> SLC29A11         4.718797e-21 0.849338428 0.912 0.457  1.415639e-17       3
#> HSPD1            6.828827e-21 0.788440300 0.791 0.451  2.048648e-17       3
#> DNPH1            8.573118e-21 0.754620795 0.865 0.488  2.571935e-17       3
#> ACTL8            1.860656e-20 0.686228224 0.385 0.102  5.581969e-17       3
#> C1orf115         2.188905e-20 0.879833862 0.588 0.240  6.566714e-17       3
#> CMBL             3.780504e-20 0.837645534 0.466 0.159  1.134151e-16       3
#> FREM2            4.162222e-20 0.807238856 0.284 0.057  1.248667e-16       3
#> PODXL21          6.798432e-20 0.856133893 0.743 0.419  2.039530e-16       3
#> CLU1             1.026682e-19 0.800621195 0.946 0.398  3.080047e-16       3
#> NANOS11          1.318577e-19 0.833452311 0.561 0.226  3.955732e-16       3
#> HIST2H2BE1       1.470355e-19 0.707216677 0.811 0.368  4.411065e-16       3
#> PTGIS            2.448485e-19 0.658932272 0.405 0.114  7.345455e-16       3
#> CPM              2.668866e-19 0.869477123 0.696 0.368  8.006599e-16       3
#> MSRB31           2.791765e-19 0.628201286 0.554 0.207  8.375294e-16       3
#> SERPINE21        4.630361e-19 0.705659849 0.723 0.355  1.389108e-15       3
#> CCT51            4.978230e-19 0.752053618 0.818 0.488  1.493469e-15       3
#> EPHX11           5.561209e-19 0.776661404 0.838 0.414  1.668363e-15       3
#> HMGB2            5.977352e-19 0.546048227 0.784 0.350  1.793205e-15       3
#> KRT182           6.312043e-19 0.806641016 0.872 0.505  1.893613e-15       3
#> TSPAN21          7.564337e-19 0.820892535 0.439 0.148  2.269301e-15       3
#> NES              7.911876e-19 0.739140335 0.682 0.296  2.373563e-15       3
#> SYCP2            8.001881e-19 0.773312637 0.534 0.209  2.400564e-15       3
#> NQO11            1.126705e-18 0.803291008 0.736 0.398  3.380114e-15       3
#> PFN2             1.691971e-18 0.776683598 0.811 0.466  5.075912e-15       3
#> VANGL11          1.699441e-18 0.755900349 0.824 0.492  5.098324e-15       3
#> BAMBI1           4.116109e-18 0.726914400 0.797 0.415  1.234833e-14       3
#> LAMB31           4.227691e-18 0.641571886 0.547 0.207  1.268307e-14       3
#> SERPINA51        4.577410e-18 0.808570417 0.642 0.296  1.373223e-14       3
#> CTHRC11          6.043015e-18 0.717973468 0.777 0.430  1.812904e-14       3
#> GOLM11           7.431732e-18 0.771624548 0.770 0.445  2.229519e-14       3
#> IDI11            9.621618e-18 0.723200640 0.770 0.417  2.886485e-14       3
#> SOD31            2.403609e-17 0.704881926 0.784 0.362  7.210827e-14       3
#> IGFBP71          3.631217e-17 0.645399420 0.709 0.375  1.089365e-13       3
#> COL11A21         4.278491e-17 0.751264678 0.507 0.187  1.283547e-13       3
#> AARD1            5.557541e-17 0.689012080 0.791 0.401  1.667262e-13       3
#> PRPS21           6.619645e-17 0.798827939 0.696 0.355  1.985893e-13       3
#> HSPA1B1          1.237259e-16 0.698462650 0.824 0.415  3.711776e-13       3
#> TSPAN121         1.327649e-16 0.568893499 0.669 0.300  3.982948e-13       3
#> PITX11           1.505650e-16 0.722154504 0.716 0.376  4.516951e-13       3
#> CRISPLD11        1.839346e-16 0.696847488 0.736 0.412  5.518037e-13       3
#> GTF3A            1.866161e-16 0.726582235 0.736 0.470  5.598484e-13       3
#> HMGN2            2.023442e-16 0.627650757 0.784 0.440  6.070326e-13       3
#> LRRC73           2.435949e-16 0.757975694 0.486 0.196  7.307847e-13       3
#> METTL7A1         2.752637e-16 0.548630395 0.669 0.325  8.257911e-13       3
#> PRR7             4.472103e-16 0.799821319 0.466 0.190  1.341631e-12       3
#> HIST1H1C1        4.531342e-16 0.681344343 0.736 0.407  1.359403e-12       3
#> HMGB11           4.934239e-16 0.647776393 0.797 0.454  1.480272e-12       3
#> XAGE2            8.561231e-16 0.567698014 0.608 0.264  2.568369e-12       3
#> GOLT1A           1.345217e-15 0.747432587 0.331 0.100  4.035652e-12       3
#> RP11-89K21.1     1.393888e-15 0.732689779 0.399 0.140  4.181665e-12       3
#> NCAN             2.214706e-15 0.636663348 0.250 0.058  6.644119e-12       3
#> TUBB1            2.330243e-15 0.648261046 0.824 0.485  6.990728e-12       3
#> RP11-25K19.11    2.414551e-15 0.569735378 0.466 0.177  7.243654e-12       3
#> CITED41          3.196570e-15 0.667135012 0.764 0.445  9.589711e-12       3
#> TUBA1B           3.323929e-15 0.606578214 0.797 0.448  9.971786e-12       3
#> ST141            4.282363e-15 0.701793219 0.770 0.454  1.284709e-11       3
#> GPM6B1           4.598248e-15 0.623931941 0.750 0.398  1.379474e-11       3
#> MDFI1            5.039475e-15 0.624027295 0.811 0.406  1.511843e-11       3
#> DCXR             9.201032e-15 0.665236652 0.736 0.407  2.760310e-11       3
#> RAMP2            9.338869e-15 0.536986998 0.520 0.221  2.801661e-11       3
#> PCBD1            1.007942e-14 0.722556851 0.716 0.479  3.023826e-11       3
#> KCNQ1OT12        1.068721e-14 0.654580570 0.730 0.459  3.206164e-11       3
#> TEKT32           1.210941e-14 0.683317018 0.770 0.477  3.632823e-11       3
#> CKS21            1.231403e-14 0.500196890 0.770 0.399  3.694208e-11       3
#> ERVMER34-1       1.619669e-14 0.600620022 0.486 0.198  4.859008e-11       3
#> SMOC11           1.829649e-14 0.673178932 0.547 0.251  5.488946e-11       3
#> NTHL1            1.982066e-14 0.723380735 0.669 0.385  5.946199e-11       3
#> PRSS82           3.174537e-14 0.676671488 0.784 0.482  9.523611e-11       3
#> KLK11            4.056030e-14 0.744588428 0.284 0.081  1.216809e-10       3
#> FXYD61           4.187404e-14 0.651575831 0.791 0.452  1.256221e-10       3
#> TUBA1A1          4.284415e-14 0.658103924 0.723 0.453  1.285324e-10       3
#> TMEM204          5.350733e-14 0.473091142 0.378 0.125  1.605220e-10       3
#> VGF              6.968297e-14 0.683323931 0.385 0.144  2.090489e-10       3
#> TNFSF13B1        9.258129e-14 0.674473127 0.730 0.405  2.777439e-10       3
#> TPM12            1.007843e-13 0.647865998 0.818 0.523  3.023530e-10       3
#> SQLE1            1.455228e-13 0.654170553 0.696 0.409  4.365683e-10       3
#> PHYH1            1.507771e-13 0.670937059 0.764 0.467  4.523313e-10       3
#> HEY21            1.762905e-13 0.589630555 0.473 0.195  5.288715e-10       3
#> COPZ2            2.081578e-13 0.518638072 0.480 0.198  6.244734e-10       3
#> PLK2             4.672363e-13 0.641506445 0.635 0.326  1.401709e-09       3
#> TOX1             8.202804e-13 0.538764219 0.527 0.231  2.460841e-09       3
#> SLBP1            1.418614e-12 0.565763957 0.743 0.466  4.255843e-09       3
#> DBNDD11          1.486123e-12 0.560854040 0.628 0.313  4.458369e-09       3
#> FBXO321          1.646435e-12 0.593338771 0.743 0.430  4.939306e-09       3
#> IL17B1           1.687359e-12 0.402090595 0.412 0.156  5.062078e-09       3
#> TOB11            1.833974e-12 0.537613538 0.764 0.418  5.501922e-09       3
#> ITIH6            1.951013e-12 0.611332995 0.291 0.092  5.853040e-09       3
#> KCNMB11          2.188432e-12 0.564359634 0.628 0.310  6.565297e-09       3
#> HN1              2.209199e-12 0.564782737 0.791 0.462  6.627596e-09       3
#> RASL12           2.460165e-12 0.541557803 0.304 0.096  7.380495e-09       3
#> RUVBL1           3.200456e-12 0.586687727 0.689 0.407  9.601369e-09       3
#> SUSD5            3.358730e-12 0.503393055 0.311 0.101  1.007619e-08       3
#> PTS1             4.177274e-12 0.609557907 0.743 0.446  1.253182e-08       3
#> COL9A32          4.215417e-12 0.555495658 0.682 0.417  1.264625e-08       3
#> MAPK13           4.688226e-12 0.611844620 0.750 0.464  1.406468e-08       3
#> PPIL1            6.014055e-12 0.604921405 0.635 0.341  1.804217e-08       3
#> P3H4             6.176799e-12 0.633728152 0.642 0.355  1.853040e-08       3
#> ELOVL5           6.370017e-12 0.541406861 0.520 0.241  1.911005e-08       3
#> LTBP1            6.661253e-12 0.490072739 0.399 0.157  1.998376e-08       3
#> SNX221           9.505422e-12 0.578228667 0.574 0.281  2.851627e-08       3
#> GAS1             1.210264e-11 0.641370030 0.581 0.301  3.630791e-08       3
#> ACTL6A1          1.442245e-11 0.576534953 0.696 0.462  4.326736e-08       3
#> CTNNB11          1.512781e-11 0.594523682 0.723 0.451  4.538344e-08       3
#> LRRCC1           1.741310e-11 0.540953682 0.662 0.380  5.223931e-08       3
#> PAQR4            1.929309e-11 0.580696358 0.595 0.300  5.787927e-08       3
#> VSNL1            2.048945e-11 0.516038144 0.291 0.096  6.146836e-08       3
#> GAMT             2.390239e-11 0.519114422 0.473 0.214  7.170717e-08       3
#> PCP4             2.430338e-11 0.363813275 0.264 0.078  7.291013e-08       3
#> TMSB15A1         2.500310e-11 0.573559859 0.399 0.169  7.500931e-08       3
#> FHL11            2.647376e-11 0.503316880 0.412 0.173  7.942128e-08       3
#> IMPA21           2.943407e-11 0.563850366 0.736 0.476  8.830220e-08       3
#> SNAI21           4.164089e-11 0.484619259 0.561 0.309  1.249227e-07       3
#> PRSS222          4.305204e-11 0.534093537 0.635 0.355  1.291561e-07       3
#> C6orf141         4.328323e-11 0.580368541 0.331 0.126  1.298497e-07       3
#> NET12            4.428409e-11 0.568813473 0.696 0.450  1.328523e-07       3
#> IFRD11           5.711265e-11 0.485599249 0.750 0.397  1.713379e-07       3
#> DNM31            5.774009e-11 0.550781732 0.500 0.243  1.732203e-07       3
#> TPD52L11         5.929429e-11 0.570239304 0.696 0.443  1.778829e-07       3
#> SCGB3A11         6.765000e-11 0.554938366 0.588 0.312  2.029500e-07       3
#> MLLT111          8.149809e-11 0.539950859 0.568 0.295  2.444943e-07       3
#> NUDT41           8.362373e-11 0.613876433 0.628 0.358  2.508712e-07       3
#> FRMD4A1          1.068372e-10 0.412686027 0.669 0.326  3.205115e-07       3
#> CRELD2           1.775894e-10 0.560608602 0.703 0.471  5.327683e-07       3
#> ADAM151          2.297256e-10 0.536787784 0.716 0.479  6.891768e-07       3
#> SNRNP25          2.345752e-10 0.536419920 0.696 0.427  7.037256e-07       3
#> HSPH11           2.509654e-10 0.524465353 0.709 0.468  7.528963e-07       3
#> BOP1             2.516461e-10 0.563415520 0.689 0.423  7.549382e-07       3
#> VASN2            2.527033e-10 0.566940453 0.642 0.440  7.581098e-07       3
#> DCTPP1           2.707600e-10 0.545799628 0.709 0.473  8.122800e-07       3
#> B3GNT71          2.727255e-10 0.571778957 0.574 0.311  8.181765e-07       3
#> PALLD            4.660873e-10 0.352284673 0.608 0.338  1.398262e-06       3
#> SENCR            5.040159e-10 0.526402239 0.257 0.086  1.512048e-06       3
#> LIPH             5.212538e-10 0.480301698 0.486 0.238  1.563761e-06       3
#> SLC43A31         5.400348e-10 0.481127082 0.736 0.417  1.620104e-06       3
#> HIST1H4E1        5.626588e-10 0.469609213 0.378 0.163  1.687976e-06       3
#> ELN1             7.241550e-10 0.391366474 0.385 0.168  2.172465e-06       3
#> KLHDC32          8.778918e-10 0.553057020 0.838 0.536  2.633676e-06       3
#> SLC2A4RG1        9.435410e-10 0.512136892 0.662 0.410  2.830623e-06       3
#> NPPC             1.064267e-09 0.502771038 0.250 0.085  3.192801e-06       3
#> PTRF1            1.082535e-09 0.344779913 0.568 0.302  3.247605e-06       3
#> CP2              1.335763e-09 0.434731182 0.595 0.378  4.007289e-06       3
#> NKD21            1.423939e-09 0.350193493 0.372 0.155  4.271817e-06       3
#> PPP1R14A         1.651659e-09 0.378039087 0.399 0.175  4.954978e-06       3
#> NREP             1.689578e-09 0.346084274 0.486 0.230  5.068733e-06       3
#> NSG11            2.090031e-09 0.476020147 0.372 0.163  6.270093e-06       3
#> PLOD22           2.221115e-09 0.501177813 0.696 0.449  6.663346e-06       3
#> CYR61            2.556749e-09 0.412505307 0.635 0.395  7.670247e-06       3
#> RP11-798M19.6    2.640609e-09 0.457139836 0.480 0.233  7.921828e-06       3
#> ODC1             2.920398e-09 0.532957398 0.635 0.361  8.761194e-06       3
#> SERPINB52        3.300196e-09 0.345285820 0.669 0.399  9.900589e-06       3
#> TMEM1391         4.774943e-09 0.491209429 0.588 0.309  1.432483e-05       3
#> SOX82            5.834425e-09 0.447748041 0.595 0.376  1.750327e-05       3
#> CRNDE1           6.534255e-09 0.491642304 0.730 0.478  1.960276e-05       3
#> DLX51            6.713408e-09 0.258029166 0.426 0.187  2.014023e-05       3
#> RAD211           7.662595e-09 0.468779471 0.770 0.505  2.298779e-05       3
#> POSTN1           9.088833e-09 0.292508662 0.432 0.219  2.726650e-05       3
#> FBLN1            9.163791e-09 0.360294649 0.493 0.239  2.749137e-05       3
#> ATP6V0E2         9.805443e-09 0.538916412 0.507 0.294  2.941633e-05       3
#> UCHL1            1.189880e-08 0.575485171 0.331 0.146  3.569639e-05       3
#> LRP2             1.267645e-08 0.458002419 0.439 0.222  3.802934e-05       3
#> C9orf40          1.357472e-08 0.474143480 0.446 0.233  4.072416e-05       3
#> BGN1             1.492935e-08 0.451342990 0.709 0.459  4.478805e-05       3
#> CA2              1.580076e-08 0.408636572 0.439 0.217  4.740229e-05       3
#> SFN1             1.799459e-08 0.292638127 0.601 0.346  5.398377e-05       3
#> MEST             2.172303e-08 0.465796364 0.284 0.114  6.516908e-05       3
#> MT1G1            2.342373e-08 0.091416473 0.419 0.209  7.027118e-05       3
#> ZG16B1           2.537455e-08 0.409743571 0.514 0.347  7.612366e-05       3
#> GJA1             2.639703e-08 0.342242115 0.365 0.161  7.919109e-05       3
#> CALD11           3.254764e-08 0.304025719 0.615 0.395  9.764292e-05       3
#> PTP4A31          3.484664e-08 0.390521309 0.709 0.450  1.045399e-04       3
#> PRKDC            3.594281e-08 0.460037800 0.689 0.487  1.078284e-04       3
#> HSPB11           4.250485e-08 0.411809913 0.682 0.473  1.275146e-04       3
#> SORBS22          4.927263e-08 0.452213531 0.743 0.457  1.478179e-04       3
#> EPHX21           6.048450e-08 0.335034850 0.473 0.238  1.814535e-04       3
#> NEXN             6.677565e-08 0.242868761 0.264 0.099  2.003269e-04       3
#> MT2A2            6.981193e-08 0.351960760 0.615 0.380  2.094358e-04       3
#> MAPK8IP21        8.556816e-08 0.443609880 0.392 0.197  2.567045e-04       3
#> PRNP             8.824841e-08 0.417681225 0.703 0.502  2.647452e-04       3
#> KRT72            9.131865e-08 0.596538296 0.872 0.583  2.739559e-04       3
#> FAM89A1          9.738908e-08 0.383939267 0.561 0.311  2.921672e-04       3
#> SLC43A1          1.069447e-07 0.354518431 0.304 0.130  3.208340e-04       3
#> CD320            1.076630e-07 0.438207656 0.500 0.276  3.229891e-04       3
#> CDH31            1.198841e-07 0.341272649 0.608 0.376  3.596523e-04       3
#> TMEM67           1.323019e-07 0.420078450 0.399 0.209  3.969057e-04       3
#> TNC              1.333932e-07 0.374484457 0.351 0.164  4.001796e-04       3
#> QPCT1            1.334106e-07 0.341365307 0.541 0.338  4.002318e-04       3
#> LY6E1            1.490848e-07 0.446625817 0.676 0.447  4.472544e-04       3
#> QPRT             1.715358e-07 0.286108600 0.500 0.252  5.146074e-04       3
#> RTN4RL21         1.796294e-07 0.456680346 0.351 0.172  5.388882e-04       3
#> UNG              1.907296e-07 0.370465608 0.459 0.247  5.721889e-04       3
#> FAM46B1          1.981969e-07 0.394384234 0.473 0.262  5.945907e-04       3
#> SLC26A7          2.297953e-07 0.306255182 0.270 0.110  6.893858e-04       3
#> SMC1B            2.378330e-07 0.385949200 0.257 0.105  7.134990e-04       3
#> HYLS1            2.536407e-07 0.350490962 0.331 0.153  7.609220e-04       3
#> RP11-19E11.11    3.157858e-07 0.417777753 0.635 0.394  9.473573e-04       3
#> HIST1H2BJ        4.044340e-07 0.373217600 0.372 0.189  1.213302e-03       3
#> CDC42EP11        4.159889e-07 0.400603250 0.608 0.436  1.247967e-03       3
#> PRELP1           4.166471e-07 0.374095836 0.399 0.213  1.249941e-03       3
#> TUBA4A1          4.207523e-07 0.377138841 0.581 0.323  1.262257e-03       3
#> GCNT11           4.343471e-07 0.351460451 0.568 0.320  1.303041e-03       3
#> SEPT4            4.891463e-07 0.220676315 0.264 0.110  1.467439e-03       3
#> TMEM1581         5.487727e-07 0.389911296 0.561 0.326  1.646318e-03       3
#> C12orf75         5.529430e-07 0.421269319 0.466 0.255  1.658829e-03       3
#> GGH              5.539524e-07 0.443545491 0.534 0.377  1.661857e-03       3
#> SCIN1            5.812880e-07 0.431930042 0.331 0.161  1.743864e-03       3
#> ASNS             5.829695e-07 0.409278122 0.405 0.223  1.748908e-03       3
#> PPDPF1           6.992031e-07 0.464170979 0.865 0.533  2.097609e-03       3
#> MT1F1            7.902930e-07 0.140315366 0.520 0.297  2.370879e-03       3
#> UACA             8.510173e-07 0.166968044 0.547 0.330  2.553052e-03       3
#> SBSPON1          9.583968e-07 0.325679356 0.473 0.295  2.875190e-03       3
#> HIST1H2BN        1.009659e-06 0.355783337 0.351 0.177  3.028978e-03       3
#> GAS61            1.020742e-06 0.360617924 0.622 0.416  3.062226e-03       3
#> HSPA51           1.304204e-06 0.444037050 0.669 0.488  3.912612e-03       3
#> SCARB11          1.385040e-06 0.424965709 0.615 0.444  4.155119e-03       3
#> C1orf1161        1.412156e-06 0.348681674 0.426 0.232  4.236469e-03       3
#> CELF41           1.415945e-06 0.371229960 0.534 0.318  4.247835e-03       3
#> KRT151           1.461592e-06 0.432123235 0.500 0.277  4.384776e-03       3
#> ELF52            1.540180e-06 0.351576715 0.622 0.405  4.620539e-03       3
#> RP1-313I6.12     1.553363e-06 0.363391446 0.426 0.243  4.660088e-03       3
#> TMEM611          1.942629e-06 0.317178816 0.568 0.341  5.827886e-03       3
#> BYSL             2.041597e-06 0.405584022 0.608 0.407  6.124792e-03       3
#> YBX2             2.061911e-06 0.445727767 0.270 0.130  6.185733e-03       3
#> SKA2             2.136945e-06 0.264683596 0.514 0.280  6.410835e-03       3
#> CDCA7L           2.189250e-06 0.338845239 0.466 0.264  6.567751e-03       3
#> TSPAN51          2.216630e-06 0.334732987 0.473 0.259  6.649890e-03       3
#> TMEM97           2.561569e-06 0.351459143 0.419 0.240  7.684706e-03       3
#> C2orf822         2.659424e-06 0.431920945 0.750 0.470  7.978273e-03       3
#> MBNL1-AS11       2.951941e-06 0.355986908 0.500 0.287  8.855823e-03       3
#> DGAT2            3.001554e-06 0.334104449 0.507 0.300  9.004662e-03       3
#> SIGMAR1          3.362065e-06 0.372242601 0.601 0.359  1.008619e-02       3
#> ETV51            3.554807e-06 0.235100725 0.351 0.169  1.066442e-02       3
#> GMPR1            3.700479e-06 0.287769767 0.453 0.249  1.110144e-02       3
#> TUBB4B1          4.223607e-06 0.378543641 0.649 0.496  1.267082e-02       3
#> TUBA1C1          4.422588e-06 0.320266994 0.628 0.449  1.326776e-02       3
#> BOC1             4.450746e-06 0.319649243 0.412 0.228  1.335224e-02       3
#> BARX11           4.541080e-06 0.350780166 0.588 0.379  1.362324e-02       3
#> PAICS            4.604955e-06 0.370991496 0.622 0.401  1.381486e-02       3
#> MESP12           4.771439e-06 0.332472096 0.615 0.431  1.431432e-02       3
#> CCND1            5.012529e-06 0.360814428 0.662 0.450  1.503759e-02       3
#> SAP30            5.544314e-06 0.360876403 0.628 0.415  1.663294e-02       3
#> MCM7             5.820942e-06 0.268045812 0.574 0.364  1.746283e-02       3
#> LOXL2            7.121446e-06 0.232189692 0.331 0.171  2.136434e-02       3
#> GADD45G          8.659685e-06 0.374127585 0.453 0.266  2.597906e-02       3
#> TMEM792          8.886719e-06 0.252617613 0.588 0.369  2.666016e-02       3
#> TSEN15           9.632384e-06 0.360488358 0.635 0.462  2.889715e-02       3
#> S100A102         9.739002e-06 0.414273692 0.723 0.527  2.921700e-02       3
#> MINCR            9.869568e-06 0.307694840 0.453 0.269  2.960870e-02       3
#> CTNNAL1          1.297330e-05 0.330407046 0.432 0.249  3.891989e-02       3
#> NUCKS11          1.441449e-05 0.331097753 0.669 0.483  4.324347e-02       3
#> NRM              1.763972e-05 0.350392140 0.358 0.205  5.291916e-02       3
#> DHTKD1           2.002350e-05 0.341530284 0.615 0.400  6.007051e-02       3
#> MTL51            2.914541e-05 0.308430472 0.554 0.387  8.743623e-02       3
#> ANXA2R1          4.708892e-05 0.309643941 0.405 0.241  1.412667e-01       3
#> CHI3L11          6.442235e-05 0.342067422 0.615 0.482  1.932671e-01       3
#> CNTNAP3B1        7.073159e-05 0.177284264 0.446 0.254  2.121948e-01       3
#> RRS1             7.163894e-05 0.352694861 0.514 0.354  2.149168e-01       3
#> FBLN22           7.383702e-05 0.345987058 0.649 0.512  2.215111e-01       3
#> TRAF3IP31        7.778741e-05 0.269241024 0.466 0.280  2.333622e-01       3
#> FKBP101          8.305121e-05 0.277064713 0.581 0.435  2.491536e-01       3
#> TMEM45A          8.313931e-05 0.150666428 0.338 0.180  2.494179e-01       3
#> SMOC2            8.593458e-05 0.139222205 0.372 0.208  2.578037e-01       3
#> MT1X1            8.839889e-05 0.136793164 0.527 0.380  2.651967e-01       3
#> KCTD13           9.941291e-05 0.338746248 0.527 0.346  2.982387e-01       3
#> DKC1             1.036218e-04 0.301100212 0.622 0.466  3.108655e-01       3
#> ID4              1.071489e-04 0.334662530 0.419 0.249  3.214466e-01       3
#> CEP70            1.697574e-04 0.210745813 0.385 0.235  5.092722e-01       3
#> PEG32            1.781365e-04 0.278462036 0.561 0.425  5.344094e-01       3
#> BACE22           1.962761e-04 0.274140569 0.588 0.455  5.888283e-01       3
#> H2AFV            2.058721e-04 0.271080825 0.622 0.477  6.176162e-01       3
#> PDZK1IP11        2.097859e-04 0.224908471 0.622 0.446  6.293576e-01       3
#> RAB251           2.267428e-04 0.274711492 0.527 0.380  6.802284e-01       3
#> ALYREF           2.320906e-04 0.233344561 0.520 0.332  6.962717e-01       3
#> VTCN12           2.350106e-04 0.274533198 0.655 0.464  7.050317e-01       3
#> GCSH1            2.452445e-04 0.235213746 0.588 0.426  7.357334e-01       3
#> THBS11           2.476496e-04 0.134985526 0.514 0.328  7.429488e-01       3
#> KRT861           2.557353e-04 0.177036830 0.412 0.248  7.672060e-01       3
#> ACTG21           2.718299e-04 0.202058763 0.541 0.372  8.154896e-01       3
#> TINAGL1          3.199513e-04 0.334755014 0.250 0.141  9.598539e-01       3
#> PLA2G4A1         3.281858e-04 0.129606914 0.453 0.275  9.845575e-01       3
#> SDC42            3.378695e-04 0.212844053 0.615 0.451  1.000000e+00       3
#> TRIB31           3.401950e-04 0.242107069 0.459 0.296  1.000000e+00       3
#> GPSM2            3.473231e-04 0.261852962 0.304 0.175  1.000000e+00       3
#> FRZB             3.509361e-04 0.414766253 0.392 0.241  1.000000e+00       3
#> KLRG21           3.593672e-04 0.252554154 0.324 0.195  1.000000e+00       3
#> C8orf461         3.734268e-04 0.208061624 0.473 0.287  1.000000e+00       3
#> KRT811           3.849369e-04 0.326258149 0.669 0.507  1.000000e+00       3
#> KRTCAP31         3.913622e-04 0.270644934 0.628 0.462  1.000000e+00       3
#> BARD11           3.924867e-04 0.213622609 0.500 0.305  1.000000e+00       3
#> TEX30            4.191586e-04 0.198783853 0.507 0.314  1.000000e+00       3
#> MEX3A1           4.873443e-04 0.242370382 0.486 0.338  1.000000e+00       3
#> TUBB6            5.293968e-04 0.259760034 0.547 0.466  1.000000e+00       3
#> CDK41            5.807671e-04 0.287843081 0.649 0.508  1.000000e+00       3
#> CDH11            6.712922e-04 0.136968710 0.622 0.423  1.000000e+00       3
#> MYBL11           8.613839e-04 0.205749817 0.392 0.249  1.000000e+00       3
#> COL6A11          1.047072e-03 0.130893972 0.554 0.434  1.000000e+00       3
#> CERCAM           1.083086e-03 0.139036269 0.358 0.229  1.000000e+00       3
#> RUVBL2           1.092829e-03 0.218105159 0.568 0.385  1.000000e+00       3
#> A2M1             1.099084e-03 0.289754893 0.581 0.475  1.000000e+00       3
#> HMGA11           1.217900e-03 0.285639541 0.662 0.502  1.000000e+00       3
#> AC005152.32      1.238116e-03 0.232149990 0.561 0.430  1.000000e+00       3
#> FGF1             1.278025e-03 0.255100779 0.270 0.165  1.000000e+00       3
#> HS3ST11          1.280606e-03 0.061973569 0.304 0.169  1.000000e+00       3
#> DNAJC9           1.289496e-03 0.215206177 0.568 0.416  1.000000e+00       3
#> PGP              1.339238e-03 0.257958979 0.608 0.487  1.000000e+00       3
#> NEDD91           1.391028e-03 0.202040466 0.372 0.242  1.000000e+00       3
#> PDLIM1           1.469251e-03 0.182359148 0.466 0.322  1.000000e+00       3
#> SELENBP12        1.838647e-03 0.231421036 0.480 0.332  1.000000e+00       3
#> SDC12            1.953923e-03 0.286111571 0.669 0.512  1.000000e+00       3
#> IRX31            2.333554e-03 0.168483139 0.574 0.434  1.000000e+00       3
#> SLC12A2          2.382669e-03 0.285525988 0.385 0.268  1.000000e+00       3
#> ANXA11           2.741793e-03 0.225639365 0.642 0.488  1.000000e+00       3
#> BNIP31           2.941606e-03 0.191508321 0.595 0.459  1.000000e+00       3
#> LBR              2.981925e-03 0.156864489 0.493 0.341  1.000000e+00       3
#> PIR1             3.398444e-03 0.219894176 0.466 0.330  1.000000e+00       3
#> WDR341           3.691311e-03 0.215378874 0.601 0.490  1.000000e+00       3
#> SMTN2            3.727306e-03 0.180392194 0.568 0.440  1.000000e+00       3
#> COL4A21          3.754244e-03 0.036387355 0.480 0.327  1.000000e+00       3
#> KIF22            3.977987e-03 0.162400621 0.547 0.425  1.000000e+00       3
#> ANP32E           4.042592e-03 0.145551350 0.561 0.428  1.000000e+00       3
#> PPA11            4.282216e-03 0.230562011 0.628 0.482  1.000000e+00       3
#> THEM61           4.305151e-03 0.224324996 0.412 0.291  1.000000e+00       3
#> ATAD2            4.436055e-03 0.056059136 0.412 0.264  1.000000e+00       3
#> TOMM40           4.575686e-03 0.220664430 0.561 0.464  1.000000e+00       3
#> FH1              5.332981e-03 0.200749547 0.588 0.455  1.000000e+00       3
#> CTD-3065J16.9    5.469172e-03 0.240731163 0.318 0.219  1.000000e+00       3
#> GAPDH1           5.818581e-03 0.219108461 0.649 0.490  1.000000e+00       3
#> MFSD31           5.915277e-03 0.199881499 0.527 0.385  1.000000e+00       3
#> TCF7L11          5.986736e-03 0.156035929 0.412 0.287  1.000000e+00       3
#> ZBTB102          6.119927e-03 0.199459069 0.541 0.427  1.000000e+00       3
#> CTNNBIP1         7.246222e-03 0.185746249 0.554 0.474  1.000000e+00       3
#> LYAR             8.165315e-03 0.112820986 0.507 0.345  1.000000e+00       3
#> MAP1B1           8.593652e-03 0.171158386 0.486 0.418  1.000000e+00       3
#> HILPDA1          9.174038e-03 0.143451005 0.514 0.389  1.000000e+00       3
#> FXYD32           9.980606e-03 0.199576778 0.635 0.481  1.000000e+00       3
#> CLDN61           7.172546e-86 1.898126552 0.860 0.145  2.151764e-82       4
#> SMIM221          1.385336e-78 1.850543957 0.860 0.165  4.156007e-75       4
#> PCAT191          2.292897e-73 1.443169163 0.676 0.079  6.878691e-70       4
#> MAGIX1           1.984059e-67 1.835862883 0.956 0.278  5.952177e-64       4
#> RAB6B1           3.165224e-63 1.653225372 0.647 0.099  9.495672e-60       4
#> PQBP11           1.131886e-60 1.681424026 0.985 0.367  3.395659e-57       4
#> ILF21            7.823751e-60 1.775014860 0.956 0.376  2.347125e-56       4
#> KIF1A            1.506799e-59 1.381771931 0.566 0.067  4.520397e-56       4
#> CT831            2.168918e-59 1.383446274 0.662 0.102  6.506753e-56       4
#> PLP21            1.963081e-58 1.706808932 0.993 0.288  5.889242e-55       4
#> ORM1             1.624452e-56 1.408884319 0.522 0.062  4.873357e-53       4
#> CRABP21          1.736358e-56 1.501811028 0.699 0.145  5.209075e-53       4
#> TIMM17B1         5.490907e-56 1.553005290 0.956 0.396  1.647272e-52       4
#> CDK161           2.927723e-54 1.537308597 0.971 0.373  8.783169e-51       4
#> UXT1             1.008661e-53 1.547465781 0.956 0.392  3.025984e-50       4
#> GAL1             1.156208e-53 1.660530310 0.801 0.210  3.468625e-50       4
#> TRPM81           1.708858e-53 1.498088982 0.691 0.145  5.126575e-50       4
#> NDUFB111         2.282581e-53 1.522044480 0.971 0.373  6.847744e-50       4
#> MESP21           1.773470e-51 1.544480854 0.772 0.223  5.320411e-48       4
#> SLC34A2          7.917985e-51 1.389243583 0.434 0.044  2.375395e-47       4
#> PRRX21           8.013079e-51 1.375969704 0.662 0.133  2.403924e-47       4
#> EBP1             4.036302e-50 1.590181835 0.919 0.391  1.210891e-46       4
#> GPKOW1           2.150522e-49 1.537791731 0.904 0.373  6.451566e-46       4
#> KCNE51           2.379720e-49 1.503559253 0.632 0.131  7.139160e-46       4
#> SLC35A21         4.909414e-49 1.487535396 0.949 0.377  1.472824e-45       4
#> CENPV1           1.748037e-48 1.383922514 0.596 0.110  5.244110e-45       4
#> WDR451           4.028124e-47 1.470217984 0.934 0.421  1.208437e-43       4
#> WDR131           8.233149e-47 1.470991472 0.912 0.373  2.469945e-43       4
#> PGBD51           1.438964e-46 1.133722531 0.551 0.087  4.316892e-43       4
#> RBM31            2.223040e-46 1.343763929 0.963 0.375  6.669120e-43       4
#> FTSJ11           9.292508e-46 1.436631720 0.868 0.322  2.787753e-42       4
#> UBA11            1.549807e-45 1.392987412 0.941 0.337  4.649422e-42       4
#> MAGEA4           1.110574e-44 1.130340888 0.419 0.046  3.331722e-41       4
#> C5orf66-AS1      1.625509e-44 1.295060965 0.485 0.073  4.876526e-41       4
#> RPS31            3.316876e-43 1.378691108 0.890 0.415  9.950627e-40       4
#> SOX181           3.425564e-43 0.951983507 0.603 0.113  1.027669e-39       4
#> VTCN13           1.828317e-42 1.328268169 0.904 0.431  5.484952e-39       4
#> HDAC61           2.913568e-41 1.223745768 0.831 0.291  8.740704e-38       4
#> SLC2A111         3.667236e-41 1.199912198 0.735 0.202  1.100171e-37       4
#> IRX41            4.643824e-41 1.014042077 0.463 0.067  1.393147e-37       4
#> SLC6A141         3.473981e-40 1.174548311 0.426 0.059  1.042194e-36       4
#> TIMM8B           2.554075e-39 1.370732433 0.860 0.427  7.662225e-36       4
#> NFASC1           3.736740e-39 1.112963506 0.412 0.056  1.121022e-35       4
#> PDXK1            3.740094e-39 1.187155202 0.934 0.431  1.122028e-35       4
#> OSR2             5.150255e-39 1.141150287 0.404 0.054  1.545077e-35       4
#> MISP1            9.847685e-39 0.969952839 0.463 0.070  2.954306e-35       4
#> DPM11            1.508372e-38 1.254218924 0.912 0.440  4.525115e-35       4
#> CALB21           2.155899e-38 1.283496609 0.750 0.258  6.467698e-35       4
#> CALML51          2.400416e-38 1.176245422 0.941 0.249  7.201247e-35       4
#> RHOV1            3.671111e-38 1.154212077 0.618 0.155  1.101333e-34       4
#> PORCN1           4.926011e-38 0.969306452 0.544 0.105  1.477803e-34       4
#> ALCAM1           5.110847e-38 1.283657630 0.787 0.280  1.533254e-34       4
#> SLC9A71          1.528571e-37 1.011032617 0.603 0.136  4.585712e-34       4
#> MORN31           1.538417e-37 1.121848643 0.544 0.111  4.615251e-34       4
#> RGN1             1.940547e-37 1.102513808 0.419 0.060  5.821641e-34       4
#> SYCE1L1          2.032248e-37 1.123828423 0.699 0.206  6.096744e-34       4
#> PTGFRN1          3.723839e-37 1.307417322 0.860 0.379  1.117152e-33       4
#> RP11-817J15.3    6.440338e-37 1.051879341 0.338 0.036  1.932101e-33       4
#> PVRL11           1.320159e-36 1.264830070 0.640 0.178  3.960477e-33       4
#> VGLL11           2.153917e-36 1.303733323 0.772 0.317  6.461752e-33       4
#> MTL52            2.426869e-36 1.282960710 0.846 0.348  7.280607e-33       4
#> OTUD51           2.718682e-36 1.201798062 0.875 0.408  8.156047e-33       4
#> ARMC11           7.173800e-36 1.104560326 0.735 0.239  2.152140e-32       4
#> CLDN11           7.560686e-36 1.010073125 0.581 0.130  2.268206e-32       4
#> HES41            1.080996e-35 1.112465130 0.838 0.347  3.242987e-32       4
#> ORM21            1.614276e-35 1.360770455 0.426 0.075  4.842828e-32       4
#> ATP9A1           4.398405e-35 1.118435248 0.779 0.291  1.319521e-31       4
#> LRP51            5.255852e-35 1.084142143 0.721 0.252  1.576755e-31       4
#> GS1-124K5.41     1.185557e-34 1.066329982 0.581 0.138  3.556671e-31       4
#> ADCYAP1          1.215450e-34 1.099527534 0.272 0.022  3.646351e-31       4
#> ELK11            1.894511e-34 1.135461450 0.713 0.234  5.683532e-31       4
#> HMGB31           1.912694e-34 1.244742878 0.882 0.440  5.738082e-31       4
#> TTC9             2.378791e-34 0.954058155 0.360 0.047  7.136372e-31       4
#> RGS201           2.840726e-34 0.992752287 0.441 0.076  8.522177e-31       4
#> NPDC11           6.221367e-34 1.001474464 0.787 0.278  1.866410e-30       4
#> CRYAB2           7.736342e-34 1.071593439 0.956 0.517  2.320903e-30       4
#> ART31            8.102727e-34 1.148930972 0.618 0.171  2.430818e-30       4
#> S100A141         8.290342e-34 1.155206479 0.853 0.371  2.487103e-30       4
#> ACTG22           1.002418e-33 1.174851002 0.809 0.336  3.007254e-30       4
#> CSF1             1.079182e-33 0.912272158 0.353 0.045  3.237547e-30       4
#> UBASH3B          3.138776e-33 0.793030840 0.375 0.051  9.416329e-30       4
#> CFB1             3.421761e-33 0.938792361 0.684 0.211  1.026528e-29       4
#> RP11-554I8.21    3.661902e-33 1.027471081 0.699 0.217  1.098571e-29       4
#> RP11-817J15.2    4.467376e-33 0.847922166 0.257 0.020  1.340213e-29       4
#> STC21            4.597617e-33 1.003538160 0.610 0.164  1.379285e-29       4
#> EFNA5            4.697823e-33 1.085283448 0.507 0.110  1.409347e-29       4
#> ARAF1            8.441476e-33 1.080284392 0.816 0.349  2.532443e-29       4
#> UBE2T            1.314790e-32 1.002174302 0.676 0.245  3.944370e-29       4
#> TM4SF1-AS12      1.510154e-32 1.195887064 0.735 0.278  4.530462e-29       4
#> PTX3             1.871851e-32 1.052649002 0.346 0.047  5.615553e-29       4
#> CLDN71           2.138563e-32 1.054672433 0.934 0.497  6.415690e-29       4
#> S100A161         3.335641e-32 1.050944115 0.721 0.241  1.000692e-28       4
#> PPP1R14B1        4.440938e-32 1.132070809 0.860 0.424  1.332281e-28       4
#> LA16c-380H5.51   5.088659e-32 1.162150667 0.625 0.195  1.526598e-28       4
#> RBBP71           7.008060e-32 1.121862172 0.875 0.472  2.102418e-28       4
#> CHST9            2.116316e-31 0.913694131 0.287 0.029  6.348947e-28       4
#> PRRT3-AS11       2.245781e-31 0.884372274 0.419 0.075  6.737343e-28       4
#> GPR271           3.574317e-31 0.990497202 0.375 0.059  1.072295e-27       4
#> CAGE1            3.767349e-31 0.939075512 0.250 0.021  1.130205e-27       4
#> C1orf61          3.829972e-31 0.988933948 0.294 0.032  1.148992e-27       4
#> RAB30-AS11       4.491545e-31 1.217980348 0.809 0.385  1.347464e-27       4
#> ARHGDIG          5.804513e-31 0.973326385 0.250 0.021  1.741354e-27       4
#> HYOU11           7.541005e-31 1.187707199 0.728 0.312  2.262301e-27       4
#> IL322            1.895629e-30 1.068412905 0.868 0.420  5.686887e-27       4
#> TBC1D251         2.117996e-30 1.081840087 0.713 0.263  6.353989e-27       4
#> ATP1B11          2.515766e-30 1.091871596 0.860 0.410  7.547297e-27       4
#> TMEM2081         2.602462e-30 1.119125735 0.853 0.457  7.807385e-27       4
#> IFITM32          2.651474e-30 1.091485729 0.963 0.377  7.954422e-27       4
#> LINC01133        3.740517e-30 1.055226926 0.353 0.055  1.122155e-26       4
#> MANEAL1          6.939547e-30 0.943664880 0.353 0.054  2.081864e-26       4
#> ZNF674-AS11      6.999698e-30 1.020224483 0.507 0.122  2.099909e-26       4
#> KRT191           9.350760e-30 1.047176041 0.919 0.471  2.805228e-26       4
#> MFSD32           1.255591e-29 1.158251982 0.772 0.352  3.766773e-26       4
#> GINS4            4.501051e-29 1.164397422 0.515 0.140  1.350315e-25       4
#> TMPRSS31         6.791597e-29 1.011438339 0.647 0.212  2.037479e-25       4
#> CDH12            7.440827e-29 1.041978539 0.816 0.398  2.232248e-25       4
#> PRAF21           7.933113e-29 1.027655221 0.787 0.332  2.379934e-25       4
#> FDX1             9.170085e-29 1.148684790 0.831 0.437  2.751025e-25       4
#> HEBP22           2.100592e-28 0.948927619 0.919 0.507  6.301775e-25       4
#> AARD2            3.219198e-28 1.112250390 0.794 0.406  9.657594e-25       4
#> SEC131           4.874985e-28 1.071619507 0.824 0.383  1.462496e-24       4
#> RRP361           5.377122e-28 1.121464892 0.809 0.431  1.613137e-24       4
#> DDR11            6.829353e-28 1.063050672 0.794 0.390  2.048806e-24       4
#> BRK11            1.140636e-27 0.970460789 0.787 0.349  3.421907e-24       4
#> NUCKS12          1.646412e-27 1.021485609 0.890 0.454  4.939237e-24       4
#> TLDC11           1.707512e-27 1.042130951 0.809 0.359  5.122537e-24       4
#> SPINT12          1.767105e-27 1.015934341 0.831 0.466  5.301316e-24       4
#> RP11-490M8.1     2.855014e-27 0.877309581 0.324 0.049  8.565041e-24       4
#> PPP1R14C1        3.396415e-27 0.932382790 0.676 0.235  1.018925e-23       4
#> CEP57            4.593031e-27 1.086927073 0.772 0.401  1.377909e-23       4
#> METRN1           6.255259e-27 0.981998088 0.824 0.426  1.876578e-23       4
#> CD242            6.425872e-27 0.890673986 0.978 0.612  1.927762e-23       4
#> RBP11            7.652511e-27 1.024337630 0.676 0.308  2.295753e-23       4
#> TCF7L12          8.299416e-27 0.952326896 0.691 0.249  2.489825e-23       4
#> PPDPF2           9.228712e-27 0.884515551 0.912 0.531  2.768613e-23       4
#> KRTCAP32         1.022823e-26 1.041344043 0.831 0.435  3.068470e-23       4
#> FSTL11           1.075824e-26 0.532977930 0.596 0.169  3.227471e-23       4
#> PSAT1            1.573417e-26 0.809832133 0.375 0.071  4.720250e-23       4
#> TM4SF12          1.673544e-26 0.971331093 0.956 0.546  5.020631e-23       4
#> LYNX1            1.740085e-26 0.921940087 0.294 0.042  5.220255e-23       4
#> LSR1             2.034434e-26 0.993380901 0.868 0.465  6.103302e-23       4
#> FBXO91           2.072487e-26 1.066056282 0.765 0.366  6.217462e-23       4
#> MLLT41           2.362881e-26 0.975710993 0.772 0.367  7.088642e-23       4
#> STARD101         2.523399e-26 1.063601964 0.787 0.381  7.570197e-23       4
#> ANXA31           3.029570e-26 0.863223832 0.375 0.072  9.088710e-23       4
#> PHGDH1           7.788014e-26 0.984979948 0.772 0.359  2.336404e-22       4
#> HSPB21           1.010190e-25 1.008002092 0.699 0.313  3.030569e-22       4
#> KRT231           2.597550e-25 0.972485414 0.838 0.462  7.792651e-22       4
#> POLD21           2.954429e-25 1.054321156 0.779 0.406  8.863286e-22       4
#> GRB101           3.274860e-25 1.000428827 0.625 0.229  9.824579e-22       4
#> KRBOX41          3.423131e-25 0.871491784 0.610 0.199  1.026939e-21       4
#> ECEL11           4.324398e-25 0.818854904 0.471 0.120  1.297319e-21       4
#> TMEM251          5.016035e-25 0.990584474 0.654 0.258  1.504810e-21       4
#> CCDC1671         5.923117e-25 1.006329812 0.846 0.459  1.776935e-21       4
#> KCNG11           1.680830e-24 0.819853521 0.632 0.235  5.042489e-21       4
#> EEF1A2           1.725039e-24 0.814284285 0.294 0.046  5.175116e-21       4
#> SEMA6A1          1.870668e-24 0.798062832 0.368 0.073  5.612003e-21       4
#> LSM51            2.419468e-24 0.999995899 0.801 0.435  7.258405e-21       4
#> APOD1            2.650692e-24 0.738470142 0.647 0.211  7.952075e-21       4
#> RANBP11          2.676641e-24 1.003126461 0.809 0.467  8.029924e-21       4
#> KLK5             3.047521e-24 0.875146330 0.250 0.032  9.142562e-21       4
#> SERPINH12        3.128555e-24 0.903218644 0.890 0.492  9.385665e-21       4
#> RPUSD31          4.682143e-24 0.954243518 0.728 0.383  1.404643e-20       4
#> CENPW            4.699974e-24 1.053852335 0.610 0.257  1.409992e-20       4
#> LCN22            7.323997e-24 0.948361761 0.743 0.338  2.197199e-20       4
#> LY6E2            7.656108e-24 0.899323250 0.853 0.425  2.296832e-20       4
#> DNAAF11          7.941555e-24 0.754291407 0.471 0.125  2.382467e-20       4
#> MESP13           7.956601e-24 0.990207722 0.794 0.408  2.386980e-20       4
#> MAD2L1BP1        9.738611e-24 0.914656590 0.809 0.392  2.921583e-20       4
#> EGLN11           1.089194e-23 0.851814689 0.699 0.261  3.267582e-20       4
#> IRX32            1.112658e-23 0.971704241 0.787 0.406  3.337975e-20       4
#> FADS21           1.281001e-23 0.878223733 0.625 0.235  3.843002e-20       4
#> KLK61            1.410883e-23 0.654059087 0.507 0.152  4.232649e-20       4
#> PDIA31           1.444194e-23 0.878170152 0.868 0.461  4.332581e-20       4
#> SLC7A51          1.600571e-23 0.890143593 0.691 0.274  4.801713e-20       4
#> DSP1             1.731135e-23 0.924951383 0.838 0.462  5.193406e-20       4
#> H2AFX            1.808647e-23 0.949872555 0.706 0.329  5.425942e-20       4
#> MYBL2            1.969246e-23 0.898361882 0.463 0.126  5.907739e-20       4
#> CAPN22           2.224577e-23 0.915126727 0.882 0.477  6.673730e-20       4
#> PITX12           3.200489e-23 0.936204245 0.750 0.376  9.601467e-20       4
#> SCD5             3.646763e-23 0.531073141 0.272 0.039  1.094029e-19       4
#> CEBPB1           4.331437e-23 0.886388973 0.860 0.436  1.299431e-19       4
#> VHL1             4.565152e-23 0.858246849 0.699 0.303  1.369546e-19       4
#> CLDN9            5.504843e-23 0.835406148 0.309 0.056  1.651453e-19       4
#> PPM1H            5.744293e-23 0.813547704 0.412 0.100  1.723288e-19       4
#> RAI141           6.502803e-23 0.806751080 0.735 0.316  1.950841e-19       4
#> SECTM11          7.739345e-23 0.678912013 0.588 0.183  2.321803e-19       4
#> DCUN1D5          1.000745e-22 1.072440881 0.750 0.440  3.002236e-19       4
#> LEMD12           1.004369e-22 0.878440159 0.875 0.448  3.013106e-19       4
#> CENPF            1.033656e-22 0.770642493 0.485 0.164  3.100968e-19       4
#> PVRL41           1.145226e-22 0.923519540 0.750 0.352  3.435678e-19       4
#> FAM92A11         1.174517e-22 0.806828322 0.669 0.245  3.523550e-19       4
#> MMP151           1.224754e-22 0.864179500 0.757 0.362  3.674261e-19       4
#> SYNCRIP1         1.258046e-22 0.951012146 0.831 0.470  3.774139e-19       4
#> FASN1            1.463820e-22 0.939994843 0.610 0.232  4.391460e-19       4
#> PPFIA12          1.828916e-22 0.868911484 0.853 0.471  5.486749e-19       4
#> DNAJC152         1.927290e-22 0.798371232 0.750 0.306  5.781871e-19       4
#> ARHGAP291        1.981651e-22 0.875934372 0.794 0.388  5.944952e-19       4
#> SIPA1L21         2.053088e-22 0.734497150 0.574 0.184  6.159264e-19       4
#> TYMS             2.067866e-22 0.787856039 0.478 0.140  6.203598e-19       4
#> ZPR1             2.107625e-22 0.993874874 0.787 0.395  6.322875e-19       4
#> TSEN151          2.360860e-22 0.965640833 0.794 0.441  7.082581e-19       4
#> KNSTRN           3.788590e-22 0.801903068 0.456 0.123  1.136577e-18       4
#> ATP2A1-AS1       3.845578e-22 0.762999878 0.294 0.051  1.153673e-18       4
#> PRICKLE11        3.887346e-22 0.781594593 0.662 0.257  1.166204e-18       4
#> CCT6A1           6.251026e-22 0.886064077 0.816 0.479  1.875308e-18       4
#> CX3CL1           8.007938e-22 0.797479430 0.324 0.065  2.402382e-18       4
#> TTF21            1.044221e-21 0.959177592 0.669 0.292  3.132664e-18       4
#> QDPR1            1.118813e-21 0.882491458 0.765 0.400  3.356438e-18       4
#> TMEM200A1        1.571217e-21 0.639140348 0.331 0.066  4.713652e-18       4
#> TMEM63A1         1.784398e-21 0.943079127 0.721 0.377  5.353193e-18       4
#> XAGE21           1.852123e-21 1.017919900 0.618 0.267  5.556368e-18       4
#> BARX12           1.939505e-21 0.863541687 0.743 0.360  5.818516e-18       4
#> CDH32            2.266755e-21 0.904092366 0.691 0.367  6.800264e-18       4
#> BCL9L1           2.527059e-21 0.798176232 0.566 0.196  7.581176e-18       4
#> SNCG             3.249866e-21 0.858434706 0.272 0.047  9.749599e-18       4
#> GALNT21          4.817527e-21 0.870688031 0.721 0.341  1.445258e-17       4
#> TM4SF18          4.998846e-21 0.789793750 0.294 0.055  1.499654e-17       4
#> HMGA12           6.121729e-21 0.864776940 0.801 0.484  1.836519e-17       4
#> TMEM132A1        7.252424e-21 0.862452999 0.801 0.415  2.175727e-17       4
#> NLRP2            8.877852e-21 0.715491352 0.279 0.049  2.663356e-17       4
#> PCSK1N1          1.088298e-20 0.709795299 0.757 0.358  3.264893e-17       4
#> SUSD21           1.092727e-20 0.654850429 0.346 0.076  3.278180e-17       4
#> RGS102           1.329875e-20 0.709271358 0.838 0.404  3.989626e-17       4
#> GCNT2            1.580787e-20 0.849219306 0.507 0.169  4.742361e-17       4
#> TNFRSF12A1       1.838591e-20 0.870482215 0.772 0.399  5.515774e-17       4
#> OAZ31            1.891812e-20 0.832251513 0.522 0.177  5.675435e-17       4
#> NRBP21           2.477724e-20 0.845908448 0.669 0.299  7.433173e-17       4
#> CACNA1A1         2.661427e-20 0.735438804 0.632 0.264  7.984280e-17       4
#> RASGRP11         3.288589e-20 0.645798694 0.353 0.080  9.865766e-17       4
#> OGG11            3.656997e-20 0.804629592 0.588 0.228  1.097099e-16       4
#> NOA11            3.813171e-20 0.860196006 0.743 0.353  1.143951e-16       4
#> SUN11            4.718501e-20 0.834747550 0.787 0.414  1.415550e-16       4
#> WNT5B            5.576143e-20 0.674498048 0.294 0.056  1.672843e-16       4
#> NUPR21           6.037937e-20 0.738095405 0.676 0.337  1.811381e-16       4
#> DUT              6.340921e-20 0.837521545 0.794 0.425  1.902276e-16       4
#> UFD1L1           6.406066e-20 0.874707172 0.757 0.462  1.921820e-16       4
#> NT5E1            7.389391e-20 0.705609682 0.529 0.176  2.216817e-16       4
#> PRRT31           8.939699e-20 0.761387473 0.529 0.175  2.681910e-16       4
#> PSORS1C11        9.166254e-20 0.755565571 0.478 0.147  2.749876e-16       4
#> RP11-273G15.2    1.123001e-19 0.676624094 0.316 0.066  3.369004e-16       4
#> AGT2             1.144623e-19 0.855279675 0.721 0.318  3.433869e-16       4
#> C11orf73         1.159242e-19 0.923241890 0.794 0.488  3.477726e-16       4
#> CMSS11           1.373589e-19 0.814851932 0.669 0.288  4.120767e-16       4
#> FERMT1           1.520976e-19 0.791992172 0.287 0.057  4.562928e-16       4
#> PDIA61           1.682910e-19 0.763545222 0.853 0.473  5.048729e-16       4
#> CUL71            2.379311e-19 0.808569330 0.691 0.301  7.137934e-16       4
#> MIS18A           2.750996e-19 0.813505863 0.640 0.258  8.252987e-16       4
#> GABRP1           3.004609e-19 0.790851792 0.809 0.415  9.013828e-16       4
#> ARHGEF12         3.010384e-19 0.877831986 0.684 0.336  9.031153e-16       4
#> CXCL10           3.770785e-19 0.414923843 0.301 0.062  1.131236e-15       4
#> SUSD4            4.942815e-19 0.679471491 0.279 0.053  1.482844e-15       4
#> SELM2            5.629395e-19 0.800221774 0.831 0.499  1.688819e-15       4
#> IRAK11           6.517210e-19 0.828163216 0.706 0.368  1.955163e-15       4
#> LINC009601       8.341241e-19 0.818383210 0.507 0.179  2.502372e-15       4
#> NMB1             9.733749e-19 0.616742890 0.640 0.245  2.920125e-15       4
#> NAT8L1           1.025930e-18 0.696199341 0.537 0.190  3.077789e-15       4
#> NTN11            1.110213e-18 0.834410738 0.566 0.221  3.330638e-15       4
#> DEK              1.147381e-18 0.898490495 0.757 0.482  3.442144e-15       4
#> APP2             1.175154e-18 0.797404369 0.801 0.463  3.525462e-15       4
#> H2AFY21          1.199626e-18 0.656656450 0.522 0.173  3.598878e-15       4
#> ELF32            1.218342e-18 0.786319650 0.831 0.456  3.655027e-15       4
#> PADI21           1.302602e-18 0.863175577 0.669 0.324  3.907806e-15       4
#> HES61            1.875496e-18 0.688951726 0.559 0.203  5.626488e-15       4
#> PLK1             2.643717e-18 0.549967087 0.382 0.103  7.931150e-15       4
#> ZNF2171          2.782593e-18 0.700930635 0.684 0.315  8.347778e-15       4
#> LSM4             3.274910e-18 0.823419909 0.765 0.463  9.824731e-15       4
#> GUCY1A3          3.557336e-18 0.345326285 0.368 0.091  1.067201e-14       4
#> EPPK11           3.896843e-18 0.699673568 0.507 0.180  1.169053e-14       4
#> RRAGD1           3.903175e-18 0.783529240 0.507 0.182  1.170953e-14       4
#> ECT2             3.982954e-18 0.878154640 0.346 0.091  1.194886e-14       4
#> LDOC11           4.108293e-18 0.776665529 0.750 0.375  1.232488e-14       4
#> WEE11            5.666601e-18 0.821487561 0.676 0.320  1.699980e-14       4
#> SMYD21           6.161925e-18 0.816136202 0.728 0.416  1.848577e-14       4
#> PSCA1            8.133416e-18 0.744447070 0.574 0.233  2.440025e-14       4
#> UBE2C            9.236940e-18 0.620291806 0.493 0.181  2.771082e-14       4
#> CDC20            9.902405e-18 0.705557707 0.426 0.135  2.970722e-14       4
#> ITGB81           1.076385e-17 0.749324826 0.669 0.296  3.229155e-14       4
#> CLDN32           1.270951e-17 0.793190997 0.912 0.480  3.812853e-14       4
#> QPRT1            1.426699e-17 0.893478293 0.551 0.248  4.280096e-14       4
#> KLHL351          1.454815e-17 0.808857128 0.669 0.337  4.364446e-14       4
#> PGP1             1.978985e-17 0.848547814 0.772 0.465  5.936954e-14       4
#> PPA12            2.012419e-17 0.749235442 0.816 0.457  6.037258e-14       4
#> TJP31            2.178822e-17 0.661193896 0.471 0.155  6.536467e-14       4
#> IRX51            2.199326e-17 0.694581911 0.566 0.213  6.597977e-14       4
#> SORL11           3.037888e-17 0.590920997 0.478 0.154  9.113665e-14       4
#> PTP4A32          3.077201e-17 0.827282896 0.765 0.445  9.231602e-14       4
#> SETD51           3.720912e-17 0.784395080 0.684 0.379  1.116274e-13       4
#> CD551            4.004568e-17 0.713078034 0.750 0.435  1.201370e-13       4
#> THUMPD3-AS11     4.260138e-17 0.796907593 0.713 0.413  1.278041e-13       4
#> KIFC1            4.475755e-17 0.742299154 0.360 0.101  1.342726e-13       4
#> PMAIP11          4.675024e-17 0.704271780 0.596 0.310  1.402507e-13       4
#> CENPA            5.308388e-17 0.674926762 0.301 0.071  1.592516e-13       4
#> THUMPD31         6.222076e-17 0.872943465 0.640 0.358  1.866623e-13       4
#> FAM64A           7.335513e-17 0.750817645 0.287 0.066  2.200654e-13       4
#> PXDN1            7.728010e-17 0.706853869 0.699 0.353  2.318403e-13       4
#> H2AFV1           7.866995e-17 0.816767567 0.772 0.458  2.360099e-13       4
#> NOTCH3           1.054527e-16 0.253343917 0.353 0.089  3.163580e-13       4
#> CDKN3            1.086727e-16 0.608816098 0.375 0.108  3.260182e-13       4
#> PLA2G16          1.142550e-16 0.715793107 0.375 0.110  3.427649e-13       4
#> RIPK2            1.374110e-16 0.841913680 0.699 0.415  4.122331e-13       4
#> DAPK2            1.390833e-16 0.543917509 0.360 0.097  4.172498e-13       4
#> FANCD2           1.540500e-16 0.890338720 0.397 0.128  4.621500e-13       4
#> PTPRF2           1.626724e-16 0.752303079 0.750 0.426  4.880172e-13       4
#> DNPH11           1.805159e-16 0.764632025 0.787 0.504  5.415477e-13       4
#> C1QBP            2.056517e-16 0.801092883 0.801 0.470  6.169551e-13       4
#> NUF2             2.522735e-16 0.712849239 0.353 0.101  7.568204e-13       4
#> KLRG22           2.761046e-16 0.738913046 0.478 0.175  8.283139e-13       4
#> SNX81            2.810423e-16 0.771861090 0.669 0.359  8.431269e-13       4
#> FH2              3.509502e-16 0.810388471 0.699 0.441  1.052851e-12       4
#> EIF3B1           3.957007e-16 0.818600252 0.750 0.450  1.187102e-12       4
#> MUC11            3.962644e-16 0.754898457 0.654 0.333  1.188793e-12       4
#> CRNDE2           4.076731e-16 0.770360368 0.779 0.475  1.223019e-12       4
#> CCND11           4.477525e-16 0.802613607 0.728 0.443  1.343258e-12       4
#> MUC5B1           4.697229e-16 0.629198356 0.529 0.204  1.409169e-12       4
#> KRT812           5.051569e-16 0.737569105 0.838 0.485  1.515471e-12       4
#> GATA6-AS11       5.334845e-16 0.608582603 0.301 0.074  1.600453e-12       4
#> LPIN11           5.852886e-16 0.716910357 0.618 0.297  1.755866e-12       4
#> PTS2             7.432408e-16 0.762679037 0.757 0.447  2.229722e-12       4
#> BCAS41           8.738785e-16 0.684901498 0.603 0.271  2.621636e-12       4
#> GBA1             9.103199e-16 0.763457461 0.750 0.475  2.730960e-12       4
#> TMPRSS131        1.055791e-15 0.708455444 0.574 0.247  3.167373e-12       4
#> BIRC5            1.193294e-15 0.526971215 0.404 0.139  3.579882e-12       4
#> SLC1A31          1.203015e-15 0.400979401 0.463 0.145  3.609045e-12       4
#> FZD81            1.324752e-15 0.605919235 0.294 0.073  3.974256e-12       4
#> MND1             1.526411e-15 0.691638416 0.309 0.081  4.579234e-12       4
#> LAD11            1.827540e-15 0.671989975 0.735 0.391  5.482621e-12       4
#> SEZ6L21          2.095850e-15 0.468474677 0.449 0.148  6.287550e-12       4
#> DIAPH3           2.141103e-15 0.639947802 0.390 0.124  6.423308e-12       4
#> CHAF1B1          2.197141e-15 0.638043026 0.441 0.151  6.591423e-12       4
#> LINC007071       2.258599e-15 0.683842274 0.331 0.094  6.775796e-12       4
#> TACSTD22         2.490200e-15 0.669659444 0.904 0.418  7.470600e-12       4
#> CDCA3            2.694052e-15 0.693437629 0.375 0.120  8.082157e-12       4
#> FGF132           2.883288e-15 0.681867438 0.691 0.333  8.649864e-12       4
#> APMAP1           2.916299e-15 0.692780871 0.794 0.470  8.748897e-12       4
#> CHORDC1          3.590291e-15 0.799332214 0.735 0.455  1.077087e-11       4
#> KRT141           3.601771e-15 0.599373527 0.265 0.061  1.080531e-11       4
#> YES11            3.654877e-15 0.669720792 0.647 0.328  1.096463e-11       4
#> BYSL1            4.602182e-15 0.707058899 0.699 0.396  1.380654e-11       4
#> GLYATL21         4.690023e-15 0.721698760 0.515 0.225  1.407007e-11       4
#> CARHSP11         4.850131e-15 0.836706601 0.721 0.467  1.455039e-11       4
#> TARBP11          5.076272e-15 0.673917608 0.596 0.278  1.522882e-11       4
#> EMC31            6.043354e-15 0.835455692 0.699 0.425  1.813006e-11       4
#> YDJC             6.186822e-15 0.751859976 0.706 0.436  1.856046e-11       4
#> MIR4435-2HG1     7.091097e-15 0.717310459 0.662 0.391  2.127329e-11       4
#> MCAM1            9.050940e-15 0.720911764 0.632 0.332  2.715282e-11       4
#> LMNB1            1.005104e-14 0.488167977 0.434 0.146  3.015312e-11       4
#> TUNAR1           1.014479e-14 0.580387985 0.265 0.062  3.043436e-11       4
#> ZNF695           1.138730e-14 0.656425195 0.309 0.083  3.416191e-11       4
#> GINS2            1.259236e-14 0.675469841 0.412 0.146  3.777709e-11       4
#> BRPF11           1.285077e-14 0.708474700 0.360 0.116  3.855230e-11       4
#> VEGFA1           1.324765e-14 0.602512950 0.632 0.335  3.974295e-11       4
#> RCCD1            1.340257e-14 0.602507726 0.515 0.204  4.020771e-11       4
#> TJP11            1.517176e-14 0.594055895 0.603 0.270  4.551527e-11       4
#> RNF81            1.552476e-14 0.702890460 0.669 0.367  4.657429e-11       4
#> TRIP13           1.615002e-14 0.601902455 0.368 0.116  4.845007e-11       4
#> PRC1             1.628529e-14 0.672034615 0.434 0.161  4.885586e-11       4
#> TIMP11           1.686248e-14 0.544462010 0.647 0.338  5.058745e-11       4
#> SDF2L11          1.807728e-14 0.714195234 0.721 0.443  5.423183e-11       4
#> IL17RC1          2.028915e-14 0.600080647 0.647 0.329  6.086746e-11       4
#> GPC11            2.476503e-14 0.574302427 0.632 0.338  7.429509e-11       4
#> RAB252           2.537488e-14 0.563763200 0.706 0.357  7.612465e-11       4
#> SEC11C1          3.302041e-14 0.714067042 0.721 0.448  9.906122e-11       4
#> ASPM             3.341792e-14 0.693214867 0.265 0.067  1.002538e-10       4
#> SDR16C51         3.756629e-14 0.493936955 0.390 0.130  1.126989e-10       4
#> GGT5             3.869321e-14 0.271234865 0.360 0.106  1.160796e-10       4
#> MTMR141          4.453414e-14 0.716257312 0.669 0.381  1.336024e-10       4
#> LTF              4.529488e-14 0.832414689 0.338 0.107  1.358846e-10       4
#> CHST71           5.038781e-14 0.477250274 0.368 0.112  1.511634e-10       4
#> S100A21          5.742583e-14 0.509764700 0.522 0.239  1.722775e-10       4
#> SAPCD2           7.185525e-14 0.633554903 0.294 0.081  2.155657e-10       4
#> ZNRF21           8.126568e-14 0.532718475 0.449 0.163  2.437970e-10       4
#> AURKB            8.475109e-14 0.611920125 0.316 0.093  2.542533e-10       4
#> PTPN141          1.010463e-13 0.570907444 0.625 0.326  3.031390e-10       4
#> OCLN1            1.053325e-13 0.617101984 0.588 0.262  3.159974e-10       4
#> NEK2             1.191163e-13 0.522672956 0.294 0.081  3.573490e-10       4
#> CP3              1.285218e-13 0.774387388 0.647 0.374  3.855653e-10       4
#> CEP55            1.387550e-13 0.591345147 0.294 0.082  4.162650e-10       4
#> CTA-293F17.11    1.430273e-13 0.644745627 0.551 0.282  4.290818e-10       4
#> ANP32E1          1.679694e-13 0.787838334 0.684 0.412  5.039083e-10       4
#> FZD7             1.687224e-13 0.651029170 0.478 0.198  5.061671e-10       4
#> RP11-410L14.2    1.715298e-13 0.578778861 0.250 0.060  5.145894e-10       4
#> RP11-783K16.51   1.870705e-13 0.505914934 0.419 0.148  5.612116e-10       4
#> FKBP102          1.923785e-13 0.633248702 0.728 0.416  5.771355e-10       4
#> SLC6A81          1.966140e-13 0.526080555 0.434 0.156  5.898419e-10       4
#> RAD51AP1         1.969820e-13 0.537842092 0.434 0.159  5.909460e-10       4
#> RIC31            1.970661e-13 0.645965343 0.654 0.361  5.911982e-10       4
#> ORC6             1.984773e-13 0.702187957 0.368 0.127  5.954318e-10       4
#> CRACR2B          2.002090e-13 0.561467950 0.471 0.182  6.006271e-10       4
#> DHCR24           2.726816e-13 0.652017652 0.397 0.142  8.180449e-10       4
#> GGCT1            3.005084e-13 0.672179669 0.713 0.457  9.015252e-10       4
#> NETO2            3.141913e-13 0.643874696 0.397 0.145  9.425739e-10       4
#> CENPN            3.145706e-13 0.782408331 0.551 0.265  9.437117e-10       4
#> IGSF23           3.150933e-13 0.593253483 0.250 0.062  9.452800e-10       4
#> DGAT21           3.154207e-13 0.712734905 0.596 0.290  9.462622e-10       4
#> KIAA0101         4.678599e-13 0.483274884 0.456 0.201  1.403580e-09       4
#> SLC12A21         5.064059e-13 0.603857790 0.537 0.248  1.519218e-09       4
#> MLPH1            5.538254e-13 0.558756034 0.360 0.121  1.661476e-09       4
#> TACC3            6.576397e-13 0.397854099 0.456 0.172  1.972919e-09       4
#> C16orf59         6.583011e-13 0.616655101 0.265 0.072  1.974903e-09       4
#> FAM129A          1.065170e-12 0.541558567 0.382 0.136  3.195511e-09       4
#> PBK              1.106203e-12 0.407644200 0.294 0.083  3.318608e-09       4
#> CENPQ1           1.240922e-12 0.765582176 0.618 0.351  3.722765e-09       4
#> FAM208B1         1.332277e-12 0.680990373 0.618 0.383  3.996830e-09       4
#> POC1A            1.347975e-12 0.557561494 0.346 0.116  4.043924e-09       4
#> LGALSL1          1.497755e-12 0.524277272 0.471 0.195  4.493265e-09       4
#> MAD2L1           1.659904e-12 0.552520296 0.404 0.163  4.979711e-09       4
#> SLC38A1          1.687041e-12 0.670816770 0.618 0.333  5.061122e-09       4
#> MTURN1           1.898869e-12 0.609441264 0.544 0.253  5.696608e-09       4
#> NPTXR1           2.082622e-12 0.440468120 0.507 0.204  6.247866e-09       4
#> TTLL31           2.137369e-12 0.504709115 0.419 0.156  6.412107e-09       4
#> RP11-19E11.12    2.408620e-12 0.646130114 0.706 0.387  7.225859e-09       4
#> APOBEC3B         2.619093e-12 0.556024214 0.375 0.133  7.857278e-09       4
#> LDLR             2.639731e-12 0.478977133 0.368 0.127  7.919192e-09       4
#> KDELC2           3.098095e-12 0.532405265 0.419 0.159  9.294285e-09       4
#> CCNA2            3.206539e-12 0.507816120 0.294 0.087  9.619618e-09       4
#> FLNB1            3.455989e-12 0.609596222 0.654 0.399  1.036797e-08       4
#> CENPU            3.686705e-12 0.587656090 0.382 0.142  1.106011e-08       4
#> MBD21            3.773690e-12 0.666714837 0.684 0.445  1.132107e-08       4
#> PROM12           4.207081e-12 0.598209380 0.750 0.462  1.262124e-08       4
#> TOP2A            4.366687e-12 0.436922690 0.346 0.119  1.310006e-08       4
#> FLNA1            4.716456e-12 0.619605886 0.757 0.430  1.414937e-08       4
#> CCDC34           5.068481e-12 0.530647490 0.441 0.181  1.520544e-08       4
#> NFIB2            5.333762e-12 0.557144300 0.735 0.422  1.600129e-08       4
#> BSPRY            5.705361e-12 0.424679653 0.559 0.252  1.711608e-08       4
#> SAC3D11          5.727368e-12 0.521368206 0.529 0.247  1.718211e-08       4
#> GRXCR11          6.128133e-12 0.456016009 0.272 0.077  1.838440e-08       4
#> C6orf151         6.637799e-12 0.436258554 0.272 0.079  1.991340e-08       4
#> PAK1IP11         9.238823e-12 0.663612037 0.662 0.385  2.771647e-08       4
#> TBC1D71          9.577128e-12 0.582466970 0.625 0.336  2.873138e-08       4
#> NR2F22           1.000748e-11 0.571438271 0.728 0.429  3.002244e-08       4
#> PIR2             1.012605e-11 0.614887856 0.574 0.316  3.037816e-08       4
#> KLHDC33          1.031445e-11 0.629293647 0.765 0.550  3.094336e-08       4
#> SMC4             1.161795e-11 0.678001807 0.581 0.351  3.485384e-08       4
#> FOXM1            1.171622e-11 0.550297403 0.316 0.104  3.514867e-08       4
#> RECQL4           1.219919e-11 0.476520417 0.331 0.110  3.659758e-08       4
#> CD3201           1.245336e-11 0.497417072 0.566 0.270  3.736007e-08       4
#> THEM62           1.290028e-11 0.518386615 0.574 0.270  3.870085e-08       4
#> E2F1             1.349070e-11 0.431125386 0.265 0.076  4.047210e-08       4
#> SMC2             1.401564e-11 0.580905945 0.544 0.276  4.204692e-08       4
#> MAP3K81          1.559895e-11 0.569021873 0.566 0.333  4.679684e-08       4
#> TPBG1            1.565607e-11 0.513875176 0.640 0.345  4.696822e-08       4
#> FKBP91           1.630589e-11 0.529601127 0.640 0.399  4.891768e-08       4
#> ACTR3B1          2.082285e-11 0.531833504 0.574 0.273  6.246855e-08       4
#> CHPF1            2.202118e-11 0.545495590 0.699 0.413  6.606354e-08       4
#> SLC5A6           2.320969e-11 0.508882702 0.390 0.148  6.962906e-08       4
#> EGLN31           2.448586e-11 0.426188853 0.316 0.104  7.345758e-08       4
#> EFNA12           2.588200e-11 0.599988524 0.750 0.484  7.764600e-08       4
#> NUSAP1           2.640235e-11 0.420728060 0.353 0.124  7.920704e-08       4
#> COTL11           3.293352e-11 0.614752443 0.676 0.478  9.880057e-08       4
#> ARPC42           3.490420e-11 0.707398993 0.647 0.441  1.047126e-07       4
#> SLC2A121         3.505364e-11 0.432396159 0.338 0.116  1.051609e-07       4
#> KRT6B1           3.528681e-11 0.451521952 0.265 0.079  1.058604e-07       4
#> EDN11            3.627912e-11 0.567971349 0.353 0.133  1.088374e-07       4
#> JAGN11           3.949974e-11 0.698830578 0.669 0.447  1.184992e-07       4
#> C21orf58         4.900476e-11 0.637186152 0.287 0.094  1.470143e-07       4
#> ADHFE11          5.595279e-11 0.687241382 0.434 0.197  1.678584e-07       4
#> NBL11            6.491244e-11 0.396920075 0.625 0.371  1.947373e-07       4
#> DHFR             6.569625e-11 0.407571324 0.390 0.147  1.970887e-07       4
#> FBXL161          7.857506e-11 0.318065721 0.375 0.133  2.357252e-07       4
#> LINC001521       8.514948e-11 0.517099288 0.676 0.395  2.554484e-07       4
#> NDC80            8.895851e-11 0.612703629 0.265 0.084  2.668755e-07       4
#> TROAP            9.061441e-11 0.529198586 0.265 0.082  2.718432e-07       4
#> TPX2             9.148973e-11 0.636834619 0.346 0.135  2.744692e-07       4
#> RDH102           1.001556e-10 0.585621826 0.618 0.381  3.004668e-07       4
#> CTD-3065J16.91   1.276757e-10 0.540115859 0.456 0.201  3.830270e-07       4
#> MYO5B1           1.276923e-10 0.478090248 0.574 0.291  3.830769e-07       4
#> TK1              1.349251e-10 0.267286448 0.404 0.191  4.047753e-07       4
#> PSIP1            1.366459e-10 0.514907722 0.581 0.320  4.099378e-07       4
#> RUVBL11          1.412427e-10 0.531117678 0.640 0.417  4.237282e-07       4
#> RMI21            1.449265e-10 0.534241600 0.471 0.223  4.347796e-07       4
#> PRR15L1          1.499948e-10 0.335945632 0.434 0.174  4.499844e-07       4
#> CENPH            1.562295e-10 0.507137478 0.382 0.152  4.686884e-07       4
#> MEX3A2           1.564805e-10 0.501727858 0.596 0.325  4.694415e-07       4
#> PCNA             1.779605e-10 0.581050221 0.588 0.390  5.338814e-07       4
#> ACOT7            1.834643e-10 0.427709931 0.507 0.226  5.503929e-07       4
#> DSCC1            1.865892e-10 0.503134332 0.279 0.092  5.597676e-07       4
#> PRSS83           1.995867e-10 0.590792149 0.743 0.491  5.987602e-07       4
#> RARRES12         2.086563e-10 0.566018143 0.426 0.255  6.259689e-07       4
#> IGSF31           2.127784e-10 0.545472824 0.662 0.422  6.383351e-07       4
#> CHMP2B1          2.160520e-10 0.693199184 0.654 0.458  6.481560e-07       4
#> MAFK1            2.216874e-10 0.443737525 0.419 0.171  6.650623e-07       4
#> MCM3             2.273758e-10 0.659810537 0.574 0.371  6.821273e-07       4
#> COL9A2           2.401458e-10 0.416753261 0.287 0.093  7.204374e-07       4
#> CCNB1            2.459881e-10 0.529801758 0.375 0.157  7.379644e-07       4
#> SLC25A51         2.601777e-10 0.513654231 0.713 0.465  7.805332e-07       4
#> IRF2BP21         2.624308e-10 0.589494738 0.706 0.460  7.872925e-07       4
#> TRIB21           2.646116e-10 0.466264933 0.485 0.217  7.938349e-07       4
#> SCCPDH2          2.922643e-10 0.587540510 0.551 0.311  8.767928e-07       4
#> ZG16B2           3.041284e-10 0.552914899 0.581 0.339  9.123853e-07       4
#> C1S2             3.185422e-10 0.516039474 0.588 0.351  9.556266e-07       4
#> KLF131           3.217550e-10 0.448847509 0.522 0.248  9.652649e-07       4
#> GMNN             3.321344e-10 0.478238048 0.478 0.248  9.964031e-07       4
#> MKI67            3.392623e-10 0.492329293 0.324 0.120  1.017787e-06       4
#> NUDT1            3.443013e-10 0.519381487 0.596 0.355  1.032904e-06       4
#> RPA3             4.363575e-10 0.593225899 0.699 0.458  1.309073e-06       4
#> MELK             4.639152e-10 0.511360917 0.301 0.107  1.391746e-06       4
#> PRR71            4.693612e-10 0.443971540 0.449 0.196  1.408084e-06       4
#> IFT1722          4.695532e-10 0.528190104 0.522 0.290  1.408660e-06       4
#> CYR611           4.811701e-10 0.512912178 0.654 0.395  1.443510e-06       4
#> RNASEH2A         4.938154e-10 0.417752716 0.544 0.250  1.481446e-06       4
#> CHAF1A           5.436734e-10 0.434122029 0.360 0.139  1.631020e-06       4
#> CYP27A11         5.513139e-10 0.258312771 0.360 0.134  1.653942e-06       4
#> PKP11            5.920912e-10 0.311956926 0.338 0.122  1.776274e-06       4
#> EGFL71           6.256990e-10 0.302603820 0.324 0.116  1.877097e-06       4
#> PARPBP           8.426461e-10 0.508088140 0.279 0.097  2.527938e-06       4
#> PRKDC1           8.467863e-10 0.557648011 0.713 0.486  2.540359e-06       4
#> PFN21            8.658471e-10 0.568880565 0.713 0.484  2.597541e-06       4
#> OBP2B            8.932846e-10 0.412207903 0.309 0.109  2.679854e-06       4
#> CGGBP11          1.055859e-09 0.683330707 0.603 0.430  3.167578e-06       4
#> ATP6V0E21        1.155793e-09 0.465785646 0.566 0.288  3.467380e-06       4
#> TADA31           1.173898e-09 0.656743954 0.632 0.442  3.521694e-06       4
#> TINAGL11         1.180480e-09 0.444819351 0.338 0.130  3.541440e-06       4
#> DCTPP11          1.236825e-09 0.545973748 0.721 0.475  3.710475e-06       4
#> C6orf1322        1.393842e-09 0.513827146 0.647 0.435  4.181527e-06       4
#> LRP21            1.493965e-09 0.478208660 0.463 0.222  4.481896e-06       4
#> TRIB12           1.949637e-09 0.630563137 0.566 0.391  5.848912e-06       4
#> GATA61           2.258342e-09 0.342744218 0.434 0.183  6.775026e-06       4
#> PKMYT1           2.272721e-09 0.451460120 0.265 0.089  6.818163e-06       4
#> CCNB2            2.281046e-09 0.649358146 0.353 0.174  6.843139e-06       4
#> KRT171           3.539176e-09 0.368626930 0.478 0.301  1.061753e-05       4
#> MAP21            3.658010e-09 0.382404568 0.434 0.190  1.097403e-05       4
#> GAN1             4.085653e-09 0.477658938 0.515 0.274  1.225696e-05       4
#> ASNS1            4.783434e-09 0.359516673 0.463 0.217  1.435030e-05       4
#> PDZD21           4.876643e-09 0.375820950 0.287 0.103  1.462993e-05       4
#> MFI21            5.045916e-09 0.462618187 0.574 0.307  1.513775e-05       4
#> RCAN2            5.453805e-09 0.362079524 0.265 0.092  1.636141e-05       4
#> HELLS            5.790630e-09 0.518437714 0.272 0.098  1.737189e-05       4
#> GGH1             5.871140e-09 0.407063588 0.581 0.373  1.761342e-05       4
#> HSPA52           6.376891e-09 0.504022496 0.735 0.481  1.913067e-05       4
#> SHISA91          7.052756e-09 0.383703180 0.397 0.169  2.115827e-05       4
#> ANLN             7.246262e-09 0.495036702 0.265 0.095  2.173879e-05       4
#> MYC              7.274208e-09 0.461386019 0.515 0.297  2.182262e-05       4
#> PRSS223          7.420206e-09 0.460744666 0.566 0.368  2.226062e-05       4
#> TUBB61           9.370556e-09 0.496423455 0.691 0.446  2.811167e-05       4
#> PSTPIP21         1.094992e-08 0.525614216 0.529 0.331  3.284975e-05       4
#> PCDHB91          1.113893e-08 0.382031543 0.449 0.207  3.341678e-05       4
#> PTTG1            1.207573e-08 0.583195142 0.368 0.212  3.622719e-05       4
#> H2AFZ            1.297573e-08 0.599453350 0.625 0.436  3.892719e-05       4
#> MTHFD1           1.306653e-08 0.438056530 0.434 0.206  3.919958e-05       4
#> SPON22           1.334907e-08 0.525190190 0.515 0.324  4.004722e-05       4
#> PODXL22          1.414824e-08 0.451653063 0.669 0.434  4.244473e-05       4
#> KRT73            1.455366e-08 0.604816531 0.824 0.593  4.366099e-05       4
#> CLMN1            1.571305e-08 0.541894213 0.603 0.424  4.713916e-05       4
#> HOMER21          1.728414e-08 0.440727265 0.426 0.202  5.185243e-05       4
#> DTYMK            1.813572e-08 0.396652781 0.566 0.345  5.440715e-05       4
#> MCM2             1.825029e-08 0.392602426 0.309 0.121  5.475086e-05       4
#> RCAN11           2.033844e-08 0.456344450 0.544 0.330  6.101532e-05       4
#> ADAM152          2.993289e-08 0.487681509 0.691 0.486  8.979866e-05       4
#> TNFRSF18         3.113750e-08 0.295487469 0.294 0.112  9.341249e-05       4
#> LINC016151       3.146482e-08 0.361315111 0.287 0.109  9.439447e-05       4
#> FANCI            3.686509e-08 0.463600458 0.265 0.100  1.105953e-04       4
#> SAA12            4.097909e-08 0.435896788 0.522 0.298  1.229373e-04       4
#> KCNN4            4.572693e-08 0.417678297 0.500 0.275  1.371808e-04       4
#> FEN1             4.764339e-08 0.319495046 0.426 0.196  1.429302e-04       4
#> SPC25            5.197571e-08 0.348549140 0.250 0.088  1.559271e-04       4
#> FADS1            6.447408e-08 0.257203931 0.382 0.168  1.934222e-04       4
#> NDUFAF61         6.474787e-08 0.448628983 0.632 0.433  1.942436e-04       4
#> RAD212           6.525812e-08 0.530106852 0.676 0.521  1.957744e-04       4
#> SCD2             6.925071e-08 0.442389546 0.581 0.385  2.077521e-04       4
#> DBNDD12          7.183482e-08 0.435792547 0.537 0.330  2.155044e-04       4
#> NANOS12          1.014883e-07 0.412170110 0.449 0.246  3.044648e-04       4
#> RPP25            1.041136e-07 0.459686974 0.566 0.355  3.123409e-04       4
#> LAPTM4B2         1.057912e-07 0.479430626 0.787 0.522  3.173736e-04       4
#> CLCN42           1.089839e-07 0.391930351 0.566 0.334  3.269516e-04       4
#> PAICS1           1.172660e-07 0.433110480 0.625 0.404  3.517979e-04       4
#> ABLIM11          1.227638e-07 0.416465627 0.257 0.098  3.682914e-04       4
#> MAL22            1.272437e-07 0.441352790 0.684 0.473  3.817312e-04       4
#> RASD11           1.298498e-07 0.487751351 0.368 0.180  3.895493e-04       4
#> MINCR1           1.337345e-07 0.461113355 0.507 0.263  4.012035e-04       4
#> PHLDA21          1.394100e-07 0.370590722 0.507 0.298  4.182299e-04       4
#> TMPO             1.411529e-07 0.417388814 0.544 0.345  4.234588e-04       4
#> CDT1             1.440011e-07 0.411608344 0.301 0.128  4.320033e-04       4
#> RAB11FIP11       1.597227e-07 0.441362346 0.507 0.303  4.791681e-04       4
#> NPR32            1.670288e-07 0.324613788 0.434 0.213  5.010865e-04       4
#> LRR1             1.737537e-07 0.348946131 0.368 0.161  5.212612e-04       4
#> PHF19            1.746332e-07 0.508735260 0.507 0.308  5.238996e-04       4
#> SCGB3A12         1.870710e-07 0.323118384 0.544 0.322  5.612130e-04       4
#> YBX21            1.894597e-07 0.431119863 0.301 0.127  5.683791e-04       4
#> CRABP11          1.988763e-07 0.408313177 0.581 0.388  5.966289e-04       4
#> ZWINT            2.141629e-07 0.335691687 0.338 0.152  6.424886e-04       4
#> MCM6             2.167414e-07 0.393468826 0.257 0.101  6.502243e-04       4
#> MDC1             2.180526e-07 0.424546293 0.265 0.106  6.541577e-04       4
#> VANGL12          2.183617e-07 0.510387696 0.721 0.511  6.550852e-04       4
#> ATAD21           2.278785e-07 0.445383333 0.471 0.258  6.836356e-04       4
#> ERBB32           2.300855e-07 0.380629174 0.625 0.428  6.902565e-04       4
#> LMNB2            2.391179e-07 0.322566510 0.353 0.157  7.173538e-04       4
#> CDCA7L1          2.392766e-07 0.424473813 0.493 0.263  7.178298e-04       4
#> C31              2.492413e-07 0.431226354 0.419 0.228  7.477238e-04       4
#> CELF42           2.776005e-07 0.402091477 0.566 0.316  8.328016e-04       4
#> CCDC64B1         2.785703e-07 0.373961558 0.647 0.417  8.357108e-04       4
#> TRIB32           3.162721e-07 0.411431877 0.515 0.290  9.488163e-04       4
#> SAA21            3.261419e-07 0.503626782 0.294 0.132  9.784258e-04       4
#> SORBS23          3.463664e-07 0.473718246 0.640 0.476  1.039099e-03       4
#> SGOL1            3.769757e-07 0.324259758 0.257 0.101  1.130927e-03       4
#> CXADR1           3.960830e-07 0.432470600 0.500 0.354  1.188249e-03       4
#> LOXL11           4.264935e-07 0.102310983 0.404 0.189  1.279480e-03       4
#> C1GALT11         4.658457e-07 0.407226769 0.559 0.394  1.397537e-03       4
#> CSTB2            4.754732e-07 0.458845313 0.640 0.468  1.426420e-03       4
#> TCF19            4.768922e-07 0.392127974 0.375 0.181  1.430677e-03       4
#> RP11-400K9.41    4.772759e-07 0.450321517 0.250 0.100  1.431828e-03       4
#> BASP11           4.826883e-07 0.190661292 0.471 0.282  1.448065e-03       4
#> CD821            5.528718e-07 0.420991175 0.574 0.411  1.658615e-03       4
#> BID1             5.547907e-07 0.348079739 0.588 0.371  1.664372e-03       4
#> TFAP2A2          5.605894e-07 0.379232848 0.603 0.418  1.681768e-03       4
#> ELN2             5.684726e-07 0.489375869 0.353 0.175  1.705418e-03       4
#> KRT183           5.815261e-07 0.496036791 0.779 0.522  1.744578e-03       4
#> CEP701           6.030625e-07 0.439107491 0.441 0.229  1.809187e-03       4
#> EHD11            6.678805e-07 0.331797536 0.493 0.262  2.003642e-03       4
#> CDK1             7.799502e-07 0.262819578 0.294 0.134  2.339850e-03       4
#> SELENBP13        8.497962e-07 0.368627148 0.537 0.326  2.549389e-03       4
#> PLXNA21          8.913549e-07 0.445999568 0.471 0.261  2.674065e-03       4
#> C19orf48         9.037026e-07 0.386299212 0.507 0.281  2.711108e-03       4
#> PHYH2            1.170099e-06 0.431472407 0.706 0.479  3.510297e-03       4
#> SOX42            1.245248e-06 0.485204951 0.787 0.550  3.735744e-03       4
#> TNF              1.462630e-06 0.322061275 0.279 0.121  4.387889e-03       4
#> CAMK1D1          1.479398e-06 0.325006437 0.404 0.211  4.438195e-03       4
#> BACE23           1.586872e-06 0.429241499 0.647 0.448  4.760617e-03       4
#> SOX92            1.670983e-06 0.321128940 0.588 0.405  5.012950e-03       4
#> PLEKHB12         1.671871e-06 0.400589264 0.676 0.467  5.015612e-03       4
#> RRS11            1.867677e-06 0.362504419 0.566 0.349  5.603031e-03       4
#> TP53BP21         1.880389e-06 0.381324235 0.588 0.425  5.641167e-03       4
#> TINCR2           1.916680e-06 0.229063549 0.500 0.305  5.750041e-03       4
#> COL4A22          1.955346e-06 0.091127333 0.537 0.320  5.866039e-03       4
#> PDGFA1           2.026895e-06 0.271875684 0.375 0.184  6.080685e-03       4
#> ARID5A1          2.067143e-06 0.352094725 0.493 0.292  6.201430e-03       4
#> IRF2BPL1         2.135051e-06 0.389030629 0.640 0.444  6.405154e-03       4
#> RTP41            2.499780e-06 0.370522749 0.456 0.249  7.499340e-03       4
#> KRT161           2.726753e-06 0.338786079 0.294 0.135  8.180258e-03       4
#> MYO62            2.977945e-06 0.381803862 0.706 0.448  8.933835e-03       4
#> FADS31           3.136586e-06 0.275501263 0.493 0.267  9.409757e-03       4
#> HIST1H2BJ1       3.168988e-06 0.443129704 0.368 0.191  9.506963e-03       4
#> ACSL12           3.183543e-06 0.195783536 0.485 0.242  9.550630e-03       4
#> MSLN1            3.325840e-06 0.290843222 0.368 0.185  9.977519e-03       4
#> KIAA05131        3.845472e-06 0.185601967 0.346 0.159  1.153642e-02       4
#> CDKN2D           3.960452e-06 0.257457802 0.456 0.233  1.188136e-02       4
#> CENPM            5.408139e-06 0.177422174 0.272 0.118  1.622442e-02       4
#> CACYBP2          6.814623e-06 0.410752302 0.669 0.519  2.044387e-02       4
#> MAP2K31          6.987221e-06 0.352956700 0.588 0.369  2.096166e-02       4
#> TMEM971          7.099025e-06 0.264816476 0.449 0.238  2.129708e-02       4
#> THBS12           7.576667e-06 0.278759031 0.522 0.329  2.273000e-02       4
#> ADAMTS91         8.685817e-06 0.269748442 0.279 0.127  2.605745e-02       4
#> PEG102           8.818414e-06 0.362028623 0.691 0.461  2.645524e-02       4
#> RUNX31           9.303674e-06 0.044347586 0.287 0.122  2.791102e-02       4
#> KIF20B           9.410181e-06 0.334576796 0.331 0.170  2.823054e-02       4
#> SIGMAR11         1.020411e-05 0.330878040 0.515 0.375  3.061233e-02       4
#> RRM1             1.032003e-05 0.272172495 0.471 0.274  3.096008e-02       4
#> ARFGEF31         1.032654e-05 0.274398558 0.331 0.159  3.097963e-02       4
#> PTPRS1           1.097659e-05 0.252306083 0.375 0.191  3.292978e-02       4
#> CTSF1            1.211647e-05 0.350228675 0.654 0.442  3.634941e-02       4
#> MTSS1L1          1.303229e-05 0.335624386 0.544 0.388  3.909688e-02       4
#> HIST1H2BN1       1.334166e-05 0.294209830 0.353 0.179  4.002499e-02       4
#> CAMK12           1.363239e-05 0.159331667 0.324 0.152  4.089716e-02       4
#> MBP1             1.456437e-05 0.347225853 0.581 0.412  4.369312e-02       4
#> MCM4             1.499955e-05 0.318766471 0.500 0.338  4.499865e-02       4
#> CDCA4            1.526920e-05 0.335888939 0.419 0.232  4.580761e-02       4
#> UBE2S            1.605859e-05 0.420504139 0.522 0.394  4.817577e-02       4
#> ST6GAL11         1.706812e-05 0.206191448 0.441 0.234  5.120436e-02       4
#> SLPI2            1.833459e-05 0.242955090 0.471 0.323  5.500378e-02       4
#> HES12            1.880771e-05 0.339415014 0.588 0.438  5.642312e-02       4
#> WDR342           1.914646e-05 0.381670491 0.625 0.488  5.743937e-02       4
#> GRB142           2.013242e-05 0.312639446 0.390 0.214  6.039726e-02       4
#> TUBB2            2.206867e-05 0.396737538 0.691 0.508  6.620600e-02       4
#> C21              2.441638e-05 0.085119301 0.382 0.220  7.324915e-02       4
#> SULF21           2.512780e-05 0.195326021 0.485 0.318  7.538341e-02       4
#> SLC2A4RG2        2.684225e-05 0.358494117 0.581 0.425  8.052674e-02       4
#> RFC3             3.011355e-05 0.276153619 0.456 0.266  9.034065e-02       4
#> CKAP5            3.365978e-05 0.334891563 0.360 0.197  1.009793e-01       4
#> SPDL1            3.956824e-05 0.441656866 0.257 0.128  1.187047e-01       4
#> TPT1-AS11        4.214063e-05 0.259073761 0.537 0.385  1.264219e-01       4
#> HLA-B2           4.501782e-05 0.414360275 0.654 0.479  1.350535e-01       4
#> SMOC21           4.914704e-05 0.179721925 0.375 0.209  1.474411e-01       4
#> TTC39A1          5.279642e-05 0.229776398 0.434 0.250  1.583893e-01       4
#> LIMD21           6.877479e-05 0.063694209 0.390 0.220  2.063244e-01       4
#> RACGAP1          6.908581e-05 0.294496506 0.301 0.159  2.072574e-01       4
#> NEBL1            6.994098e-05 0.187323197 0.463 0.289  2.098229e-01       4
#> CKB2             7.459298e-05 0.340914140 0.596 0.436  2.237789e-01       4
#> JHDM1D-AS11      7.987090e-05 0.281124411 0.485 0.301  2.396127e-01       4
#> HIST1H4C1        7.988312e-05 0.293669517 0.441 0.331  2.396493e-01       4
#> LYAR1            8.132217e-05 0.267774269 0.493 0.349  2.439665e-01       4
#> MCM71            8.678070e-05 0.186633371 0.507 0.376  2.603421e-01       4
#> OAF              9.476658e-05 0.114499649 0.375 0.199  2.842997e-01       4
#> NASP             1.043691e-04 0.279678899 0.515 0.390  3.131073e-01       4
#> CCT52            1.081620e-04 0.357256829 0.676 0.512  3.244860e-01       4
#> VASN3            1.110140e-04 0.298436624 0.574 0.453  3.330419e-01       4
#> EFHD11           1.130112e-04 0.214493746 0.596 0.409  3.390337e-01       4
#> THOP1            1.189617e-04 0.052505134 0.368 0.198  3.568850e-01       4
#> MTHFD21          1.240981e-04 0.345425177 0.537 0.405  3.722942e-01       4
#> PHACTR11         1.263437e-04 0.183563644 0.412 0.273  3.790312e-01       4
#> GSDMC2           1.264504e-04 0.201199441 0.404 0.230  3.793513e-01       4
#> BOP11            1.365451e-04 0.279007571 0.603 0.438  4.096353e-01       4
#> CRLF11           1.430162e-04 0.213785109 0.257 0.130  4.290486e-01       4
#> MAPK131          1.461620e-04 0.306815975 0.662 0.480  4.384861e-01       4
#> ZNF6541          1.513209e-04 0.325602667 0.353 0.204  4.539626e-01       4
#> S100A81          1.521621e-04 0.127792893 0.397 0.324  4.564862e-01       4
#> PSRC1            1.523133e-04 0.409421951 0.287 0.157  4.569398e-01       4
#> MB2              1.568801e-04 0.258339097 0.419 0.260  4.706402e-01       4
#> FRZB1            1.642255e-04 0.163349841 0.375 0.246  4.926765e-01       4
#> CORO1A1          1.668180e-04 0.168662416 0.441 0.303  5.004539e-01       4
#> HOPX1            1.692831e-04 0.123599149 0.250 0.122  5.078494e-01       4
#> CDK42            1.724316e-04 0.370700117 0.654 0.509  5.172948e-01       4
#> INHBB1           1.807188e-04 0.229246850 0.279 0.144  5.421564e-01       4
#> HLA-F2           2.044859e-04 0.239029612 0.471 0.339  6.134576e-01       4
#> MARC1            2.272137e-04 0.239463921 0.257 0.132  6.816411e-01       4
#> SPHK11           2.386126e-04 0.353974893 0.529 0.425  7.158378e-01       4
#> KIF221           2.552272e-04 0.236266345 0.551 0.426  7.656816e-01       4
#> NTHL11           2.935814e-04 0.230403939 0.581 0.401  8.807443e-01       4
#> TOMM401          3.060310e-04 0.276866677 0.610 0.458  9.180931e-01       4
#> IL17RE1          3.087742e-04 0.235372546 0.441 0.302  9.263227e-01       4
#> C10orf102        3.395491e-04 0.281034727 0.471 0.387  1.000000e+00       4
#> NEAT12           3.762115e-04 0.164025381 0.676 0.524  1.000000e+00       4
#> NDRG12           3.814369e-04 0.363818648 0.551 0.436  1.000000e+00       4
#> DHTKD11          3.921992e-04 0.299954383 0.529 0.415  1.000000e+00       4
#> DONSON           3.924397e-04 0.177144532 0.272 0.142  1.000000e+00       4
#> MGLL2            5.332698e-04 0.127437650 0.390 0.234  1.000000e+00       4
#> C1R2             5.930368e-04 0.060070694 0.463 0.357  1.000000e+00       4
#> SLC39A141        5.943891e-04 0.058652873 0.301 0.158  1.000000e+00       4
#> NCAPD2           6.584102e-04 0.290126561 0.257 0.143  1.000000e+00       4
#> PDIA42           6.755783e-04 0.308849133 0.632 0.494  1.000000e+00       4
#> RFC4             6.972325e-04 0.259632917 0.382 0.235  1.000000e+00       4
#> PRNP1            7.868913e-04 0.350682048 0.618 0.516  1.000000e+00       4
#> TPD52L12         8.106673e-04 0.256791337 0.596 0.460  1.000000e+00       4
#> TPM13            8.852762e-04 0.360489672 0.699 0.543  1.000000e+00       4
#> DKC11            9.575053e-04 0.274778950 0.551 0.478  1.000000e+00       4
#> FBLN23           9.743896e-04 0.301722430 0.654 0.513  1.000000e+00       4
#> GAPDH2           1.001624e-03 0.285575642 0.662 0.490  1.000000e+00       4
#> PDP11            1.005870e-03 0.185372588 0.493 0.361  1.000000e+00       4
#> PDZK1IP12        1.034149e-03 0.202757455 0.551 0.458  1.000000e+00       4
#> PHLDA31          1.256322e-03 0.282612950 0.419 0.345  1.000000e+00       4
#> HSP90AB13        1.286683e-03 0.485415971 0.824 0.651  1.000000e+00       4
#> ZNF528           1.558797e-03 0.169101217 0.316 0.182  1.000000e+00       4
#> CRELD21          1.587024e-03 0.305194034 0.566 0.493  1.000000e+00       4
#> SPEG1            1.595608e-03 0.269681981 0.346 0.222  1.000000e+00       4
#> FSCN11           1.678962e-03 0.202360001 0.441 0.308  1.000000e+00       4
#> RPL39L1          1.995661e-03 0.221695229 0.603 0.454  1.000000e+00       4
#> ALYREF1          2.022361e-03 0.194920437 0.471 0.341  1.000000e+00       4
#> PLSCR12          2.031751e-03 0.249639131 0.537 0.427  1.000000e+00       4
#> NAALADL21        2.129544e-03 0.089436607 0.412 0.252  1.000000e+00       4
#> AIF1L2           2.198051e-03 0.225065995 0.559 0.457  1.000000e+00       4
#> LBR1             2.307208e-03 0.257362420 0.419 0.354  1.000000e+00       4
#> RAB3IP2          2.570663e-03 0.144632714 0.426 0.297  1.000000e+00       4
#> GAS11            2.751679e-03 0.139585306 0.456 0.323  1.000000e+00       4
#> CXCL161          2.980621e-03 0.189368533 0.559 0.422  1.000000e+00       4
#> C3orf381         3.448373e-03 0.253050592 0.449 0.373  1.000000e+00       4
#> S100A62          3.471993e-03 0.229731245 0.588 0.441  1.000000e+00       4
#> GPSM21           3.651370e-03 0.149096246 0.294 0.178  1.000000e+00       4
#> CASP41           3.782271e-03 0.061211969 0.390 0.284  1.000000e+00       4
#> CREB51           4.002410e-03 0.178033444 0.346 0.226  1.000000e+00       4
#> EPAS1            4.026441e-03 0.004917717 0.294 0.170  1.000000e+00       4
#> MBNL1-AS12       4.361202e-03 0.237513631 0.412 0.302  1.000000e+00       4
#> PPIL11           4.732447e-03 0.195110215 0.471 0.368  1.000000e+00       4
#> B4GALT11         4.781113e-03 0.213286716 0.485 0.419  1.000000e+00       4
#> TUBB4B2          4.846298e-03 0.258198388 0.625 0.502  1.000000e+00       4
#> ITGA61           4.988437e-03 0.092619046 0.397 0.286  1.000000e+00       4
#> MARS             5.631080e-03 0.173787203 0.463 0.328  1.000000e+00       4
#> CTSV1            5.849255e-03 0.092441002 0.426 0.320  1.000000e+00       4
#> CTNNBIP11        6.001454e-03 0.162846948 0.551 0.476  1.000000e+00       4
#> SAP301           6.140305e-03 0.212395749 0.544 0.430  1.000000e+00       4
#> CHEK1            6.246778e-03 0.154637184 0.309 0.195  1.000000e+00       4
#> CKAP2            6.740986e-03 0.244806744 0.338 0.233  1.000000e+00       4
#> FGF11            7.848543e-03 0.106313046 0.272 0.166  1.000000e+00       4
#> BRCA2            8.010998e-03 0.119220953 0.250 0.152  1.000000e+00       4
#> CLDN42           8.206319e-03 0.316503835 0.699 0.556  1.000000e+00       4
#> C8orf462         9.378660e-03 0.151545333 0.397 0.300  1.000000e+00       4
#> UNG1             9.934988e-03 0.174562267 0.375 0.261  1.000000e+00       4
#> TOX2             3.680723e-46 2.080470994 0.955 0.226  1.104217e-42       5
#> SLC26A71         1.205001e-35 1.652019199 0.612 0.100  3.615003e-32       5
#> SBSPON2          9.894931e-31 1.642054373 0.940 0.279  2.968479e-27       5
#> COL11A22         3.837653e-30 1.594111928 0.761 0.195  1.151296e-26       5
#> S100P2           5.688115e-30 1.534660450 1.000 0.435  1.706435e-26       5
#> PLA2G4A2         1.483281e-29 1.550689678 0.866 0.262  4.449842e-26       5
#> CNGA1            2.154795e-29 1.561213704 0.582 0.111  6.464384e-26       5
#> DNM32            2.629300e-29 1.516525654 0.836 0.242  7.887901e-26       5
#> MATN32           6.291119e-28 1.447301390 0.910 0.285  1.887336e-24       5
#> SFN2             1.429921e-27 1.467634783 0.970 0.342  4.289764e-24       5
#> SOD32            1.641080e-27 1.497813845 1.000 0.382  4.923239e-24       5
#> FXYD62           1.906874e-26 1.428267928 0.970 0.467  5.720621e-23       5
#> RP1-27K12.2      2.324454e-26 1.197074272 0.463 0.070  6.973362e-23       5
#> EPHX12           2.349094e-26 1.449835396 0.955 0.440  7.047281e-23       5
#> C2orf80          3.216363e-26 1.519821617 0.493 0.086  9.649090e-23       5
#> CDKN2A1          3.649712e-26 1.657221852 0.851 0.319  1.094914e-22       5
#> MDFI2            1.799901e-25 1.396582011 0.985 0.426  5.399704e-22       5
#> RAB253           2.647019e-25 1.559069563 0.866 0.370  7.941057e-22       5
#> RP11-25K19.12    4.497353e-25 1.568803988 0.672 0.186  1.349206e-21       5
#> CA82             8.698705e-25 1.228660643 0.806 0.243  2.609612e-21       5
#> SLC29A12         1.114449e-24 1.240447485 1.000 0.487  3.343346e-21       5
#> BOC2             2.153510e-24 1.376203879 0.746 0.220  6.460529e-21       5
#> S100A12          3.047897e-24 1.210402644 1.000 0.483  9.143690e-21       5
#> HIST1H2AE2       1.212195e-23 1.339097953 0.896 0.382  3.636585e-20       5
#> ISLR             1.226826e-23 1.186394236 0.687 0.186  3.680477e-20       5
#> ROPN12           1.754478e-23 1.314987086 0.881 0.298  5.263434e-20       5
#> GJA11            1.800010e-22 1.130574448 0.627 0.160  5.400029e-19       5
#> CRISPLD12        1.916455e-22 1.263618990 0.955 0.423  5.749364e-19       5
#> SCRG12           1.822931e-21 1.222428854 0.970 0.351  5.468792e-18       5
#> MYOZ11           2.133262e-21 1.374134967 0.746 0.247  6.399787e-18       5
#> GLS2             2.966593e-21 1.238621863 1.000 0.432  8.899780e-18       5
#> TSPAN122         5.510522e-21 1.062775081 0.851 0.317  1.653157e-17       5
#> TMEM793          5.773542e-21 1.388977876 0.866 0.368  1.732063e-17       5
#> C1orf1162        5.996263e-21 1.271181668 0.716 0.228  1.798879e-17       5
#> DBI2             1.146797e-20 1.127263019 0.985 0.479  3.440391e-17       5
#> LIPH1            1.443679e-20 1.234825924 0.716 0.243  4.331036e-17       5
#> CCDC64B2         2.393490e-20 1.325929177 0.866 0.418  7.180469e-17       5
#> GPM6B2           2.522893e-20 1.209081865 0.940 0.414  7.568680e-17       5
#> COL11A12         2.642680e-20 1.068306383 0.866 0.356  7.928040e-17       5
#> TSPAN52          3.773651e-20 1.195699453 0.746 0.258  1.132095e-16       5
#> HIBCH2           5.064611e-20 1.268896345 0.910 0.399  1.519383e-16       5
#> MGST12           6.429089e-20 1.186382240 0.970 0.436  1.928727e-16       5
#> DMKN             9.504377e-20 0.841623611 0.627 0.172  2.851313e-16       5
#> CLDN33           1.114873e-19 1.118017839 0.955 0.506  3.344619e-16       5
#> SH3BGR2          1.299244e-19 1.240261808 0.806 0.318  3.897732e-16       5
#> MARCKSL12        1.432590e-19 0.992627386 0.985 0.554  4.297771e-16       5
#> FBXO22           1.449024e-19 1.317203052 0.851 0.338  4.347072e-16       5
#> CLU2             1.913215e-19 1.181329594 0.985 0.439  5.739644e-16       5
#> GCNT12           2.252449e-19 1.210935275 0.821 0.323  6.757346e-16       5
#> DIO2             2.871634e-19 0.793551835 0.299 0.040  8.614903e-16       5
#> HSPB12           3.664747e-19 1.094163703 0.970 0.471  1.099424e-15       5
#> TTYH12           5.992951e-19 1.163497109 0.896 0.415  1.797885e-15       5
#> AZGP12           8.319506e-19 1.205180794 0.881 0.348  2.495852e-15       5
#> SLC43A32         9.503841e-19 1.114205923 0.955 0.428  2.851152e-15       5
#> IGFBP21          1.484164e-18 1.102505106 0.940 0.381  4.452493e-15       5
#> HIST1H1C2        1.805641e-18 1.166008511 0.910 0.421  5.416922e-15       5
#> RP11-89K21.11    5.219441e-18 1.097697826 0.552 0.150  1.565832e-14       5
#> C1QL42           6.038589e-18 1.025345267 0.731 0.262  1.811577e-14       5
#> KLK12            7.266629e-18 1.040886699 0.507 0.124  2.179989e-14       5
#> NET13            1.087580e-17 1.132539847 0.940 0.453  3.262739e-14       5
#> CLDN43           1.192702e-17 1.018818041 0.985 0.547  3.578106e-14       5
#> ROPN1B2          1.493015e-17 1.091345199 0.955 0.414  4.479046e-14       5
#> PAX11            1.700398e-17 1.023619787 0.493 0.122  5.101194e-14       5
#> TUBB2B2          2.214755e-17 1.033431687 0.746 0.318  6.644266e-14       5
#> IDI12            2.734209e-17 1.102348316 0.896 0.437  8.202626e-14       5
#> HIST1H2BG2       3.933118e-17 1.014217268 0.940 0.378  1.179935e-13       5
#> GYLTL1B          4.717323e-17 1.276510533 0.478 0.129  1.415197e-13       5
#> DNM3OS2          5.213236e-17 0.978158912 0.687 0.225  1.563971e-13       5
#> SLC9A3R22        1.057028e-16 1.055771767 0.881 0.425  3.171085e-13       5
#> RP11-268P4.5     1.115620e-16 1.067419466 0.418 0.093  3.346860e-13       5
#> C1orf1861        1.119490e-16 1.070199994 0.821 0.384  3.358469e-13       5
#> METTL7A2         2.338630e-16 1.203480670 0.776 0.345  7.015889e-13       5
#> MOG              2.760754e-16 0.944509721 0.254 0.034  8.282262e-13       5
#> EXTL11           2.814167e-16 0.914114882 0.657 0.210  8.442501e-13       5
#> PEG33            1.248271e-15 1.015165816 0.881 0.415  3.744814e-12       5
#> TAGLN2           2.069709e-15 1.026412665 0.821 0.423  6.209127e-12       5
#> TIFA1            2.191886e-15 0.958275261 0.836 0.399  6.575658e-12       5
#> FXYD33           3.749756e-15 0.970144422 0.970 0.471  1.124927e-11       5
#> SNX222           5.927696e-15 0.980898023 0.716 0.295  1.778309e-11       5
#> PDZK1IP13        6.078090e-15 1.050849786 0.866 0.444  1.823427e-11       5
#> TUBB2A2          6.741770e-15 0.999123193 0.851 0.406  2.022531e-11       5
#> TUBA1A2          6.904459e-15 0.971031447 0.910 0.462  2.071338e-11       5
#> TOB12            7.176047e-15 0.965527627 0.896 0.437  2.152814e-11       5
#> GOLT1A1          1.005033e-14 0.968060650 0.433 0.112  3.015098e-11       5
#> C8orf463         1.129376e-14 0.982171183 0.716 0.285  3.388127e-11       5
#> BSPRY1           1.508937e-14 1.090596179 0.672 0.265  4.526812e-11       5
#> LIMCH12          1.875265e-14 0.934942228 0.806 0.381  5.625796e-11       5
#> ENHO             2.557034e-14 0.935197837 0.388 0.091  7.671101e-11       5
#> PLEKHB13         3.116129e-14 0.928212864 0.896 0.467  9.348386e-11       5
#> TUBB4B3          4.270779e-14 0.956878060 0.925 0.490  1.281234e-10       5
#> FHL12            7.928904e-14 0.945561919 0.552 0.183  2.378671e-10       5
#> NAAA             8.824938e-14 0.972924037 0.642 0.236  2.647481e-10       5
#> TSPAN22          9.142952e-14 0.755655274 0.537 0.164  2.742885e-10       5
#> MYL92            1.004270e-13 0.929905477 0.925 0.470  3.012809e-10       5
#> MFAP23           1.151647e-13 0.940362825 0.896 0.447  3.454942e-10       5
#> NGF              1.206329e-13 0.943510421 0.269 0.049  3.618988e-10       5
#> IDH12            1.288873e-13 0.928254516 0.896 0.384  3.866620e-10       5
#> FBN2             1.618408e-13 1.131028372 0.358 0.086  4.855223e-10       5
#> MIA2             1.685500e-13 0.896351383 0.746 0.382  5.056500e-10       5
#> PAM2             1.786880e-13 0.920308129 0.851 0.449  5.360641e-10       5
#> IGFBP51          1.812705e-13 0.659547469 0.687 0.250  5.438115e-10       5
#> SERTAD42         2.961043e-13 0.939221828 0.866 0.453  8.883128e-10       5
#> NCCRP12          3.350663e-13 1.116567991 0.612 0.238  1.005199e-09       5
#> VSNL11           4.118018e-13 0.818671044 0.403 0.104  1.235405e-09       5
#> PFN22            4.320655e-13 0.914961791 0.851 0.490  1.296196e-09       5
#> PRSS332          5.073337e-13 0.892994248 0.791 0.378  1.522001e-09       5
#> CYP26B11         5.545512e-13 0.874310961 0.448 0.130  1.663654e-09       5
#> HIST2H2BE2       6.005073e-13 0.881667634 0.851 0.400  1.801522e-09       5
#> NDRG22           6.562932e-13 0.933112641 0.836 0.427  1.968880e-09       5
#> MUC15            8.285546e-13 0.982693166 0.313 0.068  2.485664e-09       5
#> SMTN3            9.167995e-13 0.979609489 0.776 0.437  2.750398e-09       5
#> CTSV2            9.977741e-13 0.995140011 0.687 0.311  2.993322e-09       5
#> NPPC1            1.461178e-12 0.957378559 0.358 0.091  4.383535e-09       5
#> GOLM12           1.501136e-12 0.934852818 0.866 0.464  4.503409e-09       5
#> TACSTD23         1.546748e-12 0.870652947 0.881 0.452  4.640244e-09       5
#> PRR15L2          1.681106e-12 1.133634617 0.507 0.186  5.043318e-09       5
#> RIC32            1.883816e-12 0.841889469 0.746 0.375  5.651447e-09       5
#> BAMBI2           2.039325e-12 0.855415062 0.910 0.438  6.117976e-09       5
#> LINC01048        2.302518e-12 0.975324345 0.328 0.078  6.907553e-09       5
#> SDC13            2.672016e-12 0.834717003 0.896 0.510  8.016047e-09       5
#> MYLK2            2.716625e-12 0.898245920 0.791 0.422  8.149876e-09       5
#> S100A103         2.727879e-12 0.765486910 0.881 0.532  8.183637e-09       5
#> RTN4RL22         2.772223e-12 0.869746781 0.522 0.175  8.316669e-09       5
#> DTNB2            3.206385e-12 0.892795142 0.746 0.343  9.619156e-09       5
#> CAPN23           4.873088e-12 0.845350187 0.851 0.506  1.461926e-08       5
#> ERBB33           5.962831e-12 0.971106139 0.761 0.432  1.788849e-08       5
#> FAM84A1          6.686331e-12 0.786176945 0.448 0.137  2.005899e-08       5
#> LRRC731          7.738972e-12 0.740295580 0.582 0.213  2.321691e-08       5
#> PRSS84           9.207982e-12 0.856868525 0.896 0.498  2.762395e-08       5
#> ENPP52           1.143767e-11 0.737684958 0.761 0.325  3.431301e-08       5
#> SOHLH12          1.249237e-11 0.727281750 0.716 0.301  3.747712e-08       5
#> UCHL11           1.269953e-11 1.020689446 0.448 0.153  3.809858e-08       5
#> SYT82            1.284705e-11 0.727256708 0.731 0.355  3.854114e-08       5
#> KCNQ1OT13        1.664437e-11 0.824143160 0.896 0.470  4.993311e-08       5
#> IFRD12           1.886901e-11 0.782255561 0.821 0.420  5.660702e-08       5
#> HSPA1B2          2.054855e-11 0.819682945 0.881 0.444  6.164566e-08       5
#> WDR343           2.065494e-11 0.883514638 0.806 0.485  6.196483e-08       5
#> CITED42          2.558790e-11 0.805458007 0.836 0.465  7.676370e-08       5
#> PTGS21           2.670954e-11 0.522807160 0.537 0.178  8.012863e-08       5
#> ALDH1B12         3.225950e-11 0.887591256 0.701 0.332  9.677850e-08       5
#> MBNL1-AS13       4.332910e-11 0.632962967 0.701 0.290  1.299873e-07       5
#> SORBS24          4.343774e-11 0.793773988 0.851 0.473  1.303132e-07       5
#> PTGIS1           4.596538e-11 0.779066397 0.433 0.135  1.378961e-07       5
#> RBP4             4.666132e-11 0.688995844 0.284 0.064  1.399840e-07       5
#> PMP223           5.305580e-11 0.808333876 0.791 0.525  1.591674e-07       5
#> ELF53            5.718205e-11 0.899477501 0.761 0.413  1.715461e-07       5
#> CYP39A12         6.422615e-11 0.742209734 0.761 0.342  1.926784e-07       5
#> GCSH2            6.996921e-11 0.898709470 0.761 0.427  2.099076e-07       5
#> FZD51            7.543134e-11 0.715340868 0.567 0.218  2.262940e-07       5
#> CRABP12          7.799540e-11 0.803385068 0.761 0.389  2.339862e-07       5
#> HSPA1A2          8.252972e-11 0.829750426 0.985 0.476  2.475892e-07       5
#> SERTAD4-AS12     9.719413e-11 0.790350640 0.701 0.430  2.915824e-07       5
#> GATA31           9.837434e-11 0.794089471 0.537 0.209  2.951230e-07       5
#> S100B2           1.084147e-10 0.751009831 0.806 0.398  3.252441e-07       5
#> PXDNL            1.230414e-10 0.856322924 0.254 0.054  3.691242e-07       5
#> SCIN2            1.396253e-10 0.640131058 0.493 0.164  4.188758e-07       5
#> ARL6IP12         1.407179e-10 0.593999177 0.866 0.463  4.221536e-07       5
#> JHDM1D-AS12      1.492311e-10 0.768174113 0.672 0.301  4.476934e-07       5
#> TIMM102          1.657483e-10 0.797701207 0.925 0.465  4.972449e-07       5
#> TMEM1392         1.944330e-10 0.760929691 0.701 0.323  5.832991e-07       5
#> MFSD61           2.021395e-10 0.915901974 0.597 0.265  6.064186e-07       5
#> LNX12            2.667718e-10 0.719232201 0.731 0.337  8.003154e-07       5
#> SMIM5            2.899659e-10 0.918336988 0.343 0.099  8.698976e-07       5
#> FBXL162          2.987310e-10 0.966182673 0.418 0.147  8.961929e-07       5
#> FOXC12           3.259709e-10 0.757709955 0.701 0.329  9.779127e-07       5
#> SERPINB53        3.273066e-10 0.786906532 0.746 0.416  9.819199e-07       5
#> ANXA2R2          3.400711e-10 0.922162218 0.552 0.245  1.020213e-06       5
#> LMTK3            3.789007e-10 0.952394829 0.313 0.087  1.136702e-06       5
#> IL17B2           3.940064e-10 0.772631694 0.478 0.172  1.182019e-06       5
#> ADAM153          4.016102e-10 0.867483457 0.806 0.492  1.204831e-06       5
#> SUN32            5.061009e-10 0.684121251 0.537 0.207  1.518303e-06       5
#> SENCR1           5.516550e-10 0.700821824 0.343 0.094  1.654965e-06       5
#> WFDC21           6.249185e-10 0.698386660 0.388 0.121  1.874756e-06       5
#> C1orf1151        6.611849e-10 0.741605613 0.627 0.265  1.983555e-06       5
#> TMEM1582         7.744760e-10 0.626485313 0.716 0.334  2.323428e-06       5
#> PYCR12           1.098449e-09 0.733083269 0.791 0.442  3.295346e-06       5
#> AP1M22           1.573461e-09 0.769144929 0.776 0.450  4.720383e-06       5
#> SYNGR12          1.752273e-09 0.648414865 0.731 0.348  5.256818e-06       5
#> TMEM106C2        1.827983e-09 0.671042029 0.851 0.398  5.483949e-06       5
#> KRT74            1.829596e-09 0.739175229 0.955 0.600  5.488787e-06       5
#> DEFB11           1.844777e-09 0.705854677 0.493 0.185  5.534332e-06       5
#> PALLD1           1.950001e-09 0.612274570 0.731 0.351  5.850004e-06       5
#> TP53BP22         2.095778e-09 0.764815572 0.776 0.423  6.287333e-06       5
#> MYBL12           2.281556e-09 0.736165213 0.567 0.249  6.844669e-06       5
#> TEKT33           2.401866e-09 0.752841712 0.836 0.496  7.205597e-06       5
#> MB3              2.809523e-09 0.668575106 0.597 0.259  8.428570e-06       5
#> TPM22            2.824175e-09 0.612896003 0.746 0.417  8.472526e-06       5
#> RHOBTB31         2.938304e-09 0.706851943 0.746 0.439  8.814912e-06       5
#> CLIC3            3.362015e-09 0.796244640 0.299 0.083  1.008604e-05       5
#> TMX22            3.398291e-09 0.700984284 0.851 0.452  1.019487e-05       5
#> MAPK132          4.050523e-09 0.763723332 0.746 0.486  1.215157e-05       5
#> POSTN2           4.050965e-09 0.459474288 0.507 0.231  1.215289e-05       5
#> OVOL12           4.553843e-09 0.789309971 0.507 0.220  1.366153e-05       5
#> GRB143           4.935853e-09 0.607073160 0.537 0.217  1.480756e-05       5
#> KRT84            5.139678e-09 0.717348767 0.866 0.536  1.541903e-05       5
#> NUDT42           5.821375e-09 0.952266479 0.687 0.376  1.746413e-05       5
#> RDH103           6.332318e-09 0.633180561 0.701 0.391  1.899695e-05       5
#> CA11             6.549893e-09 0.603435272 0.343 0.104  1.964968e-05       5
#> CKS1B2           6.938530e-09 0.672398957 0.910 0.490  2.081559e-05       5
#> SFRP12           7.066987e-09 0.701066353 0.821 0.451  2.120096e-05       5
#> PDGFRA2          8.925358e-09 0.580553447 0.716 0.318  2.677607e-05       5
#> FBXO322          1.051643e-08 0.667295556 0.851 0.448  3.154930e-05       5
#> MAP1B2           1.416988e-08 0.742065040 0.761 0.406  4.250965e-05       5
#> ERVMER34-11      1.499008e-08 0.509465797 0.537 0.217  4.497024e-05       5
#> C2orf823         1.729655e-08 0.720876233 0.896 0.483  5.188966e-05       5
#> CLDN72           2.155565e-08 0.719965360 0.821 0.534  6.466695e-05       5
#> NUPR22           2.973987e-08 0.584016228 0.731 0.356  8.921961e-05       5
#> ID41             3.098071e-08 0.690117125 0.552 0.253  9.294212e-05       5
#> CHI3L12          3.344269e-08 0.667705425 0.776 0.482  1.003281e-04       5
#> RAB3IP3          3.368665e-08 0.703510013 0.597 0.294  1.010600e-04       5
#> HRCT12           3.670970e-08 0.690635163 0.672 0.382  1.101291e-04       5
#> CRYAB3           3.873965e-08 0.671340715 0.925 0.549  1.162189e-04       5
#> LEFTY21          3.920647e-08 0.552840729 0.463 0.183  1.176194e-04       5
#> SSRP12           4.316763e-08 0.628491153 0.910 0.474  1.295029e-04       5
#> OPRK11           4.760258e-08 0.602823366 0.463 0.189  1.428077e-04       5
#> AQP52            6.506461e-08 0.495900326 0.731 0.418  1.951938e-04       5
#> SOX83            6.771574e-08 0.611313565 0.687 0.387  2.031472e-04       5
#> MMP152           8.056731e-08 0.614786302 0.657 0.395  2.417019e-04       5
#> HEY22            8.738790e-08 0.639759915 0.493 0.216  2.621637e-04       5
#> RP11-554I8.22    1.145785e-07 0.582415734 0.552 0.259  3.437355e-04       5
#> PHLDA11          1.196565e-07 0.597416152 0.672 0.395  3.589694e-04       5
#> WNT7B            1.438968e-07 0.703202878 0.373 0.142  4.316905e-04       5
#> IFI62            1.541881e-07 0.577611971 0.716 0.415  4.625644e-04       5
#> FAM3C2           1.747058e-07 0.608186894 0.806 0.418  5.241174e-04       5
#> TUBB62           2.098892e-07 0.688935581 0.746 0.459  6.296677e-04       5
#> PTRF2            2.150849e-07 0.480181556 0.657 0.317  6.452547e-04       5
#> AC005152.33      2.535661e-07 0.569726724 0.731 0.429  7.606984e-04       5
#> TMEM612          2.839707e-07 0.517936836 0.716 0.350  8.519122e-04       5
#> C6orf1323        2.941425e-07 0.590231749 0.746 0.443  8.824275e-04       5
#> PIM11            3.009446e-07 0.556796163 0.731 0.437  9.028339e-04       5
#> FAM46B2          3.064936e-07 0.637094992 0.567 0.273  9.194807e-04       5
#> CEBPD1           3.134703e-07 0.612721479 0.851 0.481  9.404110e-04       5
#> MLLT112          3.189114e-07 0.456390387 0.672 0.310  9.567342e-04       5
#> SNHG252          3.199183e-07 0.605046460 0.672 0.363  9.597550e-04       5
#> MGP2             3.277445e-07 0.676705393 0.881 0.534  9.832336e-04       5
#> GADD45G1         3.559642e-07 0.572543988 0.582 0.272  1.067893e-03       5
#> RPL39L2          3.702637e-07 0.596974112 0.746 0.454  1.110791e-03       5
#> TNFSF13B2        3.781430e-07 0.665669148 0.716 0.431  1.134429e-03       5
#> B3GNT72          3.971848e-07 0.531438867 0.642 0.327  1.191554e-03       5
#> BGN2             4.007378e-07 0.578626288 0.806 0.473  1.202213e-03       5
#> CEP702           4.372944e-07 0.515439264 0.507 0.239  1.311883e-03       5
#> FRZB2            4.444709e-07 0.691479104 0.507 0.246  1.333413e-03       5
#> TINCR3           4.555483e-07 0.507316215 0.657 0.308  1.366645e-03       5
#> HN11             4.582287e-07 0.510285797 0.806 0.486  1.374686e-03       5
#> DSC32            4.947347e-07 0.479247740 0.642 0.309  1.484204e-03       5
#> VGF1             5.061010e-07 0.599454110 0.403 0.162  1.518303e-03       5
#> CSRP12           5.439371e-07 0.616489610 0.687 0.447  1.631811e-03       5
#> PVRL42           6.147875e-07 0.536497911 0.716 0.381  1.844362e-03       5
#> ZG16B3           6.510361e-07 0.539480394 0.642 0.351  1.953108e-03       5
#> MARC11           6.551728e-07 0.785351513 0.343 0.135  1.965518e-03       5
#> NEXN1            6.910728e-07 0.519872196 0.313 0.109  2.073218e-03       5
#> CD243            7.845320e-07 0.680892948 0.985 0.636  2.353596e-03       5
#> RP11-357H14.172  8.211295e-07 0.590007952 0.448 0.197  2.463389e-03       5
#> ARHGAP292        9.134939e-07 0.548595493 0.716 0.420  2.740482e-03       5
#> QPCT2            9.634519e-07 0.523648321 0.612 0.350  2.890356e-03       5
#> FBXL221          1.019321e-06 0.650727765 0.418 0.183  3.057963e-03       5
#> VANGL13          1.058150e-06 0.618767049 0.806 0.519  3.174449e-03       5
#> AIF1L3           1.098205e-06 0.626213821 0.731 0.452  3.294614e-03       5
#> SAT12            1.113914e-06 0.569951226 0.791 0.504  3.341741e-03       5
#> CDC42EP12        1.342722e-06 0.571334519 0.687 0.445  4.028166e-03       5
#> IMPA22           1.447120e-06 0.549649932 0.776 0.494  4.341361e-03       5
#> MMP72            1.464389e-06 0.565556240 0.672 0.354  4.393168e-03       5
#> PCBD11           1.586955e-06 0.573812151 0.776 0.494  4.760865e-03       5
#> HIST1H4E2        1.622787e-06 0.569898062 0.418 0.178  4.868360e-03       5
#> DKK11            1.627181e-06 0.358372558 0.373 0.146  4.881543e-03       5
#> PCOLCE22         1.700523e-06 0.411131046 0.627 0.318  5.101568e-03       5
#> OCLN2            1.734291e-06 0.596941138 0.552 0.286  5.202872e-03       5
#> LA16c-380H5.52   1.781748e-06 0.428291808 0.507 0.231  5.345244e-03       5
#> XAGE22           2.499338e-06 0.521383497 0.582 0.293  7.498015e-03       5
#> GMPR2            2.727019e-06 0.358704206 0.537 0.259  8.181058e-03       5
#> STMN12           2.757608e-06 0.548329424 0.791 0.517  8.272823e-03       5
#> MUC12            3.616693e-06 0.454019110 0.657 0.354  1.085008e-02       5
#> PBX12            3.863521e-06 0.532471731 0.701 0.457  1.159056e-02       5
#> RP11-798M19.61   3.915697e-06 0.687298500 0.478 0.252  1.174709e-02       5
#> NEGR12           4.128908e-06 0.476713379 0.448 0.200  1.238673e-02       5
#> MSRB32           4.249733e-06 0.457149982 0.493 0.238  1.274920e-02       5
#> CTNND22          5.088247e-06 0.505262839 0.507 0.245  1.526474e-02       5
#> KRT232           5.159096e-06 0.532862467 0.806 0.489  1.547729e-02       5
#> RAD213           5.978279e-06 0.474176315 0.851 0.520  1.793484e-02       5
#> LGALS32          6.639928e-06 0.490262590 0.806 0.515  1.991978e-02       5
#> LSR2             6.792202e-06 0.560221046 0.716 0.502  2.037661e-02       5
#> HAPLN11          7.189532e-06 0.215856093 0.403 0.161  2.156860e-02       5
#> BARD12           7.962521e-06 0.624759705 0.567 0.316  2.388756e-02       5
#> DSEL1            8.623474e-06 0.516374652 0.627 0.332  2.587042e-02       5
#> SNHG191          8.893831e-06 0.572400437 0.612 0.409  2.668149e-02       5
#> ITGB42           9.411090e-06 0.485943138 0.612 0.340  2.823327e-02       5
#> TNFRSF12A2       9.971816e-06 0.511914900 0.642 0.432  2.991545e-02       5
#> NFIB3            1.007718e-05 0.489781438 0.657 0.449  3.023153e-02       5
#> STAT11           1.163662e-05 0.430138629 0.716 0.427  3.490985e-02       5
#> KCNMB12          1.251331e-05 0.526007773 0.612 0.336  3.753992e-02       5
#> ACTG23           1.302443e-05 0.503192211 0.657 0.378  3.907329e-02       5
#> INAFM1           1.487440e-05 0.355711824 0.388 0.169  4.462320e-02       5
#> KRT813           1.784782e-05 0.545271665 0.806 0.511  5.354345e-02       5
#> DLX52            1.997018e-05 0.421040548 0.433 0.205  5.991054e-02       5
#> AARD3            2.030606e-05 0.428450847 0.672 0.440  6.091817e-02       5
#> SYCP21           2.097023e-05 0.420878769 0.478 0.238  6.291070e-02       5
#> TPD52L13         2.104951e-05 0.486412172 0.687 0.463  6.314853e-02       5
#> KRT192           2.137019e-05 0.522247493 0.866 0.505  6.411056e-02       5
#> EFHD12           2.413750e-05 0.490362759 0.687 0.416  7.241250e-02       5
#> NSG12            2.461867e-05 0.471120807 0.388 0.179  7.385601e-02       5
#> ACTA22           2.613727e-05 0.548838581 0.627 0.345  7.841181e-02       5
#> VASN4            2.631389e-05 0.487193535 0.716 0.451  7.894166e-02       5
#> MTSS1L2          2.708124e-05 0.386843030 0.687 0.389  8.124372e-02       5
#> LAYN             2.764698e-05 0.468858387 0.299 0.122  8.294094e-02       5
#> IL17RE2          3.376782e-05 0.366981514 0.567 0.303  1.013035e-01       5
#> CKS22            3.467934e-05 0.341141855 0.701 0.433  1.040380e-01       5
#> NQO12            3.488484e-05 0.474233238 0.687 0.428  1.046545e-01       5
#> TBC1D12          3.579358e-05 0.455242531 0.716 0.446  1.073808e-01       5
#> MESP14           3.679694e-05 0.482778906 0.701 0.440  1.103908e-01       5
#> MAL23            3.791363e-05 0.442799300 0.776 0.482  1.137409e-01       5
#> LIMA1            3.809947e-05 0.377050979 0.463 0.224  1.142984e-01       5
#> LINC007072       4.003035e-05 0.442414242 0.284 0.113  1.200910e-01       5
#> PHYH3            4.233541e-05 0.483418273 0.731 0.492  1.270062e-01       5
#> NAB11            4.426575e-05 0.488672205 0.507 0.268  1.327973e-01       5
#> PPP1R14A1        4.671892e-05 0.617185122 0.388 0.193  1.401568e-01       5
#> SDC43            4.982930e-05 0.411722514 0.657 0.461  1.494879e-01       5
#> LEMD13           5.167075e-05 0.494907494 0.731 0.486  1.550122e-01       5
#> FKBP42           6.378092e-05 0.469865416 0.731 0.488  1.913427e-01       5
#> HACD12           6.881185e-05 0.487349120 0.672 0.399  2.064355e-01       5
#> CPM1             7.077472e-05 0.532329866 0.612 0.399  2.123242e-01       5
#> OVOS21           7.082467e-05 0.300175533 0.597 0.356  2.124740e-01       5
#> RGS161           7.215569e-05 0.239192261 0.582 0.309  2.164671e-01       5
#> SLC52A11         7.465929e-05 0.452997566 0.373 0.178  2.239779e-01       5
#> TTC39A2          7.898564e-05 0.411966950 0.478 0.259  2.369569e-01       5
#> SPEG2            8.216745e-05 0.472352994 0.433 0.224  2.465024e-01       5
#> GABRP2           8.428608e-05 0.400044601 0.761 0.445  2.528582e-01       5
#> HILPDA2          9.039975e-05 0.341735568 0.582 0.394  2.711992e-01       5
#> MEST1            1.017272e-04 0.325277165 0.299 0.126  3.051815e-01       5
#> MEX3A3           1.185503e-04 0.404587391 0.597 0.343  3.556508e-01       5
#> ADHFE12          1.203128e-04 0.322925582 0.433 0.213  3.609385e-01       5
#> GPX11            1.232793e-04 0.437531713 0.746 0.478  3.698378e-01       5
#> NES1             1.237841e-04 0.334211812 0.627 0.330  3.713524e-01       5
#> DOK51            1.239371e-04 0.414859886 0.299 0.135  3.718113e-01       5
#> DSN1             1.266785e-04 0.450798155 0.448 0.242  3.800355e-01       5
#> CNTNAP3B2        1.288832e-04 0.313221577 0.522 0.264  3.866496e-01       5
#> VGLL12           1.448103e-04 0.333126223 0.582 0.360  4.344310e-01       5
#> SQLE2            1.470091e-04 0.401974967 0.627 0.436  4.410272e-01       5
#> CYR612           1.835546e-04 0.303403141 0.672 0.412  5.506639e-01       5
#> PMEPA1           1.948100e-04 0.393019402 0.657 0.397  5.844300e-01       5
#> ITGA102          2.046295e-04 0.445938881 0.403 0.215  6.138886e-01       5
#> IRX33            2.603701e-04 0.388721271 0.612 0.443  7.811102e-01       5
#> CD822            2.772964e-04 0.361945027 0.627 0.418  8.318892e-01       5
#> SLBP2            2.955846e-04 0.347645028 0.746 0.487  8.867538e-01       5
#> SOX43            3.150892e-04 0.488743138 0.776 0.567  9.452675e-01       5
#> NR2F23           3.212376e-04 0.409269788 0.701 0.450  9.637127e-01       5
#> LRP22            3.473506e-04 0.347540745 0.448 0.239  1.000000e+00       5
#> GAS62            3.596913e-04 0.399739648 0.597 0.434  1.000000e+00       5
#> NME11            3.762513e-04 0.394989012 0.657 0.489  1.000000e+00       5
#> EFNA51           3.931595e-04 0.364284149 0.313 0.150  1.000000e+00       5
#> RAMP21           4.141845e-04 0.336845582 0.448 0.250  1.000000e+00       5
#> RAC31            4.150656e-04 0.198859456 0.433 0.208  1.000000e+00       5
#> ELF33            4.225858e-04 0.388998640 0.657 0.492  1.000000e+00       5
#> RP11-19E11.13    4.262496e-04 0.362760612 0.657 0.412  1.000000e+00       5
#> DGAT22           4.523765e-04 0.294726306 0.552 0.314  1.000000e+00       5
#> FAM89A2          5.188629e-04 0.279089322 0.567 0.330  1.000000e+00       5
#> SLPI3            5.717963e-04 0.405289442 0.478 0.332  1.000000e+00       5
#> TUBA4A2          6.420528e-04 0.405960323 0.582 0.344  1.000000e+00       5
#> SHISA92          6.423888e-04 0.397055679 0.358 0.186  1.000000e+00       5
#> PAQR41           6.958458e-04 0.431867135 0.522 0.328  1.000000e+00       5
#> RP11-161M6.21    7.356579e-04 0.387496785 0.433 0.251  1.000000e+00       5
#> EPHX22           7.419670e-04 0.385169119 0.448 0.258  1.000000e+00       5
#> GGCT2            7.594745e-04 0.431014168 0.701 0.475  1.000000e+00       5
#> AC022007.51      9.538582e-04 0.277060743 0.403 0.219  1.000000e+00       5
#> HLA-C2           9.639796e-04 0.356623441 0.731 0.516  1.000000e+00       5
#> P3H41            1.004071e-03 0.352707729 0.612 0.380  1.000000e+00       5
#> LTBP11           1.022583e-03 0.338630835 0.343 0.180  1.000000e+00       5
#> LAPTM4B3         1.086510e-03 0.414014693 0.791 0.540  1.000000e+00       5
#> NANOS13          1.132923e-03 0.376306446 0.448 0.259  1.000000e+00       5
#> HSP90AB14        1.136829e-03 0.568510713 0.985 0.652  1.000000e+00       5
#> ABCA51           1.309097e-03 0.396494820 0.403 0.239  1.000000e+00       5
#> FKBP103          1.393494e-03 0.324468821 0.597 0.446  1.000000e+00       5
#> HES13            1.458014e-03 0.335992582 0.567 0.450  1.000000e+00       5
#> CACYBP3          1.575659e-03 0.391678018 0.761 0.523  1.000000e+00       5
#> ACTL6A2          1.603715e-03 0.357457664 0.687 0.481  1.000000e+00       5
#> NEBL2            1.672324e-03 0.334080103 0.507 0.298  1.000000e+00       5
#> DSP2             1.870574e-03 0.328996605 0.642 0.500  1.000000e+00       5
#> CHI3L22          1.902427e-03 0.314526067 0.552 0.385  1.000000e+00       5
#> SLC2A122         1.908582e-03 0.357112922 0.269 0.135  1.000000e+00       5
#> CLUL11           1.967158e-03 0.359591468 0.299 0.159  1.000000e+00       5
#> ZNF5281          2.001701e-03 0.358068856 0.343 0.189  1.000000e+00       5
#> CLMN2            2.488962e-03 0.305740559 0.642 0.433  1.000000e+00       5
#> CIART1           2.511068e-03 0.432747970 0.448 0.282  1.000000e+00       5
#> TFAP2A3          2.569072e-03 0.290872686 0.657 0.427  1.000000e+00       5
#> ST142            2.576388e-03 0.326526510 0.687 0.484  1.000000e+00       5
#> NDRG13           2.646737e-03 0.192672279 0.642 0.438  1.000000e+00       5
#> CARHSP12         2.858296e-03 0.355907490 0.687 0.486  1.000000e+00       5
#> CALD12           2.869797e-03 0.221763261 0.612 0.413  1.000000e+00       5
#> TUBA1C2          2.891051e-03 0.297481253 0.642 0.462  1.000000e+00       5
#> KLRG23           2.960597e-03 0.362809119 0.343 0.204  1.000000e+00       5
#> LINC014361       3.039867e-03 0.269081383 0.493 0.311  1.000000e+00       5
#> TUBB3            3.158204e-03 0.371801006 0.612 0.525  1.000000e+00       5
#> CKB3             3.388644e-03 0.312446117 0.552 0.450  1.000000e+00       5
#> C1orf562         4.179149e-03 0.289326933 0.612 0.435  1.000000e+00       5
#> KRT184           4.575145e-03 0.401513081 0.791 0.539  1.000000e+00       5
#> CRNDE3           4.632410e-03 0.356470915 0.657 0.503  1.000000e+00       5
#> MAP22            4.790219e-03 0.312580078 0.358 0.212  1.000000e+00       5
#> NKD22            4.807649e-03 0.294030620 0.313 0.176  1.000000e+00       5
#> PODXL23          4.990078e-03 0.287484571 0.597 0.454  1.000000e+00       5
#> QPRT2            5.367459e-03 0.352279997 0.448 0.275  1.000000e+00       5
#> RGCC1            5.539087e-03 0.295926730 0.552 0.386  1.000000e+00       5
#> PRELP2           5.649211e-03 0.163003788 0.403 0.227  1.000000e+00       5
#> PROM13           5.678529e-03 0.289770668 0.687 0.485  1.000000e+00       5
#> CDH111           7.150563e-03 0.082229501 0.328 0.186  1.000000e+00       5
#> CLPSL11          7.238339e-03 0.122098981 0.493 0.325  1.000000e+00       5
#> IGFBP72          7.405392e-03 0.134630052 0.597 0.409  1.000000e+00       5
#> IGKC1            8.169671e-03 0.109122948 0.507 0.321  1.000000e+00       5
#> ATP1B12          8.379147e-03 0.264392159 0.642 0.454  1.000000e+00       5
#> C6orf1411        8.437130e-03 0.271153589 0.269 0.147  1.000000e+00       5
#> LAMB11           8.786008e-03 0.199004496 0.313 0.180  1.000000e+00       5
#> APOBEC3B1        8.825660e-03 0.369724807 0.269 0.156  1.000000e+00       5
#> SOX93            9.291404e-03 0.242875829 0.582 0.417  1.000000e+00       5
#> MITF1            9.490490e-03 0.175542184 0.358 0.207  1.000000e+00       5
#> CFH             1.549368e-162 3.855307455 0.852 0.012 4.648103e-159       6
#> NTM             1.774126e-160 3.644233889 0.741 0.004 5.322378e-157       6
#> PDGFRB          6.010230e-159 3.979637019 0.907 0.018 1.803069e-155       6
#> FBN1            1.553514e-150 3.901691667 0.852 0.016 4.660541e-147       6
#> ECM2            5.919022e-149 3.562661774 0.722 0.006 1.775707e-145       6
#> CHN1            7.750026e-148 3.859078837 0.852 0.017 2.325008e-144       6
#> FAP             8.784725e-148 3.923384250 0.889 0.021 2.635418e-144       6
#> COL5A1          5.076205e-147 3.929573918 0.944 0.028 1.522862e-143       6
#> ITGA1           2.103140e-145 3.615170773 0.815 0.014 6.309420e-142       6
#> COL5A2          1.481966e-139 4.300739035 1.000 0.041 4.445898e-136       6
#> LAMA4           9.820803e-139 3.702339957 0.852 0.021 2.946241e-135       6
#> EMILIN1         3.227461e-138 4.022240065 0.944 0.034 9.682384e-135       6
#> SULF1           1.619303e-137 3.692688603 0.778 0.014 4.857910e-134       6
#> SGIP1           4.219697e-137 3.332624886 0.630 0.003 1.265909e-133       6
#> WNT2            4.703665e-136 3.284501570 0.593 0.001 1.411099e-132       6
#> COL10A1         3.948281e-134 3.681758801 0.685 0.008 1.184484e-130       6
#> VCAN            1.277715e-130 4.006838285 0.944 0.038 3.833144e-127       6
#> VCAM1           2.827827e-122 3.233753732 0.759 0.018 8.483480e-119       6
#> NDN             2.175687e-121 3.462044571 0.759 0.019 6.527061e-118       6
#> COL6A3          3.518074e-121 4.273125377 0.981 0.053 1.055422e-117       6
#> PLAC9           1.724417e-118 3.514381787 0.685 0.013 5.173251e-115       6
#> DCN             1.943554e-118 4.193990586 0.981 0.052 5.830661e-115       6
#> DKK3            3.422451e-118 3.532767444 0.722 0.017 1.026735e-114       6
#> PODN            3.941256e-118 3.036853102 0.519 0.001 1.182377e-114       6
#> LRRC32          1.056351e-114 2.926810900 0.537 0.003 3.169052e-111       6
#> CAV1            5.728788e-114 3.563914367 0.778 0.026 1.718636e-110       6
#> NID2            1.240036e-113 3.278147193 0.611 0.009 3.720107e-110       6
#> ITGBL1          1.305380e-113 3.086810999 0.574 0.006 3.916140e-110       6
#> ITGA11          1.411623e-110 3.153285727 0.574 0.007 4.234868e-107       6
#> COX7A1          2.674726e-110 3.620613669 0.796 0.031 8.024178e-107       6
#> THY1            1.351084e-109 4.355899004 1.000 0.065 4.053252e-106       6
#> GNG11           1.801073e-108 3.478746378 0.852 0.040 5.403219e-105       6
#> EDNRA           3.506232e-108 3.145028642 0.537 0.005 1.051870e-104       6
#> GPX8            1.146550e-106 3.695438654 0.907 0.052 3.439651e-103       6
#> MFAP5           4.955760e-106 3.244239584 0.611 0.012 1.486728e-102       6
#> TFPI            2.456858e-104 3.383827838 0.741 0.027 7.370573e-101       6
#> SPARCL1         2.741904e-104 3.674670723 0.778 0.033 8.225712e-101       6
#> PLAT            5.321875e-104 3.300609306 0.667 0.018 1.596563e-100       6
#> HEG1            1.794661e-103 2.975430052 0.593 0.011 5.383982e-100       6
#> ANGPTL2         5.517836e-102 3.023772207 0.537 0.007  1.655351e-98       6
#> RCN3            6.550294e-102 3.845294315 0.852 0.047  1.965088e-98       6
#> ADAMTS12        9.291128e-102 3.047509338 0.537 0.007  2.787339e-98       6
#> PRRX1           1.353940e-100 3.701573017 0.852 0.048  4.061820e-97       6
#> CPXM1            3.893542e-97 3.192361154 0.611 0.015  1.168063e-93       6
#> TIMP3            7.015956e-97 3.270606603 0.722 0.029  2.104787e-93       6
#> ASPN             1.671827e-95 3.206777048 0.593 0.014  5.015481e-92       6
#> AEBP1            1.489803e-93 4.215691168 0.963 0.078  4.469408e-90       6
#> C11orf96         3.111576e-93 3.583988352 0.907 0.067  9.334728e-90       6
#> ENPEP            2.747213e-92 2.732441355 0.426 0.002  8.241638e-89       6
#> EBF1             1.405558e-90 2.777108374 0.500 0.008  4.216674e-87       6
#> OLFML3           1.978961e-89 3.409625232 0.796 0.048  5.936884e-86       6
#> DPYSL3           4.282383e-89 3.044082947 0.574 0.015  1.284715e-85       6
#> WISP2            2.207113e-85 2.674943252 0.426 0.004  6.621340e-82       6
#> CD248            6.493572e-85 3.138947879 0.611 0.022  1.948072e-81       6
#> CREB3L1          2.180259e-84 2.660454436 0.407 0.003  6.540776e-81       6
#> CYGB             2.625753e-84 3.067124088 0.630 0.025  7.877260e-81       6
#> KCNE4            9.919817e-84 2.668818627 0.481 0.009  2.975945e-80       6
#> RARRES2          6.349696e-82 4.069281096 0.926 0.085  1.904909e-78       6
#> SFRP2            1.381656e-81 3.635066469 0.704 0.038  4.144967e-78       6
#> SFRP4            2.234503e-81 2.630929013 0.500 0.012  6.703509e-78       6
#> COL15A1          2.438789e-81 2.783005401 0.556 0.017  7.316367e-78       6
#> COL3A1           1.728165e-79 4.403334607 0.981 0.091  5.184495e-76       6
#> PLXDC1           2.352659e-79 3.013079638 0.667 0.034  7.057978e-76       6
#> SPRY1            4.319900e-78 3.098032413 0.685 0.038  1.295970e-74       6
#> COL5A3           1.613450e-77 2.796800278 0.463 0.010  4.840350e-74       6
#> ADAMTS2          3.512070e-77 2.784100137 0.519 0.015  1.053621e-73       6
#> FIBIN            9.530667e-77 3.004930169 0.593 0.025  2.859200e-73       6
#> EDIL3            4.073693e-75 2.956581826 0.556 0.021  1.222108e-71       6
#> SEMA5A           9.471473e-75 3.127449967 0.648 0.035  2.841442e-71       6
#> COL12A1          1.156752e-74 3.772094625 0.907 0.103  3.470256e-71       6
#> LHFP             1.536408e-74 3.449352022 0.852 0.083  4.609223e-71       6
#> ZFHX4            1.978204e-74 2.304353540 0.315 0.000  5.934611e-71       6
#> MEG3             2.715689e-72 3.107066580 0.574 0.026  8.147068e-69       6
#> TNFAIP6          3.410466e-72 3.251091756 0.870 0.090  1.023140e-68       6
#> FOXS1            1.736151e-71 2.426747405 0.352 0.003  5.208453e-68       6
#> HEYL             2.252031e-71 2.357610907 0.352 0.003  6.756093e-68       6
#> IFI27            2.265267e-71 3.311371137 0.981 0.107  6.795801e-68       6
#> ADGRF5           1.319111e-70 2.436239507 0.389 0.006  3.957334e-67       6
#> SSPN             2.376422e-69 3.206293366 0.722 0.057  7.129265e-66       6
#> HIC1             3.166835e-69 2.492801923 0.426 0.010  9.500504e-66       6
#> NOX4             3.611067e-68 2.472819283 0.389 0.007  1.083320e-64       6
#> LUM              5.967029e-68 4.352056556 0.981 0.107  1.790109e-64       6
#> NID1             1.015152e-65 3.155527921 0.870 0.105  3.045457e-62       6
#> BICC1            6.607654e-65 2.838452372 0.556 0.029  1.982296e-61       6
#> OLFML2B1         1.028486e-64 3.214906552 0.944 0.131  3.085457e-61       6
#> FILIP1           2.636101e-63 2.603036287 0.426 0.012  7.908302e-60       6
#> COL1A1           4.124377e-63 4.374868515 0.981 0.110  1.237313e-59       6
#> PRSS231          5.979980e-63 4.121581893 0.981 0.163  1.793994e-59       6
#> FRMD6            9.264836e-63 2.888937159 0.537 0.028  2.779451e-59       6
#> IGFBP4           7.322112e-62 3.759895130 0.963 0.160  2.196634e-58       6
#> FILIP1L1         2.893712e-61 2.979298613 0.815 0.094  8.681135e-58       6
#> ECM11            4.193633e-61 2.896948566 0.796 0.087  1.258090e-57       6
#> COL8A1           2.294526e-60 3.261537107 0.704 0.066  6.883578e-57       6
#> COL1A2           2.406506e-59 4.387524014 0.981 0.133  7.219517e-56       6
#> AKAP12           4.590554e-59 2.648789564 0.481 0.022  1.377166e-55       6
#> PCOLCE           2.544099e-58 4.097761096 0.963 0.172  7.632298e-55       6
#> MXRA5            3.059226e-58 3.526517877 0.889 0.142  9.177679e-55       6
#> MAP1A            8.245436e-58 2.214531878 0.444 0.017  2.473631e-54       6
#> ANGPT2           8.931419e-58 2.775161218 0.537 0.033  2.679426e-54       6
#> GJA4             3.272047e-57 1.943844671 0.259 0.001  9.816142e-54       6
#> CCDC102B         4.146109e-57 2.700545833 0.500 0.027  1.243833e-53       6
#> COL18A1          4.351316e-57 3.612756820 0.852 0.128  1.305395e-53       6
#> INHBA            1.068712e-56 3.270034714 0.685 0.068  3.206137e-53       6
#> GGT51            2.197746e-56 3.167473844 0.796 0.104  6.593238e-53       6
#> ABCC9            4.228113e-56 2.394824261 0.333 0.007  1.268434e-52       6
#> ANTXR1           6.137325e-56 3.122040899 0.796 0.105  1.841197e-52       6
#> VGLL3            7.251079e-56 2.261497489 0.333 0.007  2.175324e-52       6
#> LRRC15           1.187604e-55 2.273616835 0.296 0.004  3.562812e-52       6
#> EGFLAM           1.187729e-55 2.040458154 0.296 0.004  3.563186e-52       6
#> P3H3             1.212165e-55 2.498660285 0.389 0.012  3.636494e-52       6
#> EFEMP1           2.195035e-54 3.053240931 0.574 0.044  6.585104e-51       6
#> LGALS3BP1        2.262354e-54 3.408887711 0.963 0.184  6.787061e-51       6
#> NOTCH31          2.364195e-54 3.156026790 0.741 0.090  7.092585e-51       6
#> GUCY1A2          7.597100e-54 2.350661355 0.333 0.008  2.279130e-50       6
#> CTSK1            3.077780e-53 2.391542252 0.944 0.139  9.233339e-50       6
#> GUCY1A31         5.526880e-53 3.057487125 0.741 0.093  1.658064e-49       6
#> AXL              3.733208e-52 2.240195980 0.648 0.058  1.119962e-48       6
#> C3orf80          5.271434e-52 2.143210005 0.315 0.007  1.581430e-48       6
#> SERPING1         7.677982e-52 3.739683692 0.981 0.213  2.303395e-48       6
#> CEMIP            1.727984e-51 2.113310934 0.278 0.004  5.183952e-48       6
#> SYTL2            2.706193e-51 2.592367753 0.537 0.039  8.118579e-48       6
#> COL4A1           2.833397e-51 3.452377482 0.907 0.165  8.500190e-48       6
#> ID31             2.895137e-51 2.767912964 0.852 0.129  8.685411e-48       6
#> EID11            6.932214e-51 2.884078901 1.000 0.230  2.079664e-47       6
#> C1QTNF6          2.228459e-50 2.808366330 0.556 0.045  6.685378e-47       6
#> SPARC1           1.442667e-49 3.956257220 1.000 0.194  4.328001e-46       6
#> CRISPLD2         2.133991e-49 2.692328937 0.500 0.035  6.401974e-46       6
#> POSTN3           3.116264e-49 3.247799457 1.000 0.209  9.348791e-46       6
#> SERPINF11        3.231887e-49 3.421957077 0.963 0.204  9.695660e-46       6
#> PDE1A            9.246000e-49 1.984813168 0.278 0.005  2.773800e-45       6
#> CPE              1.462139e-48 3.082895109 0.741 0.105  4.386416e-45       6
#> CDH112           2.716680e-48 3.312717064 0.852 0.161  8.150039e-45       6
#> UNC5B            1.592034e-47 2.719609865 0.630 0.068  4.776101e-44       6
#> FMO1             2.723721e-47 1.834487312 0.259 0.004  8.171164e-44       6
#> ACKR4            3.858856e-47 1.571848098 0.259 0.004  1.157657e-43       6
#> CERCAM1          9.928304e-46 3.067806145 0.907 0.212  2.978491e-42       6
#> FN11             2.890942e-45 2.982994020 1.000 0.198  8.672827e-42       6
#> PDGFRL           8.151654e-45 2.246246712 0.370 0.017  2.445496e-41       6
#> ZEB1             1.327001e-44 2.170471011 0.463 0.032  3.981004e-41       6
#> FSTL12           4.171585e-44 3.015344695 0.870 0.188  1.251475e-40       6
#> KCNJ8            4.837711e-44 1.970417622 0.315 0.011  1.451313e-40       6
#> GUCY1B3          1.047869e-43 2.563056854 0.444 0.031  3.143608e-40       6
#> PODNL1           2.454537e-43 2.415877609 0.444 0.031  7.363611e-40       6
#> RUNX1T1          5.205177e-43 2.154488065 0.333 0.013  1.561553e-39       6
#> HTRA11           1.706611e-42 2.471329131 0.944 0.196  5.119834e-39       6
#> PLPP4            3.937270e-42 2.078306296 0.370 0.019  1.181181e-38       6
#> IFI161           7.853713e-42 2.303824222 1.000 0.248  2.356114e-38       6
#> CLEC11A          1.119015e-41 2.660794449 0.981 0.268  3.357044e-38       6
#> SRPX2            3.880280e-41 2.004458640 0.352 0.017  1.164084e-37       6
#> CCDC801          3.092401e-40 2.938501338 0.889 0.220  9.277204e-37       6
#> RBMS3            4.639429e-40 2.241260156 0.407 0.028  1.391829e-36       6
#> PTRF3            8.312303e-40 2.762354242 0.981 0.305  2.493691e-36       6
#> UACA1            1.220179e-39 3.255650779 0.981 0.327  3.660537e-36       6
#> CAV2             2.738855e-39 2.747474705 0.648 0.098  8.216566e-36       6
#> COMP             4.062476e-39 1.957107163 0.296 0.012  1.218743e-35       6
#> ADAM12           5.522637e-39 2.730110654 0.704 0.120  1.656791e-35       6
#> GJC1             5.544925e-39 1.995292521 0.370 0.022  1.663477e-35       6
#> C5orf46          7.426844e-39 2.071594519 0.481 0.043  2.228053e-35       6
#> LIMA11           1.839954e-38 2.522398606 0.870 0.206  5.519862e-35       6
#> SMOC22           2.001557e-38 2.854320838 0.833 0.198  6.004670e-35       6
#> MMP111           5.363156e-38 3.565295780 0.870 0.225  1.608947e-34       6
#> TGFB3            6.603677e-38 2.066965955 0.333 0.017  1.981103e-34       6
#> EHD2             7.811890e-38 2.368903950 0.593 0.080  2.343567e-34       6
#> XIST1            9.486135e-37 1.961861012 0.963 0.225  2.845840e-33       6
#> CARMN            1.780994e-36 2.201617590 0.519 0.058  5.342981e-33       6
#> TSHZ2            1.929637e-36 2.140266902 0.426 0.035  5.788912e-33       6
#> S1PR3            3.463840e-36 2.118694369 0.352 0.022  1.039152e-32       6
#> ADAMTS4          6.076527e-36 2.208737000 0.426 0.036  1.822958e-32       6
#> CALD13           1.671206e-35 2.862095947 1.000 0.395  5.013617e-32       6
#> SYDE1            1.698418e-35 2.211983006 0.389 0.030  5.095253e-32       6
#> COL6A21          2.177854e-35 3.058490497 1.000 0.408  6.533562e-32       6
#> RARRES3          2.202785e-35 2.078761780 0.593 0.083  6.608356e-32       6
#> IGFBP73          4.340433e-35 2.591044484 1.000 0.390  1.302130e-31       6
#> FGF7             1.568165e-34 2.872210438 0.611 0.101  4.704495e-31       6
#> INAFM11          2.345663e-34 2.558462246 0.722 0.154  7.036990e-31       6
#> BMP1             2.995461e-34 2.506716068 0.593 0.091  8.986384e-31       6
#> COL6A12          2.364211e-33 2.825155721 1.000 0.422  7.092634e-30       6
#> VAMP51           2.479480e-33 2.144325397 0.852 0.217  7.438439e-30       6
#> TCF41            7.404277e-33 2.106556229 0.889 0.243  2.221283e-29       6
#> SPON1            8.117802e-33 2.013711361 0.315 0.019  2.435340e-29       6
#> LZTS1            9.753210e-33 2.010634217 0.296 0.016  2.925963e-29       6
#> MTUS1            1.482786e-32 2.122005194 0.463 0.052  4.448357e-29       6
#> LAMB12           2.721836e-32 2.450634948 0.722 0.160  8.165509e-29       6
#> SALRNA2          5.913131e-32 1.595792884 0.259 0.012  1.773939e-28       6
#> GLT8D2           1.132811e-31 2.234796109 0.556 0.083  3.398434e-28       6
#> ABCA6            1.473994e-31 1.902445864 0.278 0.014  4.421981e-28       6
#> TIMP12           2.309606e-31 2.587034439 0.963 0.346  6.928818e-28       6
#> THBS21           3.135950e-31 2.969945635 0.833 0.256  9.407849e-28       6
#> CFI              3.171906e-31 1.547996581 0.315 0.020  9.515718e-28       6
#> SRPX             7.232811e-31 1.883436680 0.519 0.071  2.169843e-27       6
#> FOXP11           7.627390e-31 2.010470335 0.833 0.212  2.288217e-27       6
#> TMEM119          1.588419e-30 1.794099361 0.352 0.028  4.765257e-27       6
#> TPM23            2.441376e-30 2.360417535 0.963 0.410  7.324127e-27       6
#> TBX2             4.810537e-30 1.976776428 0.333 0.026  1.443161e-26       6
#> ANKRD28          6.191132e-30 2.188495064 0.741 0.178  1.857339e-26       6
#> COL4A23          1.337187e-29 3.039818357 0.870 0.320  4.011560e-26       6
#> BGN3             2.187730e-29 1.817992857 0.981 0.468  6.563189e-26       6
#> HTRA3            2.249645e-29 2.034086110 0.500 0.070  6.748935e-26       6
#> C1R3             2.801950e-29 2.503450683 0.926 0.341  8.405851e-26       6
#> F2R              7.614713e-29 2.563457856 0.574 0.105  2.284414e-25       6
#> SCG5             1.313098e-28 1.847995973 0.407 0.044  3.939295e-25       6
#> IFITM21          1.793353e-28 2.090551860 0.963 0.365  5.380059e-25       6
#> PSMB91           2.180323e-28 1.588869186 0.833 0.210  6.540968e-25       6
#> UBB1             2.234752e-28 1.775088376 1.000 0.260  6.704255e-25       6
#> HLA-A2           7.825838e-28 1.693311218 1.000 0.401  2.347751e-24       6
#> IFITM12          3.584521e-27 1.987790526 0.889 0.299  1.075356e-23       6
#> GRP              3.597217e-27 1.651742393 0.463 0.061  1.079165e-23       6
#> SMIM3            1.015054e-26 1.802022430 0.389 0.043  3.045162e-23       6
#> HEY1             1.072741e-26 1.783429902 0.278 0.019  3.218223e-23       6
#> IFITM33          1.558777e-26 1.626829333 1.000 0.421  4.676331e-23       6
#> MIR4435-2HG2     2.458896e-26 1.938035018 0.963 0.397  7.376687e-23       6
#> RAMP1            2.977966e-26 1.704133147 0.315 0.027  8.933897e-23       6
#> APOL1            3.415209e-26 1.792114932 0.315 0.027  1.024563e-22       6
#> HSPG2            4.148310e-26 2.100407278 0.574 0.115  1.244493e-22       6
#> NRN1             5.334173e-26 1.975892728 0.426 0.057  1.600252e-22       6
#> IGFBP6           8.260370e-26 1.678285360 0.407 0.050  2.478111e-22       6
#> EPAS11           1.067094e-25 2.402559959 0.648 0.161  3.201282e-22       6
#> BST21            1.945856e-25 1.932898811 0.815 0.240  5.837568e-22       6
#> PCDH7            3.357824e-25 2.104021254 0.407 0.054  1.007347e-21       6
#> B2M1             3.946753e-25 1.532875679 1.000 0.417  1.184026e-21       6
#> HOXB2            4.172805e-25 1.808135570 0.315 0.029  1.251842e-21       6
#> PLAU1            6.382766e-25 1.778741231 0.722 0.179  1.914830e-21       6
#> HIGD1B           8.146749e-25 1.754208067 0.259 0.018  2.444025e-21       6
#> STEAP1           9.235987e-25 1.613883398 0.426 0.057  2.770796e-21       6
#> FNDC1            1.072439e-24 2.360026190 0.519 0.097  3.217317e-21       6
#> TNFSF101         1.141252e-24 1.836774130 0.556 0.109  3.423757e-21       6
#> NBL12            3.132193e-24 2.584093248 0.889 0.378  9.396579e-21       6
#> LMCD11           3.689428e-24 2.429918321 0.704 0.213  1.106829e-20       6
#> PALLD2           8.099558e-24 2.091952025 0.870 0.349  2.429867e-20       6
#> TMEM45A1         9.491132e-24 1.863180946 0.667 0.177  2.847340e-20       6
#> JAG1             1.488408e-23 1.918253060 0.389 0.052  4.465223e-20       6
#> NREP1            2.572383e-23 2.146310034 0.722 0.241  7.717150e-20       6
#> OAF1             2.990027e-23 1.933457821 0.685 0.197  8.970082e-20       6
#> COPZ21           3.480416e-23 2.024900981 0.704 0.212  1.044125e-19       6
#> EVA1A            3.605502e-23 1.829916802 0.519 0.102  1.081651e-19       6
#> HOPX2            1.510481e-22 2.166445678 0.537 0.117  4.531444e-19       6
#> IFI44L           3.060843e-22 1.484750970 0.463 0.077  9.182529e-19       6
#> LXN1             3.309060e-22 1.601010657 0.519 0.100  9.927180e-19       6
#> NNMT2            3.546704e-22 1.674162855 0.907 0.405  1.064011e-18       6
#> TMEM2041         4.399771e-22 1.953447852 0.574 0.138  1.319931e-18       6
#> NTRK2            8.244491e-22 1.788086930 0.259 0.022  2.473347e-18       6
#> FMOD             9.307510e-22 1.354857547 0.444 0.071  2.792253e-18       6
#> CYP1B11          1.961917e-21 1.909814334 0.630 0.165  5.885752e-18       6
#> PTN              1.991222e-21 1.687613752 0.370 0.051  5.973665e-18       6
#> NEXN2            2.672316e-21 1.865034504 0.500 0.102  8.016949e-18       6
#> LOXL12           2.898317e-21 2.642345809 0.630 0.195  8.694952e-18       6
#> ENPP2            4.295200e-21 1.169329638 0.278 0.026  1.288560e-17       6
#> PLEKHA4          4.957426e-21 1.710145221 0.315 0.036  1.487228e-17       6
#> FKBP11           6.448794e-21 2.140238234 0.667 0.222  1.934638e-17       6
#> RAB311           1.205569e-20 1.634137841 0.870 0.314  3.616707e-17       6
#> PDPN1            2.867833e-20 1.384012200 0.389 0.058  8.603499e-17       6
#> MARCKS1          3.712467e-20 1.379520240 0.981 0.334  1.113740e-16       6
#> SEPT41           3.765173e-20 2.086999189 0.500 0.111  1.129552e-16       6
#> GEM              5.628582e-20 1.608705593 0.481 0.096  1.688575e-16       6
#> LINC001522       2.784185e-19 1.673832976 0.889 0.407  8.352556e-16       6
#> SPOCK1           6.056598e-19 1.731757204 0.389 0.066  1.816979e-15       6
#> SLC2A31          8.156849e-19 1.366382853 0.704 0.206  2.447055e-15       6
#> IFI27L21         1.315453e-18 1.393060146 0.759 0.255  3.946359e-15       6
#> OLFML2A          1.356270e-18 1.968126963 0.315 0.043  4.068810e-15       6
#> CLEC2B1          3.236067e-18 1.418499691 0.741 0.236  9.708200e-15       6
#> TRIM221          4.265295e-18 1.336407038 0.426 0.079  1.279588e-14       6
#> XAF11            5.272729e-18 1.269726988 0.463 0.093  1.581819e-14       6
#> PTGIR            5.570882e-18 1.416538744 0.278 0.033  1.671265e-14       6
#> SPON23           8.879002e-18 1.560786021 0.833 0.322  2.663701e-14       6
#> C1S3             9.025835e-18 1.690019414 0.815 0.358  2.707750e-14       6
#> ISG151           1.333629e-17 1.527718058 0.889 0.350  4.000887e-14       6
#> MMP21            1.691475e-17 1.814293797 0.685 0.230  5.074426e-14       6
#> RFTN1            1.963525e-17 1.237819742 0.444 0.090  5.890576e-14       6
#> TPPP3            3.339400e-17 2.155014568 0.352 0.061  1.001820e-13       6
#> CXCL121          4.571498e-17 1.814791665 0.611 0.206  1.371450e-13       6
#> TMEM47           4.685345e-17 1.484603955 0.278 0.035  1.405604e-13       6
#> CD200            4.776845e-17 1.650787436 0.426 0.090  1.433053e-13       6
#> GBP11            4.998983e-17 1.058252443 0.444 0.089  1.499695e-13       6
#> HLA-B3           6.942665e-17 1.087924873 1.000 0.475  2.082800e-13       6
#> HSD17B111        8.070967e-17 1.036249341 0.722 0.207  2.421290e-13       6
#> PSMB81           1.269695e-16 1.271698272 0.815 0.315  3.809085e-13       6
#> NR2F1-AS1        1.435802e-16 1.671310870 0.315 0.049  4.307407e-13       6
#> RP11-394O4.5     3.306230e-16 1.448066417 0.315 0.049  9.918689e-13       6
#> FRZB3            3.470024e-16 1.124010847 0.722 0.238  1.041007e-12       6
#> FAM198B1         5.071721e-16 1.230884469 0.519 0.132  1.521516e-12       6
#> MSRB33           9.784643e-16 1.702737951 0.630 0.234  2.935393e-12       6
#> CST31            1.083348e-15 1.100291101 0.963 0.416  3.250045e-12       6
#> LOXL21           1.105505e-15 1.988513526 0.537 0.174  3.316515e-12       6
#> DAB21            2.087425e-15 1.004069862 0.741 0.247  6.262276e-12       6
#> MYL93            7.558678e-15 1.115856639 0.963 0.474  2.267603e-11       6
#> MARVELD1         7.872436e-15 1.537155361 0.704 0.311  2.361731e-11       6
#> PLA2G161         9.742872e-15 1.205204529 0.481 0.126  2.922862e-11       6
#> SELM3            1.012885e-14 1.076253295 0.889 0.523  3.038654e-11       6
#> TPM14            1.232894e-14 1.087767452 0.926 0.544  3.698682e-11       6
#> PCDH18           1.278487e-14 1.744023342 0.352 0.072  3.835461e-11       6
#> MX2              1.595388e-14 1.228291617 0.407 0.091  4.786163e-11       6
#> EPSTI11          2.032391e-14 0.918753990 0.481 0.118  6.097172e-11       6
#> SAMD91           2.172163e-14 1.168413686 0.444 0.109  6.516490e-11       6
#> COLEC121         3.057480e-14 1.286416607 0.444 0.111  9.172440e-11       6
#> DIO21            3.883797e-14 1.589845591 0.278 0.044  1.165139e-10       6
#> EPB41L21         5.733336e-14 1.350572687 0.648 0.250  1.720001e-10       6
#> HES42            6.744518e-14 1.503438241 0.796 0.387  2.023356e-10       6
#> ZEB21            2.100255e-13 0.949885902 0.630 0.202  6.300764e-10       6
#> ISLR1            2.312754e-13 1.589402499 0.556 0.199  6.938262e-10       6
#> FBLN5            2.339015e-13 1.340265427 0.259 0.040  7.017045e-10       6
#> C1orf541         2.382232e-13 0.896096643 0.685 0.221  7.146696e-10       6
#> CAMK2N11         2.644735e-13 1.561816957 0.463 0.140  7.934204e-10       6
#> ADA              2.835176e-13 1.235594939 0.315 0.060  8.505529e-10       6
#> SLC40A11         2.858027e-13 0.969590628 0.500 0.143  8.574081e-10       6
#> VIM1             3.219810e-13 0.957602620 0.963 0.521  9.659429e-10       6
#> ACTA23           3.342290e-13 1.128099783 0.833 0.337  1.002687e-09       6
#> PPP1R14A2        3.603820e-13 1.662069054 0.537 0.188  1.081146e-09       6
#> IFIT3            3.618800e-13 1.159551525 0.407 0.102  1.085640e-09       6
#> HNMT1            4.050702e-13 0.864433614 0.778 0.255  1.215211e-09       6
#> DUSP11           4.062703e-13 1.058495796 0.907 0.356  1.218811e-09       6
#> FBLN11           4.244702e-13 1.820490478 0.611 0.256  1.273411e-09       6
#> GAS63            1.344062e-12 1.237786935 0.852 0.423  4.032185e-09       6
#> SERPINH13        1.373557e-12 1.023715562 0.963 0.520  4.120670e-09       6
#> MVP1             2.285379e-12 1.134307473 0.685 0.282  6.856138e-09       6
#> CTHRC12          2.871430e-12 1.295695701 0.759 0.462  8.614290e-09       6
#> LPAR61           2.882574e-12 0.833683473 0.481 0.131  8.647722e-09       6
#> FKBP104          3.699118e-12 1.132646348 0.815 0.436  1.109735e-08       6
#> EPS81            4.350853e-12 1.548158612 0.537 0.209  1.305256e-08       6
#> CD401            7.265599e-12 0.933673557 0.389 0.095  2.179680e-08       6
#> FHL13            9.147104e-12 0.902428636 0.556 0.187  2.744131e-08       6
#> EGR1             1.134431e-11 1.329035781 0.759 0.414  3.403293e-08       6
#> BASP12           1.291085e-11 1.080269115 0.722 0.284  3.873254e-08       6
#> GJA12            1.717364e-11 1.479981900 0.481 0.174  5.152092e-08       6
#> PMEPA11          2.123321e-11 1.385522215 0.741 0.396  6.369962e-08       6
#> PLPP31           3.902525e-11 1.473325899 0.444 0.149  1.170757e-07       6
#> LBH              4.042381e-11 1.453764492 0.407 0.123  1.212714e-07       6
#> IFIT11           6.147413e-11 0.971310395 0.426 0.127  1.844224e-07       6
#> IFI63            8.078196e-11 1.189410797 0.778 0.415  2.423459e-07       6
#> ARID5B1          8.921355e-11 1.113748330 0.778 0.398  2.676406e-07       6
#> COL11A13         9.359773e-11 1.413753959 0.722 0.370  2.807932e-07       6
#> GREM1            1.886583e-10 1.528859524 0.259 0.054  5.659749e-07       6
#> IFI441           2.523591e-10 0.911500819 0.389 0.110  7.570773e-07       6
#> EMP31            2.647675e-10 0.708628173 0.796 0.271  7.943026e-07       6
#> TAGLN3           3.164228e-10 0.932101151 0.796 0.430  9.492685e-07       6
#> DUT1             3.403318e-10 0.850999423 0.852 0.451  1.020996e-06       6
#> IFI351           4.333120e-10 0.990683827 0.444 0.145  1.299936e-06       6
#> HSD17B41         5.440115e-10 0.773228789 0.407 0.115  1.632035e-06       6
#> TGFB11           5.551486e-10 0.988357866 0.537 0.198  1.665446e-06       6
#> PTPRE1           5.823962e-10 0.876276418 0.648 0.266  1.747188e-06       6
#> PDLIM11          6.335358e-10 1.483743792 0.611 0.328  1.900607e-06       6
#> OAS2             7.405032e-10 0.822700205 0.278 0.060  2.221510e-06       6
#> HLA-C3           8.723895e-10 0.711390486 0.963 0.506  2.617168e-06       6
#> CCDC85B1         8.751372e-10 0.998571841 0.759 0.406  2.625412e-06       6
#> GJB2             2.010671e-09 0.903836699 0.259 0.057  6.032014e-06       6
#> TGM2             3.379568e-09 1.345540676 0.296 0.079  1.013870e-05       6
#> NEURL1B          3.427788e-09 1.214846386 0.259 0.060  1.028337e-05       6
#> S100A162         3.486227e-09 0.714628398 0.722 0.279  1.045868e-05       6
#> BHLHE411         4.813008e-09 0.712062213 0.333 0.088  1.443902e-05       6
#> GADD45A          5.863320e-09 1.217827795 0.500 0.216  1.758996e-05       6
#> ANO11            6.258178e-09 1.166373076 0.407 0.145  1.877453e-05       6
#> SAMD9L1          6.768022e-09 0.680943345 0.333 0.089  2.030406e-05       6
#> CCL21            7.062202e-09 0.833694214 0.407 0.132  2.118661e-05       6
#> SULF22           7.166528e-09 1.173305094 0.611 0.325  2.149958e-05       6
#> HLA-F3           7.485530e-09 0.952290327 0.667 0.339  2.245659e-05       6
#> CHPF2            9.570727e-09 1.024340464 0.741 0.433  2.871218e-05       6
#> NDUFA4L2         1.404299e-08 1.997350518 0.296 0.086  4.212898e-05       6
#> PARP141          1.702727e-08 0.893830172 0.500 0.205  5.108180e-05       6
#> CRIP11           1.866462e-08 1.193500108 0.407 0.151  5.599385e-05       6
#> BST1             2.378789e-08 0.901602759 0.278 0.072  7.136368e-05       6
#> SPHK12           3.482011e-08 0.990015827 0.704 0.424  1.044603e-04       6
#> TSC22D31         4.443804e-08 0.976154675 0.519 0.230  1.333141e-04       6
#> RGS5             5.817301e-08 1.831941510 0.333 0.116  1.745190e-04       6
#> CTD-3193K9.4     6.541690e-08 1.339946744 0.259 0.070  1.962507e-04       6
#> CTGF2            8.748887e-08 1.085457442 0.667 0.363  2.624666e-04       6
#> MYO1B2           1.462006e-07 0.936028617 0.704 0.401  4.386017e-04       6
#> BICD1            1.978048e-07 1.278351534 0.278 0.084  5.934144e-04       6
#> FLNA2            2.243022e-07 0.726124296 0.778 0.454  6.729067e-04       6
#> TGFBI1           5.388809e-07 0.308131597 0.648 0.230  1.616643e-03       6
#> IRF71            5.777857e-07 0.668831710 0.278 0.080  1.733357e-03       6
#> C12orf751        5.793340e-07 0.900215451 0.519 0.271  1.738002e-03       6
#> CTSB1            1.001164e-06 0.432209637 0.704 0.314  3.003493e-03       6
#> NR2F24           1.207814e-06 0.998493399 0.704 0.453  3.623443e-03       6
#> NRP11            1.588253e-06 0.726624904 0.519 0.236  4.764758e-03       6
#> TRIB22           1.793136e-06 0.923741004 0.463 0.240  5.379407e-03       6
#> MX12             1.883350e-06 0.799610702 0.352 0.140  5.650050e-03       6
#> EMP11            2.521334e-06 0.686267647 0.722 0.453  7.564003e-03       6
#> ELOVL51          3.551816e-06 1.008948724 0.481 0.268  1.065545e-02       6
#> TAP11            3.615639e-06 0.908879573 0.519 0.281  1.084692e-02       6
#> SKA21            4.510344e-06 0.689402061 0.556 0.299  1.353103e-02       6
#> SDC22            4.956561e-06 0.596416648 0.722 0.465  1.486968e-02       6
#> PDIA32           5.398897e-06 0.550460438 0.778 0.498  1.619669e-02       6
#> CYR613           5.940546e-06 1.187540030 0.648 0.416  1.782164e-02       6
#> TJP12            6.168254e-06 1.006372088 0.500 0.301  1.850476e-02       6
#> THBS13           7.796138e-06 1.154738316 0.556 0.342  2.338841e-02       6
#> LY961            7.983465e-06 0.366537103 0.611 0.245  2.395040e-02       6
#> OAS11            1.053944e-05 0.764193599 0.352 0.144  3.161831e-02       6
#> SLC39A142        1.087784e-05 1.172311093 0.352 0.167  3.263353e-02       6
#> RHOBTB32         1.122315e-05 0.584705493 0.796 0.440  3.366944e-02       6
#> PHLDA32          1.169427e-05 0.802963652 0.556 0.344  3.508282e-02       6
#> SNAI22           1.194789e-05 0.646549959 0.648 0.327  3.584367e-02       6
#> ID42             1.338946e-05 0.684584084 0.500 0.260  4.016838e-02       6
#> LY6E3            1.441752e-05 0.607648464 0.704 0.466  4.325255e-02       6
#> CTSL1            1.678500e-05 0.407792214 0.630 0.376  5.035499e-02       6
#> FAM46A           1.879787e-05 0.946068782 0.519 0.306  5.639362e-02       6
#> RGS162           1.905112e-05 0.634021823 0.611 0.311  5.715336e-02       6
#> KDELC21          2.018497e-05 0.972121519 0.370 0.182  6.055492e-02       6
#> STK17A           2.474577e-05 0.852460461 0.463 0.261  7.423730e-02       6
#> CD552            2.667440e-05 0.692919978 0.704 0.462  8.002321e-02       6
#> CD1091           2.931006e-05 0.509787865 0.296 0.110  8.793018e-02       6
#> MYOF1            4.881756e-05 0.808694594 0.481 0.287  1.464527e-01       6
#> CYBA1            5.636199e-05 0.368879065 0.722 0.387  1.690860e-01       6
#> AKR1B11          5.711109e-05 0.412865692 0.611 0.294  1.713333e-01       6
#> WLS              5.966935e-05 0.664351857 0.296 0.124  1.790081e-01       6
#> SUSD22           6.914588e-05 0.952532709 0.259 0.102  2.074377e-01       6
#> PGAM1            8.441232e-05 0.445275913 0.685 0.469  2.532370e-01       6
#> SEPT61           1.048648e-04 0.617476206 0.352 0.162  3.145944e-01       6
#> MEST2            1.109400e-04 0.703776992 0.296 0.128  3.328199e-01       6
#> FCGRT1           1.892428e-04 0.427323563 0.685 0.355  5.677284e-01       6
#> CASP42           2.093475e-04 0.499511912 0.519 0.286  6.280424e-01       6
#> PRDM11           2.094907e-04 0.402508748 0.352 0.153  6.284720e-01       6
#> MT2A3            3.446939e-04 0.386327964 0.685 0.398  1.000000e+00       6
#> TPBG2            5.127947e-04 0.630434310 0.537 0.374  1.000000e+00       6
#> GBP21            5.601882e-04 0.449816030 0.426 0.226  1.000000e+00       6
#> PXDN2            6.292087e-04 0.584998129 0.574 0.386  1.000000e+00       6
#> MEF2C1           7.551519e-04 0.750608155 0.407 0.237  1.000000e+00       6
#> LIMD22           8.201571e-04 0.597191762 0.407 0.232  1.000000e+00       6
#> LDLR1            8.227439e-04 0.779834260 0.296 0.150  1.000000e+00       6
#> PSIP11           8.993141e-04 0.607240570 0.519 0.344  1.000000e+00       6
#> PRRX22           9.573065e-04 0.536965690 0.352 0.191  1.000000e+00       6
#> P3H42            1.084100e-03 0.566999984 0.556 0.385  1.000000e+00       6
#> SCCPDH3          1.138040e-03 0.551852009 0.500 0.333  1.000000e+00       6
#> DSE1             1.212548e-03 0.519791515 0.444 0.255  1.000000e+00       6
#> PAPSS2           1.517310e-03 0.580732984 0.333 0.186  1.000000e+00       6
#> PDGFA2           1.566663e-03 0.955930177 0.333 0.201  1.000000e+00       6
#> PRAF22           1.639216e-03 0.565679695 0.537 0.381  1.000000e+00       6
#> MAF1             1.789172e-03 0.515452078 0.352 0.202  1.000000e+00       6
#> GALM1            2.002115e-03 0.391612994 0.333 0.170  1.000000e+00       6
#> NAAA1            2.031431e-03 0.651361255 0.407 0.253  1.000000e+00       6
#> BCAT11           3.114260e-03 0.343305489 0.259 0.123  1.000000e+00       6
#> ARHGEF121        3.128792e-03 0.591693367 0.500 0.373  1.000000e+00       6
#> RARRES13         3.154821e-03 0.405968645 0.463 0.267  1.000000e+00       6
#> ARID5A2          3.229490e-03 0.605011097 0.444 0.311  1.000000e+00       6
#> CCDC28B          4.302149e-03 0.591186079 0.333 0.205  1.000000e+00       6
#> ARHGAP181        5.028938e-03 0.249385297 0.426 0.230  1.000000e+00       6
#> LMO42            5.080365e-03 0.564902401 0.463 0.336  1.000000e+00       6
#> ARL4C1           5.169520e-03 0.371413998 0.574 0.341  1.000000e+00       6
#> UBE2L62          5.812692e-03 0.450269494 0.537 0.362  1.000000e+00       6
#> MMP91            6.029580e-03 0.475650495 0.481 0.221  1.000000e+00       6
#> MAP7D31          6.453011e-03 0.571973218 0.407 0.288  1.000000e+00       6
#> MAP1B3           6.709486e-03 0.536710558 0.574 0.420  1.000000e+00       6
#> PHLDA12          6.710439e-03 0.372741667 0.630 0.401  1.000000e+00       6
#> PPP1R12A1        7.952746e-03 0.363010030 0.648 0.453  1.000000e+00       6
#> LDHB1            8.130576e-03 0.342745487 0.722 0.489  1.000000e+00       6
#> ZFP361           9.546519e-03 0.385401797 0.630 0.405  1.000000e+00       6
#> FCER1A           6.304196e-77 2.772432858 0.408 0.005  1.891259e-73       7
#> CD1E             1.332495e-69 2.541608377 0.347 0.003  3.997484e-66       7
#> CD1A             1.046284e-63 2.135410986 0.306 0.002  3.138853e-60       7
#> CLEC10A          2.176725e-54 1.811671161 0.265 0.002  6.530175e-51       7
#> CD1C             6.947607e-51 1.854845307 0.265 0.003  2.084282e-47       7
#> ASGR2            4.367738e-49 1.843507708 0.286 0.005  1.310321e-45       7
#> GAPT             1.975682e-46 1.956022211 0.347 0.011  5.927047e-43       7
#> ABI31            3.693230e-45 2.183826173 0.837 0.120  1.107969e-41       7
#> JAML             4.859610e-45 2.701840169 0.551 0.047  1.457883e-41       7
#> IL1RN1           1.338100e-43 2.126285963 0.592 0.055  4.014301e-40       7
#> CXorf211         1.214437e-39 2.034325629 0.653 0.078  3.643310e-36       7
#> MNDA1            6.373534e-39 2.173458682 0.755 0.121  1.912060e-35       7
#> NCKAP1L1         2.767745e-38 1.640400878 0.857 0.136  8.303235e-35       7
#> ALOX51           1.245185e-37 2.146396644 0.735 0.115  3.735555e-34       7
#> HLA-DQA21        2.138415e-37 2.105182053 0.959 0.191  6.415246e-34       7
#> CRIP12           7.041378e-37 1.787263198 0.796 0.134  2.112414e-33       7
#> IL1B1            1.179397e-36 1.768536096 0.673 0.091  3.538192e-33       7
#> ADAM281          2.650686e-36 1.788314406 0.714 0.102  7.952057e-33       7
#> HLA-DRB51        2.824591e-36 2.220519181 0.980 0.216  8.473772e-33       7
#> PKIB1            8.676267e-36 2.352657089 0.857 0.173  2.602880e-32       7
#> CLEC4A1          4.082116e-35 1.842530647 0.735 0.115  1.224635e-31       7
#> TNFRSF9          8.421169e-35 1.614782901 0.469 0.041  2.526351e-31       7
#> HLA-DQA11        1.265890e-34 2.043370348 0.959 0.190  3.797669e-31       7
#> SAMHD11          1.640156e-34 1.845126869 0.939 0.198  4.920469e-31       7
#> HCK1             2.503260e-34 1.670827532 0.837 0.151  7.509779e-31       7
#> CD531            3.461599e-34 1.724828414 0.939 0.189  1.038480e-30       7
#> LST11            3.498771e-34 1.999895815 0.918 0.198  1.049631e-30       7
#> PYCARD1          6.479092e-34 1.958130429 0.959 0.217  1.943728e-30       7
#> LCP11            9.871747e-34 2.036889020 0.878 0.182  2.961524e-30       7
#> NRG1             1.009182e-33 1.878425271 0.306 0.015  3.027545e-30       7
#> HCLS11           1.096219e-33 1.809479713 0.898 0.179  3.288658e-30       7
#> CD481            2.335869e-33 1.608340828 0.837 0.151  7.007608e-30       7
#> FGL21            3.781942e-33 1.870400993 0.592 0.076  1.134583e-29       7
#> RNASE61          4.727031e-33 1.799540459 0.796 0.142  1.418109e-29       7
#> CSF1R1           1.049515e-32 1.790079805 0.898 0.189  3.148546e-29       7
#> FCGBP            1.121467e-32 1.993359143 0.408 0.033  3.364402e-29       7
#> PTPRC1           1.765415e-32 1.813247318 0.918 0.198  5.296245e-29       7
#> CYTIP1           1.851352e-32 1.754765906 0.694 0.110  5.554055e-29       7
#> PRKCB            2.020216e-32 1.749102348 0.306 0.016  6.060649e-29       7
#> ALOX5AP1         3.130179e-32 2.142154847 0.816 0.160  9.390536e-29       7
#> WAS1             3.938740e-32 1.890252778 0.714 0.123  1.181622e-28       7
#> SPI11            8.555842e-32 1.879881323 0.918 0.199  2.566752e-28       7
#> MS4A6A1          1.300231e-31 1.604106115 0.837 0.169  3.900694e-28       7
#> CD861            2.202858e-31 1.539044790 0.776 0.136  6.608574e-28       7
#> DOK21            1.149321e-30 1.632162691 0.612 0.088  3.447964e-27       7
#> RGS11            1.301983e-30 1.954224532 0.918 0.185  3.905948e-27       7
#> FMNL11           2.394067e-30 1.713698001 0.816 0.172  7.182201e-27       7
#> SNX20            7.322836e-30 1.835310441 0.469 0.052  2.196851e-26       7
#> SRGN1            9.488879e-30 1.802858782 1.000 0.210  2.846664e-26       7
#> BTK1             1.097891e-29 1.557682324 0.592 0.083  3.293673e-26       7
#> HLA-DRB11        1.580055e-29 2.142820671 0.980 0.240  4.740164e-26       7
#> HLA-DQB11        2.470837e-29 2.177736823 0.959 0.248  7.412512e-26       7
#> FCN1             4.887970e-29 2.285211336 0.429 0.046  1.466391e-25       7
#> GMFG1            7.710209e-29 1.610057240 0.939 0.229  2.313063e-25       7
#> FCGR2B1          1.109941e-28 1.847445964 0.673 0.122  3.329822e-25       7
#> LILRB21          1.483555e-28 1.376097220 0.612 0.090  4.450666e-25       7
#> CSTA1            1.537698e-28 1.571575793 0.673 0.115  4.613094e-25       7
#> LSP11            2.134683e-28 1.998680909 0.918 0.225  6.404050e-25       7
#> CYBB1            2.603268e-28 1.516919442 0.857 0.182  7.809805e-25       7
#> GPR1831          2.738290e-28 2.063626238 0.878 0.245  8.214871e-25       7
#> MYO1F1           2.860625e-28 1.761614657 0.898 0.220  8.581875e-25       7
#> SLA1             2.949690e-28 1.433998665 0.694 0.123  8.849070e-25       7
#> ITGB21           3.089955e-28 1.970436451 0.959 0.237  9.269866e-25       7
#> TNFAIP81         3.409248e-28 1.652213691 0.837 0.197  1.022774e-24       7
#> HLA-DPA11        9.881073e-28 2.081808903 0.980 0.211  2.964322e-24       7
#> PLEK21           1.334307e-27 1.730545454 0.551 0.081  4.002922e-24       7
#> AIF11            1.457619e-27 1.776192125 0.959 0.200  4.372858e-24       7
#> CD831            2.394830e-27 1.812554324 0.816 0.193  7.184490e-24       7
#> HLA-DPB11        4.012684e-27 2.112278029 0.980 0.205  1.203805e-23       7
#> NAIP1            4.033842e-27 1.328945742 0.653 0.113  1.210152e-23       7
#> GNA151           4.051968e-27 1.730058701 0.694 0.136  1.215590e-23       7
#> NFKBID1          1.139455e-26 1.728161592 0.673 0.134  3.418364e-23       7
#> PLEK1            1.141237e-26 1.347203458 0.816 0.172  3.423710e-23       7
#> IL7R1            1.339247e-26 1.843254560 0.735 0.158  4.017741e-23       7
#> EVI2A1           2.287896e-26 1.358937537 0.735 0.142  6.863687e-23       7
#> CLECL11          2.622495e-26 1.459673405 0.531 0.073  7.867486e-23       7
#> IFI162           3.223525e-26 1.591073786 0.939 0.255  9.670576e-23       7
#> PARVG1           3.654475e-26 1.446946800 0.571 0.089  1.096342e-22       7
#> LAPTM51          3.963008e-26 1.696296798 1.000 0.215  1.188903e-22       7
#> IL2RG1           3.968411e-26 1.739872534 0.571 0.094  1.190523e-22       7
#> EID12            5.032481e-26 1.328153485 0.959 0.236  1.509744e-22       7
#> EVI2B1           5.095817e-26 1.461551011 0.755 0.158  1.528745e-22       7
#> IL10RA1          5.819887e-26 1.473917382 0.694 0.136  1.745966e-22       7
#> RRM2             6.446666e-26 2.289455371 0.429 0.053  1.934000e-22       7
#> DUSP21           7.265896e-26 2.412491354 0.673 0.159  2.179769e-22       7
#> RUNX32           9.071955e-26 2.151566534 0.612 0.120  2.721587e-22       7
#> HLA-DMA1         1.420433e-25 1.838661920 0.939 0.256  4.261298e-22       7
#> LPXN1            1.427430e-25 1.462763730 0.714 0.144  4.282290e-22       7
#> RTN1             1.460334e-25 1.553003738 0.429 0.050  4.381003e-22       7
#> KYNU1            1.718813e-25 1.416941714 0.735 0.151  5.156438e-22       7
#> CST7             1.862516e-25 2.134885944 0.388 0.043  5.587547e-22       7
#> TGFBI2           1.895996e-25 1.738854581 0.857 0.222  5.687989e-22       7
#> CHST2            2.058498e-25 1.847134668 0.429 0.052  6.175495e-22       7
#> HLA-DMB1         2.194951e-25 1.638400267 0.796 0.197  6.584854e-22       7
#> EMP32            2.688111e-25 1.798767628 0.959 0.266  8.064332e-22       7
#> ESCO2            3.085003e-25 2.007851436 0.367 0.038  9.255008e-22       7
#> HSD17B112        4.150338e-25 1.471589372 0.837 0.204  1.245101e-21       7
#> C10orf541        5.724965e-25 1.762309757 0.673 0.140  1.717490e-21       7
#> FERMT31          6.127298e-25 1.376065245 0.816 0.184  1.838190e-21       7
#> IFI301           6.206257e-25 1.389939663 0.755 0.162  1.861877e-21       7
#> CORO1A2          6.709981e-25 2.068753230 0.918 0.292  2.012994e-21       7
#> XIST2            1.170507e-24 1.349243087 0.918 0.231  3.511521e-21       7
#> IGSF61           2.001016e-24 1.242321920 0.673 0.130  6.003048e-21       7
#> UBB2             3.056268e-24 1.708625859 0.980 0.264  9.168805e-21       7
#> CASP11           3.193919e-24 1.661439191 0.735 0.187  9.581758e-21       7
#> ADAM81           4.101835e-24 1.493876427 0.633 0.123  1.230550e-20       7
#> HIST1H1D         5.293704e-24 2.070061577 0.388 0.046  1.588111e-20       7
#> CD300C           5.450017e-24 1.254032196 0.408 0.047  1.635005e-20       7
#> HLA-DRA1         7.941748e-24 2.009423078 0.980 0.225  2.382524e-20       7
#> FAM111A          8.038891e-24 1.950175069 0.531 0.093  2.411667e-20       7
#> CD331            1.180240e-23 1.421630065 0.490 0.071  3.540720e-20       7
#> MFNG             1.779811e-23 1.410274238 0.408 0.049  5.339432e-20       7
#> POU2F21          2.337859e-23 1.203493838 0.510 0.077  7.013576e-20       7
#> SERPINB91        4.415015e-23 1.274424177 0.490 0.073  1.324505e-19       7
#> LYZ1             4.993280e-23 1.491788458 0.837 0.206  1.497984e-19       7
#> CD741            7.492115e-23 1.946075303 0.980 0.230  2.247634e-19       7
#> PALD1            1.093091e-22 1.272432309 0.306 0.027  3.279274e-19       7
#> CCL3L31          1.312379e-22 1.555442005 0.592 0.115  3.937138e-19       7
#> CD371            1.493966e-22 1.335103047 0.735 0.164  4.481897e-19       7
#> CD41             1.969369e-22 1.238628700 0.796 0.188  5.908107e-19       7
#> ARHGDIB1         3.006486e-22 1.649741419 0.980 0.434  9.019457e-19       7
#> NCF41            3.510762e-22 1.261181230 0.612 0.116  1.053228e-18       7
#> SLC8A11          4.794800e-22 1.206292610 0.592 0.108  1.438440e-18       7
#> TNFAIP31         4.874339e-22 1.239808118 0.796 0.209  1.462302e-18       7
#> AOAH1            5.977603e-22 1.373827230 0.469 0.072  1.793281e-18       7
#> RAC21            6.547139e-22 1.288986772 0.592 0.116  1.964142e-18       7
#> BIN2             7.287542e-22 1.496896589 0.469 0.073  2.186263e-18       7
#> HSD17B42         1.251597e-21 1.143424384 0.592 0.108  3.754790e-18       7
#> HLA-DQB2         1.378151e-21 1.729939121 0.367 0.045  4.134452e-18       7
#> CARD161          1.395544e-21 1.297085528 0.796 0.205  4.186631e-18       7
#> CTSS1            2.180348e-21 1.512257090 0.939 0.284  6.541044e-18       7
#> CD300A1          2.829120e-21 1.284613566 0.490 0.078  8.487361e-18       7
#> PSMB92           3.250485e-21 1.598088716 0.776 0.216  9.751454e-18       7
#> MARCH11          3.757421e-21 1.095733063 0.694 0.150  1.127226e-17       7
#> EFHD21           8.795010e-21 1.923094146 0.898 0.318  2.638503e-17       7
#> UCP21            9.976542e-21 1.693175472 0.878 0.336  2.992963e-17       7
#> RGS21            1.052441e-20 1.577612606 0.878 0.249  3.157323e-17       7
#> SLC2A32          1.472981e-20 1.268953755 0.796 0.204  4.418943e-17       7
#> HCST1            1.982094e-20 1.337724995 0.735 0.188  5.946281e-17       7
#> GPR841           3.480285e-20 1.099673442 0.469 0.074  1.044085e-16       7
#> HLA-DOA1         3.956911e-20 0.995790447 0.490 0.079  1.187073e-16       7
#> MEF2C2           4.253134e-20 1.310702405 0.755 0.221  1.275940e-16       7
#> LY962            4.512995e-20 1.103292797 0.857 0.236  1.353899e-16       7
#> CSF2RA1          5.483990e-20 1.755153307 0.633 0.160  1.645197e-16       7
#> SEPT62           7.427580e-20 1.665495908 0.612 0.151  2.228274e-16       7
#> LY861            8.225272e-20 1.300654932 0.633 0.145  2.467581e-16       7
#> CD69             9.381191e-20 1.661498213 0.388 0.057  2.814357e-16       7
#> PSMB82           1.257063e-19 1.604231272 0.857 0.316  3.771189e-16       7
#> DUSP12           1.527151e-19 1.569798031 0.939 0.357  4.581454e-16       7
#> CASC5            1.657640e-19 1.346975551 0.327 0.038  4.972919e-16       7
#> PTGS11           1.739439e-19 1.072526637 0.490 0.083  5.218318e-16       7
#> DSE2             2.267313e-19 1.391497377 0.796 0.240  6.801939e-16       7
#> FCER1G1          2.671146e-19 1.431373511 0.918 0.227  8.013437e-16       7
#> PTAFR1           3.202765e-19 1.257330950 0.571 0.121  9.608295e-16       7
#> SLC7A71          3.382804e-19 1.057012408 0.673 0.152  1.014841e-15       7
#> CYBA2            3.395717e-19 1.521981908 1.000 0.376  1.018715e-15       7
#> CLEC5A1          5.499001e-19 0.926079046 0.449 0.071  1.649700e-15       7
#> CENPE            8.410395e-19 1.253837109 0.429 0.070  2.523118e-15       7
#> BST22            9.383276e-19 1.043934268 0.857 0.240  2.814983e-15       7
#> CD300LF          1.027816e-18 1.214432853 0.265 0.025  3.083449e-15       7
#> GPR651           1.066417e-18 0.980860046 0.592 0.123  3.199250e-15       7
#> TBXAS11          1.127007e-18 1.191927421 0.653 0.159  3.381021e-15       7
#> PARM1            1.133978e-18 1.390751722 0.408 0.065  3.401933e-15       7
#> ADAP21           1.292062e-18 1.033653599 0.694 0.163  3.876185e-15       7
#> CD841            1.974955e-18 0.749140547 0.735 0.169  5.924866e-15       7
#> CASP43           2.100625e-18 1.374966409 0.837 0.272  6.301874e-15       7
#> RASAL3           3.498722e-18 1.183339012 0.286 0.031  1.049616e-14       7
#> TYROBP1          4.162600e-18 1.547435154 0.939 0.210  1.248780e-14       7
#> RNASET22         4.299799e-18 1.467739942 0.918 0.365  1.289940e-14       7
#> SAMSN11          4.684932e-18 1.071853238 0.653 0.152  1.405480e-14       7
#> FPR31            6.198077e-18 0.959565910 0.714 0.173  1.859423e-14       7
#> LGALS91          6.238494e-18 1.214059065 0.755 0.212  1.871548e-14       7
#> CENPK            7.752456e-18 1.647264122 0.449 0.088  2.325737e-14       7
#> CD141            7.988430e-18 1.428883438 0.857 0.222  2.396529e-14       7
#> IL4I11           8.963955e-18 1.077749656 0.653 0.156  2.689186e-14       7
#> MS4A71           1.073438e-17 0.929031220 0.735 0.189  3.220315e-14       7
#> CKLF1            1.223673e-17 1.664159343 0.898 0.413  3.671018e-14       7
#> GPR341           1.302391e-17 0.890509825 0.612 0.136  3.907174e-14       7
#> IFI352           2.311163e-17 1.106903963 0.592 0.139  6.933490e-14       7
#> CLSPN            2.754935e-17 1.608595011 0.408 0.073  8.264805e-14       7
#> PRDM12           3.077109e-17 1.159357833 0.592 0.143  9.231327e-14       7
#> FGR1             3.150007e-17 1.199304897 0.388 0.061  9.450022e-14       7
#> GALM2            5.206191e-17 1.241064062 0.612 0.157  1.561857e-13       7
#> PPM1N1           5.503066e-17 1.492682712 0.388 0.066  1.650920e-13       7
#> CST6             6.433762e-17 1.334219654 0.286 0.034  1.930129e-13       7
#> SCARF1           1.025920e-16 1.081632198 0.347 0.050  3.077759e-13       7
#> MTHFD22          1.088992e-16 1.381734488 0.898 0.399  3.266976e-13       7
#> KIF11            1.145002e-16 1.549422869 0.367 0.060  3.435006e-13       7
#> SLC16A101        1.321278e-16 0.871818086 0.510 0.102  3.963833e-13       7
#> RP11-1143G9.41   1.427741e-16 1.313043048 0.653 0.174  4.283223e-13       7
#> ZFP362           1.525471e-16 1.466204490 0.898 0.393  4.576413e-13       7
#> SH2B31           1.589483e-16 1.072093920 0.510 0.106  4.768450e-13       7
#> KCTD121          1.904566e-16 1.022108547 0.714 0.193  5.713697e-13       7
#> GTSE1            1.936754e-16 1.582867409 0.429 0.086  5.810263e-13       7
#> HMMR             2.251014e-16 1.348170945 0.449 0.092  6.753043e-13       7
#> HLA-A3           2.380959e-16 1.276640329 1.000 0.404  7.142877e-13       7
#> FYB1             2.530024e-16 1.087723255 0.837 0.263  7.590072e-13       7
#> CXCR41           2.749163e-16 1.298484465 0.673 0.191  8.247488e-13       7
#> SAMD9L2          2.909750e-16 1.103297176 0.449 0.085  8.729251e-13       7
#> C1orf1621        3.183121e-16 1.060456789 0.612 0.153  9.549362e-13       7
#> SGK11            3.636983e-16 1.469883757 0.796 0.313  1.091095e-12       7
#> CLEC2B2          3.902704e-16 1.101792037 0.776 0.237  1.170811e-12       7
#> TFEC1            4.411828e-16 0.914022679 0.531 0.113  1.323548e-12       7
#> TGFB12           4.900439e-16 1.096427732 0.694 0.193  1.470132e-12       7
#> CCR11            5.598172e-16 1.094558400 0.531 0.119  1.679451e-12       7
#> ZEB22            6.251839e-16 1.045234180 0.714 0.200  1.875552e-12       7
#> FAM26F1          6.600717e-16 1.287616070 0.653 0.192  1.980215e-12       7
#> PLAUR1           6.777894e-16 1.218296235 0.857 0.312  2.033368e-12       7
#> CENPM1           7.150523e-16 1.711454640 0.490 0.120  2.145157e-12       7
#> NR4A21           7.228734e-16 1.346483600 0.878 0.376  2.168620e-12       7
#> SERPINF12        7.694009e-16 0.914575561 0.735 0.219  2.308203e-12       7
#> DOCK101          7.708154e-16 1.012098747 0.469 0.094  2.312446e-12       7
#> F13A1            8.759641e-16 0.890118302 0.327 0.047  2.627892e-12       7
#> PLAU2            9.206496e-16 0.981617939 0.673 0.184  2.761949e-12       7
#> OSM1             9.744112e-16 1.072825387 0.429 0.082  2.923234e-12       7
#> CTSC1            1.308570e-15 1.068734079 0.898 0.299  3.925709e-12       7
#> CD681            1.413207e-15 0.996233441 0.878 0.226  4.239620e-12       7
#> FOXP12           1.440150e-15 1.147531230 0.714 0.220  4.320450e-12       7
#> HNMT2            1.562659e-15 1.108083202 0.816 0.256  4.687976e-12       7
#> CYTH41           1.990651e-15 0.836251074 0.551 0.124  5.971952e-12       7
#> LTB              2.087668e-15 1.395872477 0.469 0.106  6.263005e-12       7
#> FBP11            2.305138e-15 1.138235315 0.367 0.062  6.915414e-12       7
#> CDKN2C           2.353652e-15 1.577382362 0.490 0.124  7.060957e-12       7
#> RPPH1            3.168236e-15 1.026533686 0.367 0.063  9.504707e-12       7
#> PDE4B1           3.415840e-15 0.989662240 0.490 0.108  1.024752e-11       7
#> B2M2             3.830358e-15 1.225322831 1.000 0.420  1.149107e-11       7
#> LAIR11           3.845290e-15 0.957663555 0.633 0.168  1.153587e-11       7
#> LIMS11           4.002341e-15 1.272241996 0.878 0.342  1.200702e-11       7
#> C3AR11           4.102268e-15 1.044417149 0.571 0.142  1.230680e-11       7
#> ARHGAP182        4.520746e-15 1.156270456 0.714 0.218  1.356224e-11       7
#> SNX101           4.541829e-15 0.943719644 0.735 0.218  1.362549e-11       7
#> ALDH21           4.908314e-15 1.422692749 0.735 0.278  1.472494e-11       7
#> GPAT3            6.097022e-15 1.560938850 0.327 0.053  1.829107e-11       7
#> PTPRE2           6.124117e-15 1.085266415 0.776 0.261  1.837235e-11       7
#> PTPN221          9.112401e-15 0.795507416 0.388 0.069  2.733720e-11       7
#> HJURP            9.529640e-15 1.165812578 0.347 0.059  2.858892e-11       7
#> TCF42            1.206237e-14 0.808169478 0.796 0.250  3.618711e-11       7
#> ADAM19           1.298759e-14 1.457943262 0.265 0.035  3.896277e-11       7
#> NCAPG            1.330037e-14 1.242975036 0.388 0.074  3.990112e-11       7
#> LRRC251          1.359681e-14 0.935708240 0.490 0.109  4.079044e-11       7
#> HAVCR21          1.375739e-14 0.952949015 0.531 0.123  4.127217e-11       7
#> LILRB11          1.516801e-14 0.792961649 0.469 0.095  4.550403e-11       7
#> HIST1H1E         1.775280e-14 1.302424731 0.286 0.041  5.325839e-11       7
#> ARRB21           1.813523e-14 1.098863254 0.735 0.260  5.440570e-11       7
#> MKI671           2.293606e-14 1.447131309 0.490 0.129  6.880818e-11       7
#> IFI27L22         2.920498e-14 0.946369403 0.776 0.257  8.761493e-11       7
#> KIF15            2.976799e-14 1.241648512 0.347 0.061  8.930397e-11       7
#> ARHGAP11B        3.516348e-14 1.223633427 0.327 0.055  1.054904e-10       7
#> CPVL1            3.542902e-14 1.312346783 0.653 0.206  1.062871e-10       7
#> PTPN61           3.615205e-14 1.317487421 0.776 0.329  1.084561e-10       7
#> MPEG11           3.996595e-14 0.825912279 0.571 0.142  1.198978e-10       7
#> BCL2A11          4.641213e-14 1.044960462 0.388 0.073  1.392364e-10       7
#> AC092484.11      4.682844e-14 0.982246436 0.367 0.068  1.404853e-10       7
#> CCL4L21          5.209006e-14 1.183688915 0.571 0.155  1.562702e-10       7
#> SPC251           5.441495e-14 1.084324329 0.429 0.094  1.632449e-10       7
#> CDK11            5.893896e-14 1.750235463 0.490 0.138  1.768169e-10       7
#> HLX              6.261539e-14 1.031014636 0.306 0.048  1.878462e-10       7
#> PID1             7.693506e-14 1.285796779 0.306 0.050  2.308052e-10       7
#> C1QC1            7.777750e-14 1.186535646 0.755 0.191  2.333325e-10       7
#> CCNA21           7.825853e-14 1.213278878 0.429 0.098  2.347756e-10       7
#> CST32            7.875346e-14 1.106491875 0.918 0.421  2.362604e-10       7
#> PRR11            8.383958e-14 1.269033227 0.367 0.073  2.515187e-10       7
#> CKAP2L           8.879140e-14 1.515113593 0.265 0.038  2.663742e-10       7
#> COTL12           1.011963e-13 1.247894657 0.918 0.483  3.035888e-10       7
#> CD402            1.031486e-13 0.766555700 0.449 0.094  3.094459e-10       7
#> FCGR2A1          1.255917e-13 0.781694403 0.714 0.202  3.767750e-10       7
#> CAPG2            2.216943e-13 1.096868586 0.898 0.372  6.650830e-10       7
#> SLC25A52         2.472247e-13 1.242704972 0.878 0.478  7.416740e-10       7
#> ANPEP2           2.490678e-13 1.156698121 0.735 0.273  7.472033e-10       7
#> TOP2A1           2.551145e-13 1.714822359 0.469 0.132  7.653435e-10       7
#> CXCL81           2.637957e-13 0.844516050 0.551 0.150  7.913870e-10       7
#> NCAPH            2.782982e-13 1.259411202 0.265 0.039  8.348945e-10       7
#> OLFML2B2         2.912418e-13 0.591121882 0.592 0.152  8.737255e-10       7
#> IRF72            3.289058e-13 0.791782213 0.388 0.075  9.867173e-10       7
#> LMNB11           3.433407e-13 1.160771761 0.551 0.164  1.030022e-09       7
#> BMP2K1           3.699477e-13 0.943930358 0.673 0.208  1.109843e-09       7
#> C15orf481        3.734765e-13 1.359683667 0.694 0.251  1.120430e-09       7
#> ASF1B            4.289030e-13 1.256064580 0.408 0.094  1.286709e-09       7
#> STX111           4.389885e-13 0.726067033 0.347 0.060  1.316965e-09       7
#> RHOF             4.591901e-13 1.360914897 0.408 0.094  1.377570e-09       7
#> FAM105A1         5.305569e-13 1.291379995 0.714 0.286  1.591671e-09       7
#> RASGRP3          5.903274e-13 0.757645222 0.347 0.064  1.770982e-09       7
#> AURKB1           6.108513e-13 1.270712052 0.429 0.106  1.832554e-09       7
#> NFKBIA2          6.428698e-13 1.253905564 0.837 0.399  1.928610e-09       7
#> ID21             6.762353e-13 1.223429358 0.857 0.353  2.028706e-09       7
#> MS4A4A1          7.200098e-13 0.742041269 0.531 0.135  2.160029e-09       7
#> CDCA8            7.493974e-13 1.102896651 0.367 0.076  2.248192e-09       7
#> CELF22           7.503302e-13 0.982434007 0.776 0.282  2.250991e-09       7
#> EBI31            7.878192e-13 0.872062880 0.408 0.087  2.363458e-09       7
#> NUSAP11          9.066744e-13 1.556394980 0.469 0.137  2.720023e-09       7
#> RAB312           1.416618e-12 1.086725496 0.796 0.321  4.249853e-09       7
#> CD931            1.496844e-12 0.903760592 0.367 0.074  4.490532e-09       7
#> GPX12            1.668066e-12 1.050281534 0.878 0.476  5.004198e-09       7
#> PLXNC11          1.700513e-12 0.710156571 0.429 0.094  5.101540e-09       7
#> FCGRT2           1.811874e-12 1.040357049 0.878 0.347  5.435623e-09       7
#> LIMD23           1.917912e-12 1.397627681 0.633 0.222  5.753735e-09       7
#> IL181            1.935268e-12 0.656698568 0.653 0.190  5.805803e-09       7
#> PTGER41          1.990017e-12 1.032773574 0.469 0.122  5.970051e-09       7
#> SLC16A31         3.584574e-12 1.181952068 0.816 0.409  1.075372e-08       7
#> SERPINA11        3.590946e-12 0.870673991 0.837 0.312  1.077284e-08       7
#> SELPLG1          3.610755e-12 0.768558630 0.408 0.090  1.083226e-08       7
#> HBEGF1           3.770374e-12 1.013676182 0.408 0.095  1.131112e-08       7
#> CACNA2D41        4.883562e-12 0.629320288 0.429 0.094  1.465069e-08       7
#> STEAP11          5.341330e-12 0.977681963 0.327 0.063  1.602399e-08       7
#> SLAMF81          6.014543e-12 0.690151599 0.469 0.113  1.804363e-08       7
#> GIMAP41          6.360038e-12 0.827096770 0.551 0.154  1.908012e-08       7
#> SHCBP1           6.980094e-12 1.086722636 0.367 0.083  2.094028e-08       7
#> KIAA01011        8.224493e-12 1.331520393 0.592 0.216  2.467348e-08       7
#> GK1              1.010759e-11 0.975659163 0.531 0.161  3.032276e-08       7
#> FILIP1L2         1.106955e-11 0.488131240 0.469 0.114  3.320865e-08       7
#> MMP92            1.111622e-11 1.068192993 0.653 0.215  3.334865e-08       7
#> MYO1G1           1.390584e-11 1.063259902 0.327 0.067  4.171753e-08       7
#> MAFB1            1.574014e-11 1.055512054 0.755 0.292  4.722043e-08       7
#> CCL31            1.692571e-11 1.166839059 0.694 0.187  5.077712e-08       7
#> TYMS1            1.759636e-11 1.639869702 0.490 0.168  5.278908e-08       7
#> VSIG41           1.797299e-11 0.874937470 0.408 0.097  5.391897e-08       7
#> C1QA1            1.989962e-11 0.999510692 0.755 0.191  5.969886e-08       7
#> DLGAP5           2.127694e-11 1.056153363 0.327 0.068  6.383083e-08       7
#> MCOLN21          2.972241e-11 0.838053465 0.286 0.051  8.916724e-08       7
#> CDKN31           3.060609e-11 1.042711710 0.449 0.127  9.181826e-08       7
#> NCF11            3.133700e-11 0.760021018 0.571 0.186  9.401101e-08       7
#> CEP551           3.466625e-11 1.056240391 0.388 0.095  1.039987e-07       7
#> FTL1             3.827289e-11 0.857120538 0.878 0.296  1.148187e-07       7
#> FBXO5            4.424248e-11 1.307907567 0.306 0.064  1.327274e-07       7
#> CD300E1          4.519169e-11 0.876529882 0.306 0.060  1.355751e-07       7
#> LILRB41          4.623817e-11 0.535697508 0.571 0.158  1.387145e-07       7
#> APOBEC3C1        5.951870e-11 0.867343090 0.408 0.100  1.785561e-07       7
#> MRC1             6.030120e-11 0.930531467 0.265 0.046  1.809036e-07       7
#> RGS103           6.341555e-11 0.878758930 0.857 0.439  1.902467e-07       7
#> RAB7B            6.496041e-11 0.718240418 0.286 0.052  1.948812e-07       7
#> ASPM1            6.994539e-11 0.825721166 0.347 0.079  2.098362e-07       7
#> IER51            8.514182e-11 1.036458612 0.857 0.409  2.554255e-07       7
#> CENPA1           9.240121e-11 0.891822259 0.367 0.087  2.772036e-07       7
#> TSC22D32         9.523627e-11 0.928193405 0.633 0.226  2.857088e-07       7
#> GNG2             9.681613e-11 0.777327836 0.327 0.067  2.904484e-07       7
#> VIM2             1.026395e-10 0.899328636 0.898 0.526  3.079185e-07       7
#> RASSF41          1.200834e-10 1.159177027 0.714 0.329  3.602502e-07       7
#> ARL4C2           1.425823e-10 1.053791528 0.755 0.334  4.277468e-07       7
#> SPC24            1.498753e-10 1.081135921 0.265 0.050  4.496260e-07       7
#> OSCAR1           1.886852e-10 0.791246474 0.327 0.070  5.660555e-07       7
#> CD361            2.085720e-10 0.553267288 0.327 0.069  6.257159e-07       7
#> EPSTI12          2.286999e-10 0.773194831 0.449 0.121  6.860998e-07       7
#> OAS12            2.652612e-10 0.666862571 0.490 0.138  7.957835e-07       7
#> ANLN1            2.999321e-10 0.946650310 0.388 0.103  8.997964e-07       7
#> EGR2             3.037841e-10 0.999785539 0.367 0.093  9.113523e-07       7
#> NR4A3            3.199025e-10 0.869149575 0.265 0.050  9.597076e-07       7
#> DEPDC1           3.212409e-10 0.973025373 0.306 0.066  9.637227e-07       7
#> HIST1H4C2        3.883841e-10 1.675557416 0.714 0.327  1.165152e-06       7
#> MVP2             4.006215e-10 0.800177978 0.714 0.282  1.201865e-06       7
#> C1orf211         4.570156e-10 1.051873226 0.510 0.177  1.371047e-06       7
#> DHFR1            5.497155e-10 1.274392354 0.469 0.163  1.649146e-06       7
#> OLR11            5.536232e-10 0.831992549 0.327 0.074  1.660870e-06       7
#> PLTP1            6.441470e-10 0.757898846 0.612 0.202  1.932441e-06       7
#> GATM1            6.730380e-10 0.875998823 0.367 0.091  2.019114e-06       7
#> KIF20B1          7.658802e-10 1.076887979 0.490 0.176  2.297641e-06       7
#> IL27RA           8.292321e-10 0.835859989 0.265 0.051  2.487696e-06       7
#> ADAMDEC11        8.313882e-10 0.747569850 0.429 0.121  2.494164e-06       7
#> BRCA1            8.323569e-10 1.233137389 0.347 0.092  2.497071e-06       7
#> BIRC3            8.574909e-10 0.757454390 0.388 0.100  2.572473e-06       7
#> SLC39A81         8.639204e-10 0.720361222 0.592 0.207  2.591761e-06       7
#> KIF2C            1.074382e-09 1.005345831 0.306 0.071  3.223147e-06       7
#> AKR1B12          1.094718e-09 0.804505955 0.735 0.290  3.284155e-06       7
#> NCF21            1.150036e-09 0.754746154 0.408 0.111  3.450107e-06       7
#> APOBEC3G         1.193916e-09 0.793317746 0.286 0.059  3.581749e-06       7
#> ARHGAP11A        1.323049e-09 0.786950241 0.286 0.060  3.969147e-06       7
#> GHRL             1.336000e-09 1.148831505 0.286 0.063  4.007999e-06       7
#> BIRC51           1.375210e-09 1.107441219 0.469 0.158  4.125630e-06       7
#> RFTN11           1.498616e-09 0.853642717 0.367 0.095  4.495848e-06       7
#> AXL1             1.502176e-09 0.772135783 0.327 0.076  4.506529e-06       7
#> HIVEP3           1.528934e-09 0.927669793 0.286 0.062  4.586801e-06       7
#> ACP51            1.782204e-09 0.404004645 0.673 0.223  5.346611e-06       7
#> AP1S21           1.853642e-09 1.008900992 0.755 0.358  5.560927e-06       7
#> IGSF211          2.131997e-09 0.642208108 0.408 0.110  6.395991e-06       7
#> CDCA5            2.518081e-09 0.975837239 0.327 0.081  7.554242e-06       7
#> VASH11           2.650948e-09 0.798373740 0.367 0.095  7.952843e-06       7
#> LPAR62           3.993233e-09 0.691487265 0.449 0.135  1.197970e-05       7
#> FADS11           4.425552e-09 1.021877849 0.490 0.180  1.327665e-05       7
#> NDC801           4.902763e-09 0.982631344 0.347 0.095  1.470829e-05       7
#> TRPV21           5.307179e-09 0.612408837 0.408 0.111  1.592154e-05       7
#> MND11            5.947056e-09 1.078851856 0.347 0.098  1.784117e-05       7
#> RPS32            6.343515e-09 0.791691936 0.837 0.457  1.903054e-05       7
#> C1QB1            6.787399e-09 0.952352059 0.653 0.194  2.036220e-05       7
#> SGOL2            6.806733e-09 0.849935956 0.388 0.116  2.042020e-05       7
#> TUBA1B1          7.574358e-09 1.094140232 0.755 0.483  2.272307e-05       7
#> PILRA1           8.873800e-09 0.548553550 0.531 0.170  2.662140e-05       7
#> NUF21            9.483204e-09 0.544942662 0.408 0.119  2.844961e-05       7
#> UBE2C1           1.105200e-08 1.125544851 0.531 0.205  3.315599e-05       7
#> SMC21            1.416969e-08 1.321658422 0.571 0.297  4.250907e-05       7
#> CAMK13           1.593066e-08 0.958803865 0.449 0.160  4.779198e-05       7
#> FOXM11           1.702622e-08 0.872884965 0.388 0.118  5.107867e-05       7
#> SLC25A191        1.725983e-08 0.783381376 0.429 0.142  5.177949e-05       7
#> AURKA            2.074524e-08 0.920793044 0.408 0.137  6.223573e-05       7
#> CDCA2            2.541198e-08 1.023256692 0.327 0.092  7.623595e-05       7
#> FAM198B2         2.628255e-08 0.546215717 0.449 0.137  7.884765e-05       7
#> BUB1             2.803238e-08 0.776904541 0.286 0.069  8.409714e-05       7
#> CFD1             3.235970e-08 0.555452732 0.408 0.124  9.707911e-05       7
#> MXD3             3.652581e-08 0.840681814 0.327 0.091  1.095774e-04       7
#> HIST2H2AC        4.076665e-08 0.885310899 0.347 0.104  1.223000e-04       7
#> PRC11            4.683325e-08 0.890977984 0.469 0.182  1.404997e-04       7
#> ADGRE21          4.904433e-08 0.608726373 0.347 0.094  1.471330e-04       7
#> PSTPIP11         5.311524e-08 0.615819507 0.429 0.139  1.593457e-04       7
#> HLA-B4           5.439776e-08 0.729997463 0.898 0.482  1.631933e-04       7
#> CCL41            6.481431e-08 0.858013372 0.571 0.193  1.944429e-04       7
#> CCDC85B2         6.547479e-08 0.945482723 0.714 0.409  1.964244e-04       7
#> EPB41L31         6.730802e-08 0.673051212 0.490 0.176  2.019241e-04       7
#> COQ21            7.264782e-08 0.930871588 0.449 0.171  2.179435e-04       7
#> TAP12            8.184248e-08 0.809950690 0.612 0.278  2.455274e-04       7
#> C19orf481        9.010219e-08 0.843839121 0.612 0.295  2.703066e-04       7
#> KIF4A            9.557612e-08 0.689780453 0.286 0.073  2.867284e-04       7
#> TACC31           1.002782e-07 0.776247088 0.490 0.194  3.008345e-04       7
#> PTPN7            1.019860e-07 0.809445787 0.327 0.093  3.059579e-04       7
#> STAB11           1.116701e-07 0.632424680 0.490 0.174  3.350103e-04       7
#> PLXND11          1.143737e-07 0.485598492 0.633 0.243  3.431210e-04       7
#> CKAP21           1.215156e-07 0.771662436 0.531 0.233  3.645469e-04       7
#> ZWINT1           1.280589e-07 1.338098593 0.408 0.164  3.841768e-04       7
#> PLK11            1.523340e-07 0.506594221 0.388 0.126  4.570020e-04       7
#> NUDT11           1.693161e-07 0.902543695 0.694 0.370  5.079483e-04       7
#> CCNB11           1.922241e-07 0.697014077 0.449 0.172  5.766723e-04       7
#> CCL51            2.038747e-07 0.592976220 0.551 0.217  6.116240e-04       7
#> BCAT12           2.329562e-07 0.451703373 0.388 0.117  6.988686e-04       7
#> LAT22            2.447176e-07 0.766796649 0.653 0.306  7.341528e-04       7
#> SKA1             2.529488e-07 0.722728322 0.265 0.068  7.588464e-04       7
#> CENPF1           3.074876e-07 0.722654938 0.510 0.190  9.224627e-04       7
#> ANKRD281         3.206387e-07 0.769285785 0.469 0.194  9.619161e-04       7
#> SQRDL2           3.527744e-07 0.701943054 0.796 0.406  1.058323e-03       7
#> ATF51            4.017962e-07 0.713750261 0.551 0.232  1.205389e-03       7
#> GAS2L3           4.038576e-07 0.751577993 0.286 0.078  1.211573e-03       7
#> DNMT1            4.380802e-07 1.040340251 0.633 0.373  1.314241e-03       7
#> KIAA1524         5.111310e-07 0.885852381 0.286 0.083  1.533393e-03       7
#> SLCO2B11         5.633033e-07 0.370774239 0.449 0.146  1.689910e-03       7
#> PTTG11           5.853578e-07 0.953644668 0.510 0.219  1.756073e-03       7
#> GBP22            6.047235e-07 0.704601168 0.510 0.223  1.814170e-03       7
#> MAD2L11          7.148753e-07 0.758649791 0.449 0.181  2.144626e-03       7
#> LDHB2            7.571266e-07 0.726061686 0.837 0.485  2.271380e-03       7
#> SGOL11           8.502316e-07 0.571967537 0.347 0.110  2.550695e-03       7
#> KIF23            8.632842e-07 0.549083695 0.306 0.090  2.589853e-03       7
#> GIMAP71          8.952011e-07 0.611311246 0.306 0.090  2.685603e-03       7
#> UPP11            9.963086e-07 1.059001518 0.592 0.312  2.988926e-03       7
#> IQGAP21          1.020010e-06 0.602749780 0.367 0.121  3.060031e-03       7
#> CTSB2            1.198342e-06 0.570413446 0.714 0.316  3.595026e-03       7
#> MCM5             1.255010e-06 0.973739475 0.571 0.307  3.765031e-03       7
#> LINC009361       1.491208e-06 0.934750899 0.449 0.200  4.473624e-03       7
#> CD1631           1.698935e-06 0.487776351 0.449 0.158  5.096805e-03       7
#> TK11             1.717530e-06 0.889796563 0.469 0.206  5.152589e-03       7
#> UBE2S1           1.768900e-06 1.087010272 0.653 0.399  5.306700e-03       7
#> PGAM11           1.865623e-06 0.804305889 0.735 0.468  5.596870e-03       7
#> S100A42          1.894189e-06 0.632066866 0.857 0.459  5.682567e-03       7
#> FEN11            2.083393e-06 1.033017551 0.449 0.214  6.250180e-03       7
#> CD521            2.116506e-06 1.229024699 0.490 0.240  6.349518e-03       7
#> DTYMK1           2.249346e-06 1.050383924 0.612 0.362  6.748038e-03       7
#> SKA22            2.695674e-06 0.751468156 0.592 0.299  8.087022e-03       7
#> FCGR1A1          2.751376e-06 0.597338601 0.388 0.138  8.254128e-03       7
#> SOCS32           2.783522e-06 0.689710412 0.776 0.427  8.350567e-03       7
#> TROAP1           2.822057e-06 0.520378935 0.306 0.095  8.466170e-03       7
#> GBP12            3.089678e-06 0.656216959 0.306 0.097  9.269034e-03       7
#> CCNB21           3.759840e-06 0.582405292 0.449 0.184  1.127952e-02       7
#> LXN2             3.841876e-06 0.396872041 0.347 0.110  1.152563e-02       7
#> LINC010941       4.187261e-06 0.380103964 0.388 0.130  1.256178e-02       7
#> CDCA31           4.306783e-06 0.426561611 0.388 0.140  1.292035e-02       7
#> NASP1            4.333384e-06 0.891999753 0.673 0.393  1.300015e-02       7
#> ETV52            4.512770e-06 0.539952266 0.449 0.181  1.353831e-02       7
#> SAMD92           4.633379e-06 0.483019984 0.347 0.115  1.390014e-02       7
#> H2AFZ1           5.803412e-06 1.093557314 0.673 0.449  1.741024e-02       7
#> PRIM1            6.098941e-06 0.924815123 0.347 0.136  1.829682e-02       7
#> OAS3             6.125472e-06 0.577490241 0.265 0.078  1.837642e-02       7
#> TPX21            6.126735e-06 0.575002227 0.388 0.151  1.838021e-02       7
#> ID32             6.612087e-06 0.516108005 0.388 0.155  1.983626e-02       7
#> TRIM222          7.395982e-06 0.514214848 0.286 0.087  2.218795e-02       7
#> IFI442           7.600249e-06 0.664250496 0.327 0.115  2.280075e-02       7
#> CXCL21           8.684837e-06 0.491064727 0.388 0.144  2.605451e-02       7
#> IRF81            8.806099e-06 0.258801427 0.327 0.101  2.641830e-02       7
#> HMGB21           9.503107e-06 0.850407397 0.694 0.395  2.850932e-02       7
#> HMGN21           1.040996e-05 0.953041010 0.673 0.478  3.122988e-02       7
#> CDC201           1.221524e-05 0.716753338 0.388 0.161  3.664573e-02       7
#> HS3ST12          1.390215e-05 0.742513144 0.408 0.177  4.170645e-02       7
#> LMO43            1.533947e-05 0.668794821 0.633 0.328  4.601841e-02       7
#> XAF12            1.733804e-05 0.431647389 0.306 0.102  5.201413e-02       7
#> FANCI1           1.813832e-05 0.774589800 0.306 0.112  5.441496e-02       7
#> ITGAX1           1.880149e-05 0.488680875 0.306 0.103  5.640448e-02       7
#> RAD51AP11        1.997748e-05 0.746709840 0.408 0.183  5.993244e-02       7
#> LRR11            2.111198e-05 0.780592209 0.388 0.177  6.333593e-02       7
#> IER31            2.241104e-05 0.758177606 0.673 0.387  6.723311e-02       7
#> PARPBP1          2.306380e-05 0.622562569 0.306 0.111  6.919141e-02       7
#> KIFC11           2.347728e-05 0.638709212 0.327 0.124  7.043183e-02       7
#> THBD1            2.591411e-05 0.449491776 0.265 0.084  7.774233e-02       7
#> G0S21            2.916756e-05 0.561248236 0.408 0.172  8.750267e-02       7
#> APOC11           3.189479e-05 0.634420290 0.551 0.202  9.568437e-02       7
#> TNFSF102         3.476159e-05 0.491503910 0.327 0.122  1.042848e-01       7
#> ZNF3311          3.656824e-05 0.595677070 0.449 0.206  1.097047e-01       7
#> PHACTR12         4.100613e-05 0.630416253 0.531 0.279  1.230184e-01       7
#> FABP51           4.307844e-05 0.597130316 0.694 0.427  1.292353e-01       7
#> DAB22            4.633690e-05 0.364542399 0.551 0.259  1.390107e-01       7
#> SMC41            5.116313e-05 0.876856050 0.592 0.369  1.534894e-01       7
#> PPIF1            5.659865e-05 0.823096414 0.612 0.401  1.697959e-01       7
#> TCIRG11          5.720659e-05 0.487010187 0.633 0.364  1.716198e-01       7
#> TRAC             5.881071e-05 0.523033778 0.265 0.092  1.764321e-01       7
#> HLA-F4           6.801842e-05 0.582268031 0.592 0.344  2.040553e-01       7
#> WARS             7.004913e-05 0.532222226 0.469 0.241  2.101474e-01       7
#> TMPO1            7.603060e-05 1.013273644 0.551 0.362  2.280918e-01       7
#> FAM129A1         8.683465e-05 0.588610880 0.367 0.157  2.605039e-01       7
#> DUSP61           9.345776e-05 0.733164248 0.612 0.361  2.803733e-01       7
#> ICAM12           1.052067e-04 0.615892094 0.592 0.364  3.156202e-01       7
#> DNAJC153         1.057892e-04 0.472180756 0.653 0.347  3.173676e-01       7
#> DMXL21           1.071868e-04 0.408026236 0.367 0.149  3.215604e-01       7
#> DUT2             1.082958e-04 0.786990197 0.694 0.460  3.248873e-01       7
#> LINC001523       1.162152e-04 0.466916104 0.714 0.417  3.486457e-01       7
#> MPP11            1.200488e-04 0.323749059 0.429 0.184  3.601464e-01       7
#> LIMA12           1.228428e-04 0.487341227 0.469 0.228  3.685283e-01       7
#> CENPH1           1.257974e-04 0.668750508 0.367 0.172  3.773921e-01       7
#> GPX31            1.292471e-04 0.459887473 0.449 0.205  3.877413e-01       7
#> PARVB1           1.317551e-04 0.647558747 0.612 0.379  3.952654e-01       7
#> ERO1A1           1.330521e-04 0.494139628 0.449 0.212  3.991564e-01       7
#> SAC3D12          1.342342e-04 0.932083182 0.469 0.273  4.027025e-01       7
#> CDT11            1.366880e-04 0.641832163 0.327 0.141  4.100639e-01       7
#> C1orf542         1.448315e-04 0.566366281 0.490 0.233  4.344944e-01       7
#> BTG22            1.760083e-04 0.669405376 0.673 0.441  5.280250e-01       7
#> CENPU1           1.795103e-04 0.774398897 0.347 0.163  5.385310e-01       7
#> EZH2             2.001152e-04 0.776726873 0.429 0.250  6.003455e-01       7
#> HMOX11           2.119600e-04 0.278966872 0.571 0.302  6.358799e-01       7
#> SPP11            2.217626e-04 0.651658727 0.490 0.187  6.652877e-01       7
#> KCNN41           2.595751e-04 0.674705210 0.510 0.293  7.787252e-01       7
#> DAPP12           2.669152e-04 0.782614493 0.408 0.215  8.007455e-01       7
#> THOP11           3.101837e-04 0.353019518 0.408 0.210  9.305510e-01       7
#> IRF1             3.153546e-04 0.768619551 0.449 0.259  9.460638e-01       7
#> ORC61            3.488121e-04 0.528760895 0.327 0.149  1.000000e+00       7
#> NCAPD21          3.718547e-04 0.596725188 0.327 0.149  1.000000e+00       7
#> SORL12           3.795010e-04 0.695140676 0.367 0.186  1.000000e+00       7
#> BRCA21           4.168729e-04 0.650841283 0.327 0.156  1.000000e+00       7
#> MASTL            4.189487e-04 0.586055629 0.265 0.109  1.000000e+00       7
#> PARP142          4.443162e-04 0.458464010 0.429 0.210  1.000000e+00       7
#> CXCL31           4.532653e-04 0.500766375 0.327 0.137  1.000000e+00       7
#> RACGAP11         4.602254e-04 0.583413515 0.347 0.169  1.000000e+00       7
#> ARPC43           4.800506e-04 0.389478364 0.735 0.454  1.000000e+00       7
#> PKMYT11          5.026391e-04 0.426290212 0.265 0.104  1.000000e+00       7
#> MSR11            5.049890e-04 0.190828720 0.408 0.176  1.000000e+00       7
#> CLEC4E1          5.252845e-04 0.319927053 0.265 0.101  1.000000e+00       7
#> KNSTRN1          5.377588e-04 0.642842983 0.327 0.156  1.000000e+00       7
#> MELK1            5.398141e-04 0.606247599 0.286 0.124  1.000000e+00       7
#> MARCKS2          5.416471e-04 0.530135502 0.571 0.356  1.000000e+00       7
#> NEK21            5.980410e-04 0.163264573 0.265 0.100  1.000000e+00       7
#> NRP12            8.172964e-04 0.401293228 0.469 0.240  1.000000e+00       7
#> EPB41L22         9.832876e-04 0.495192288 0.469 0.260  1.000000e+00       7
#> OTOA1            9.991745e-04 0.159845984 0.286 0.114  1.000000e+00       7
#> LSM41            1.111359e-03 0.504144123 0.633 0.494  1.000000e+00       7
#> NETO21           1.129210e-03 0.622359330 0.327 0.169  1.000000e+00       7
#> TUBA1C3          1.151631e-03 0.683836028 0.612 0.467  1.000000e+00       7
#> H2AFX1           1.208996e-03 0.833010749 0.510 0.369  1.000000e+00       7
#> DCXR1            1.283530e-03 0.568389378 0.592 0.445  1.000000e+00       7
#> C5AR11           1.315142e-03 0.319872347 0.327 0.147  1.000000e+00       7
#> LGALS3BP2        1.450818e-03 0.171026216 0.449 0.212  1.000000e+00       7
#> HTRA12           1.463403e-03 0.601264873 0.429 0.223  1.000000e+00       7
#> MYBL21           1.734818e-03 0.391379882 0.327 0.160  1.000000e+00       7
#> C1QBP1           1.981964e-03 0.558985312 0.653 0.505  1.000000e+00       7
#> CD1092           2.027682e-03 0.302099379 0.265 0.113  1.000000e+00       7
#> GPSM22           2.126310e-03 0.504369298 0.347 0.185  1.000000e+00       7
#> GMNN1            2.160893e-03 0.679459586 0.429 0.269  1.000000e+00       7
#> RRM11            2.178270e-03 0.903901980 0.408 0.293  1.000000e+00       7
#> ACOT71           2.473787e-03 0.567551875 0.408 0.254  1.000000e+00       7
#> RPA31            2.682592e-03 0.641493323 0.612 0.482  1.000000e+00       7
#> ANP32E2          2.768581e-03 0.355670256 0.653 0.436  1.000000e+00       7
#> BID2             2.923996e-03 0.584357323 0.551 0.391  1.000000e+00       7
#> CAMK1D2          3.047971e-03 0.385859872 0.408 0.227  1.000000e+00       7
#> PPA13            3.537710e-03 0.563077886 0.653 0.494  1.000000e+00       7
#> CKAP51           3.663371e-03 0.600343675 0.347 0.211  1.000000e+00       7
#> RFC41            4.306279e-03 0.612093506 0.388 0.247  1.000000e+00       7
#> PAPSS21          4.408720e-03 0.398691386 0.347 0.186  1.000000e+00       7
#> MX13             4.875639e-03 0.420347992 0.286 0.144  1.000000e+00       7
#> CMSS12           4.957748e-03 0.470952214 0.490 0.328  1.000000e+00       7
#> CLDN12           4.979674e-03 0.468802101 0.327 0.179  1.000000e+00       7
#> TREM21           5.425556e-03 0.173591995 0.265 0.117  1.000000e+00       7
#> VAMP52           5.759391e-03 0.312915230 0.429 0.240  1.000000e+00       7
#> CBR3             6.040726e-03 0.509036548 0.306 0.174  1.000000e+00       7
#> CDKN1A1          6.568586e-03 0.259398733 0.490 0.311  1.000000e+00       7
#> BASP13           6.764971e-03 0.387246860 0.510 0.296  1.000000e+00       7
#> NKG71            6.890730e-03 0.367726749 0.306 0.165  1.000000e+00       7
#> AHCY             7.698775e-03 0.514417419 0.510 0.385  1.000000e+00       7
#> C32              9.206649e-03 0.458642569 0.388 0.245  1.000000e+00       7
#> SDS1             9.892722e-03 0.376534503 0.286 0.155  1.000000e+00       7
#> NEK22            1.803522e-61 4.430332166 1.000 0.084  5.410566e-58       8
#> PBK1             2.064247e-60 4.235198148 1.000 0.086  6.192741e-57       8
#> CENPA2           4.710395e-58 3.448688563 0.964 0.077  1.413118e-54       8
#> CDC25C           3.077076e-57 2.845775146 0.607 0.023  9.231229e-54       8
#> NMU              1.576128e-54 3.614908923 0.786 0.051  4.728384e-51       8
#> KIF20A           6.266232e-54 3.619764637 0.786 0.051  1.879870e-50       8
#> KIF231           1.437770e-53 3.369947825 0.929 0.078  4.313309e-50       8
#> HMMR1            2.231641e-53 3.233181459 0.964 0.085  6.694922e-50       8
#> LGR6             3.641513e-53 3.085744347 0.679 0.036  1.092454e-49       8
#> BUB11            1.303934e-52 3.132917061 0.821 0.059  3.911802e-49       8
#> CCNA22           6.603385e-52 3.382151542 0.964 0.091  1.981016e-48       8
#> PRSS3            1.281809e-51 2.747182887 0.786 0.051  3.845427e-48       8
#> DLGAP51          3.277439e-51 2.969607680 0.821 0.060  9.832316e-48       8
#> CENPE1           1.016810e-50 2.899038608 0.857 0.065  3.050430e-47       8
#> NCAPG1           1.950198e-50 3.127182249 0.857 0.068  5.850595e-47       8
#> PLK12            4.736306e-50 4.507881747 1.000 0.115  1.420892e-46       8
#> NUF22            4.916198e-50 3.703560066 1.000 0.109  1.474859e-46       8
#> FAM72D           3.816459e-48 2.036294457 0.536 0.022  1.144938e-44       8
#> FAM83D           1.809405e-47 3.813681343 0.857 0.078  5.428216e-44       8
#> FAM72C           1.900262e-45 2.654638301 0.536 0.024  5.700786e-42       8
#> KIF151           5.287403e-45 2.607999418 0.750 0.056  1.586221e-41       8
#> FAM64A1          5.140728e-44 2.826585686 0.821 0.074  1.542218e-40       8
#> KIF2C1           5.067235e-42 3.087486197 0.750 0.064  1.520171e-38       8
#> DEPDC11          5.198403e-42 2.325326193 0.750 0.059  1.559521e-38       8
#> SGOL12           7.101720e-42 3.392791698 0.893 0.100  2.130516e-38       8
#> GTSE11           8.819126e-42 2.205328602 0.857 0.081  2.645738e-38       8
#> KIF4A1           2.443348e-41 3.066245497 0.750 0.065  7.330045e-38       8
#> CDCA32           4.135499e-41 3.494643968 0.964 0.130  1.240650e-37       8
#> SHCBP11          1.434848e-40 2.495811654 0.821 0.077  4.304545e-37       8
#> SPC252           4.287665e-40 2.554571338 0.857 0.089  1.286299e-36       8
#> CDKN32           4.955203e-40 3.419642283 0.929 0.121  1.486561e-36       8
#> BIRC52           1.895528e-39 3.449883753 1.000 0.151  5.686585e-36       8
#> TPX22            3.340557e-39 3.572372762 0.964 0.140  1.002167e-35       8
#> RACGAP12         9.051207e-39 3.424851202 1.000 0.155  2.715362e-35       8
#> AURKB2           3.346313e-38 2.354231290 0.893 0.100  1.003894e-34       8
#> ANLN2            1.062881e-37 2.470009297 0.857 0.096  3.188644e-34       8
#> TROAP2           1.325424e-37 2.485193289 0.821 0.086  3.976273e-34       8
#> KIF14            3.759502e-37 2.304805094 0.714 0.063  1.127851e-33       8
#> CDCA81           2.389675e-36 2.287461852 0.750 0.072  7.169024e-33       8
#> CDC202           6.142882e-36 3.241593735 0.964 0.151  1.842865e-32       8
#> NUSAP12          1.073409e-35 2.802451593 0.929 0.132  3.220227e-32       8
#> CCNB12           1.187648e-35 3.900791432 0.964 0.164  3.562944e-32       8
#> KIFC12           1.254778e-35 2.493021640 0.893 0.113  3.764335e-32       8
#> AURKA1           1.881015e-35 3.598159474 0.893 0.130  5.643046e-32       8
#> TOP2A2           5.541599e-35 2.545345657 0.929 0.126  1.662480e-31       8
#> KIAA15241        6.867450e-35 2.404399656 0.750 0.075  2.060235e-31       8
#> MKI672           1.579632e-34 2.800900294 0.893 0.125  4.738895e-31       8
#> CEP552           6.132816e-34 2.625163490 0.786 0.091  1.839845e-30       8
#> MAD2L12          1.139878e-33 3.414996026 0.964 0.173  3.419635e-30       8
#> ECT21            4.470679e-32 2.082660986 0.821 0.104  1.341204e-28       8
#> CENPF2           4.623303e-32 3.565250419 1.000 0.183  1.386991e-28       8
#> CCNB22           5.448151e-32 3.273434720 0.964 0.176  1.634445e-28       8
#> SPAG5            1.506086e-31 1.955484709 0.607 0.052  4.518258e-28       8
#> ASPM2            3.133447e-31 1.926191408 0.714 0.075  9.400340e-28       8
#> KIF20B2          5.753860e-31 3.209996456 0.929 0.170  1.726158e-27       8
#> PKMYT12          9.225366e-31 1.944974874 0.786 0.094  2.767610e-27       8
#> PARPBP2          9.396401e-31 2.566824576 0.786 0.102  2.818920e-27       8
#> FOXM12           2.757751e-30 2.321306141 0.821 0.112  8.273253e-27       8
#> ASF1B1           4.015531e-30 1.665662169 0.786 0.091  1.204659e-26       8
#> PRC12            4.945195e-30 3.026342115 0.929 0.176  1.483559e-26       8
#> TACC32           1.147902e-29 3.467778966 0.929 0.188  3.443707e-26       8
#> SGOL21           2.401564e-29 2.655223697 0.786 0.111  7.204693e-26       8
#> FAM72B           3.744942e-29 1.652812891 0.321 0.013  1.123483e-25       8
#> TK12             6.900653e-29 2.568843168 1.000 0.197  2.070196e-25       8
#> KIF111           7.367877e-28 1.672451035 0.607 0.060  2.210363e-24       8
#> CKAP22           2.496115e-27 3.453144612 0.964 0.227  7.488344e-24       8
#> UBE2C2           6.203923e-27 3.033580266 0.964 0.200  1.861177e-23       8
#> SKA11            7.539861e-27 1.593660672 0.607 0.063  2.261958e-23       8
#> PTTG12           9.523225e-27 3.017704868 0.964 0.212  2.856968e-23       8
#> CDCA21           1.861318e-26 2.365291962 0.679 0.087  5.583955e-23       8
#> LMNB12           3.945769e-26 2.613898013 0.857 0.164  1.183731e-22       8
#> NDC802           6.287462e-26 1.812526041 0.714 0.091  1.886239e-22       8
#> TTK              4.176608e-25 1.958681710 0.607 0.069  1.252982e-21       8
#> SKA23            4.298513e-24 2.914424611 1.000 0.294  1.289554e-20       8
#> MELK2            5.714449e-24 1.902513775 0.750 0.115  1.714335e-20       8
#> DBF4             8.187302e-24 3.109278021 0.964 0.266  2.456191e-20       8
#> KIF18A           8.600896e-24 1.525324565 0.393 0.028  2.580269e-20       8
#> MYBL22           4.813579e-23 1.917244348 0.821 0.151  1.444074e-19       8
#> GPSM23           6.998735e-23 2.313026186 0.857 0.175  2.099621e-19       8
#> POC1A1           8.347274e-23 2.262606681 0.750 0.128  2.504182e-19       8
#> TRIP131          1.331363e-22 2.180311348 0.750 0.131  3.994088e-19       8
#> SKA3             3.307556e-22 1.659494580 0.500 0.051  9.922669e-19       8
#> ZMYND10          2.010233e-21 1.859399737 0.286 0.016  6.030698e-18       8
#> IQGAP3           2.246053e-21 1.262933838 0.464 0.045  6.738158e-18       8
#> LBR2             5.432993e-21 2.680227662 1.000 0.345  1.629898e-17       8
#> ATAD22           7.927352e-21 2.088291674 0.964 0.267  2.378206e-17       8
#> TRAF3IP32        8.399913e-21 2.227532615 1.000 0.287  2.519974e-17       8
#> DEPDC1B          2.089899e-20 1.403251020 0.464 0.048  6.269696e-17       8
#> UBE2S2           2.270068e-20 3.034309897 1.000 0.395  6.810204e-17       8
#> ACTL81           9.062259e-20 2.419644224 0.679 0.126  2.718678e-16       8
#> CDKN2C1          1.030860e-19 1.723770431 0.714 0.125  3.092581e-16       8
#> UBE2T1           1.223378e-19 2.458390506 0.929 0.282  3.670134e-16       8
#> DIAPH31          1.410941e-19 1.628510674 0.750 0.141  4.232823e-16       8
#> CKS23            1.697196e-19 3.236675776 1.000 0.435  5.091587e-16       8
#> SMC42            1.860831e-19 2.432019674 1.000 0.363  5.582493e-16       8
#> NCAPD22          1.915171e-19 1.617676267 0.750 0.141  5.745513e-16       8
#> HMGB22           2.720921e-19 2.895305250 1.000 0.393  8.162763e-16       8
#> PSRC11           2.891880e-19 1.752731872 0.786 0.157  8.675639e-16       8
#> TEX301           2.980324e-19 2.338047534 0.964 0.324  8.940973e-16       8
#> CDKN2D1          3.076773e-19 1.985509136 0.893 0.244  9.230320e-16       8
#> RIBC2            3.540562e-19 1.214515585 0.357 0.029  1.062168e-15       8
#> REEP41           3.694799e-19 2.521341406 0.893 0.286  1.108440e-15       8
#> HN12             3.824979e-19 2.516794745 1.000 0.493  1.147494e-15       8
#> CDCA51           4.736131e-19 1.545518600 0.571 0.080  1.420839e-15       8
#> RAD214           6.379547e-19 2.294202501 1.000 0.529  1.913864e-15       8
#> MXD31            8.036584e-19 1.485658080 0.607 0.088  2.410975e-15       8
#> ARL6IP13         8.561475e-19 2.599305864 1.000 0.474  2.568443e-15       8
#> ANP32E3          1.117912e-18 2.572911578 1.000 0.431  3.353737e-15       8
#> CENPW1           1.269794e-18 2.047866999 0.929 0.284  3.809382e-15       8
#> KPNA2            1.429887e-18 2.826143474 1.000 0.453  4.289662e-15       8
#> PRPH             1.605055e-18 1.134604339 0.286 0.019  4.815166e-15       8
#> KIF222           1.612286e-18 2.533934195 1.000 0.427  4.836859e-15       8
#> MCM72            1.651495e-18 2.323266113 1.000 0.376  4.954484e-15       8
#> AC004381.6       1.891856e-18 1.562553799 0.464 0.054  5.675569e-15       8
#> WFDC11           2.535609e-18 1.502330115 0.607 0.094  7.606826e-15       8
#> GGH2             3.070649e-18 2.383081638 1.000 0.383  9.211946e-15       8
#> HMGB12           3.620566e-18 2.154658300 1.000 0.487  1.086170e-14       8
#> ZWINT2           3.804495e-18 1.346260844 0.786 0.159  1.141349e-14       8
#> CKAP2L1          4.453100e-18 1.549121978 0.393 0.039  1.335930e-14       8
#> GAS2L31          5.423733e-18 1.671876756 0.536 0.076  1.627120e-14       8
#> CA83             6.512660e-18 2.293910766 0.857 0.262  1.953798e-14       8
#> CKS1B3           1.156581e-17 1.864874279 1.000 0.503  3.469742e-14       8
#> HSPD11           1.601911e-17 2.196429244 1.000 0.484  4.805732e-14       8
#> SAPCD21          1.682900e-17 1.986424543 0.571 0.095  5.048700e-14       8
#> CENPM2           1.888662e-17 1.386060036 0.679 0.123  5.665986e-14       8
#> BRCA22           2.026462e-17 1.477654764 0.750 0.149  6.079386e-14       8
#> SSRP13           3.207905e-17 1.850453701 1.000 0.487  9.623714e-14       8
#> CENPU2           3.680847e-17 1.447681985 0.750 0.156  1.104254e-13       8
#> CDT12            3.964921e-17 1.353411585 0.714 0.135  1.189476e-13       8
#> KIAA01012        4.139124e-17 1.355898578 0.893 0.215  1.241737e-13       8
#> RAC32            4.473490e-17 2.075580038 0.786 0.207  1.342047e-13       8
#> STMN13           5.086392e-17 1.722287494 1.000 0.522  1.525918e-13       8
#> CDK12            6.959271e-17 1.657933623 0.679 0.140  2.087781e-13       8
#> EME1             1.116102e-16 1.896986658 0.321 0.029  3.348305e-13       8
#> YBX22            1.193058e-16 1.680864762 0.679 0.135  3.579175e-13       8
#> MT1E1            1.232629e-16 2.157093252 0.964 0.383  3.697886e-13       8
#> CCNF             1.246926e-16 1.712179475 0.607 0.108  3.740777e-13       8
#> TMX23            1.555526e-16 1.928905770 1.000 0.463  4.666579e-13       8
#> DTYMK2           1.720249e-16 2.034385417 0.964 0.357  5.160748e-13       8
#> HSP90AB15        2.261826e-16 1.224779882 1.000 0.664  6.785477e-13       8
#> TUBA1B2          4.095243e-16 1.835393298 1.000 0.482  1.228573e-12       8
#> RAD51AP12        6.703256e-16 1.271009324 0.786 0.178  2.010977e-12       8
#> CHEK11           9.194587e-16 1.565236592 0.786 0.194  2.758376e-12       8
#> CENPN1           1.242240e-15 1.580037349 0.893 0.285  3.726719e-12       8
#> CKAP52           1.303732e-15 1.770147942 0.786 0.202  3.911195e-12       8
#> ACOT72           1.774176e-15 1.881998436 0.821 0.246  5.322528e-12       8
#> FANCI2           2.230836e-15 1.269666242 0.607 0.108  6.692509e-12       8
#> CCT53            2.268468e-15 1.831255536 0.964 0.521  6.805405e-12       8
#> CACYBP4          2.446879e-15 1.612422150 1.000 0.526  7.340637e-12       8
#> RTKN2            2.676069e-15 1.758882798 0.429 0.059  8.028206e-12       8
#> ODC11            2.775440e-15 1.983159995 0.929 0.384  8.326321e-12       8
#> LSM52            3.989762e-15 1.737227886 1.000 0.467  1.196929e-11       8
#> SNRNP251         4.652624e-15 1.775092634 1.000 0.449  1.395787e-11       8
#> TYMS2            5.615050e-15 0.830934023 0.786 0.167  1.684515e-11       8
#> TMPO2            6.060792e-15 1.615080804 0.964 0.355  1.818238e-11       8
#> NME12            6.061101e-15 1.763647918 1.000 0.486  1.818330e-11       8
#> DKC12            8.702092e-15 2.067443142 0.929 0.475  2.610628e-11       8
#> HYLS11           1.152988e-14 1.874698556 0.679 0.164  3.458964e-11       8
#> KRT85            1.326733e-14 1.522363320 1.000 0.544  3.980198e-11       8
#> SMC22            1.661923e-14 1.418964642 0.929 0.293  4.985768e-11       8
#> KNSTRN2          1.897023e-14 1.406010760 0.679 0.151  5.691069e-11       8
#> MT2A4            2.285150e-14 1.814410546 1.000 0.397  6.855450e-11       8
#> H2AFV2           2.420404e-14 1.608387994 0.964 0.485  7.261213e-11       8
#> ACTL6A3          2.618500e-14 1.783248842 1.000 0.480  7.855500e-11       8
#> NUCKS13          2.765234e-14 1.498865621 1.000 0.495  8.295702e-11       8
#> RUVBL21          2.882282e-14 1.919946871 0.893 0.397  8.646847e-11       8
#> HSPH12           3.797229e-14 1.852210411 0.964 0.488  1.139169e-10       8
#> CKB4             3.966510e-14 1.617031968 1.000 0.442  1.189953e-10       8
#> RUVBL12          5.781531e-14 1.849559652 0.929 0.432  1.734459e-10       8
#> PIF1             5.950342e-14 1.084212777 0.393 0.051  1.785102e-10       8
#> TIMM103          6.012914e-14 1.543378563 1.000 0.480  1.803874e-10       8
#> H2AFX2           6.249385e-14 1.578815528 0.929 0.361  1.874816e-10       8
#> CEP703           8.527292e-14 1.707936716 0.786 0.241  2.558188e-10       8
#> DNMT11           8.576610e-14 1.445765945 1.000 0.369  2.572983e-10       8
#> RAD51            8.861021e-14 0.837730063 0.500 0.077  2.658306e-10       8
#> DSC33            1.087473e-13 1.558743595 0.893 0.314  3.262419e-10       8
#> MDC11            1.189266e-13 1.481789741 0.571 0.114  3.567797e-10       8
#> GMNN2            1.193623e-13 1.364553260 0.857 0.261  3.580868e-10       8
#> PCBD12           1.226400e-13 1.589258854 0.964 0.500  3.679201e-10       8
#> LYAR2            1.283711e-13 2.209599419 0.857 0.354  3.851133e-10       8
#> AZGP13           1.384792e-13 1.670784482 1.000 0.364  4.154376e-10       8
#> H2AFZ2           1.600007e-13 1.549612706 1.000 0.445  4.800020e-10       8
#> SLBP3            2.616515e-13 1.783477125 0.929 0.492  7.849545e-10       8
#> CASC51           2.984136e-13 0.808085287 0.357 0.043  8.952408e-10       8
#> OBP2B1           3.930708e-13 1.651647274 0.571 0.123  1.179212e-09       8
#> FKBP43           5.704696e-13 1.592132475 0.929 0.492  1.711409e-09       8
#> RRS12            1.400851e-12 1.638408260 0.893 0.362  4.202553e-09       8
#> HMGN22           1.408672e-12 1.480846974 1.000 0.473  4.226015e-09       8
#> TOMM402          1.778490e-12 1.604947200 0.929 0.465  5.335469e-09       8
#> DCXR2            2.756254e-12 1.574553496 0.929 0.439  8.268762e-09       8
#> CLSPN1           2.805378e-12 0.938599408 0.464 0.078  8.416134e-09       8
#> CA61             4.081175e-12 1.233436156 0.536 0.109  1.224352e-08       8
#> BSPRY2           7.598390e-12 1.569840924 0.821 0.276  2.279517e-08       8
#> HIST1H2BN2       8.965383e-12 1.324386441 0.679 0.188  2.689615e-08       8
#> EFHD13           9.896041e-12 1.538232457 0.929 0.419  2.968812e-08       8
#> HSPA1A3          1.088611e-11 1.326878319 0.929 0.496  3.265832e-08       8
#> ZNF6951          1.199001e-11 1.254239563 0.500 0.101  3.597002e-08       8
#> NCAPH1           1.414905e-11 1.180795921 0.321 0.042  4.244715e-08       8
#> MTHFD23          1.432351e-11 1.395988559 0.929 0.408  4.297052e-08       8
#> SPDL11           1.587576e-11 1.259331254 0.571 0.133  4.762727e-08       8
#> MT1G2            1.781159e-11 0.829908908 0.786 0.223  5.343478e-08       8
#> NTHL12           1.911022e-11 1.411918981 0.893 0.411  5.733067e-08       8
#> RFC42            2.656903e-11 0.999681485 0.821 0.239  7.970710e-08       8
#> TMEM106C3        2.696336e-11 1.377900840 1.000 0.411  8.089007e-08       8
#> PCNA1            2.899125e-11 1.171016784 0.964 0.400  8.697375e-08       8
#> LRRCC11          2.959298e-11 1.594317055 0.857 0.407  8.877894e-08       8
#> HJURP1           2.964885e-11 1.191366371 0.393 0.064  8.894656e-08       8
#> TUBA1C4          3.405148e-11 1.600380147 0.929 0.461  1.021545e-07       8
#> PPA14            3.643914e-11 1.209527056 1.000 0.488  1.093174e-07       8
#> PHGDH2           3.870324e-11 1.361259100 0.929 0.397  1.161097e-07       8
#> SYCP22           4.202752e-11 1.369567077 0.750 0.239  1.260825e-07       8
#> NES2             4.307003e-11 1.284637012 0.857 0.335  1.292101e-07       8
#> SLC2A4RG3        5.897636e-11 1.244735720 0.964 0.430  1.769291e-07       8
#> DNAJC91          6.552496e-11 1.340740475 0.893 0.425  1.965749e-07       8
#> MCM21            7.711925e-11 1.060124497 0.571 0.133  2.313578e-07       8
#> RPL39L3          1.242703e-10 1.315080184 0.964 0.459  3.728110e-07       8
#> PBX13            1.347167e-10 1.344134619 0.929 0.460  4.041500e-07       8
#> ALYREF2          1.446851e-10 1.293464500 0.857 0.344  4.340553e-07       8
#> TTF22            1.472145e-10 1.544501524 0.821 0.326  4.416434e-07       8
#> ASNS2            1.509833e-10 1.140285752 0.750 0.235  4.529500e-07       8
#> EZH21            1.902533e-10 1.212398567 0.750 0.245  5.707600e-07       8
#> CENPH2           1.972550e-10 1.011723212 0.643 0.168  5.917651e-07       8
#> CDC45            2.071566e-10 0.827265307 0.393 0.065  6.214697e-07       8
#> ORC62            2.360128e-10 0.777023897 0.607 0.145  7.080385e-07       8
#> IGFBP22          2.468032e-10 1.253699061 0.964 0.400  7.404096e-07       8
#> PDGFRA3          4.225113e-10 1.065596182 0.821 0.330  1.267534e-06       8
#> PRNP2            4.263977e-10 1.280506915 0.929 0.518  1.279193e-06       8
#> SMC1B1           4.582433e-10 1.215006790 0.500 0.116  1.374730e-06       8
#> ZNF367           4.750775e-10 0.921923053 0.286 0.038  1.425232e-06       8
#> SLC9A3R23        5.122910e-10 1.219772503 0.964 0.440  1.536873e-06       8
#> TUBB4B4          5.133284e-10 1.362712205 0.929 0.506  1.539985e-06       8
#> CHAF1A1          5.862057e-10 0.969560158 0.607 0.155  1.758617e-06       8
#> CDCA7            6.100555e-10 1.242795069 0.393 0.072  1.830166e-06       8
#> CCT6A2           6.821719e-10 1.181692426 0.929 0.510  2.046516e-06       8
#> DSN11            7.338979e-10 1.068921335 0.750 0.241  2.201694e-06       8
#> VANGL14          8.730755e-10 1.208684480 1.000 0.525  2.619227e-06       8
#> C9orf401         9.897646e-10 1.088659630 0.750 0.249  2.969294e-06       8
#> UHRF1            1.183774e-09 0.934748187 0.464 0.098  3.551323e-06       8
#> BOP12            1.236755e-09 1.446768521 0.857 0.448  3.710266e-06       8
#> MSX1             1.245801e-09 1.209262089 0.429 0.090  3.737404e-06       8
#> GAL2             1.311695e-09 1.001277407 0.750 0.271  3.935085e-06       8
#> RAD54B           1.384997e-09 1.060126727 0.357 0.063  4.154991e-06       8
#> PLOD32           1.457330e-09 1.272911864 0.929 0.418  4.371990e-06       8
#> PAICS2           1.517309e-09 1.405624735 0.893 0.419  4.551927e-06       8
#> GINS21           1.690925e-09 0.683312344 0.643 0.167  5.072775e-06       8
#> CHORDC11         1.777893e-09 1.381986364 0.893 0.479  5.333680e-06       8
#> RASL11B          1.972584e-09 0.656835010 0.536 0.123  5.917751e-06       8
#> HIST1H4C3        2.011531e-09 0.877165158 0.893 0.330  6.034593e-06       8
#> TUBB2A3          2.142326e-09 1.160834766 0.929 0.420  6.426978e-06       8
#> VGF2             2.533373e-09 0.986572415 0.607 0.166  7.600119e-06       8
#> C12orf752        2.543085e-09 1.137567803 0.750 0.271  7.629255e-06       8
#> CTHRC13          2.554947e-09 1.174256895 0.929 0.465  7.664842e-06       8
#> MIS18A1          2.654580e-09 1.255299040 0.750 0.294  7.963740e-06       8
#> ARHGAP11A1       3.074010e-09 0.805552549 0.357 0.063  9.222029e-06       8
#> PRR111           3.184991e-09 1.128700649 0.393 0.078  9.554972e-06       8
#> FBLN12           3.341908e-09 1.055335787 0.750 0.261  1.002573e-05       8
#> LAPTM4B4         3.635727e-09 1.100599167 0.929 0.545  1.090718e-05       8
#> CD3202           4.084952e-09 1.524054168 0.714 0.296  1.225486e-05       8
#> TMEM972          6.794633e-09 1.128908224 0.714 0.253  2.038390e-05       8
#> HSPA1B3          8.302120e-09 1.127976257 0.929 0.458  2.490636e-05       8
#> RFC31            8.682388e-09 0.944382108 0.786 0.277  2.604716e-05       8
#> HIST1H2BJ2       9.013035e-09 0.957280588 0.643 0.202  2.703910e-05       8
#> RMI22            1.014682e-08 0.869487419 0.714 0.241  3.044045e-05       8
#> ST143            1.162041e-08 1.081851227 0.893 0.486  3.486122e-05       8
#> FRMD4A2          1.164688e-08 1.131522498 0.821 0.360  3.494065e-05       8
#> GTF3A1           1.171702e-08 1.185916270 0.893 0.496  3.515105e-05       8
#> NT5DC21          1.172919e-08 1.145346090 0.893 0.463  3.518756e-05       8
#> PPIL12           1.190488e-08 1.156622463 0.857 0.369  3.571465e-05       8
#> RNASEH2A1        1.296565e-08 0.965157617 0.750 0.274  3.889695e-05       8
#> PPP1R12A2        1.352263e-08 1.165723269 0.929 0.450  4.056788e-05       8
#> UCHL12           1.444606e-08 0.802932125 0.571 0.161  4.333818e-05       8
#> CDC42EP13        1.588895e-08 1.188352133 0.857 0.449  4.766685e-05       8
#> TNFRSF181        1.669038e-08 0.920567081 0.500 0.125  5.007113e-05       8
#> SH3BGR3          1.967817e-08 1.130098027 0.750 0.338  5.903451e-05       8
#> DCTPP12          2.029965e-08 1.214256910 0.857 0.496  6.089894e-05       8
#> SLC43A33         2.061934e-08 1.086898879 0.964 0.447  6.185801e-05       8
#> PYCR13           2.144649e-08 1.206222037 0.893 0.452  6.433948e-05       8
#> SOX84            2.712456e-08 1.362156306 0.821 0.395  8.137369e-05       8
#> CCDC341          3.545055e-08 1.084175480 0.607 0.203  1.063516e-04       8
#> SIGMAR12         3.792314e-08 1.166628080 0.857 0.380  1.137694e-04       8
#> CCDC85B3         4.183500e-08 1.240060329 0.821 0.413  1.255050e-04       8
#> EXOSC8           4.671878e-08 1.175808284 0.821 0.433  1.401563e-04       8
#> BYSL2            4.871858e-08 1.140316152 0.786 0.425  1.461557e-04       8
#> RCCD11           4.877046e-08 0.992314731 0.679 0.231  1.463114e-04       8
#> SERTAD43         5.431754e-08 1.004975010 0.964 0.466  1.629526e-04       8
#> MCM10            6.515399e-08 0.701899996 0.357 0.072  1.954620e-04       8
#> MGST13           6.608563e-08 1.070308027 0.964 0.456  1.982569e-04       8
#> MND12            7.117909e-08 0.779320685 0.429 0.101  2.135373e-04       8
#> SQLE3            7.307646e-08 0.884705787 0.929 0.435  2.192294e-04       8
#> KRT185           8.507044e-08 1.033932851 0.929 0.544  2.552113e-04       8
#> NAV2             9.996339e-08 1.148585718 0.357 0.079  2.998902e-04       8
#> TIMELESS         1.153064e-07 0.803514444 0.500 0.135  3.459193e-04       8
#> ADAMTS41         1.215263e-07 0.477793789 0.286 0.050  3.645789e-04       8
#> PRKDC2           1.462233e-07 1.063617698 0.821 0.506  4.386700e-04       8
#> DUSP9            1.511922e-07 0.768571133 0.250 0.040  4.535765e-04       8
#> CENPQ2           1.542959e-07 0.855976228 0.857 0.371  4.628877e-04       8
#> C6orf1412        1.946800e-07 0.912656614 0.500 0.145  5.840399e-04       8
#> SUSD51           2.050492e-07 0.855178662 0.464 0.121  6.151477e-04       8
#> LMNB21           2.160072e-07 0.831984336 0.571 0.171  6.480216e-04       8
#> PARVB2           2.570586e-07 1.000351615 0.786 0.379  7.711757e-04       8
#> DSCC11           2.667335e-07 0.747519576 0.429 0.107  8.002005e-04       8
#> IL17B3           3.125555e-07 0.724854379 0.571 0.181  9.376666e-04       8
#> CTNNAL11         3.389944e-07 1.098156029 0.643 0.264  1.016983e-03       8
#> CMSS13           3.816266e-07 0.858805317 0.786 0.324  1.144880e-03       8
#> HIBCH3           4.230080e-07 0.973447201 0.929 0.417  1.269024e-03       8
#> GAMT1            4.352272e-07 1.018190693 0.643 0.239  1.305682e-03       8
#> MYC1             4.427584e-07 0.893361795 0.750 0.312  1.328275e-03       8
#> GOLM13           4.464887e-07 0.989111613 0.821 0.480  1.339466e-03       8
#> HACD13           5.794534e-07 1.049440796 0.821 0.405  1.738360e-03       8
#> CENPK1           5.863425e-07 0.688237019 0.393 0.096  1.759028e-03       8
#> PITX13           6.030778e-07 1.016970698 0.821 0.412  1.809233e-03       8
#> PSIP12           6.070374e-07 0.826023603 0.821 0.341  1.821112e-03       8
#> FRMD31           6.461785e-07 0.847482761 0.607 0.216  1.938536e-03       8
#> APOLD1           6.593680e-07 1.005855862 0.536 0.180  1.978104e-03       8
#> BARD13           7.040025e-07 0.856784804 0.714 0.321  2.112007e-03       8
#> SDC23            7.151394e-07 0.989084036 0.893 0.467  2.145418e-03       8
#> APOBEC3B2        7.650246e-07 0.930073946 0.500 0.154  2.295074e-03       8
#> MCM31            8.258977e-07 0.832812941 0.821 0.385  2.477693e-03       8
#> GINS1            8.393624e-07 0.595826436 0.357 0.081  2.518087e-03       8
#> YDJC1            9.075771e-07 0.951501494 0.821 0.460  2.722731e-03       8
#> PSAT11           9.393658e-07 0.687952351 0.393 0.101  2.818097e-03       8
#> MCM41            1.051367e-06 0.789060231 0.750 0.348  3.154101e-03       8
#> LSM42            1.153974e-06 0.950940549 0.893 0.490  3.461923e-03       8
#> CDC6             1.283628e-06 0.823467804 0.286 0.059  3.850884e-03       8
#> COL2A12          1.449203e-06 0.670770636 0.714 0.312  4.347610e-03       8
#> NUDT12           1.554535e-06 0.938264679 0.786 0.374  4.663606e-03       8
#> LY6E4            1.576824e-06 0.914059379 0.893 0.467  4.730473e-03       8
#> GINS41           1.703729e-06 0.474352682 0.571 0.177  5.111186e-03       8
#> RANBP12          1.866332e-06 0.767357924 0.929 0.499  5.598997e-03       8
#> DKK12            2.241414e-06 0.559940741 0.464 0.152  6.724241e-03       8
#> DEK1             3.071460e-06 0.799379930 0.857 0.507  9.214380e-03       8
#> SAP302           3.172089e-06 1.030099340 0.750 0.436  9.516266e-03       8
#> THOP12           3.185716e-06 0.344445498 0.607 0.209  9.557149e-03       8
#> ADAM154          3.198678e-06 0.877230536 0.857 0.502  9.596035e-03       8
#> RBBP72           3.201799e-06 0.867993576 0.857 0.514  9.605398e-03       8
#> CTSV3            3.210396e-06 0.653566075 0.679 0.325  9.631189e-03       8
#> ELOVL52          3.524022e-06 0.867073919 0.643 0.269  1.057207e-02       8
#> LRR12            4.456066e-06 0.683907772 0.536 0.178  1.336820e-02       8
#> ACTA24           5.293022e-06 0.856037600 0.821 0.350  1.587907e-02       8
#> C21orf581        5.369788e-06 0.617829142 0.393 0.110  1.610936e-02       8
#> SERPINE22        5.424636e-06 0.735899018 0.857 0.393  1.627391e-02       8
#> TCF191           6.205823e-06 0.543037257 0.571 0.196  1.861747e-02       8
#> TM4SF13          6.263501e-06 0.790152650 0.929 0.588  1.879050e-02       8
#> AHCY1            6.595451e-06 1.003140035 0.786 0.381  1.978635e-02       8
#> DCUN1D51         7.896449e-06 0.824381807 0.750 0.471  2.368935e-02       8
#> NEMP2            8.310511e-06 0.580581992 0.357 0.094  2.493153e-02       8
#> RECQL41          9.051925e-06 0.640874678 0.429 0.130  2.715578e-02       8
#> MYBL13           9.648024e-06 0.749103007 0.643 0.258  2.894407e-02       8
#> TUBB2B3          9.982815e-06 0.669580519 0.750 0.334  2.994844e-02       8
#> DACT3            1.015188e-05 0.613003529 0.286 0.065  3.045565e-02       8
#> ERVMER34-12      1.072226e-05 0.859873222 0.571 0.228  3.216677e-02       8
#> MT1X2            1.158507e-05 0.675875333 0.786 0.390  3.475522e-02       8
#> ITGA62           1.165398e-05 0.758800123 0.679 0.290  3.496195e-02       8
#> TBC1D13          1.364400e-05 0.830606043 0.857 0.452  4.093200e-02       8
#> NQO13            1.739576e-05 0.837824178 0.786 0.435  5.218727e-02       8
#> FEN12            1.810121e-05 0.605642418 0.571 0.215  5.430362e-02       8
#> CDCA41           1.997810e-05 0.646968053 0.607 0.246  5.993431e-02       8
#> NTF3             2.195058e-05 0.500202866 0.321 0.081  6.585173e-02       8
#> CPNE7            2.253050e-05 0.465395813 0.250 0.054  6.759151e-02       8
#> FABP71           2.254675e-05 0.823722384 0.250 0.056  6.764026e-02       8
#> POLD22           2.579159e-05 0.824692770 0.786 0.443  7.737476e-02       8
#> RRM12            2.689529e-05 0.696697571 0.679 0.288  8.068588e-02       8
#> SEPT42           2.852932e-05 0.605076425 0.393 0.123  8.558795e-02       8
#> EXO1             2.867428e-05 0.461689138 0.250 0.055  8.602284e-02       8
#> HAPLN12          3.250909e-05 0.459452719 0.464 0.168  9.752728e-02       8
#> LEFTY22          3.963709e-05 0.806218440 0.500 0.192  1.189113e-01       8
#> TOB13            3.998772e-05 0.708087605 0.857 0.455  1.199632e-01       8
#> PCOLCE23         4.350640e-05 0.563183812 0.679 0.328  1.305192e-01       8
#> PALLD3           4.630546e-05 0.600276727 0.786 0.364  1.389164e-01       8
#> LDHB3            4.839460e-05 0.752170107 0.857 0.491  1.451838e-01       8
#> MAPK133          5.147408e-05 0.888958164 0.821 0.494  1.544223e-01       8
#> YES12            5.306117e-05 0.739100996 0.750 0.357  1.591835e-01       8
#> PTS3             5.310996e-05 0.806246019 0.786 0.478  1.593299e-01       8
#> DSP3             5.641913e-05 0.724061662 0.786 0.501  1.692574e-01       8
#> TPM24            6.407968e-05 0.640276755 0.857 0.427  1.922390e-01       8
#> CRABP13          6.480574e-05 0.738839375 0.714 0.404  1.944172e-01       8
#> ART32            6.758091e-05 0.788241070 0.536 0.218  2.027427e-01       8
#> ACHE             7.117167e-05 0.611971264 0.286 0.075  2.135150e-01       8
#> SMYD22           7.233304e-05 0.712284916 0.786 0.446  2.169991e-01       8
#> TNFRSF212        7.667270e-05 0.625431471 0.821 0.417  2.300181e-01       8
#> KLHL352          8.171627e-05 0.633507193 0.750 0.369  2.451488e-01       8
#> KCTD14           1.116508e-04 0.808682573 0.679 0.362  3.349525e-01       8
#> SYNCRIP2         1.125139e-04 0.710376659 0.750 0.509  3.375416e-01       8
#> SLC38A11         1.132845e-04 0.708227867 0.679 0.360  3.398536e-01       8
#> CDCA7L2          1.379669e-04 0.694401987 0.607 0.283  4.139007e-01       8
#> ANO12            1.468512e-04 0.530444353 0.429 0.151  4.405536e-01       8
#> E2F11            1.512011e-04 0.362802027 0.321 0.094  4.536034e-01       8
#> NRM1             1.569622e-04 0.614519521 0.536 0.218  4.708865e-01       8
#> BARX13           1.732390e-04 0.776108195 0.714 0.399  5.197170e-01       8
#> MARVELD11        1.795207e-04 0.636371967 0.679 0.321  5.385622e-01       8
#> KLK111           1.841707e-04 0.558152851 0.321 0.103  5.525122e-01       8
#> RP11-357H14.173  2.010168e-04 0.469375850 0.536 0.204  6.030505e-01       8
#> PRSS333          2.010534e-04 0.486248677 0.821 0.392  6.031603e-01       8
#> NEURL1B1         2.108196e-04 0.426358463 0.250 0.065  6.324587e-01       8
#> DHTKD12          2.121363e-04 0.797792979 0.643 0.424  6.364089e-01       8
#> IDH13            2.154710e-04 0.627661724 0.821 0.405  6.464130e-01       8
#> HMGB32           2.166270e-04 0.641377105 0.786 0.487  6.498809e-01       8
#> RPA32            2.255720e-04 0.656318469 0.750 0.481  6.767160e-01       8
#> RRM21            2.410594e-04 0.275426742 0.250 0.065  7.231782e-01       8
#> CYP39A13         2.597428e-04 0.584045756 0.679 0.359  7.792284e-01       8
#> ACTA11           2.717195e-04 0.572040855 0.393 0.138  8.151585e-01       8
#> PRIM11           2.730891e-04 0.495769955 0.393 0.139  8.192672e-01       8
#> AP1M23           2.830391e-04 0.658599185 0.821 0.460  8.491173e-01       8
#> UACA2            2.831056e-04 0.461753844 0.607 0.353  8.493167e-01       8
#> PHF191           2.883180e-04 0.644870783 0.679 0.324  8.649539e-01       8
#> AIF1L4           3.103858e-04 0.635030424 0.821 0.460  9.311575e-01       8
#> PDIA43           3.610585e-04 0.611438726 0.857 0.502  1.000000e+00       8
#> EPHX13           4.079766e-04 0.629392942 0.857 0.461  1.000000e+00       8
#> LNX13            4.429242e-04 0.558467538 0.607 0.355  1.000000e+00       8
#> CPED12           4.895459e-04 0.420521999 0.607 0.278  1.000000e+00       8
#> YEATS4           5.296060e-04 0.605910146 0.679 0.327  1.000000e+00       8
#> MEX3A4           5.385129e-04 0.526925926 0.643 0.351  1.000000e+00       8
#> DTL              5.558262e-04 0.224309427 0.250 0.068  1.000000e+00       8
#> TUBB63           5.699673e-04 0.668560950 0.714 0.471  1.000000e+00       8
#> PPIF2            6.419447e-04 0.518443567 0.750 0.401  1.000000e+00       8
#> DNPH12           6.535283e-04 0.614336168 0.821 0.531  1.000000e+00       8
#> ANXA2R3          6.649978e-04 0.328608484 0.536 0.256  1.000000e+00       8
#> SCARB12          7.160219e-04 0.559018386 0.786 0.458  1.000000e+00       8
#> CRISPLD13        7.306391e-04 0.583501841 0.821 0.446  1.000000e+00       8
#> COQ22            7.877390e-04 0.241014880 0.464 0.176  1.000000e+00       8
#> SLC25A192        7.909724e-04 0.420034053 0.393 0.149  1.000000e+00       8
#> P3H43            7.968217e-04 0.573468783 0.714 0.385  1.000000e+00       8
#> IDI13            8.064929e-04 0.523220410 0.786 0.457  1.000000e+00       8
#> NASP2            8.919990e-04 0.579595562 0.679 0.399  1.000000e+00       8
#> GLS3             9.406762e-04 0.607415298 0.893 0.456  1.000000e+00       8
#> MTHFD11          9.569770e-04 0.554036922 0.500 0.227  1.000000e+00       8
#> HAUS1            9.741793e-04 0.496113068 0.643 0.417  1.000000e+00       8
#> PIM12            1.186154e-03 0.536079256 0.786 0.446  1.000000e+00       8
#> TMPRSS32         1.192560e-03 0.520875219 0.536 0.259  1.000000e+00       8
#> INSIG11          1.292775e-03 0.432436299 0.643 0.385  1.000000e+00       8
#> SERPINH14        1.404334e-03 0.579807300 0.857 0.533  1.000000e+00       8
#> LRP23            1.457621e-03 0.429292593 0.536 0.244  1.000000e+00       8
#> EBP2             1.682318e-03 0.437875895 0.786 0.448  1.000000e+00       8
#> MARCKSL13        1.857197e-03 0.618429579 0.929 0.572  1.000000e+00       8
#> EIF3B2           1.947545e-03 0.493475495 0.821 0.478  1.000000e+00       8
#> FAM3C3           2.007551e-03 0.535630144 0.786 0.433  1.000000e+00       8
#> KDELC22          2.015845e-03 0.486099820 0.429 0.185  1.000000e+00       8
#> RASL121          2.095596e-03 0.384598941 0.321 0.119  1.000000e+00       8
#> WEE12            2.108971e-03 0.528681044 0.643 0.357  1.000000e+00       8
#> CITED43          2.121453e-03 0.548718764 0.750 0.481  1.000000e+00       8
#> EPB41L23         2.137424e-03 0.535704981 0.536 0.263  1.000000e+00       8
#> SLC7A52          2.544513e-03 0.507906119 0.607 0.318  1.000000e+00       8
#> HSPB13           2.648104e-03 0.499786100 0.750 0.495  1.000000e+00       8
#> IMPA23           2.671994e-03 0.577196892 0.821 0.503  1.000000e+00       8
#> FAM208B2         2.771237e-03 0.483101628 0.643 0.406  1.000000e+00       8
#> MTL53            2.818941e-03 0.462759793 0.679 0.402  1.000000e+00       8
#> AARD4            2.868421e-03 0.521057189 0.679 0.448  1.000000e+00       8
#> PAQR42           2.879469e-03 0.556529594 0.607 0.333  1.000000e+00       8
#> DONSON1          2.909386e-03 0.551424635 0.357 0.152  1.000000e+00       8
#> MCM51            2.941177e-03 0.388593813 0.643 0.311  1.000000e+00       8
#> CSRP13           3.073259e-03 0.463995272 0.821 0.452  1.000000e+00       8
#> SNHG192          3.239726e-03 0.386495798 0.714 0.413  1.000000e+00       8
#> RTN4RL23         3.294322e-03 0.398213928 0.429 0.190  1.000000e+00       8
#> NKD23            3.561602e-03 0.239216110 0.429 0.178  1.000000e+00       8
#> CTNNBIP12        3.721850e-03 0.486012089 0.750 0.478  1.000000e+00       8
#> SFN3             3.757620e-03 0.352228643 0.607 0.374  1.000000e+00       8
#> TNFRSF12A3       3.932220e-03 0.487582047 0.643 0.440  1.000000e+00       8
#> RAMP22           3.967533e-03 0.620173492 0.464 0.256  1.000000e+00       8
#> TTC39A3          4.083034e-03 0.671297738 0.500 0.267  1.000000e+00       8
#> CDK43            4.404635e-03 0.561992255 0.786 0.520  1.000000e+00       8
#> ENPP53           4.425091e-03 0.459884658 0.607 0.345  1.000000e+00       8
#> STAT12           4.794205e-03 0.408762243 0.714 0.438  1.000000e+00       8
#> TUBA1A3          4.806509e-03 0.469795366 0.750 0.483  1.000000e+00       8
#> NDUFAF62         4.891988e-03 0.442144760 0.679 0.452  1.000000e+00       8
#> COL9A21          5.021269e-03 0.417596527 0.286 0.112  1.000000e+00       8
#> FABP52           5.113971e-03 0.476244835 0.679 0.433  1.000000e+00       8
#> PPP1R1B2         5.212279e-03 0.443522095 0.679 0.385  1.000000e+00       8
#> C1QBP2           5.291949e-03 0.479016276 0.714 0.506  1.000000e+00       8
#> PHYH4            5.704040e-03 0.522381686 0.750 0.500  1.000000e+00       8
#> UCP22            5.743782e-03 0.300420781 0.571 0.355  1.000000e+00       8
#> COL11A14         5.955778e-03 0.344899244 0.571 0.383  1.000000e+00       8
#> C1QL43           6.054592e-03 0.119171768 0.607 0.283  1.000000e+00       8
#> HMGA13           6.125720e-03 0.477155813 0.786 0.516  1.000000e+00       8
#> MARS1            6.833537e-03 0.378255172 0.607 0.338  1.000000e+00       8
#> SYT83            6.907490e-03 0.244946052 0.643 0.371  1.000000e+00       8
#> BNIP32           7.125312e-03 0.439439026 0.679 0.472  1.000000e+00       8
#> MASTL1           7.375371e-03 0.257490045 0.286 0.111  1.000000e+00       8
#> DDIT41           7.682623e-03 0.352761994 0.679 0.436  1.000000e+00       8
#> LDHA1            7.935389e-03 0.411126331 0.821 0.488  1.000000e+00       8
#> PGP2             8.155310e-03 0.439400548 0.786 0.496  1.000000e+00       8
#> CLPSL12          8.576266e-03 0.279210578 0.571 0.329  1.000000e+00       8
#> ATP6V0E22        8.604068e-03 0.404826215 0.571 0.316  1.000000e+00       8
#> HELLS1           8.748685e-03 0.278874134 0.286 0.115  1.000000e+00       8
#> IGSF32           9.545220e-03 0.380103418 0.643 0.447  1.000000e+00       8
#> GGCT3            9.818306e-03 0.399167847 0.786 0.481  1.000000e+00       8
#> RIBC21           9.166234e-77 3.465912741 0.704 0.021  2.749870e-73       9
#> MCM101           7.017606e-59 3.548696129 0.889 0.059  2.105282e-55       9
#> DTL1             2.417974e-58 3.473248196 0.852 0.053  7.253921e-55       9
#> EXO11            9.729034e-52 2.967511156 0.741 0.043  2.918710e-48       9
#> CDC451           5.391993e-48 3.232693019 0.778 0.056  1.617598e-44       9
#> CDC61            4.520882e-47 2.414207727 0.741 0.048  1.356265e-43       9
#> PKMYT13          1.519505e-46 3.364424037 0.926 0.091  4.558515e-43       9
#> DSCC12           1.267003e-45 3.487349283 0.926 0.094  3.801008e-42       9
#> MCM22            2.236965e-43 3.376771234 1.000 0.122  6.710896e-40       9
#> CLSPN2           7.681818e-43 2.462834830 0.815 0.069  2.304545e-39       9
#> CCNE2            7.943647e-43 3.119607167 0.741 0.058  2.383094e-39       9
#> UHRF11           1.010436e-42 2.804629558 0.889 0.088  3.031309e-39       9
#> ASF1B2           1.759741e-41 2.741484796 0.889 0.089  5.279223e-38       9
#> EME11            4.364728e-41 1.948592324 0.519 0.024  1.309418e-37       9
#> RRM22            4.631820e-41 1.249190050 0.741 0.053  1.389546e-37       9
#> ESCO21           1.511404e-40 1.595191518 0.630 0.038  4.534211e-37       9
#> SMC1B2           3.053335e-40 3.030311040 0.926 0.106  9.160005e-37       9
#> ZNF3671          5.144910e-39 2.230125643 0.556 0.032  1.543473e-35       9
#> SGOL13           5.170660e-38 2.030634001 0.926 0.100  1.551198e-34       9
#> GINS11           6.418684e-38 2.359460491 0.778 0.071  1.925605e-34       9
#> SPC253           1.452824e-37 2.545574877 0.852 0.090  4.358473e-34       9
#> NCAPG2           6.345294e-36 1.678363658 0.778 0.071  1.903588e-32       9
#> CDCA52           1.180833e-35 2.147218129 0.778 0.075  3.542499e-32       9
#> HELLS2           2.801929e-35 1.964024983 0.889 0.100  8.405787e-32       9
#> E2F12            6.387118e-35 2.770756324 0.778 0.082  1.916136e-31       9
#> CDCA71           1.572401e-34 2.296275304 0.704 0.064  4.717204e-31       9
#> C16orf591        7.120523e-34 1.927295469 0.778 0.079  2.136157e-30       9
#> FAM111B          8.079333e-34 1.726359227 0.630 0.050  2.423800e-30       9
#> CDT13            1.093574e-33 2.653298036 0.926 0.130  3.280722e-30       9
#> CDK13            1.509456e-33 2.099861614 0.963 0.134  4.528368e-30       9
#> DHFR2            6.150298e-33 2.460171612 1.000 0.156  1.845089e-29       9
#> KIFC13           7.183056e-33 2.187583437 0.889 0.114  2.154917e-29       9
#> PBK2             1.724039e-32 2.009263059 0.815 0.092  5.172118e-29       9
#> ZWINT3           9.270253e-31 2.826191192 0.926 0.156  2.781076e-27       9
#> MYBL23           1.087114e-30 2.272429425 0.963 0.148  3.261342e-27       9
#> RAD511           1.782316e-30 2.085117722 0.704 0.072  5.346948e-27       9
#> RAD54L           2.208169e-30 1.818834741 0.407 0.021  6.624507e-27       9
#> CHEK12           2.698001e-30 2.735785550 1.000 0.189  8.094002e-27       9
#> RAD51AP13        6.829190e-30 2.741190113 0.963 0.174  2.048757e-26       9
#> NUF23            2.027681e-29 1.236543562 0.889 0.113  6.083043e-26       9
#> TK13             2.137193e-29 3.160035768 1.000 0.198  6.411578e-26       9
#> RMI23            1.270013e-28 3.146704734 1.000 0.235  3.810040e-25       9
#> TCF192           1.718982e-28 2.768682674 0.963 0.186  5.156947e-25       9
#> MELK3            2.101476e-28 2.282138920 0.815 0.114  6.304428e-25       9
#> RECQL42          5.720463e-28 2.449929628 0.815 0.121  1.716139e-24       9
#> CENPM3           9.774113e-28 2.318113206 0.815 0.120  2.932234e-24       9
#> TRIP132          1.934553e-27 2.187271938 0.852 0.129  5.803660e-24       9
#> RAD54B1          3.581113e-27 2.004572819 0.593 0.057  1.074334e-23       9
#> TPX23            3.983754e-27 1.383512657 0.926 0.142  1.195126e-23       9
#> CENPU3           6.806646e-27 2.245841025 0.889 0.153  2.041994e-23       9
#> ELFN1-AS1        1.416834e-26 2.119905429 0.519 0.044  4.250502e-23       9
#> GINS22           1.753755e-26 2.665538919 0.889 0.161  5.261266e-23       9
#> KIAA01013        2.706021e-26 2.971554120 0.963 0.214  8.118062e-23       9
#> AC004381.61      9.450787e-26 1.859849077 0.556 0.052  2.835236e-22       9
#> BRCA11           2.922982e-25 1.595359413 0.704 0.088  8.768947e-22       9
#> BRCA23           3.063126e-25 1.870376517 0.889 0.146  9.189379e-22       9
#> NCAPH2           3.216041e-25 1.353484408 0.481 0.038  9.648123e-22       9
#> PDZK1            4.030576e-25 1.547578061 0.259 0.009  1.209173e-21       9
#> GMNN3            4.124506e-25 2.741204200 1.000 0.258  1.237352e-21       9
#> DONSON2          5.042186e-25 2.072990018 0.852 0.140  1.512656e-21       9
#> CNIH2            5.695340e-25 1.868817741 0.370 0.022  1.708602e-21       9
#> LMNB13           6.688514e-25 1.821781508 0.926 0.163  2.006554e-21       9
#> SPC241           9.388367e-25 1.795450430 0.519 0.048  2.816510e-21       9
#> KIF152           9.963041e-25 1.455061582 0.593 0.061  2.988912e-21       9
#> ATAD23           1.056149e-24 2.749857695 1.000 0.266  3.168447e-21       9
#> FBXO51           2.128992e-24 1.383716357 0.593 0.062  6.386977e-21       9
#> TYMS3            2.983465e-24 1.522535313 0.926 0.164  8.950394e-21       9
#> MAD2L13          3.657825e-24 1.871259364 0.926 0.175  1.097348e-20       9
#> BIRC53           5.939174e-24 1.518954394 0.889 0.154  1.781752e-20       9
#> UBE2T2           1.029231e-23 2.907979924 1.000 0.280  3.087694e-20       9
#> TOP2A3           1.221693e-23 1.163509270 0.852 0.129  3.665079e-20       9
#> RFC32            1.774542e-23 2.776249292 0.963 0.273  5.323627e-20       9
#> PRC13            2.269421e-23 1.418732323 0.963 0.176  6.808262e-20       9
#> NUSAP13          3.682341e-23 1.267519038 0.852 0.135  1.104702e-19       9
#> CCNA23           4.772797e-23 1.014475056 0.741 0.097  1.431839e-19       9
#> RTKN21           4.888657e-23 1.152410564 0.556 0.056  1.466597e-19       9
#> FEN13            7.845412e-23 2.224960149 0.926 0.207  2.353624e-19       9
#> RFC43            1.148152e-22 2.225107450 0.963 0.236  3.444455e-19       9
#> CDCA22           3.667421e-22 1.520144425 0.667 0.088  1.100226e-18       9
#> MASTL2           8.340919e-22 1.523928625 0.704 0.101  2.502276e-18       9
#> TIMELESS1        9.170494e-22 2.319670482 0.741 0.129  2.751148e-18       9
#> MCM42            1.408568e-21 2.927691777 1.000 0.342  4.225705e-18       9
#> RNASEH2A2        2.171944e-21 2.421032159 0.963 0.269  6.515831e-18       9
#> C9orf135         2.900761e-21 1.660842823 0.444 0.040  8.702283e-18       9
#> CENPK2           2.985702e-21 1.259533934 0.667 0.090  8.957107e-18       9
#> C21orf582        5.387879e-21 1.255472109 0.704 0.103  1.616364e-17       9
#> CENPH3           1.175923e-20 1.699492832 0.852 0.164  3.527768e-17       9
#> RASL11B1         1.631659e-20 1.449251086 0.741 0.119  4.894978e-17       9
#> ORC63            2.545641e-20 1.760363469 0.778 0.141  7.636924e-17       9
#> UBE2C3           2.784761e-20 1.811240696 0.926 0.202  8.354283e-17       9
#> MCM73            3.805459e-20 2.978868039 1.000 0.377  1.141638e-16       9
#> ERVMER34-13      6.803712e-20 2.223854652 0.889 0.221  2.041113e-16       9
#> RRM13            1.272294e-19 1.957963613 1.000 0.280  3.816881e-16       9
#> ACTL82           1.562397e-19 2.066717983 0.704 0.126  4.687191e-16       9
#> DSN12            2.031976e-19 2.165414636 0.889 0.238  6.095927e-16       9
#> DIAPH32          2.058714e-19 1.536869166 0.778 0.141  6.176142e-16       9
#> PCNA2            2.108634e-19 2.832037929 1.000 0.400  6.325902e-16       9
#> TMEM973          3.052217e-19 2.126112031 0.926 0.248  9.156652e-16       9
#> DNAJC92          3.811630e-19 2.543583294 1.000 0.422  1.143489e-15       9
#> GINS42           7.072866e-19 1.382411870 0.852 0.170  2.121860e-15       9
#> MND13            1.114783e-18 1.500049893 0.630 0.096  3.344350e-15       9
#> FANCI3           1.317671e-18 1.452725121 0.667 0.107  3.953014e-15       9
#> CENPF3           1.426396e-18 1.255956137 0.889 0.187  4.279188e-15       9
#> UNG2             1.682890e-18 2.489546312 0.889 0.260  5.048670e-15       9
#> WFDC12           2.257555e-18 1.273949937 0.630 0.093  6.772664e-15       9
#> RAC33            2.565151e-18 1.705815956 0.889 0.205  7.695453e-15       9
#> SLBP4            4.457261e-18 2.301472689 1.000 0.491  1.337178e-14       9
#> TEX302           4.804308e-18 2.025751135 1.000 0.323  1.441292e-14       9
#> CHAF1A2          4.831824e-18 1.567585706 0.778 0.151  1.449547e-14       9
#> RACGAP13         6.309707e-18 1.324861696 0.815 0.161  1.892912e-14       9
#> MCM32            7.515459e-18 2.494802860 0.963 0.382  2.254638e-14       9
#> CENPW2           7.537442e-18 1.560128325 1.000 0.283  2.261233e-14       9
#> TTF23            1.081723e-17 1.854126518 0.963 0.323  3.245170e-14       9
#> DHTKD13          1.095677e-17 2.080069543 1.000 0.415  3.287031e-14       9
#> PRIM12           1.199732e-17 1.560916741 0.704 0.132  3.599196e-14       9
#> POC1A2           1.304398e-17 1.595408629 0.704 0.130  3.913195e-14       9
#> SHCBP12          1.590072e-17 1.017518518 0.593 0.083  4.770216e-14       9
#> FAM83D1          1.662393e-17 0.875564300 0.593 0.085  4.987179e-14       9
#> WDR76            1.737765e-17 1.026148432 0.444 0.049  5.213295e-14       9
#> GTSE12           2.170320e-17 1.205357728 0.593 0.089  6.510961e-14       9
#> MKI673           2.256279e-17 0.992676408 0.741 0.130  6.768837e-14       9
#> ATP2A1-AS11      2.433793e-17 1.428307944 0.519 0.070  7.301379e-14       9
#> CHAF1B2          2.504701e-17 1.851702319 0.778 0.172  7.514103e-14       9
#> KIF20B3          6.814889e-17 1.004570927 0.852 0.173  2.044467e-13       9
#> DBF41            6.868241e-17 1.614874814 0.926 0.267  2.060472e-13       9
#> CCNB13           1.578454e-16 0.801671230 0.815 0.168  4.735362e-13       9
#> NRM2             1.696412e-16 1.760837199 0.852 0.210  5.089237e-13       9
#> SSRP14           1.933573e-16 1.827335490 1.000 0.488  5.800720e-13       9
#> SYCP23           2.121647e-16 1.532198032 0.926 0.236  6.364940e-13       9
#> MTHFD12          2.639778e-16 1.574950311 0.889 0.218  7.919333e-13       9
#> HMMR2            2.759882e-16 0.493613769 0.630 0.094  8.279647e-13       9
#> TACC33           3.116612e-16 1.398838035 0.815 0.192  9.349837e-13       9
#> BARD14           3.437661e-16 1.617316099 1.000 0.314  1.031298e-12       9
#> PAICS3           3.662947e-16 1.968634220 1.000 0.417  1.098884e-12       9
#> AURKB3           4.126195e-16 1.271349064 0.630 0.107  1.237858e-12       9
#> ANLN3            5.019238e-16 0.879257325 0.630 0.103  1.505771e-12       9
#> CDCA82           5.556525e-16 0.665063981 0.556 0.078  1.666957e-12       9
#> CDC203           6.183629e-16 0.818195154 0.778 0.156  1.855089e-12       9
#> KIF18A1          6.326007e-16 0.781321029 0.333 0.030  1.897802e-12       9
#> CKS1B4           6.573768e-16 1.657139667 1.000 0.504  1.972131e-12       9
#> GGH3             7.170243e-16 2.058467850 1.000 0.383  2.151073e-12       9
#> CENPN2           7.222360e-16 1.530575889 0.963 0.284  2.166708e-12       9
#> SPSB4            7.662908e-16 0.965506628 0.444 0.052  2.298872e-12       9
#> C6orf1413        7.822953e-16 1.045098943 0.741 0.139  2.346886e-12       9
#> LRRCC12          1.050538e-15 1.920061179 0.963 0.405  3.151613e-12       9
#> CDCA42           1.085494e-15 1.633754816 0.889 0.239  3.256481e-12       9
#> KIF4A2           1.805127e-15 0.697633064 0.519 0.071  5.415382e-12       9
#> DNMT12           2.110754e-15 1.790671007 1.000 0.369  6.332263e-12       9
#> NEMP21           3.537282e-15 1.254605869 0.556 0.089  1.061184e-11       9
#> SKA12            3.781388e-15 1.072311872 0.481 0.066  1.134416e-11       9
#> SLC29A2          4.367912e-15 1.376678825 0.407 0.050  1.310374e-11       9
#> PRKDC3           7.405053e-15 1.831824185 0.963 0.503  2.221516e-11       9
#> TMEM106C4        7.641888e-15 1.783887262 1.000 0.411  2.292566e-11       9
#> MCM52            8.155837e-15 1.715530959 0.926 0.304  2.446751e-11       9
#> KIF232           8.385236e-15 0.942147717 0.556 0.088  2.515571e-11       9
#> EZH22            8.899497e-15 1.695735172 0.852 0.243  2.669849e-11       9
#> YEATS41          8.962110e-15 1.556637006 0.963 0.321  2.688633e-11       9
#> DTYMK3           1.258269e-14 1.645391348 0.963 0.358  3.774806e-11       9
#> CDCA7L3          1.728744e-14 1.639026467 0.889 0.277  5.186233e-11       9
#> CA62             2.624855e-14 1.358621993 0.593 0.108  7.874565e-11       9
#> TUBA1B3          2.914023e-14 1.662281888 1.000 0.482  8.742070e-11       9
#> CKS24            3.209426e-14 1.570534656 1.000 0.436  9.628277e-11       9
#> C16orf741        3.970719e-14 1.074052637 0.852 0.207  1.191216e-10       9
#> THEM63           4.198294e-14 1.721815048 0.852 0.293  1.259488e-10       9
#> SGOL22           4.273749e-14 0.856183076 0.630 0.116  1.282125e-10       9
#> H2AFZ3           4.569315e-14 1.678084342 1.000 0.446  1.370794e-10       9
#> CCNB23           4.967240e-14 0.888898094 0.778 0.181  1.490172e-10       9
#> HAPLN13          5.007823e-14 1.451828190 0.704 0.163  1.502347e-10       9
#> DKC13            6.034377e-14 1.583542816 1.000 0.474  1.810313e-10       9
#> KIAA15242        7.203661e-14 0.878141772 0.519 0.081  2.161098e-10       9
#> DEK2             8.132998e-14 1.443337109 1.000 0.504  2.439899e-10       9
#> PRPS22           8.259471e-14 1.727729869 0.963 0.387  2.477841e-10       9
#> SKA31            8.692511e-14 1.273920574 0.407 0.054  2.607753e-10       9
#> RPA33            8.806345e-14 1.655643897 0.963 0.476  2.641904e-10       9
#> TMSB15A2         9.087244e-14 1.312644961 0.778 0.185  2.726173e-10       9
#> HMGB23           1.097416e-13 1.632246138 1.000 0.393  3.292247e-10       9
#> PARPBP3          1.656330e-13 0.901965885 0.593 0.107  4.968989e-10       9
#> HACD14           1.796857e-13 1.523284993 0.963 0.402  5.390572e-10       9
#> KIF2C2           2.261313e-13 0.586252404 0.481 0.071  6.783940e-10       9
#> FOXM13           2.385944e-13 1.333226887 0.593 0.119  7.157831e-10       9
#> GTF3A2           2.693962e-13 1.598738469 1.000 0.493  8.081885e-10       9
#> DCTPP13          3.087801e-13 1.722740257 0.963 0.493  9.263404e-10       9
#> LSM43            3.161538e-13 1.399112477 1.000 0.488  9.484614e-10       9
#> ALYREF3          3.550501e-13 1.781137501 0.926 0.343  1.065150e-09       9
#> CDK44            3.551914e-13 1.629470280 1.000 0.515  1.065574e-09       9
#> TROAP3           3.592012e-13 0.650827799 0.556 0.093  1.077604e-09       9
#> HMGB13           3.771417e-13 1.519886458 0.963 0.489  1.131425e-09       9
#> HIST1H1A         3.988777e-13 0.735449567 0.259 0.022  1.196633e-09       9
#> C9orf402         4.970814e-13 1.599136585 0.815 0.248  1.491244e-09       9
#> RTN4RL24         5.000524e-13 1.297467864 0.741 0.182  1.500157e-09       9
#> MIS18A2          5.493662e-13 1.302713519 0.926 0.290  1.648099e-09       9
#> PRR72            5.591622e-13 1.514786162 0.741 0.214  1.677487e-09       9
#> VGF3             5.793519e-13 1.248285201 0.704 0.164  1.738056e-09       9
#> TOMM403          6.098609e-13 1.570798386 0.963 0.464  1.829583e-09       9
#> SNRNP252         6.978879e-13 1.610090619 0.963 0.450  2.093664e-09       9
#> NME13            7.448747e-13 1.637779061 0.963 0.488  2.234624e-09       9
#> CTNNAL12         8.251438e-13 1.366359875 0.852 0.259  2.475431e-09       9
#> AURKA2           1.021344e-12 0.593929684 0.667 0.136  3.064033e-09       9
#> GPR12            1.150413e-12 0.908665610 0.259 0.023  3.451240e-09       9
#> SMC23            1.298558e-12 1.340188273 0.889 0.294  3.895674e-09       9
#> HSPD12           1.328937e-12 1.409836657 1.000 0.484  3.986810e-09       9
#> HYLS12           1.516218e-12 1.165051127 0.704 0.164  4.548654e-09       9
#> CACYBP5          1.936290e-12 1.284667183 1.000 0.526  5.808871e-09       9
#> CCNE1            2.072121e-12 1.250982796 0.333 0.041  6.216363e-09       9
#> WDR344           2.170053e-12 1.602867877 0.963 0.493  6.510159e-09       9
#> CDCA33           2.359292e-12 1.107570158 0.630 0.139  7.077877e-09       9
#> TIMM104          2.659094e-12 1.341113650 1.000 0.480  7.977281e-09       9
#> APOBEC3B3        2.698842e-12 1.168049624 0.667 0.150  8.096526e-09       9
#> CENPA3           3.577728e-12 0.637097353 0.519 0.089  1.073318e-08       9
#> KIF141           4.945259e-12 0.771360948 0.444 0.070  1.483578e-08       9
#> CCT54            4.973928e-12 1.296378997 1.000 0.521  1.492178e-08       9
#> HIST1H1D1        5.103096e-12 1.290853874 0.370 0.053  1.530929e-08       9
#> MCM61            5.725958e-12 0.962331827 0.556 0.109  1.717787e-08       9
#> LAPTM4B5         6.391362e-12 1.363360388 0.963 0.545  1.917409e-08       9
#> EXOSC81          7.821428e-12 1.608192757 0.963 0.430  2.346428e-08       9
#> HMGN23           8.620857e-12 1.335973961 0.963 0.475  2.586257e-08       9
#> NES3             8.984338e-12 1.288634179 0.852 0.336  2.695301e-08       9
#> BUB12            1.226988e-11 0.493685211 0.444 0.069  3.680965e-08       9
#> PTTG13           1.270467e-11 0.710629032 0.778 0.218  3.811400e-08       9
#> STMN14           1.337384e-11 1.343670327 1.000 0.522  4.012153e-08       9
#> YBX23            1.434866e-11 0.837702186 0.630 0.136  4.304597e-08       9
#> NETO22           1.706001e-11 0.917138468 0.704 0.163  5.118003e-08       9
#> HIST1H4C4        1.762749e-11 1.863948765 0.778 0.334  5.288246e-08       9
#> GYLTL1B1         2.149723e-11 1.378676528 0.593 0.139  6.449170e-08       9
#> SLC43A34         2.380298e-11 1.372473884 1.000 0.447  7.140893e-08       9
#> SMC43            2.495977e-11 1.322360383 0.889 0.366  7.487932e-08       9
#> KLK13            2.924524e-11 0.892366384 0.630 0.136  8.773571e-08       9
#> KIF223           3.222651e-11 1.605358112 0.926 0.429  9.667954e-08       9
#> SFN4             3.274100e-11 1.251591667 0.926 0.366  9.822300e-08       9
#> LDHB4            3.632934e-11 1.485327972 0.926 0.490  1.089880e-07       9
#> ODC12            3.715354e-11 1.371036059 0.926 0.385  1.114606e-07       9
#> DUT3             3.810060e-11 1.356526722 0.926 0.459  1.143018e-07       9
#> PDGFRA4          4.004877e-11 1.176944182 0.963 0.327  1.201463e-07       9
#> GGCT4            4.413390e-11 1.244736728 0.926 0.478  1.324017e-07       9
#> HAUS11           4.692991e-11 1.405432050 0.926 0.410  1.407897e-07       9
#> LYAR3            5.245505e-11 1.206645841 0.889 0.353  1.573652e-07       9
#> PLEKHB14         5.681631e-11 1.314222854 1.000 0.480  1.704489e-07       9
#> RAD215           6.088764e-11 1.185484732 1.000 0.529  1.826629e-07       9
#> NASP3            6.131431e-11 1.512862947 0.889 0.393  1.839429e-07       9
#> SOSTDC1          6.243707e-11 0.756777682 0.370 0.054  1.873112e-07       9
#> TBC1D14          6.349775e-11 1.296017846 0.926 0.450  1.904932e-07       9
#> NTHL13           6.534594e-11 1.384291087 0.926 0.410  1.960378e-07       9
#> GDPD21           6.559062e-11 0.716432225 0.593 0.128  1.967719e-07       9
#> NDC803           7.796356e-11 0.540049595 0.519 0.096  2.338907e-07       9
#> SPDEF1           8.323027e-11 0.804594471 0.407 0.065  2.496908e-07       9
#> CRABP14          8.359972e-11 1.285708118 1.000 0.397  2.507992e-07       9
#> LRR13            8.851409e-11 1.171551907 0.667 0.175  2.655423e-07       9
#> MXD32            9.619384e-11 0.864238930 0.481 0.092  2.885815e-07       9
#> AZGP14           9.767599e-11 1.436561300 1.000 0.364  2.930280e-07       9
#> FANCD21          1.001286e-10 0.812724762 0.630 0.150  3.003858e-07       9
#> IGFBP23          1.023354e-10 1.335372258 0.926 0.402  3.070062e-07       9
#> FKBP44           1.070474e-10 1.254038376 0.963 0.492  3.211421e-07       9
#> CA84             1.300934e-10 1.206624011 0.815 0.264  3.902801e-07       9
#> CKB5             1.855600e-10 1.275219016 0.926 0.444  5.566799e-07       9
#> RAMP23           2.180628e-10 0.858257079 0.815 0.248  6.541885e-07       9
#> PAQR43           2.196774e-10 1.443943021 0.852 0.327  6.590323e-07       9
#> DNPH13           2.203103e-10 1.223612197 1.000 0.527  6.609309e-07       9
#> P3H44            2.356606e-10 1.266776837 0.926 0.380  7.069819e-07       9
#> RANBP13          2.525875e-10 1.245333907 0.963 0.498  7.577624e-07       9
#> BYSL3            2.528311e-10 1.299320634 0.963 0.421  7.584933e-07       9
#> PHGDH3           2.553577e-10 1.342830769 0.889 0.398  7.660732e-07       9
#> LEFTY23          2.984506e-10 1.252742679 0.630 0.189  8.953518e-07       9
#> GAMT2            3.377131e-10 1.218719209 0.741 0.236  1.013139e-06       9
#> CCNF1            3.481412e-10 0.882890517 0.519 0.110  1.044424e-06       9
#> NUCKS14          3.528500e-10 1.183148237 0.963 0.496  1.058550e-06       9
#> NXPH4            3.769474e-10 0.818823940 0.407 0.071  1.130842e-06       9
#> CLPSL13          3.778081e-10 1.459995647 0.778 0.324  1.133424e-06       9
#> PYCR14           4.286946e-10 1.376407433 0.926 0.451  1.286084e-06       9
#> TUBB4            4.301130e-10 1.181630893 0.963 0.520  1.290339e-06       9
#> LNX14            4.330760e-10 1.142846703 0.926 0.347  1.299228e-06       9
#> RP11-357H14.174  4.508914e-10 0.914905753 0.741 0.199  1.352674e-06       9
#> DUSP91           4.515517e-10 0.833295624 0.296 0.039  1.354655e-06       9
#> DEPDC12          4.610663e-10 0.560485025 0.407 0.068  1.383199e-06       9
#> ANP32E4          6.803729e-10 1.157787752 0.963 0.433  2.041119e-06       9
#> FAM111A1         7.183530e-10 0.415441504 0.519 0.102  2.155059e-06       9
#> RCCD12           7.188947e-10 0.915094145 0.778 0.229  2.156684e-06       9
#> VSNL12           8.453905e-10 0.898387219 0.519 0.112  2.536172e-06       9
#> MYBL14           9.482435e-10 1.041870088 0.778 0.255  2.844730e-06       9
#> GOLT1A2          1.036971e-09 0.705310526 0.556 0.121  3.110914e-06       9
#> PLK13            1.084837e-09 0.440194032 0.556 0.127  3.254511e-06       9
#> PPIL13           1.145209e-09 1.248802218 0.889 0.368  3.435628e-06       9
#> HIST2H2AC1       1.451854e-09 1.248905096 0.481 0.106  4.355561e-06       9
#> SUN33            1.525316e-09 1.033205586 0.704 0.215  4.575948e-06       9
#> NCAPD23          1.532463e-09 0.843659956 0.593 0.146  4.597388e-06       9
#> AIF1L5           2.051578e-09 1.161519899 0.963 0.457  6.154734e-06       9
#> MYEOV            2.807658e-09 0.560895478 0.444 0.085  8.422974e-06       9
#> TMPO3            2.928167e-09 1.096189532 0.926 0.356  8.784501e-06       9
#> SYT84            2.939794e-09 1.184644811 0.889 0.365  8.819382e-06       9
#> SUSD52           2.995906e-09 0.960731621 0.519 0.120  8.987719e-06       9
#> B4GALNT1         3.102052e-09 0.902437126 0.259 0.034  9.306157e-06       9
#> PSRC12           3.133714e-09 0.751943586 0.630 0.162  9.401143e-06       9
#> SLC2A4RG4        3.215325e-09 1.369365100 0.889 0.433  9.645974e-06       9
#> AQP53            3.261328e-09 1.252538842 0.889 0.426  9.783983e-06       9
#> SLC25A193        3.888703e-09 0.629006709 0.593 0.144  1.166611e-05       9
#> RPL39L4          4.174218e-09 1.124835214 0.926 0.461  1.252265e-05       9
#> FAM19A3          4.605675e-09 0.906833121 0.333 0.055  1.381702e-05       9
#> CACNB4           4.715123e-09 0.563657135 0.333 0.052  1.414537e-05       9
#> C12orf753        4.810866e-09 1.037757093 0.778 0.271  1.443260e-05       9
#> CDKN2C2          4.830695e-09 0.655126259 0.556 0.130  1.449209e-05       9
#> CMBL1            4.837201e-09 0.861438311 0.667 0.189  1.451160e-05       9
#> AHCY2            5.183419e-09 1.158870789 0.852 0.379  1.555026e-05       9
#> HIST1H2BN3       5.263942e-09 0.802151477 0.667 0.189  1.579183e-05       9
#> FBLN13           5.511318e-09 0.719333088 0.741 0.262  1.653395e-05       9
#> SYNGR13          6.164912e-09 1.139145325 0.815 0.360  1.849474e-05       9
#> SHISA2           6.691727e-09 0.667366568 0.296 0.044  2.007518e-05       9
#> C19orf482        6.912441e-09 1.226243940 0.778 0.297  2.073732e-05       9
#> HJURP2           6.924006e-09 0.435614413 0.370 0.064  2.077202e-05       9
#> FH3              7.081868e-09 1.195475322 0.889 0.463  2.124560e-05       9
#> PITX14           7.253910e-09 1.086241359 0.889 0.410  2.176173e-05       9
#> AP1M24           7.910937e-09 1.198227196 0.852 0.460  2.373281e-05       9
#> COL2A13          8.236325e-09 1.113897748 0.778 0.310  2.470897e-05       9
#> SCRG13           8.372779e-09 1.110844549 0.963 0.375  2.511834e-05       9
#> TMX24            8.471129e-09 1.095746874 1.000 0.464  2.541339e-05       9
#> CKAP23           8.641358e-09 0.654421343 0.741 0.234  2.592407e-05       9
#> VDR1             8.880524e-09 0.814647583 0.593 0.152  2.664157e-05       9
#> CBR31            9.441250e-09 0.723889278 0.630 0.168  2.832375e-05       9
#> HMGA14           1.094615e-08 1.064568472 0.963 0.512  3.283845e-05       9
#> SKA24            1.195889e-08 0.995546562 0.815 0.299  3.587667e-05       9
#> HIST1H1E1        1.224119e-08 0.723832887 0.296 0.046  3.672356e-05       9
#> ADAM155          1.257917e-08 1.166982336 0.889 0.502  3.773751e-05       9
#> KPNA21           1.299689e-08 1.233489114 0.852 0.457  3.899067e-05       9
#> PRSS334          1.338963e-08 1.128378676 0.926 0.390  4.016890e-05       9
#> NT5DC22          1.656443e-08 1.230464111 0.852 0.464  4.969329e-05       9
#> HSD11B2          1.787605e-08 0.601063189 0.296 0.047  5.362814e-05       9
#> CENPQ3           1.814277e-08 0.945565339 0.815 0.373  5.442832e-05       9
#> MARC12           1.815424e-08 0.692323801 0.556 0.137  5.446273e-05       9
#> POLD23           1.816529e-08 1.075007267 0.926 0.440  5.449588e-05       9
#> EFHD14           1.927750e-08 1.132859399 0.889 0.421  5.783249e-05       9
#> PSAT12           2.157514e-08 0.873200252 0.444 0.100  6.472541e-05       9
#> BOP13            2.250495e-08 1.167691934 0.852 0.449  6.751486e-05       9
#> RUVBL13          2.271507e-08 1.173261654 0.889 0.434  6.814520e-05       9
#> H2AFX3           2.603240e-08 0.958911331 0.852 0.364  7.809719e-05       9
#> PHLDA22          3.370558e-08 0.983160722 0.778 0.312  1.011167e-04       9
#> SLC43A11         3.948636e-08 1.079246489 0.519 0.144  1.184591e-04       9
#> FAM64A2          4.085970e-08 0.563824381 0.407 0.085  1.225791e-04       9
#> DHCR241          4.110490e-08 0.737480320 0.593 0.163  1.233147e-04       9
#> ACOT73           4.150699e-08 1.100055150 0.704 0.250  1.245210e-04       9
#> MGST14           4.435919e-08 1.103888288 0.963 0.456  1.330776e-04       9
#> ABCA7            4.576206e-08 0.962095570 0.444 0.105  1.372862e-04       9
#> MDC12            4.989845e-08 0.692946285 0.481 0.117  1.496954e-04       9
#> HSP90AB16        5.027599e-08 0.831753464 1.000 0.664  1.508280e-04       9
#> HN13             5.190391e-08 0.974869988 0.963 0.494  1.557117e-04       9
#> NUDT13           5.312705e-08 1.009199660 0.778 0.375  1.593812e-04       9
#> MAPK134          5.786838e-08 1.143930560 0.889 0.493  1.736051e-04       9
#> LMNB22           6.771048e-08 0.810693584 0.593 0.171  2.031315e-04       9
#> SIGMAR13         6.787577e-08 1.112253128 0.778 0.382  2.036273e-04       9
#> MARCKSL14        6.791322e-08 0.918374760 0.963 0.571  2.037397e-04       9
#> PODXL24          7.590436e-08 1.097983054 0.889 0.452  2.277131e-04       9
#> SERTAD44         8.028865e-08 1.045465880 0.926 0.467  2.408659e-04       9
#> RUVBL22          8.513826e-08 0.983329112 0.815 0.399  2.554148e-04       9
#> KRT186           9.349074e-08 0.988361150 0.963 0.544  2.804722e-04       9
#> TTYH13           9.775408e-08 1.044376495 0.926 0.432  2.932623e-04       9
#> SAPCD22          9.884580e-08 0.456686401 0.444 0.099  2.965374e-04       9
#> SEPT43           1.050338e-07 0.668270719 0.481 0.121  3.151015e-04       9
#> CMSS14           1.065065e-07 0.940285029 0.778 0.324  3.195195e-04       9
#> CEP553           1.120243e-07 0.463831167 0.444 0.100  3.360728e-04       9
#> CTHRC14          1.123453e-07 1.018733021 0.963 0.464  3.370358e-04       9
#> ITGA103          1.190169e-07 0.809663088 0.630 0.216  3.570508e-04       9
#> CDKN33           1.286419e-07 0.464348017 0.519 0.132  3.859257e-04       9
#> PRSS31           1.335652e-07 0.602293215 0.333 0.064  4.006956e-04       9
#> PGP3             1.371747e-07 1.074221386 0.963 0.492  4.115242e-04       9
#> MSX11            1.511181e-07 0.640143879 0.407 0.091  4.533544e-04       9
#> PSIP13           1.648930e-07 0.960221873 0.852 0.340  4.946790e-04       9
#> NKD24            1.838076e-07 0.673387634 0.593 0.174  5.514229e-04       9
#> SERPINA52        1.870041e-07 1.125521570 0.778 0.332  5.610123e-04       9
#> H2AFV3           1.881187e-07 0.972290607 0.889 0.487  5.643562e-04       9
#> KLHDC34          1.906763e-07 0.893569693 0.963 0.567  5.720290e-04       9
#> TRAF3IP33        1.944731e-07 0.876540865 0.704 0.295  5.834193e-04       9
#> SPDL12           2.068582e-07 0.499862863 0.519 0.135  6.205746e-04       9
#> DEPDC1B1         2.124186e-07 0.547217770 0.296 0.052  6.372557e-04       9
#> CD3203           2.177868e-07 0.876792383 0.741 0.295  6.533603e-04       9
#> HRCT13           2.305732e-07 0.908504220 0.889 0.387  6.917196e-04       9
#> PCBD13           2.339505e-07 0.972981563 0.926 0.501  7.018516e-04       9
#> RBBP73           2.443207e-07 0.952725454 0.926 0.512  7.329621e-04       9
#> SERPINE1         2.506323e-07 0.364272378 0.333 0.064  7.518969e-04       9
#> THOP13           2.525054e-07 0.375719115 0.667 0.207  7.575163e-04       9
#> IDI14            2.817476e-07 0.916292205 0.926 0.453  8.452428e-04       9
#> CYP39A14         3.348755e-07 0.790692211 0.889 0.354  1.004626e-03       9
#> YDJC2            3.762628e-07 0.971066083 0.852 0.460  1.128788e-03       9
#> TTK1             3.867703e-07 0.261410590 0.370 0.076  1.160311e-03       9
#> WIF11            3.905685e-07 0.362395602 0.481 0.121  1.171706e-03       9
#> CCDC342          3.965889e-07 0.508904990 0.667 0.202  1.189767e-03       9
#> MYOZ12           4.036439e-07 0.905171073 0.667 0.267  1.210932e-03       9
#> CPNE71           4.072635e-07 0.407427721 0.296 0.053  1.221791e-03       9
#> HIBCH4           4.266531e-07 0.981198955 0.889 0.419  1.279959e-03       9
#> NANOS14          4.421919e-07 0.742757572 0.741 0.259  1.326576e-03       9
#> CKAP2L2          4.502001e-07 0.498174133 0.259 0.043  1.350600e-03       9
#> RRS13            4.965720e-07 0.998719915 0.778 0.365  1.489716e-03       9
#> DGAT23           5.115685e-07 0.851409412 0.778 0.317  1.534706e-03       9
#> EBP3             5.503563e-07 0.856059531 0.926 0.445  1.651069e-03       9
#> PPA15            5.613702e-07 0.959986999 0.852 0.493  1.684111e-03       9
#> SLC9A3R24        5.766926e-07 0.956498657 0.926 0.441  1.730078e-03       9
#> SAC3D13          6.216097e-07 0.942859901 0.704 0.271  1.864829e-03       9
#> MT1G3            6.298352e-07 0.248521080 0.556 0.229  1.889505e-03       9
#> UCHL13           6.402543e-07 0.905458073 0.519 0.163  1.920763e-03       9
#> DCXR3            6.462340e-07 1.029123116 0.852 0.441  1.938702e-03       9
#> GOLM14           7.123303e-07 0.963823503 0.889 0.479  2.136991e-03       9
#> ACTL6A4          7.236460e-07 1.034237938 0.889 0.483  2.170938e-03       9
#> ROPN1B3          7.284008e-07 0.949909445 0.926 0.435  2.185203e-03       9
#> HSPH13           7.984788e-07 0.866200180 0.963 0.489  2.395436e-03       9
#> TUBB2B4          8.035236e-07 0.878534364 0.852 0.332  2.410571e-03       9
#> UBE2S3           8.617558e-07 0.809889699 0.741 0.402  2.585267e-03       9
#> TUBB2A4          8.642803e-07 0.892194583 0.889 0.421  2.592841e-03       9
#> MT1E2            9.062413e-07 0.914001663 0.778 0.388  2.718724e-03       9
#> TUBA4A3          9.231212e-07 1.008054906 0.741 0.349  2.769364e-03       9
#> PPP1R12A3        1.088425e-06 0.914853546 0.852 0.452  3.265275e-03       9
#> ZNF6952          1.117252e-06 0.755208266 0.407 0.104  3.351755e-03       9
#> CLDN111          1.124566e-06 1.048390503 0.259 0.048  3.373699e-03       9
#> VANGL15          1.145680e-06 0.868686319 0.963 0.526  3.437040e-03       9
#> GAPDH3           1.342012e-06 0.745105173 1.000 0.499  4.026035e-03       9
#> ANO13            1.376200e-06 0.549996634 0.519 0.149  4.128600e-03       9
#> CA21             1.522394e-06 0.591285598 0.667 0.236  4.567183e-03       9
#> PARVB3           1.763219e-06 0.850137811 0.852 0.378  5.289656e-03       9
#> ROPN13           1.817606e-06 0.704814434 0.667 0.325  5.452819e-03       9
#> ENPP54           1.851093e-06 0.758971974 0.815 0.340  5.553278e-03       9
#> LGR61            1.982437e-06 0.440296510 0.259 0.047  5.947312e-03       9
#> CCDC1672         2.018426e-06 0.858621118 0.926 0.496  6.055277e-03       9
#> LINC010481       2.024886e-06 0.511104544 0.370 0.086  6.074658e-03       9
#> GAL3             2.037556e-06 0.641428180 0.630 0.275  6.112669e-03       9
#> LY6E5            2.166930e-06 0.898289546 0.852 0.468  6.500789e-03       9
#> MAPK8IP22        2.191068e-06 0.650519229 0.630 0.213  6.573205e-03       9
#> S100B3           2.209759e-06 0.906568758 0.852 0.412  6.629278e-03       9
#> DBI3             2.442200e-06 0.861413747 0.926 0.499  7.326600e-03       9
#> PDLIM12          2.546253e-06 0.638052200 0.704 0.333  7.638758e-03       9
#> ASPM3            2.569485e-06 0.165208341 0.370 0.084  7.708454e-03       9
#> CSRP14           2.934145e-06 0.835897757 0.889 0.450  8.802436e-03       9
#> PFN23            3.189078e-06 0.878889382 0.889 0.503  9.567234e-03       9
#> HIST1H2BJ3       3.239570e-06 0.824746168 0.593 0.204  9.718709e-03       9
#> NSG13            3.379191e-06 0.733674683 0.556 0.182  1.013757e-02       9
#> TPM25            3.426497e-06 0.704823249 0.926 0.425  1.027949e-02       9
#> PBX14            3.686173e-06 0.914098701 0.889 0.462  1.105852e-02       9
#> LSM53            3.693425e-06 0.836304372 0.852 0.471  1.108027e-02       9
#> HILPDA3          3.932103e-06 0.647363399 0.815 0.395  1.179631e-02       9
#> NOTCH4           3.999264e-06 0.303203353 0.370 0.088  1.199779e-02       9
#> CKLF2            4.328819e-06 0.577149177 0.963 0.421  1.298646e-02       9
#> XAGE23           4.472415e-06 0.732130259 0.630 0.303  1.341725e-02       9
#> ACTA25           4.542113e-06 0.874885463 0.852 0.350  1.362634e-02       9
#> CTNND23          4.579262e-06 0.779713302 0.630 0.251  1.373779e-02       9
#> FXYD63           4.954669e-06 0.820427692 0.963 0.486  1.486401e-02       9
#> KLK112           5.038386e-06 0.934973475 0.370 0.102  1.511516e-02       9
#> SOHLH13          5.271372e-06 0.576094962 0.704 0.317  1.581412e-02       9
#> ELOVL53          6.438825e-06 0.688655521 0.704 0.268  1.931648e-02       9
#> CCND12           6.660923e-06 0.798352445 0.889 0.468  1.998277e-02       9
#> ATP6V0E23        7.050958e-06 1.026177211 0.667 0.314  2.115288e-02       9
#> S100A13          7.129762e-06 0.959512931 1.000 0.503  2.138929e-02       9
#> TUBB4B5          7.169308e-06 0.885634142 0.889 0.507  2.150792e-02       9
#> GPSM24           7.680341e-06 0.524053257 0.556 0.183  2.304102e-02       9
#> PCOLCE24         8.646961e-06 0.577561609 0.778 0.326  2.594088e-02       9
#> KIF112           8.657786e-06 0.663841664 0.296 0.068  2.597336e-02       9
#> KLHL353          9.384778e-06 0.684408625 0.852 0.366  2.815433e-02       9
#> ARL6IP14         9.552032e-06 0.692880501 1.000 0.475  2.865610e-02       9
#> PDIA44           1.022140e-05 0.752264947 0.889 0.502  3.066421e-02       9
#> FBXL222          1.035989e-05 0.480037249 0.556 0.188  3.107966e-02       9
#> TUBA1C5          1.049696e-05 0.884678814 0.815 0.464  3.149089e-02       9
#> PTS4             1.098245e-05 0.813828310 0.889 0.476  3.294734e-02       9
#> PRR112           1.103879e-05 0.364428257 0.333 0.079  3.311638e-02       9
#> MIA3             1.151974e-05 0.736197440 0.889 0.392  3.455923e-02       9
#> KDELC23          1.255589e-05 0.433893587 0.556 0.182  3.766767e-02       9
#> FKBP111          1.262682e-05 0.330911890 0.667 0.234  3.788046e-02       9
#> NEK23            1.320192e-05 0.485552447 0.370 0.101  3.960575e-02       9
#> PPIF3            1.393386e-05 0.783588966 0.778 0.401  4.180157e-02       9
#> CCT6A3           1.411809e-05 0.748574541 0.926 0.510  4.235426e-02       9
#> MTHFD24          1.485137e-05 0.847374872 0.852 0.410  4.455411e-02       9
#> COL11A23         1.574388e-05 0.740366021 0.556 0.221  4.723165e-02       9
#> CENPE2           1.792127e-05 0.194263432 0.333 0.079  5.376380e-02       9
#> COPZ22           1.884602e-05 0.474082533 0.630 0.226  5.653805e-02       9
#> KRT87            1.885672e-05 0.767790092 0.963 0.546  5.657017e-02       9
#> ST144            1.940557e-05 0.816844383 0.778 0.490  5.821670e-02       9
#> EXTL12           1.942723e-05 0.652354715 0.593 0.228  5.828169e-02       9
#> DKK13            2.041794e-05 0.597310957 0.444 0.152  6.125382e-02       9
#> DMKN1            2.055991e-05 0.425520952 0.556 0.191  6.167973e-02       9
#> IDH14            2.062480e-05 0.759865238 0.852 0.405  6.187439e-02       9
#> UCP23            2.322053e-05 0.548640012 0.778 0.350  6.966159e-02       9
#> COL4A24          2.346856e-05 0.312819899 0.667 0.339  7.040569e-02       9
#> PGAM12           2.348424e-05 0.723740054 0.815 0.471  7.045273e-02       9
#> TEKT34           2.419966e-05 0.819247632 0.889 0.507  7.259897e-02       9
#> MYC2             2.541067e-05 0.847307707 0.667 0.315  7.623201e-02       9
#> MYL94            2.588717e-05 0.812454592 0.926 0.487  7.766151e-02       9
#> AARD5            2.610869e-05 0.720341233 0.852 0.444  7.832606e-02       9
#> KCNN42           2.613182e-05 0.578093620 0.630 0.294  7.839547e-02       9
#> FABP72           2.629326e-05 0.273693329 0.259 0.056  7.887978e-02       9
#> SNHG253          2.721485e-05 0.657888114 0.815 0.371  8.164454e-02       9
#> CITED44          2.857225e-05 0.793489842 0.926 0.477  8.571675e-02       9
#> RIPPLY3          2.984146e-05 0.622769526 0.296 0.073  8.952439e-02       9
#> PTPRS2           3.211881e-05 0.643002666 0.556 0.206  9.635644e-02       9
#> LBR3             3.490180e-05 0.572883531 0.778 0.351  1.047054e-01       9
#> TRIB33           3.668741e-05 0.627956536 0.741 0.307  1.100622e-01       9
#> PALLD4           3.731050e-05 0.514490630 0.852 0.363  1.119315e-01       9
#> HES62            3.756304e-05 0.558699550 0.630 0.237  1.126891e-01       9
#> SERPINE23        4.263608e-05 0.766810872 0.704 0.397  1.279082e-01       9
#> TSPAN23          4.276235e-05 0.521998079 0.519 0.179  1.282870e-01       9
#> PPP1R14A3        4.859561e-05 0.237621624 0.481 0.198  1.457868e-01       9
#> HOXB9            5.166546e-05 0.367398964 0.259 0.058  1.549964e-01       9
#> PLOD33           5.458263e-05 0.767672959 0.778 0.422  1.637479e-01       9
#> IGFBP74          5.866825e-05 0.636576865 0.815 0.410  1.760048e-01       9
#> SFRP13           6.008739e-05 0.764401134 0.852 0.464  1.802622e-01       9
#> PPP1R1B3         6.045208e-05 0.759050809 0.778 0.383  1.813562e-01       9
#> NDUFAF63         6.307920e-05 0.788752797 0.778 0.450  1.892376e-01       9
#> GAS12            6.362020e-05 0.651004769 0.704 0.330  1.908606e-01       9
#> LAMB32           6.431579e-05 0.328591520 0.519 0.246  1.929474e-01       9
#> LRP24            6.528167e-05 0.477334258 0.630 0.242  1.958450e-01       9
#> SH3BGR4          6.743539e-05 0.482616068 0.704 0.339  2.023062e-01       9
#> DLGAP52          7.487120e-05 0.245512513 0.296 0.074  2.246136e-01       9
#> LIMCH13          7.641169e-05 0.690813924 0.778 0.397  2.292351e-01       9
#> CTNNB12          7.738574e-05 0.898528761 0.704 0.482  2.321572e-01       9
#> CDKN2D2          7.952598e-05 0.856765982 0.556 0.253  2.385780e-01       9
#> NQO14            9.062579e-05 0.695224784 0.741 0.436  2.718774e-01       9
#> TNFRSF213        9.077787e-05 0.667120876 0.815 0.418  2.723336e-01       9
#> UACA3            9.271790e-05 0.336129821 0.815 0.348  2.781537e-01       9
#> COL11A15         9.532456e-05 0.461886718 0.741 0.379  2.859737e-01       9
#> TTC39A4          1.057083e-04 0.751961428 0.593 0.264  3.171248e-01       9
#> C1QL44           1.081763e-04 0.315214359 0.519 0.285  3.245289e-01       9
#> LTBP12           1.092283e-04 0.365090701 0.519 0.181  3.276850e-01       9
#> RPP251           1.097444e-04 0.783171212 0.741 0.372  3.292332e-01       9
#> QDPR2            1.130683e-04 0.619915916 0.778 0.436  3.392049e-01       9
#> ACTA12           1.302447e-04 0.500866562 0.407 0.138  3.907342e-01       9
#> BARX14           1.307756e-04 0.731462378 0.778 0.398  3.923268e-01       9
#> DNM33            1.318507e-04 0.484129051 0.593 0.270  3.955520e-01       9
#> MT2A5            1.381658e-04 0.636389355 0.741 0.404  4.144974e-01       9
#> RP11-294J22.6    1.512003e-04 0.753817810 0.519 0.220  4.536010e-01       9
#> IMPA24           1.532575e-04 0.678146688 0.963 0.500  4.597726e-01       9
#> FREM21           1.672080e-04 0.517163967 0.296 0.082  5.016241e-01       9
#> SNHG193          1.688808e-04 0.602635749 0.778 0.412  5.066423e-01       9
#> TINAGL12         1.768283e-04 0.403647141 0.444 0.149  5.304849e-01       9
#> DSC34            1.953640e-04 0.607091195 0.630 0.321  5.860921e-01       9
#> PEG103           2.043300e-04 0.685891845 0.852 0.480  6.129901e-01       9
#> NTF31            2.051708e-04 0.442070991 0.296 0.082  6.155125e-01       9
#> DACT31           2.504959e-04 0.325314511 0.259 0.066  7.514878e-01       9
#> FAM3C4           2.620669e-04 0.648879618 0.815 0.433  7.862008e-01       9
#> RAB254           2.721734e-04 0.663409163 0.704 0.393  8.165202e-01       9
#> CPED13           3.067116e-04 0.469719740 0.593 0.279  9.201349e-01       9
#> BAMBI3           3.383234e-04 0.628493891 0.852 0.457  1.000000e+00       9
#> NUDT43           3.388201e-04 0.516208347 0.704 0.387  1.000000e+00       9
#> SELENBP14        3.434407e-04 0.408949594 0.667 0.344  1.000000e+00       9
#> CENPV2           3.437783e-04 0.397671058 0.444 0.164  1.000000e+00       9
#> FHL14            3.449273e-04 0.127094945 0.481 0.198  1.000000e+00       9
#> IL17B4           3.516740e-04 0.259646793 0.444 0.184  1.000000e+00       9
#> MTL54            3.545282e-04 0.463398257 0.778 0.400  1.000000e+00       9
#> KIF20A1          3.571794e-04 0.062325251 0.259 0.065  1.000000e+00       9
#> SCARB13          3.622858e-04 0.649734215 0.778 0.459  1.000000e+00       9
#> FABP53           3.875819e-04 0.728450540 0.704 0.433  1.000000e+00       9
#> KRTCAP33         3.968300e-04 0.621898128 0.741 0.478  1.000000e+00       9
#> RGCC2            3.968306e-04 0.609688951 0.667 0.390  1.000000e+00       9
#> B3GNT73          3.996431e-04 0.497285467 0.667 0.338  1.000000e+00       9
#> CCDC28B1         4.436445e-04 0.411412917 0.519 0.204  1.000000e+00       9
#> N4BP2            4.458815e-04 0.431356486 0.444 0.164  1.000000e+00       9
#> SERTAD4-AS13     4.494826e-04 0.550066316 0.741 0.439  1.000000e+00       9
#> STRADB           4.542464e-04 0.208388951 0.519 0.195  1.000000e+00       9
#> MYO1B3           4.764955e-04 0.507132219 0.741 0.407  1.000000e+00       9
#> C2orf801         4.821453e-04 0.435957285 0.333 0.106  1.000000e+00       9
#> OBP2B2           4.859231e-04 0.408172322 0.370 0.128  1.000000e+00       9
#> MTSS11           4.909733e-04 0.585338940 0.556 0.301  1.000000e+00       9
#> LOXL13           4.918150e-04 0.264397384 0.519 0.208  1.000000e+00       9
#> PHF192           4.989982e-04 0.484484695 0.630 0.325  1.000000e+00       9
#> RP11-273G15.21   5.033838e-04 0.435550276 0.296 0.092  1.000000e+00       9
#> NEGR13           5.341213e-04 0.440961779 0.519 0.207  1.000000e+00       9
#> MYLK3            5.484035e-04 0.585235484 0.852 0.435  1.000000e+00       9
#> SLC29A13         5.720734e-04 0.720942093 0.963 0.507  1.000000e+00       9
#> GCHFR1           6.023476e-04 0.201896561 0.593 0.234  1.000000e+00       9
#> ITGA63           6.199220e-04 0.455773361 0.556 0.293  1.000000e+00       9
#> SMOC12           6.489097e-04 0.455327714 0.630 0.282  1.000000e+00       9
#> LIPH2            6.677814e-04 0.341793439 0.630 0.263  1.000000e+00       9
#> CDC42EP14        6.730262e-04 0.634371584 0.778 0.451  1.000000e+00       9
#> SERPINH15        6.839581e-04 0.644860124 0.852 0.534  1.000000e+00       9
#> IDH1-AS11        7.064106e-04 0.239489343 0.370 0.121  1.000000e+00       9
#> ETV53            7.081175e-04 0.398138155 0.444 0.187  1.000000e+00       9
#> NREP2            7.454207e-04 0.483774905 0.556 0.257  1.000000e+00       9
#> SAP303           7.969470e-04 0.601062154 0.778 0.436  1.000000e+00       9
#> FKBP105          8.015859e-04 0.546116958 0.741 0.448  1.000000e+00       9
#> GMPR3            8.191847e-04 0.526700615 0.593 0.268  1.000000e+00       9
#> CCDC85B4         8.235513e-04 0.502819737 0.741 0.415  1.000000e+00       9
#> TNFRSF182        8.536403e-04 0.279905584 0.370 0.129  1.000000e+00       9
#> ZG16B4           8.614297e-04 0.445404377 0.704 0.361  1.000000e+00       9
#> ANGPT11          8.721577e-04 0.443030182 0.444 0.178  1.000000e+00       9
#> TSEN152          8.778175e-04 0.595478390 0.778 0.478  1.000000e+00       9
#> HSPA61           9.164056e-04 0.161389665 0.519 0.220  1.000000e+00       9
#> COL4A11          9.652490e-04 0.073952270 0.333 0.198  1.000000e+00       9
#> USP181           1.041191e-03 0.332289833 0.444 0.176  1.000000e+00       9
#> GLS4             1.047775e-03 0.616282523 0.963 0.454  1.000000e+00       9
#> WEE13            1.085703e-03 0.484771200 0.630 0.358  1.000000e+00       9
#> TPD52L14         1.101173e-03 0.600595617 0.741 0.470  1.000000e+00       9
#> SDC24            1.173963e-03 0.639298730 0.741 0.471  1.000000e+00       9
#> FASN2            1.187344e-03 0.325800985 0.630 0.270  1.000000e+00       9
#> HSPA1A4          1.201700e-03 0.651008536 0.852 0.498  1.000000e+00       9
#> FRMD4A3          1.220072e-03 0.538410978 0.667 0.364  1.000000e+00       9
#> TNFSF13B3        1.302626e-03 0.598863165 0.704 0.442  1.000000e+00       9
#> ARRB22           1.358684e-03 0.099078356 0.593 0.274  1.000000e+00       9
#> HES43            1.410475e-03 0.359001450 0.667 0.401  1.000000e+00       9
#> NAB12            1.566243e-03 0.466371840 0.519 0.277  1.000000e+00       9
#> QPCT3            1.598091e-03 0.605323754 0.667 0.358  1.000000e+00       9
#> MESP22           1.698443e-03 0.235773892 0.667 0.281  1.000000e+00       9
#> TMEM794          1.719494e-03 0.350947088 0.778 0.389  1.000000e+00       9
#> TMPRSS132        1.723131e-03 0.406022687 0.593 0.279  1.000000e+00       9
#> NPPC2            1.784847e-03 0.483124303 0.296 0.103  1.000000e+00       9
#> IGF2             1.840848e-03 0.540445149 0.259 0.085  1.000000e+00       9
#> GJA13            1.863452e-03 0.328098056 0.444 0.182  1.000000e+00       9
#> FMOD1            1.872121e-03 0.524813993 0.259 0.085  1.000000e+00       9
#> EPHX14           1.975373e-03 0.560908176 0.741 0.464  1.000000e+00       9
#> PIR3             1.984637e-03 0.397010972 0.630 0.341  1.000000e+00       9
#> CEP704           1.997912e-03 0.286229463 0.556 0.248  1.000000e+00       9
#> SNX223           2.024363e-03 0.264429114 0.519 0.316  1.000000e+00       9
#> BSPRY3           2.112183e-03 0.497419249 0.593 0.282  1.000000e+00       9
#> NFIB4            2.158434e-03 0.557006946 0.778 0.453  1.000000e+00       9
#> PIM13            2.215346e-03 0.568284706 0.778 0.447  1.000000e+00       9
#> TAGLN4           2.428437e-03 0.535545746 0.852 0.437  1.000000e+00       9
#> DBNDD13          2.599711e-03 0.472168449 0.667 0.348  1.000000e+00       9
#> APOLD11          2.633557e-03 0.215924505 0.444 0.182  1.000000e+00       9
#> MDFI3            2.713371e-03 0.511025377 0.889 0.450  1.000000e+00       9
#> C11orf731        2.784745e-03 0.521338922 0.704 0.521  1.000000e+00       9
#> SMYD23           3.001297e-03 0.559606968 0.741 0.448  1.000000e+00       9
#> ASNS3            3.106560e-03 0.332024745 0.519 0.241  1.000000e+00       9
#> SOX85            3.115160e-03 0.380008425 0.667 0.399  1.000000e+00       9
#> KLRG24           3.322763e-03 0.417161859 0.444 0.207  1.000000e+00       9
#> SMOC23           3.422037e-03 0.117100334 0.481 0.223  1.000000e+00       9
#> PRSS85           3.490584e-03 0.554060605 0.889 0.513  1.000000e+00       9
#> INSR             3.518333e-03 0.301202910 0.481 0.215  1.000000e+00       9
#> CHORDC12         3.546720e-03 0.491174807 0.815 0.481  1.000000e+00       9
#> CLU3             3.655542e-03 0.552172181 0.889 0.462  1.000000e+00       9
#> INHBB2           3.766427e-03 0.383756590 0.370 0.155  1.000000e+00       9
#> DUSP22           3.795634e-03 0.020659997 0.444 0.176  1.000000e+00       9
#> UPP12            3.939376e-03 0.236562528 0.667 0.316  1.000000e+00       9
#> CELF43           3.968132e-03 0.331998267 0.630 0.340  1.000000e+00       9
#> QPRT3            3.969646e-03 0.278795700 0.556 0.279  1.000000e+00       9
#> BARX2            4.045772e-03 0.453989764 0.296 0.112  1.000000e+00       9
#> PPDPF3           4.484528e-03 0.525992060 0.852 0.571  1.000000e+00       9
#> FRMD32           4.648954e-03 0.298035937 0.481 0.220  1.000000e+00       9
#> REEP42           4.670534e-03 0.180238743 0.519 0.296  1.000000e+00       9
#> TMEM2042         4.723356e-03 0.295400562 0.370 0.154  1.000000e+00       9
#> CRNDE4           4.868891e-03 0.541292382 0.778 0.506  1.000000e+00       9
#> LAD12            4.928062e-03 0.387816196 0.630 0.429  1.000000e+00       9
#> CRELD22          5.078094e-03 0.562474110 0.704 0.497  1.000000e+00       9
#> SLC25A53         5.189080e-03 0.490098225 0.778 0.489  1.000000e+00       9
#> CTC-425F1.4      5.221191e-03 0.377830029 0.259 0.094  1.000000e+00       9
#> C1QBP3           5.335898e-03 0.510068057 0.667 0.507  1.000000e+00       9
#> SCGB3A13         5.556365e-03 0.257845660 0.556 0.344  1.000000e+00       9
#> HSPA1B4          5.666927e-03 0.483998510 0.741 0.464  1.000000e+00       9
#> PMAIP12          5.804867e-03 0.289549811 0.519 0.341  1.000000e+00       9
#> PRNP3            5.807904e-03 0.449293555 0.815 0.521  1.000000e+00       9
#> TPM15            5.880306e-03 0.595231055 0.741 0.558  1.000000e+00       9
#> GRB102           6.234249e-03 0.200713456 0.556 0.271  1.000000e+00       9
#> DTNB3            6.348253e-03 0.310279785 0.556 0.363  1.000000e+00       9
#> FADS22           6.447687e-03 0.220036390 0.556 0.277  1.000000e+00       9
#> NAAA2            6.447726e-03 0.238399569 0.519 0.254  1.000000e+00       9
#> SQLE4            6.530820e-03 0.432779858 0.667 0.442  1.000000e+00       9
#> FAM105A2         6.536220e-03 0.283134724 0.556 0.299  1.000000e+00       9
#> MATN33           6.597889e-03 0.241847010 0.481 0.320  1.000000e+00       9
#> ELN3             6.829039e-03 0.324206302 0.407 0.192  1.000000e+00       9
#> ATP1B13          6.977766e-03 0.507179447 0.778 0.458  1.000000e+00       9
#> MARVELD12        7.359403e-03 0.281523198 0.630 0.322  1.000000e+00       9
#> LMTK31           7.435479e-03 0.348077113 0.259 0.097  1.000000e+00       9
#> PTP4A33          7.544664e-03 0.410166956 0.630 0.481  1.000000e+00       9
#> TUBA1A4          7.778580e-03 0.468503300 0.704 0.484  1.000000e+00       9
#> EPPK12           7.947295e-03 0.139445221 0.481 0.214  1.000000e+00       9
#> ERBB34           7.985922e-03 0.346814473 0.778 0.444  1.000000e+00       9
#> ITIH61           8.283235e-03 0.313420595 0.296 0.114  1.000000e+00       9
#> COL9A33          8.431743e-03 0.336106802 0.741 0.446  1.000000e+00       9
#> MSRB34           8.702711e-03 0.291456084 0.519 0.247  1.000000e+00       9
#> BACE24           9.097363e-03 0.395017546 0.704 0.467  1.000000e+00       9
#> TUBB64           9.381150e-03 0.458329755 0.704 0.471  1.000000e+00       9
#> CTGF3            9.386277e-03 0.312657568 0.519 0.375  1.000000e+00       9
#> MAL24            9.722727e-03 0.403756791 0.741 0.493  1.000000e+00       9
#> ATP6V0D2        7.748189e-119 5.681203083 0.810 0.010 2.324457e-115      10
#> RUFY4           4.422049e-117 5.317588306 0.714 0.007 1.326615e-113      10
#> SIGLEC15        2.321180e-106 6.131331569 0.952 0.022 6.963540e-103      10
#> STRIP2          4.587842e-103 4.320152506 0.571 0.004  1.376353e-99      10
#> DPP4            3.257920e-101 5.366486087 0.810 0.015  9.773761e-98      10
#> SLC9B2           2.241746e-89 6.456338908 1.000 0.035  6.725238e-86      10
#> PAGE2            4.110456e-69 2.472692785 0.286 0.000  1.233137e-65      10
#> AKAP6            2.072038e-68 4.124261557 0.571 0.011  6.216113e-65      10
#> AK5              1.258559e-66 5.421229038 0.762 0.027  3.775676e-63      10
#> SUCNR1           3.848134e-66 4.309167379 0.810 0.031  1.154440e-62      10
#> SH3GL2           1.441955e-62 2.881685650 0.333 0.002  4.325865e-59      10
#> F5               1.667475e-48 3.303440194 0.429 0.009  5.002424e-45      10
#> SLC37A2          7.140009e-45 4.193933397 0.857 0.062  2.142003e-41      10
#> ITGB3            7.393785e-45 3.913469484 0.524 0.019  2.218136e-41      10
#> DGKI             2.449118e-40 3.728609837 0.810 0.062  7.347354e-37      10
#> TM4SF19          3.174119e-40 3.137789835 0.333 0.007  9.522356e-37      10
#> MATK             5.559526e-39 3.859890723 0.810 0.066  1.667858e-35      10
#> HTRA31           1.599696e-38 3.545050124 0.857 0.076  4.799089e-35      10
#> RP11-702B10.2    2.416775e-37 2.702098227 0.333 0.007  7.250326e-34      10
#> ANGPTL6          6.560098e-37 3.220949193 0.429 0.015  1.968029e-33      10
#> CD1093           1.292615e-35 3.346606293 0.952 0.103  3.877846e-32      10
#> SCD51            2.073112e-35 3.913856440 0.714 0.055  6.219337e-32      10
#> ECSCR.1          2.336864e-34 1.742458083 0.286 0.006  7.010592e-31      10
#> PTPN222          8.867736e-33 3.712413514 0.762 0.070  2.660321e-29      10
#> LAT              5.693420e-32 3.482101065 0.857 0.095  1.708026e-28      10
#> MAP1A1           9.310605e-32 2.491951545 0.524 0.029  2.793182e-28      10
#> RAC22            1.596322e-30 4.105109661 0.905 0.123  4.788965e-27      10
#> GNG21            8.072777e-30 3.087414648 0.714 0.066  2.421833e-26      10
#> RFTN2            1.452505e-29 2.441913901 0.286 0.007  4.357516e-26      10
#> NRIP3            8.055118e-29 4.091171435 0.667 0.061  2.416536e-25      10
#> RAB38            2.173269e-27 2.754301662 0.476 0.029  6.519807e-24      10
#> CTSK2            6.755835e-25 4.467459201 0.952 0.164  2.026751e-21      10
#> HPGD             1.524828e-23 2.930244757 0.333 0.016  4.574485e-20      10
#> ITGA2            1.889590e-23 3.317777525 0.571 0.055  5.668769e-20      10
#> CD842            2.980050e-23 3.060606084 0.952 0.179  8.940150e-20      10
#> STRADB1          8.503142e-23 4.536609845 0.905 0.190  2.550943e-19      10
#> HSD3B7           8.675411e-23 3.884804331 0.810 0.132  2.602623e-19      10
#> SCARF11          9.204904e-23 2.230983715 0.571 0.053  2.761471e-19      10
#> PARM11           6.306276e-22 3.001272921 0.619 0.070  1.891883e-18      10
#> THOP14           1.049946e-21 5.052604497 0.905 0.205  3.149838e-18      10
#> TNFRSF11A        2.867584e-21 2.580622710 0.619 0.071  8.602753e-18      10
#> DUSP4            3.399264e-21 2.658528610 0.476 0.040  1.019779e-17      10
#> CCR12            8.262580e-21 2.191141894 0.810 0.125  2.478774e-17      10
#> GAPLINC1         1.050923e-20 1.799700859 0.667 0.077  3.152769e-17      10
#> CFI1             4.828767e-19 2.371265183 0.381 0.028  1.448630e-15      10
#> ACP52            5.059605e-19 4.461856857 0.952 0.230  1.517881e-15      10
#> FJX1             1.356744e-18 2.383556950 0.667 0.099  4.070231e-15      10
#> SELPLG2          3.085654e-18 1.961230380 0.667 0.093  9.256963e-15      10
#> OSCAR2           4.135239e-18 2.397513174 0.571 0.072  1.240572e-14      10
#> SPINK2           4.199733e-18 2.003888695 0.286 0.016  1.259920e-14      10
#> IGFBP61          4.236161e-18 2.363623476 0.524 0.059  1.270848e-14      10
#> SNX102           4.644536e-18 2.386610360 0.952 0.227  1.393361e-14      10
#> SPN              6.718800e-18 2.562447930 0.429 0.039  2.015640e-14      10
#> ANPEP3           1.568349e-17 3.017726812 0.952 0.281  4.705047e-14      10
#> F2R1             3.431543e-17 1.602550690 0.714 0.117  1.029463e-13      10
#> EVI2A2           7.938223e-17 2.128074160 0.810 0.156  2.381467e-13      10
#> CYTH42           1.068251e-16 1.971274213 0.762 0.131  3.204754e-13      10
#> CLEC4A2          2.106537e-16 1.656189752 0.762 0.131  6.319610e-13      10
#> GNPTAB1          2.982788e-16 3.628951264 0.952 0.336  8.948365e-13      10
#> EID13            4.278025e-16 1.903955737 1.000 0.254  1.283407e-12      10
#> MMP93            5.152383e-16 3.101592637 0.952 0.220  1.545715e-12      10
#> CD332            7.398546e-16 2.027566239 0.571 0.080  2.219564e-12      10
#> FAM43A           9.626192e-16 1.851208431 0.286 0.019  2.887858e-12      10
#> TFEC2            1.520863e-15 1.532177114 0.714 0.120  4.562588e-12      10
#> CSF1R2           2.251389e-15 2.161374112 0.857 0.208  6.754167e-12      10
#> JDP2             2.521322e-15 3.622313050 0.714 0.160  7.563965e-12      10
#> EVI2B2           2.825147e-15 1.973611563 0.810 0.173  8.475441e-12      10
#> TCIRG12          3.107965e-15 4.072135424 0.952 0.364  9.323896e-12      10
#> COL5A11          3.491587e-15 1.112364210 0.524 0.064  1.047476e-11      10
#> VASH12           1.459403e-14 1.567951239 0.619 0.098  4.378208e-11      10
#> CAV11            2.906808e-14 0.695346691 0.476 0.055  8.720425e-11      10
#> ID33             2.966821e-14 1.644313638 0.762 0.153  8.900464e-11      10
#> MITF2            3.121435e-14 2.106087992 0.810 0.204  9.364306e-11      10
#> CKLF3            3.491657e-14 2.855044511 1.000 0.424  1.047497e-10      10
#> SPP12            6.206627e-14 2.430969874 1.000 0.185  1.861988e-10      10
#> FAM134B          6.488208e-14 2.875618671 0.571 0.100  1.946462e-10      10
#> CD682            9.349320e-14 2.089158683 0.952 0.242  2.804796e-10      10
#> KIAA00401        1.764537e-13 2.593612231 0.857 0.278  5.293612e-10      10
#> RGS104           4.779377e-13 3.246720729 0.952 0.448  1.433813e-09      10
#> EDIL31           6.580523e-13 1.136820213 0.381 0.041  1.974157e-09      10
#> F13A11           7.414418e-13 1.242673413 0.429 0.052  2.224325e-09      10
#> RNASE62          1.357218e-12 1.781847188 0.714 0.161  4.071655e-09      10
#> CST33            1.668191e-12 2.684635033 0.952 0.433  5.004572e-09      10
#> CLEC11A1         2.363229e-12 2.058958193 0.905 0.292  7.089688e-09      10
#> PRDM13           2.863706e-12 1.454290126 0.714 0.152  8.591118e-09      10
#> CACNA2D42        3.692577e-12 1.472162676 0.571 0.100  1.107773e-08      10
#> IFI302           5.027368e-12 1.452564590 0.762 0.178  1.508210e-08      10
#> SPI12            5.080005e-12 1.373672794 0.952 0.217  1.524002e-08      10
#> SEPT63           6.143766e-12 1.419608864 0.714 0.161  1.843130e-08      10
#> ERO1A2           9.528137e-12 1.462230549 0.810 0.211  2.858441e-08      10
#> GPX13            1.686945e-11 1.617253431 1.000 0.484  5.060834e-08      10
#> WAS2             2.734915e-11 1.250451443 0.667 0.139  8.204744e-08      10
#> SPARC2           3.000965e-11 1.103536567 0.857 0.221  9.002896e-08      10
#> LAPTM52          3.481841e-11 1.594954525 1.000 0.235  1.044552e-07      10
#> EHD21            3.702479e-11 1.716768761 0.524 0.097  1.110744e-07      10
#> UBB3             6.428492e-11 1.562709078 1.000 0.283  1.928548e-07      10
#> PSTPIP12         6.551890e-11 1.851809016 0.619 0.143  1.965567e-07      10
#> FERMT32          7.641208e-11 1.227892221 0.810 0.201  2.292363e-07      10
#> TGFB13           9.404430e-11 1.669360301 0.762 0.204  2.821329e-07      10
#> MVP3             1.039866e-10 1.794585322 0.857 0.291  3.119597e-07      10
#> LY963            1.157458e-10 1.154161023 0.952 0.250  3.472374e-07      10
#> DOCK102          1.610776e-10 1.038050582 0.571 0.102  4.832328e-07      10
#> APOBEC3C2        1.758218e-10 1.004526736 0.571 0.105  5.274653e-07      10
#> LSP12            4.439082e-10 1.575909742 0.857 0.244  1.331725e-06      10
#> LIMS12           5.837627e-10 1.582682751 0.952 0.354  1.751288e-06      10
#> TSC22D33         6.642522e-10 1.711200165 0.762 0.234  1.992756e-06      10
#> RCAN12           9.163225e-10 2.045062921 0.857 0.347  2.748967e-06      10
#> GPR1832          1.035457e-09 1.290193202 0.857 0.262  3.106371e-06      10
#> XIST3            1.266385e-09 1.066670128 0.905 0.249  3.799154e-06      10
#> EMP33            1.273704e-09 1.534734660 0.905 0.285  3.821111e-06      10
#> PTPRC2           1.580753e-09 1.207187842 0.810 0.218  4.742258e-06      10
#> FILIP1L3         1.591411e-09 1.078525148 0.571 0.121  4.774234e-06      10
#> CLECL12          1.645950e-09 1.118249117 0.476 0.086  4.937849e-06      10
#> LCP12            2.614131e-09 1.196144219 0.762 0.203  7.842393e-06      10
#> CD300C1          3.156137e-09 1.054986649 0.381 0.057  9.468412e-06      10
#> CYTIP2           3.460945e-09 1.219585238 0.571 0.127  1.038284e-05      10
#> PGAM13           4.108446e-09 1.920797291 0.905 0.471  1.232534e-05      10
#> ITGB22           4.684899e-09 1.362448397 0.905 0.257  1.405470e-05      10
#> CSTB3            4.960508e-09 1.855176976 0.905 0.481  1.488152e-05      10
#> N4BP21           5.534183e-09 2.640611628 0.571 0.163  1.660255e-05      10
#> CD300A2          6.771747e-09 1.063302180 0.476 0.089  2.031524e-05      10
#> IGSF212          7.690696e-09 1.307360579 0.524 0.115  2.307209e-05      10
#> BST23            8.262667e-09 0.938014295 0.905 0.256  2.478800e-05      10
#> OAF2             1.131190e-08 1.655533694 0.667 0.212  3.393569e-05      10
#> FAM162A1         1.234732e-08 1.475739116 0.952 0.455  3.704197e-05      10
#> SYCE1L2          1.324002e-08 1.699757654 0.714 0.258  3.972006e-05      10
#> TBXAS12          2.500067e-08 1.273866393 0.619 0.173  7.500200e-05      10
#> NFATC12          4.925291e-08 1.915402885 0.810 0.392  1.477587e-04      10
#> FCER1G2          5.106664e-08 1.322083056 0.905 0.245  1.531999e-04      10
#> TYROBP2          5.226273e-08 1.534136572 0.952 0.229  1.567882e-04      10
#> MAP2K32          5.826572e-08 1.330987660 0.905 0.387  1.747972e-04      10
#> HLX1             6.770453e-08 1.158915271 0.333 0.054  2.031136e-04      10
#> NCKAP1L2         7.414758e-08 0.995189947 0.619 0.160  2.224427e-04      10
#> CTNNBIP13        1.084256e-07 1.955028117 0.905 0.477  3.252767e-04      10
#> B2M3             1.099200e-07 1.232505749 1.000 0.435  3.297600e-04      10
#> DUSP62           1.120167e-07 1.309529793 0.857 0.362  3.360501e-04      10
#> ADM1             1.189660e-07 0.799584100 0.714 0.204  3.568981e-04      10
#> CITED21          1.232543e-07 1.763361990 0.810 0.398  3.697628e-04      10
#> SLC25A54         1.436996e-07 1.515363845 0.810 0.490  4.310988e-04      10
#> CYBA3            1.537205e-07 1.330110836 1.000 0.392  4.611614e-04      10
#> STEAP12          1.802166e-07 0.860349780 0.381 0.069  5.406497e-04      10
#> VIM3             1.880580e-07 1.095746416 1.000 0.533  5.641740e-04      10
#> FCGRT3           1.967954e-07 1.261121475 0.905 0.361  5.903861e-04      10
#> ADAM282          2.177396e-07 0.877591223 0.524 0.122  6.532187e-04      10
#> CLEC2B3          2.467207e-07 1.074785331 0.762 0.251  7.401622e-04      10
#> S100A43          2.811467e-07 1.183684169 0.952 0.467  8.434401e-04      10
#> APOE1            3.586972e-07 1.382793557 0.905 0.226  1.076092e-03      10
#> RHOF1            4.175788e-07 1.338445819 0.429 0.101  1.252736e-03      10
#> LILRB42          4.196254e-07 0.757281604 0.619 0.168  1.258876e-03      10
#> LAT23            4.219813e-07 1.311417430 0.762 0.313  1.265944e-03      10
#> TRPV22           6.349494e-07 1.127470899 0.476 0.117  1.904848e-03      10
#> BTK2             8.501790e-07 1.143877804 0.429 0.099  2.550537e-03      10
#> ZEB11            1.313003e-06 1.032521725 0.286 0.048  3.939010e-03      10
#> ABI32            1.489241e-06 0.955166341 0.524 0.145  4.467723e-03      10
#> RPS33            1.656312e-06 0.920923845 1.000 0.464  4.968935e-03      10
#> PDLIM42          1.710607e-06 1.302934524 0.667 0.257  5.131820e-03      10
#> LDHB5            1.743804e-06 1.025086816 0.952 0.492  5.231411e-03      10
#> TGM21            1.973090e-06 0.965321240 0.381 0.084  5.919271e-03      10
#> GNA152           2.107644e-06 1.184915033 0.524 0.154  6.322931e-03      10
#> APOC12           2.181645e-06 0.703538006 0.714 0.208  6.544936e-03      10
#> SQRDL3           2.676581e-06 1.174124392 0.810 0.416  8.029744e-03      10
#> HNMT3            3.002589e-06 0.712781971 0.857 0.270  9.007766e-03      10
#> HLA-A4           5.962281e-06 1.041031906 1.000 0.419  1.788684e-02      10
#> SLC16A32         7.085284e-06 0.972044283 0.905 0.418  2.125585e-02      10
#> HSD17B113        8.583190e-06 0.649445352 0.714 0.223  2.574957e-02      10
#> SLC39A82         1.106169e-05 0.821658302 0.619 0.217  3.318507e-02      10
#> PALLD5           1.302815e-05 0.969058980 0.762 0.367  3.908444e-02      10
#> GAPDH4           1.303563e-05 0.930253277 0.905 0.504  3.910688e-02      10
#> C1orf212         1.402597e-05 1.160546658 0.524 0.186  4.207792e-02      10
#> CD532            1.490515e-05 0.710091778 0.667 0.214  4.471546e-02      10
#> IL27RA1          1.507753e-05 0.751946117 0.286 0.056  4.523260e-02      10
#> PTTG14           1.701069e-05 0.607343429 0.667 0.223  5.103207e-02      10
#> FTL2             1.710322e-05 0.605608305 0.952 0.309  5.130967e-02      10
#> CKB6             1.771764e-05 1.520395129 0.762 0.450  5.315292e-02      10
#> FAM46A1          2.248436e-05 1.118174034 0.714 0.309  6.745308e-02      10
#> AIF12            2.655514e-05 0.425278415 0.619 0.227  7.966543e-02      10
#> ANKRD282         2.862943e-05 1.150292286 0.524 0.200  8.588830e-02      10
#> PSMB83           3.177812e-05 0.831076574 0.810 0.331  9.533435e-02      10
#> CD42             3.199003e-05 0.732065345 0.619 0.207  9.597008e-02      10
#> HAVCR22          3.980910e-05 0.572941998 0.476 0.135  1.194273e-01      10
#> ID22             4.035389e-05 0.934075533 0.762 0.368  1.210617e-01      10
#> ZEB23            4.190705e-05 0.514022762 0.667 0.215  1.257212e-01      10
#> ADAM82           4.244305e-05 0.555951675 0.476 0.139  1.273291e-01      10
#> C10orf542        4.276916e-05 0.889686850 0.476 0.158  1.283075e-01      10
#> COL18A11         4.352584e-05 0.402173061 0.476 0.157  1.305775e-01      10
#> CRIP13           5.988477e-05 0.819506133 0.476 0.157  1.796543e-01      10
#> CA22             6.130431e-05 1.059220629 0.571 0.241  1.839129e-01      10
#> BCAT13           6.946700e-05 0.638970021 0.429 0.124  2.084010e-01      10
#> MPP12            7.735605e-05 0.823797844 0.524 0.189  2.320681e-01      10
#> SLC2A33          7.756708e-05 0.450570548 0.667 0.222  2.327012e-01      10
#> HTRA13           9.587806e-05 0.561500832 0.667 0.224  2.876342e-01      10
#> COX7A11          1.049998e-04 0.543801665 0.286 0.064  3.149993e-01      10
#> IFI163           1.054541e-04 0.559095968 0.714 0.277  3.163623e-01      10
#> RAI142           1.058003e-04 0.832275896 0.762 0.361  3.174009e-01      10
#> HIVEP31          1.122964e-04 0.784675576 0.286 0.068  3.368892e-01      10
#> FAM26F2          1.418337e-04 0.604122377 0.571 0.205  4.255011e-01      10
#> VWA5A1           1.453730e-04 1.129996358 0.571 0.269  4.361191e-01      10
#> RNASET23         1.504122e-04 0.703148123 0.857 0.380  4.512366e-01      10
#> EFHD22           1.583700e-04 0.708212717 0.762 0.336  4.751101e-01      10
#> FAM198B3         1.634887e-04 0.473357693 0.476 0.145  4.904662e-01      10
#> MYO1F2           1.927109e-04 0.714828516 0.619 0.243  5.781328e-01      10
#> TRAC1            2.068170e-04 0.877420535 0.333 0.095  6.204511e-01      10
#> GALM3            2.222653e-04 0.763499950 0.476 0.172  6.667959e-01      10
#> ENO2             2.253730e-04 0.697035140 0.286 0.072  6.761191e-01      10
#> LST12            2.867963e-04 0.387482291 0.619 0.223  8.603890e-01      10
#> OTOA2            3.044504e-04 0.572148231 0.381 0.116  9.133511e-01      10
#> PTPN71           3.267139e-04 0.831943779 0.333 0.099  9.801417e-01      10
#> CASP44           3.578057e-04 0.694855492 0.667 0.290  1.000000e+00      10
#> MYO1B4           3.795664e-04 0.929667398 0.762 0.409  1.000000e+00      10
#> OAS13            4.315851e-04 0.714392626 0.429 0.149  1.000000e+00      10
#> PLXND12          4.610171e-04 0.549562785 0.619 0.254  1.000000e+00      10
#> FKBP112          5.051461e-04 0.732377578 0.524 0.239  1.000000e+00      10
#> GPNMB2           5.877544e-04 0.652944728 0.762 0.446  1.000000e+00      10
#> SAMHD12          5.979495e-04 0.628381410 0.571 0.224  1.000000e+00      10
#> ABCC32           6.163855e-04 0.582880623 0.524 0.203  1.000000e+00      10
#> MTSS12           6.253020e-04 0.923437900 0.619 0.301  1.000000e+00      10
#> GBP23            6.339003e-04 0.810457509 0.524 0.230  1.000000e+00      10
#> LYZ2             6.670236e-04 0.351574075 0.619 0.227  1.000000e+00      10
#> GCHFR2           7.146469e-04 0.749451124 0.524 0.237  1.000000e+00      10
#> TCF43            7.316523e-04 0.489122893 0.667 0.267  1.000000e+00      10
#> LGALS92          9.169784e-04 0.600057643 0.571 0.230  1.000000e+00      10
#> BMP2K2           9.258704e-04 0.511637121 0.571 0.222  1.000000e+00      10
#> ABL21            9.835308e-04 0.826238578 0.667 0.365  1.000000e+00      10
#> CTSC2            1.005930e-03 0.420857602 0.762 0.317  1.000000e+00      10
#> NCF12            1.076094e-03 0.399626792 0.524 0.197  1.000000e+00      10
#> PLEK3            1.092326e-03 0.440942924 0.524 0.194  1.000000e+00      10
#> FUCA11           1.118580e-03 0.510042888 0.667 0.319  1.000000e+00      10
#> SRGAP31          1.141734e-03 0.811159941 0.381 0.142  1.000000e+00      10
#> ARFGEF32         1.145812e-03 0.950211485 0.429 0.176  1.000000e+00      10
#> RAB313           1.181928e-03 0.663713984 0.714 0.335  1.000000e+00      10
#> RNASE11          1.208753e-03 0.419296998 0.619 0.250  1.000000e+00      10
#> PFKP2            1.336759e-03 0.719897212 0.762 0.430  1.000000e+00      10
#> IFI27L23         1.629943e-03 0.575358768 0.619 0.273  1.000000e+00      10
#> DHTKD14          1.823544e-03 0.731458135 0.714 0.424  1.000000e+00      10
#> DOK22            1.841042e-03 0.453668413 0.333 0.107  1.000000e+00      10
#> ARL4C3           2.266601e-03 0.421809551 0.762 0.345  1.000000e+00      10
#> FOXP13           2.381156e-03 0.617534517 0.524 0.237  1.000000e+00      10
#> NCF42            2.401690e-03 0.455704133 0.381 0.134  1.000000e+00      10
#> PYCARD2          2.407143e-03 0.175456871 0.667 0.242  1.000000e+00      10
#> RRAGD2           3.423367e-03 0.588743327 0.476 0.217  1.000000e+00      10
#> GMFG2            3.465659e-03 0.464203643 0.619 0.254  1.000000e+00      10
#> DCXR4            5.303776e-03 0.708856370 0.714 0.446  1.000000e+00      10
#> CTSB3            5.315329e-03 0.242152772 0.810 0.324  1.000000e+00      10
#> SAMD9L3          5.592145e-03 0.444926520 0.286 0.098  1.000000e+00      10
#> LGALS33          5.594234e-03 0.558849982 0.762 0.528  1.000000e+00      10
#> ECM12            5.950262e-03 0.265761588 0.333 0.118  1.000000e+00      10
#> HSD17B43         6.002370e-03 0.704214692 0.333 0.125  1.000000e+00      10
#> COL1A21          6.065609e-03 0.069584436 0.429 0.170  1.000000e+00      10
#> MBP2             6.731421e-03 0.894038997 0.667 0.428  1.000000e+00      10
#> NRP22            7.211044e-03 0.458249082 0.810 0.413  1.000000e+00      10
#> DUT4             7.689111e-03 0.518906532 0.714 0.466  1.000000e+00      10
#> LRRC252          8.914669e-03 0.333464314 0.333 0.122  1.000000e+00      10
#> ARHGAP183        9.021814e-03 0.457351188 0.524 0.234  1.000000e+00      10
#> DNMT13           9.329277e-03 0.559201674 0.667 0.379  1.000000e+00      10
#> SLC9A72          9.453143e-03 0.818687217 0.381 0.191  1.000000e+00      10
#> PSMB93           9.500636e-03 0.310366484 0.524 0.235  1.000000e+00      10
#>                           gene
#> IFRD1                    IFRD1
#> CEBPD                    CEBPD
#> SAT1                      SAT1
#> CLU                        CLU
#> MMP7                      MMP7
#> SBSPON                  SBSPON
#> RGS16                    RGS16
#> MDFI                      MDFI
#> MAOB                      MAOB
#> SOD3                      SOD3
#> CNTNAP3B              CNTNAP3B
#> C1orf186              C1orf186
#> PHLDA1                  PHLDA1
#> MGP                        MGP
#> GLS                        GLS
#> HIST2H2BE            HIST2H2BE
#> FBXO32                  FBXO32
#> SLC43A3                SLC43A3
#> PGM2L1                  PGM2L1
#> CHI3L1                  CHI3L1
#> TSPAN12                TSPAN12
#> HSPA1B                  HSPA1B
#> TOB1                      TOB1
#> CHI3L2                  CHI3L2
#> TMEM61                  TMEM61
#> BAMBI                    BAMBI
#> SLC29A1                SLC29A1
#> CLDN4                    CLDN4
#> SFRP1                    SFRP1
#> PLA2G4A                PLA2G4A
#> S100P                    S100P
#> TIFA                      TIFA
#> NNMT                      NNMT
#> HSPA1A                  HSPA1A
#> CDC25B                  CDC25B
#> TUBB2B                  TUBB2B
#> SFN                        SFN
#> NCOA7                    NCOA7
#> S100A1                  S100A1
#> EFNA1                    EFNA1
#> BGN                        BGN
#> CRISPLD1              CRISPLD1
#> OVOS2                    OVOS2
#> TIMM10                  TIMM10
#> FAM3C                    FAM3C
#> GCSH                      GCSH
#> QPCT                      QPCT
#> SOX4                      SOX4
#> SNAI2                    SNAI2
#> NRP2                      NRP2
#> RHOBTB3                RHOBTB3
#> MLLT11                  MLLT11
#> TUBB2A                  TUBB2A
#> HIST1H2AE            HIST1H2AE
#> PCOLCE2                PCOLCE2
#> EFNB2                    EFNB2
#> EPHX2                    EPHX2
#> FZD5                      FZD5
#> GPM6B                    GPM6B
#> MGST1                    MGST1
#> CDKN1A                  CDKN1A
#> FAM89A                  FAM89A
#> SDC4                      SDC4
#> PDZK1IP1              PDZK1IP1
#> TUBA1A                  TUBA1A
#> HSPB1                    HSPB1
#> PIM1                      PIM1
#> ARID5B                  ARID5B
#> GMPR                      GMPR
#> PAM                        PAM
#> IDH1                      IDH1
#> SLC6A15                SLC6A15
#> MATN3                    MATN3
#> TNFSF13B              TNFSF13B
#> PRELP                    PRELP
#> DLX5                      DLX5
#> TMX2                      TMX2
#> FXYD3                    FXYD3
#> TM4SF1                  TM4SF1
#> IGFBP7                  IGFBP7
#> PTGS2                    PTGS2
#> DEFB1                    DEFB1
#> SSRP1                    SSRP1
#> ACOX2                    ACOX2
#> ISG20                    ISG20
#> SLC25A37              SLC25A37
#> SLPI                      SLPI
#> TNFRSF21              TNFRSF21
#> HILPDA                  HILPDA
#> SOHLH1                  SOHLH1
#> IMPA2                    IMPA2
#> HIST1H1C              HIST1H1C
#> KLF10                    KLF10
#> SNX22                    SNX22
#> MIA                        MIA
#> GAS6                      GAS6
#> C10orf10              C10orf10
#> PROM1                    PROM1
#> TPD52L1                TPD52L1
#> PLOD2                    PLOD2
#> METTL7A                METTL7A
#> S100B                    S100B
#> MYLK                      MYLK
#> LGALS3                  LGALS3
#> RPL39L                  RPL39L
#> CITED4                  CITED4
#> CKS2                      CKS2
#> SERPINE2              SERPINE2
#> EPHX1                    EPHX1
#> DBI                        DBI
#> DSC3                      DSC3
#> C2orf82                C2orf82
#> NDRG2                    NDRG2
#> MYO1B                    MYO1B
#> CTHRC1                  CTHRC1
#> ZBTB10                  ZBTB10
#> HIST1H2BG            HIST1H2BG
#> MMP2                      MMP2
#> LCN2                      LCN2
#> ROPN1B                  ROPN1B
#> GATA3                    GATA3
#> TINCR                    TINCR
#> TMEM158                TMEM158
#> AC005152.3          AC005152.3
#> B3GNT7                  B3GNT7
#> SOCS3                    SOCS3
#> BTG2                      BTG2
#> PPP1R1B                PPP1R1B
#> NDRG1                    NDRG1
#> NAB1                      NAB1
#> CSRP1                    CSRP1
#> KLK1                      KLK1
#> MAP1B                    MAP1B
#> IFITM1                  IFITM1
#> STAT1                    STAT1
#> GPRC5A                  GPRC5A
#> MYL9                      MYL9
#> IDI1                      IDI1
#> SCIN                      SCIN
#> MAFF                      MAFF
#> ELF3                      ELF3
#> C1R                        C1R
#> AGT                        AGT
#> AQP3                      AQP3
#> A2M                        A2M
#> SERPINA5              SERPINA5
#> SAA1                      SAA1
#> PMP22                    PMP22
#> SCRG1                    SCRG1
#> DDIT3                    DDIT3
#> KCNMB1                  KCNMB1
#> TOX                        TOX
#> TMEM79                  TMEM79
#> C1QL4                    C1QL4
#> SDC2                      SDC2
#> ST3GAL1                ST3GAL1
#> PROX1                    PROX1
#> MAL2                      MAL2
#> TBC1D1                  TBC1D1
#> CKS1B                    CKS1B
#> HES1                      HES1
#> SOX8                      SOX8
#> RAD21                    RAD21
#> C1S                        C1S
#> CDKN2A                  CDKN2A
#> VWA5A                    VWA5A
#> SERTAD4                SERTAD4
#> COL11A1                COL11A1
#> CTGF                      CTGF
#> CYP26B1                CYP26B1
#> DTNB                      DTNB
#> SLC9A3R2              SLC9A3R2
#> HIBCH                    HIBCH
#> MARCKSL1              MARCKSL1
#> MAPK8IP2              MAPK8IP2
#> CACYBP                  CACYBP
#> NEGR1                    NEGR1
#> HMGB1                    HMGB1
#> SH3BGR                  SH3BGR
#> SDC1                      SDC1
#> FMO2                      FMO2
#> FXYD6                    FXYD6
#> ANXA1                    ANXA1
#> ACAN                      ACAN
#> C8orf46                C8orf46
#> ABCA5                    ABCA5
#> MFAP2                    MFAP2
#> CPED1                    CPED1
#> ROPN1                    ROPN1
#> BOC                        BOC
#> CCDC80                  CCDC80
#> NEAT1                    NEAT1
#> ENPP5                    ENPP5
#> ARL6IP1                ARL6IP1
#> TACSTD2                TACSTD2
#> GLRX                      GLRX
#> C1orf56                C1orf56
#> CP                          CP
#> TEKT3                    TEKT3
#> COL9A3                  COL9A3
#> SERPINB5              SERPINB5
#> NOTUM                    NOTUM
#> IFI6                      IFI6
#> CRYAB                    CRYAB
#> PEG10                    PEG10
#> KRT15                    KRT15
#> PLPP3                    PLPP3
#> GOLM1                    GOLM1
#> TAGLN                    TAGLN
#> FKBP4                    FKBP4
#> TMEM106C              TMEM106C
#> FAM84A                  FAM84A
#> HRCT1                    HRCT1
#> KRT7                      KRT7
#> PLOD3                    PLOD3
#> RP11-357H14.17  RP11-357H14.17
#> CD24                      CD24
#> SMTN                      SMTN
#> ERRFI1                  ERRFI1
#> SLC26A2                SLC26A2
#> NET1                      NET1
#> ITGB4                    ITGB4
#> SLBP                      SLBP
#> GSDMC                    GSDMC
#> KLF2                      KLF2
#> VTCN1                    VTCN1
#> SERTAD4-AS1        SERTAD4-AS1
#> LNX1                      LNX1
#> CLDN3                    CLDN3
#> ACTL6A                  ACTL6A
#> DNM3OS                  DNM3OS
#> SOX9                      SOX9
#> CYP1B1                  CYP1B1
#> PPFIA1                  PPFIA1
#> IFIH1                    IFIH1
#> LEMD1                    LEMD1
#> CARHSP1                CARHSP1
#> CLCN4                    CLCN4
#> IER5                      IER5
#> LINC01436            LINC01436
#> CTNND2                  CTNND2
#> SELM                      SELM
#> HSPH1                    HSPH1
#> PTRF                      PTRF
#> SUN3                      SUN3
#> SQLE                      SQLE
#> KLF5                      KLF5
#> JHDM1D-AS1          JHDM1D-AS1
#> AZGP1                    AZGP1
#> ATF3                      ATF3
#> ACTA2                    ACTA2
#> SPON2                    SPON2
#> TSPAN2                  TSPAN2
#> LAT2                      LAT2
#> NFKBIE                  NFKBIE
#> VANGL1                  VANGL1
#> AP1M2                    AP1M2
#> KRT18                    KRT18
#> TMSB15A                TMSB15A
#> TRIB1                    TRIB1
#> MBNL1-AS1            MBNL1-AS1
#> TTC39A                  TTC39A
#> IFT172                  IFT172
#> CYP39A1                CYP39A1
#> RP11-25K19.1      RP11-25K19.1
#> NCCRP1                  NCCRP1
#> HEY2                      HEY2
#> MALL                      MALL
#> CCT5                      CCT5
#> HIST1H4E              HIST1H4E
#> ALDH1B1                ALDH1B1
#> NR2F2                    NR2F2
#> ICAM1                    ICAM1
#> DNM3                      DNM3
#> KLHDC3                  KLHDC3
#> RAB3IP                  RAB3IP
#> CIART                    CIART
#> NKD2                      NKD2
#> ITGA10                  ITGA10
#> ANXA2R                  ANXA2R
#> MYO6                      MYO6
#> PRSS33                  PRSS33
#> KCTD1                    KCTD1
#> RDH10                    RDH10
#> PRPS2                    PRPS2
#> LAMB3                    LAMB3
#> ELF5                      ELF5
#> FBXO2                    FBXO2
#> NFKBIA                  NFKBIA
#> S100A6                  S100A6
#> TPM2                      TPM2
#> PRSS8                    PRSS8
#> TFAP2A                  TFAP2A
#> HACD1                    HACD1
#> NFATC1                  NFATC1
#> TPM1                      TPM1
#> FHL1                      FHL1
#> PFKP                      PFKP
#> MT2A                      MT2A
#> SORBS2                  SORBS2
#> FGF13                    FGF13
#> CALD1                    CALD1
#> GGCT                      GGCT
#> SYT8                      SYT8
#> KRT8                      KRT8
#> RP11-19E11.1      RP11-19E11.1
#> PEG3                      PEG3
#> TUBB4B                  TUBB4B
#> MSRB3                    MSRB3
#> CDC42EP1              CDC42EP1
#> S100A10                S100A10
#> FBLN2                    FBLN2
#> POSTN                    POSTN
#> SLC40A1                SLC40A1
#> FAM46B                  FAM46B
#> WDR34                    WDR34
#> BACE2                    BACE2
#> LMCD1                    LMCD1
#> NPR3                      NPR3
#> IL32                      IL32
#> PBX1                      PBX1
#> TUBB                      TUBB
#> PRSS22                  PRSS22
#> SLC39A14              SLC39A14
#> CA8                        CA8
#> CTNNB1                  CTNNB1
#> OPRK1                    OPRK1
#> IL17B                    IL17B
#> CSF3R                    CSF3R
#> LAPTM4B                LAPTM4B
#> SOD2                      SOD2
#> PDGFRA                  PDGFRA
#> HEBP2                    HEBP2
#> PYCR1                    PYCR1
#> TSPAN5                  TSPAN5
#> MYBL1                    MYBL1
#> NQO1                      NQO1
#> SPINT1                  SPINT1
#> C6orf132              C6orf132
#> RTN4RL2                RTN4RL2
#> PLSCR1                  PLSCR1
#> CDK4                      CDK4
#> STMN1                    STMN1
#> HSP90AB1              HSP90AB1
#> MESP1                    MESP1
#> PTS                        PTS
#> LIMCH1                  LIMCH1
#> DNTTIP1                DNTTIP1
#> PIR                        PIR
#> SYNGR1                  SYNGR1
#> APP                        APP
#> ABL2                      ABL2
#> TM4SF1-AS1          TM4SF1-AS1
#> TRAF3IP3              TRAF3IP3
#> GPNMB                    GPNMB
#> MFSD6                    MFSD6
#> CAPN2                    CAPN2
#> OVOL1                    OVOL1
#> HMGA1                    HMGA1
#> MB                          MB
#> S100A9                  S100A9
#> EPS8                      EPS8
#> VASN                      VASN
#> NUDT4                    NUDT4
#> NSG1                      NSG1
#> COL11A2                COL11A2
#> UBE2L6                  UBE2L6
#> PDLIM4                  PDLIM4
#> LAMB1                    LAMB1
#> ABCC3                    ABCC3
#> FOXC1                    FOXC1
#> KIAA0040              KIAA0040
#> INSIG1                  INSIG1
#> BARD1                    BARD1
#> KCNQ1OT1              KCNQ1OT1
#> SNHG25                  SNHG25
#> SELENBP1              SELENBP1
#> PTPRF                    PTPRF
#> GRB14                    GRB14
#> ERBB3                    ERBB3
#> CD163                    CD163
#> MS4A7                    MS4A7
#> CYBB                      CYBB
#> FPR3                      FPR3
#> SLCO2B1                SLCO2B1
#> MSR1                      MSR1
#> SPI1                      SPI1
#> SRGN                      SRGN
#> C1QC                      C1QC
#> PLTP                      PLTP
#> PILRA                    PILRA
#> C3AR1                    C3AR1
#> LYZ                        LYZ
#> LY96                      LY96
#> CD53                      CD53
#> LST1                      LST1
#> FCGR1A                  FCGR1A
#> LILRB4                  LILRB4
#> CD4                        CD4
#> PYCARD                  PYCARD
#> C5AR1                    C5AR1
#> SAMSN1                  SAMSN1
#> CD68                      CD68
#> HNMT                      HNMT
#> PLEK                      PLEK
#> MS4A6A                  MS4A6A
#> FCGR2A                  FCGR2A
#> LAPTM5                  LAPTM5
#> MARCH1                  MARCH1
#> RGS1                      RGS1
#> LINC01094            LINC01094
#> TGFBI                    TGFBI
#> AIF1                      AIF1
#> HLA-DQA1              HLA-DQA1
#> MPEG1                    MPEG1
#> ADAP2                    ADAP2
#> HLA-DQA2              HLA-DQA2
#> NPL                        NPL
#> FYB                        FYB
#> GMFG                      GMFG
#> FCER1G                  FCER1G
#> SLC7A7                  SLC7A7
#> C1QA                      C1QA
#> TYROBP                  TYROBP
#> FCGR3A                  FCGR3A
#> ACP5                      ACP5
#> HLA-DRB5              HLA-DRB5
#> TREM2                    TREM2
#> KCTD12                  KCTD12
#> PTPRC                    PTPRC
#> ARHGAP18              ARHGAP18
#> HLA-DMB                HLA-DMB
#> IL18                      IL18
#> GPR34                    GPR34
#> C1QB                      C1QB
#> DMXL2                    DMXL2
#> HCLS1                    HCLS1
#> KYNU                      KYNU
#> MS4A4A                  MS4A4A
#> CTSS                      CTSS
#> STAB1                    STAB1
#> CD86                      CD86
#> FERMT3                  FERMT3
#> CD14                      CD14
#> C1orf162              C1orf162
#> HLA-DPA1              HLA-DPA1
#> LY86                      LY86
#> ITGB2                    ITGB2
#> APOC1                    APOC1
#> RP11-1143G9.4    RP11-1143G9.4
#> GIMAP4                  GIMAP4
#> HLA-DMA                HLA-DMA
#> CD84                      CD84
#> CCL3                      CCL3
#> PLA2G7                  PLA2G7
#> ZEB2                      ZEB2
#> MPP1                      MPP1
#> FTL                        FTL
#> IGSF6                    IGSF6
#> SAMHD1                  SAMHD1
#> IL7R                      IL7R
#> LAIR1                    LAIR1
#> SLAMF8                  SLAMF8
#> CXCL3                    CXCL3
#> IFI30                    IFI30
#> HCK                        HCK
#> CTSB                      CTSB
#> CARD16                  CARD16
#> BMP2K                    BMP2K
#> HLA-DRB1              HLA-DRB1
#> HLA-DQB1              HLA-DQB1
#> CYBA                      CYBA
#> APOE                      APOE
#> RNASE1                  RNASE1
#> CXCL2                    CXCL2
#> LCP1                      LCP1
#> EVI2B                    EVI2B
#> CSF1R                    CSF1R
#> DAB2                      DAB2
#> CLEC2B                  CLEC2B
#> CLEC4E                  CLEC4E
#> CD48                      CD48
#> CD74                      CD74
#> MAFB                      MAFB
#> HLA-DRA                HLA-DRA
#> IL4I1                    IL4I1
#> CD37                      CD37
#> HSD17B11              HSD17B11
#> BST2                      BST2
#> HLA-DPB1              HLA-DPB1
#> CAPG                      CAPG
#> CTSC                      CTSC
#> IRF8                      IRF8
#> AKR1B1                  AKR1B1
#> HAVCR2                  HAVCR2
#> IQGAP2                  IQGAP2
#> ADAMDEC1              ADAMDEC1
#> GPX1                      GPX1
#> SLA                        SLA
#> ALOX5AP                ALOX5AP
#> NCKAP1L                NCKAP1L
#> PTAFR                    PTAFR
#> GPR65                    GPR65
#> FCGRT                    FCGRT
#> EMP3                      EMP3
#> HPSE                      HPSE
#> UBB                        UBB
#> HMOX1                    HMOX1
#> MYO1F                    MYO1F
#> RGS2                      RGS2
#> SLC8A1                  SLC8A1
#> LRRC25                  LRRC25
#> PLAUR                    PLAUR
#> EPB41L3                EPB41L3
#> NCF2                      NCF2
#> XIST                      XIST
#> LGALS9                  LGALS9
#> MARCKS                  MARCKS
#> SPP1                      SPP1
#> C1orf54                C1orf54
#> C15orf48              C15orf48
#> RNASE6                  RNASE6
#> CXCR4                    CXCR4
#> B2M                        B2M
#> FCGR2B                  FCGR2B
#> NCF4                      NCF4
#> SNX10                    SNX10
#> AP1S2                    AP1S2
#> EBI3                      EBI3
#> PLXND1                  PLXND1
#> VSIG4                    VSIG4
#> IGSF21                  IGSF21
#> LSP1                      LSP1
#> CXCL8                    CXCL8
#> CCR1                      CCR1
#> ITGAX                    ITGAX
#> ARRB2                    ARRB2
#> CFD                        CFD
#> CTSL                      CTSL
#> RNASET2                RNASET2
#> DUSP1                    DUSP1
#> CD83                      CD83
#> LPAR6                    LPAR6
#> LILRB1                  LILRB1
#> TFEC                      TFEC
#> NR1H3                    NR1H3
#> KCNMA1                  KCNMA1
#> OTOA                      OTOA
#> HLA-A                    HLA-A
#> ATF5                      ATF5
#> TCF4                      TCF4
#> SLC2A3                  SLC2A3
#> CSTA                      CSTA
#> GATM                      GATM
#> ARHGDIB                ARHGDIB
#> EVI2A                    EVI2A
#> NCF1                      NCF1
#> EFHD2                    EFHD2
#> GPR183                  GPR183
#> CCL4                      CCL4
#> IFI27L2                IFI27L2
#> SLC16A10              SLC16A10
#> LIMS1                    LIMS1
#> SDS                        SDS
#> FOLR2                    FOLR2
#> CYTH4                    CYTH4
#> LPXN                      LPXN
#> LILRB2                  LILRB2
#> HTRA1                    HTRA1
#> PLAU                      PLAU
#> MNDA                      MNDA
#> CMKLR1                  CMKLR1
#> ADAM8                    ADAM8
#> GPNMB1                   GPNMB
#> CST3                      CST3
#> SOD21                     SOD2
#> ARL4C                    ARL4C
#> APOC2                    APOC2
#> ID2                        ID2
#> TGFB1                    TGFB1
#> GPR84                    GPR84
#> PLXNC1                  PLXNC1
#> GPX3                      GPX3
#> HLA-DOA                HLA-DOA
#> TNFAIP8                TNFAIP8
#> OSM                        OSM
#> ABCA1                    ABCA1
#> EID1                      EID1
#> IL10RA                  IL10RA
#> SLC16A3                SLC16A3
#> SGK1                      SGK1
#> PKIB                      PKIB
#> NCEH1                    NCEH1
#> ALOX5                    ALOX5
#> CCL3L3                  CCL3L3
#> SH2B3                    SH2B3
#> SERPINF1              SERPINF1
#> CCL5                      CCL5
#> MT1H                      MT1H
#> NAIP                      NAIP
#> ALDH2                    ALDH2
#> RAB31                    RAB31
#> VMO1                      VMO1
#> RASSF4                  RASSF4
#> OLR1                      OLR1
#> PARVG                    PARVG
#> GK                          GK
#> TRPV2                    TRPV2
#> FMNL1                    FMNL1
#> CPVL                      CPVL
#> EPSTI1                  EPSTI1
#> CLEC5A                  CLEC5A
#> HCST                      HCST
#> ADAM28                  ADAM28
#> CCL18                    CCL18
#> FAM198B                FAM198B
#> C2                          C2
#> SLC11A1                SLC11A1
#> BHLHE41                BHLHE41
#> MEF2C                    MEF2C
#> BCAT1                    BCAT1
#> TBXAS1                  TBXAS1
#> IFI16                    IFI16
#> STX11                    STX11
#> SLC39A8                SLC39A8
#> DOCK10                  DOCK10
#> ABI3                      ABI3
#> PRDM1                    PRDM1
#> DSE                        DSE
#> AOAH                      AOAH
#> CYTIP                    CYTIP
#> BTK                        BTK
#> PTGS1                    PTGS1
#> CACNA2D4              CACNA2D4
#> ADGRE2                  ADGRE2
#> NRP1                      NRP1
#> PECAM1                  PECAM1
#> CLEC4A                  CLEC4A
#> LGMN                      LGMN
#> AQP9                      AQP9
#> DOK2                      DOK2
#> FAM26F                  FAM26F
#> PTPRE                    PTPRE
#> CASP1                    CASP1
#> IL10                      IL10
#> MCOLN2                  MCOLN2
#> SEPP1                    SEPP1
#> SLC1A3                  SLC1A3
#> MMP9                      MMP9
#> LXN                        LXN
#> GNA15                    GNA15
#> BCL2A1                  BCL2A1
#> PDE4B                    PDE4B
#> PARP14                  PARP14
#> CXorf21                CXorf21
#> CD93                      CD93
#> SERPINA1              SERPINA1
#> CD40                      CD40
#> NR4A2                    NR4A2
#> DLEU7                    DLEU7
#> IL1B                      IL1B
#> CCL4L2                  CCL4L2
#> TLR4                      TLR4
#> CCL2                      CCL2
#> MAF                        MAF
#> GAPLINC                GAPLINC
#> ZNF331                  ZNF331
#> VIM                        VIM
#> AC092484.1          AC092484.1
#> MT1M                      MT1M
#> CD300A                  CD300A
#> SIGLEC10              SIGLEC10
#> LIPA                      LIPA
#> MVP                        MVP
#> TNFAIP3                TNFAIP3
#> NFKBIA1                 NFKBIA
#> OAS1                      OAS1
#> ZFP36                    ZFP36
#> SLC2A5                  SLC2A5
#> TNFRSF4                TNFRSF4
#> SERPINB9              SERPINB9
#> CXCL16                  CXCL16
#> CD300E                  CD300E
#> CASP4                    CASP4
#> IL2RG                    IL2RG
#> G0S2                      G0S2
#> HSD17B4                HSD17B4
#> WAS                        WAS
#> VASH1                    VASH1
#> PTGER4                  PTGER4
#> IFI35                    IFI35
#> OLFML2B                OLFML2B
#> APOBEC3C              APOBEC3C
#> FGL2                      FGL2
#> GIMAP7                  GIMAP7
#> HS3ST1                  HS3ST1
#> CKLF                      CKLF
#> FBP1                      FBP1
#> GALM                      GALM
#> HLA-B                    HLA-B
#> CLECL1                  CLECL1
#> OSCAR                    OSCAR
#> TREM1                    TREM1
#> SGPP1                    SGPP1
#> FN1                        FN1
#> BASP1                    BASP1
#> SPARC                    SPARC
#> TNFAIP2                TNFAIP2
#> CXCL1                    CXCL1
#> PSMB9                    PSMB9
#> RGS10                    RGS10
#> POU2F2                  POU2F2
#> GBP1                      GBP1
#> SELPLG                  SELPLG
#> ETV5                      ETV5
#> THBD                      THBD
#> SAMD9L                  SAMD9L
#> C10orf54              C10orf54
#> HBEGF                    HBEGF
#> UCP2                      UCP2
#> ABCC31                   ABCC3
#> CD36                      CD36
#> IER3                      IER3
#> FGR                        FGR
#> COQ2                      COQ2
#> PPM1N                    PPM1N
#> IRF7                      IRF7
#> IL1RN                    IL1RN
#> MT1G                      MT1G
#> SQRDL                    SQRDL
#> PTPN6                    PTPN6
#> PDPN                      PDPN
#> GBP2                      GBP2
#> CELF2                    CELF2
#> PSTPIP1                PSTPIP1
#> ERO1A                    ERO1A
#> CD33                      CD33
#> FABP5                    FABP5
#> PSMB8                    PSMB8
#> MYO1G                    MYO1G
#> CORO1A                  CORO1A
#> PTPN22                  PTPN22
#> BID                        BID
#> CSF2RA                  CSF2RA
#> MITF                      MITF
#> SAMD9                    SAMD9
#> FUCA1                    FUCA1
#> PLEK2                    PLEK2
#> RAC2                      RAC2
#> MGLL                      MGLL
#> LGALS3BP              LGALS3BP
#> IFI44                    IFI44
#> TRIM22                  TRIM22
#> COLEC12                COLEC12
#> LINC00936            LINC00936
#> S100A91                 S100A9
#> ICAM11                   ICAM1
#> TSC22D3                TSC22D3
#> MT1F                      MT1F
#> SEPT6                    SEPT6
#> ANPEP                    ANPEP
#> LIMD2                    LIMD2
#> FOXP1                    FOXP1
#> CSTB                      CSTB
#> CYP27A1                CYP27A1
#> RHOB                      RHOB
#> ECM1                      ECM1
#> PLSCR11                 PLSCR1
#> VAMP5                    VAMP5
#> CD109                    CD109
#> XAF1                      XAF1
#> FABP3                    FABP3
#> NFKBID                  NFKBID
#> FRMD4A                  FRMD4A
#> BTG21                     BTG2
#> SLC25A5                SLC25A5
#> DUSP6                    DUSP6
#> CRIP1                    CRIP1
#> FILIP1L                FILIP1L
#> TCIRG1                  TCIRG1
#> LAT21                     LAT2
#> ID3                        ID3
#> SCD                        SCD
#> C3                          C3
#> FAM105A                FAM105A
#> NKG7                      NKG7
#> TNFSF10                TNFSF10
#> CTSK                      CTSK
#> HLA-F                    HLA-F
#> CXCL12                  CXCL12
#> CITED2                  CITED2
#> TAP1                      TAP1
#> DUSP2                    DUSP2
#> MX1                        MX1
#> IFI61                     IFI6
#> SORL1                    SORL1
#> CD82                      CD82
#> GNPTAB                  GNPTAB
#> RARRES1                RARRES1
#> PARVB                    PARVB
#> ACSL1                    ACSL1
#> SAT11                     SAT1
#> RUNX3                    RUNX3
#> S100A4                  S100A4
#> UBE2L61                 UBE2L6
#> COTL1                    COTL1
#> SLC25A19              SLC25A19
#> FAM162A                FAM162A
#> ISG15                    ISG15
#> IFITM2                  IFITM2
#> GCHFR                    GCHFR
#> EPB41L2                EPB41L2
#> DAPP1                    DAPP1
#> MT1X                      MT1X
#> CAMK1                    CAMK1
#> NRP21                     NRP2
#> MTSS1                    MTSS1
#> SCCPDH                  SCCPDH
#> PPIF                      PPIF
#> CCDC85B                CCDC85B
#> IFIH11                   IFIH1
#> HLA-C                    HLA-C
#> ARPC4                    ARPC4
#> HIST1H4C              HIST1H4C
#> LMO4                      LMO4
#> ANKRD37                ANKRD37
#> CD52                      CD52
#> DNAJC15                DNAJC15
#> TUBA1C                  TUBA1C
#> REEP4                    REEP4
#> C1orf21                C1orf21
#> IGKC                      IGKC
#> ARID5A                  ARID5A
#> MT2A1                     MT2A
#> LGALS31                 LGALS3
#> UPP1                      UPP1
#> MTHFD2                  MTHFD2
#> SOCS31                   SOCS3
#> HBB                        HBB
#> IFITM3                  IFITM3
#> ST3GAL11               ST3GAL1
#> PCSK1N                  PCSK1N
#> CALML5                  CALML5
#> CFB                        CFB
#> SEZ6L2                  SEZ6L2
#> S100A61                 S100A6
#> APOD                      APOD
#> S100A16                S100A16
#> CRLF1                    CRLF1
#> RBM3                      RBM3
#> PLP2                      PLP2
#> UBA1                      UBA1
#> BRK1                      BRK1
#> IFITM31                 IFITM3
#> CDK16                    CDK16
#> KLK6                      KLK6
#> NDUFB11                NDUFB11
#> MAGIX                    MAGIX
#> KRT19                    KRT19
#> KRT17                    KRT17
#> APP1                       APP
#> PKP1                      PKP1
#> PPDPF                    PPDPF
#> FMO21                     FMO2
#> S100A14                S100A14
#> PDXK                      PDXK
#> ATP9A                    ATP9A
#> KLK8                      KLK8
#> WDR13                    WDR13
#> CDH1                      CDH1
#> GABRP                    GABRP
#> PDIA6                    PDIA6
#> EPPK1                    EPPK1
#> CD241                     CD24
#> BPIFB1                  BPIFB1
#> LINC01615            LINC01615
#> SUN1                      SUN1
#> SEC13                    SEC13
#> RPUSD3                  RPUSD3
#> UXT                        UXT
#> DPM1                      DPM1
#> PQBP1                    PQBP1
#> CUL7                      CUL7
#> PPP1R14C              PPP1R14C
#> GRXCR1                  GRXCR1
#> NUPR2                    NUPR2
#> THUMPD3-AS1        THUMPD3-AS1
#> SETD5                    SETD5
#> HEBP21                   HEBP2
#> MEIS3                    MEIS3
#> NPDC1                    NPDC1
#> TIMP1                    TIMP1
#> SMIM22                  SMIM22
#> SOX18                    SOX18
#> FAM92A1                FAM92A1
#> TACSTD21               TACSTD2
#> METRN                    METRN
#> TIMM17B                TIMM17B
#> FNDC4                    FNDC4
#> PDIA3                    PDIA3
#> HDAC6                    HDAC6
#> KRT81                    KRT81
#> NEAT11                   NEAT1
#> DNAJC151               DNAJC15
#> S100A41                 S100A4
#> KCNG1                    KCNG1
#> FADS3                    FADS3
#> FSCN1                    FSCN1
#> C6orf15                C6orf15
#> KLK7                      KLK7
#> ST6GAL1                ST6GAL1
#> PPP1R14B              PPP1R14B
#> IRF2BP2                IRF2BP2
#> PTPRF1                   PTPRF
#> KANK4                    KANK4
#> LRP5                      LRP5
#> FBLN21                   FBLN2
#> SPINT11                 SPINT1
#> TADA3                    TADA3
#> MAD2L1BP              MAD2L1BP
#> KIAA0513              KIAA0513
#> CT83                      CT83
#> PRICKLE1              PRICKLE1
#> APMAP                    APMAP
#> C6orf1321             C6orf132
#> TMEM132A              TMEM132A
#> C1GALT1                C1GALT1
#> ECEL1                    ECEL1
#> PTGFRN                  PTGFRN
#> KRT6B                    KRT6B
#> LINC00176            LINC00176
#> GPC1                      GPC1
#> TPT1-AS1              TPT1-AS1
#> SLC35A2                SLC35A2
#> C2orf821               C2orf82
#> KRBOX4                  KRBOX4
#> S100A2                  S100A2
#> IL17RC                  IL17RC
#> PXDN                      PXDN
#> CHPF                      CHPF
#> TUNAR                    TUNAR
#> DNAAF1                  DNAAF1
#> TRIB2                    TRIB2
#> IL17RE                  IL17RE
#> TRIM58                  TRIM58
#> CACNA1A                CACNA1A
#> PHACTR1                PHACTR1
#> DSP                        DSP
#> SHISA9                  SHISA9
#> GPKOW                    GPKOW
#> NAALADL2              NAALADL2
#> VHL                        VHL
#> FTSJ1                    FTSJ1
#> HOMER2                  HOMER2
#> WDR45                    WDR45
#> PCDHB9                  PCDHB9
#> NOA1                      NOA1
#> SECTM1                  SECTM1
#> EMC3                      EMC3
#> TRDC                      TRDC
#> PVRL4                    PVRL4
#> TTLL3                    TTLL3
#> PRRT3                    PRRT3
#> IRF2BPL                IRF2BPL
#> OTUD5                    OTUD5
#> CALB2                    CALB2
#> JAGN1                    JAGN1
#> CHMP2B                  CHMP2B
#> TARBP1                  TARBP1
#> LAMC2                    LAMC2
#> SLC2A11                SLC2A11
#> CTSF                      CTSF
#> KRT23                    KRT23
#> ITGB8                    ITGB8
#> FUT2                      FUT2
#> FLNA                      FLNA
#> PTPN14                  PTPN14
#> CEBPB                    CEBPB
#> MISP                      MISP
#> CAMK1D                  CAMK1D
#> FGF131                   FGF13
#> TPBG                      TPBG
#> MYO5B                    MYO5B
#> LDOC1                    LDOC1
#> SDR16C5                SDR16C5
#> PRRX2                    PRRX2
#> SEC11C                  SEC11C
#> FKBP9                    FKBP9
#> C1S1                       C1S
#> NPTXR                    NPTXR
#> LINC00960            LINC00960
#> ARAF                      ARAF
#> AGT1                       AGT
#> TRPM8                    TRPM8
#> VEGFA                    VEGFA
#> KCNK12                  KCNK12
#> MUC5B                    MUC5B
#> ZNF674-AS1          ZNF674-AS1
#> NT5E                      NT5E
#> LAD1                      LAD1
#> ARMC1                    ARMC1
#> SIPA1L2                SIPA1L2
#> SDC11                     SDC1
#> ACTR3B                  ACTR3B
#> LINC00518            LINC00518
#> AC022007.5          AC022007.5
#> MYO61                     MYO6
#> CLDN6                    CLDN6
#> WFDC2                    WFDC2
#> SAA11                     SAA1
#> ARPC41                   ARPC4
#> SEMA6A                  SEMA6A
#> GPR27                    GPR27
#> MFI2                      MFI2
#> CRYAB1                   CRYAB
#> KRT16                    KRT16
#> LSR                        LSR
#> CLCN41                   CLCN4
#> PORCN                    PORCN
#> THUMPD3                THUMPD3
#> PRAF2                    PRAF2
#> MTSS1L                  MTSS1L
#> EGLN1                    EGLN1
#> IGSF3                    IGSF3
#> SYNCRIP                SYNCRIP
#> NRBP2                    NRBP2
#> IRX4                      IRX4
#> CHST7                    CHST7
#> MAL21                     MAL2
#> MLLT4                    MLLT4
#> IFIT2                    IFIT2
#> FBXO9                    FBXO9
#> DDR1                      DDR1
#> SLC15A1                SLC15A1
#> NDUFAF6                NDUFAF6
#> THBS1                    THBS1
#> FLNB                      FLNB
#> CLDN7                    CLDN7
#> SLC52A1                SLC52A1
#> PEG31                     PEG3
#> SLC2A12                SLC2A12
#> KCNQ1OT11             KCNQ1OT1
#> NMB                        NMB
#> IL321                     IL32
#> RASD1                    RASD1
#> PADI2                    PADI2
#> FZD8                      FZD8
#> KLF51                     KLF5
#> MMP15                    MMP15
#> MAP2                      MAP2
#> FBXO21                   FBXO2
#> EFNA11                   EFNA1
#> BCAS4                    BCAS4
#> FXYD31                   FXYD3
#> CD55                      CD55
#> KRT71                     KRT7
#> RRP36                    RRP36
#> ALDOC                    ALDOC
#> RGN                        RGN
#> ANPEP1                   ANPEP
#> TBC1D25                TBC1D25
#> KRT86                    KRT86
#> ATP1B1                  ATP1B1
#> LEMD11                   LEMD1
#> OGG1                      OGG1
#> MORN3                    MORN3
#> MAFK                      MAFK
#> GBA                        GBA
#> EGLN3                    EGLN3
#> MANEAL                  MANEAL
#> MUC1                      MUC1
#> TMEM63A                TMEM63A
#> SLC2A1                  SLC2A1
#> TNNC2                    TNNC2
#> STARD10                STARD10
#> KRT14                    KRT14
#> C1R1                       C1R
#> TBC1D7                  TBC1D7
#> PVRL1                    PVRL1
#> GAN                        GAN
#> COL6A2                  COL6A2
#> TLDC1                    TLDC1
#> ZNF654                  ZNF654
#> SULF2                    SULF2
#> AHNAK2                  AHNAK2
#> SAA2                      SAA2
#> GATA6                    GATA6
#> TFAP2A1                 TFAP2A
#> TJP1                      TJP1
#> SLC26A21               SLC26A2
#> FOXC11                   FOXC1
#> ILF2                      ILF2
#> COL9A31                 COL9A3
#> KLK10                    KLK10
#> PSORS1C1              PSORS1C1
#> AQP31                     AQP3
#> PRRT3-AS1            PRRT3-AS1
#> GS1-124K5.4        GS1-124K5.4
#> CGGBP1                  CGGBP1
#> RBP1                      RBP1
#> GALNT2                  GALNT2
#> RGS20                    RGS20
#> MYOF                      MYOF
#> ITGB41                   ITGB4
#> CCDC167                CCDC167
#> EMP1                      EMP1
#> ELK1                      ELK1
#> ANXA3                    ANXA3
#> HSPB2                    HSPB2
#> TMEM208                TMEM208
#> DSEL                      DSEL
#> SORBS21                 SORBS2
#> NTN1                      NTN1
#> GSDMC1                   GSDMC
#> TMEM200A              TMEM200A
#> CRABP2                  CRABP2
#> OVOL11                   OVOL1
#> HYOU1                    HYOU1
#> H2AFY2                  H2AFY2
#> FBXL16                  FBXL16
#> RNF8                      RNF8
#> SOX91                     SOX9
#> RTP4                      RTP4
#> CCDC64B                CCDC64B
#> SLC6A14                SLC6A14
#> LPIN1                    LPIN1
#> PAX1                      PAX1
#> GATA6-AS1            GATA6-AS1
#> ERRFI11                 ERRFI1
#> PGBD5                    PGBD5
#> SERPINB51             SERPINB5
#> BRPF1                    BRPF1
#> MFSD3                    MFSD3
#> ABHD11-AS1          ABHD11-AS1
#> HSPA5                    HSPA5
#> RP11-161M6.2      RP11-161M6.2
#> RIC3                      RIC3
#> MTMR14                  MTMR14
#> CLDN31                   CLDN3
#> ZNRF2                    ZNRF2
#> ARFGEF3                ARFGEF3
#> DDIT4                    DDIT4
#> LSM5                      LSM5
#> NBL1                      NBL1
#> EDN1                      EDN1
#> ACSL11                   ACSL1
#> RHOV                      RHOV
#> ZNF217                  ZNF217
#> ARHGAP29              ARHGAP29
#> FAM208B                FAM208B
#> NEBL                      NEBL
#> CLMN                      CLMN
#> RP11-783K16.5    RP11-783K16.5
#> SCGB3A1                SCGB3A1
#> LYPD3                    LYPD3
#> LDHA                      LDHA
#> MMP11                    MMP11
#> KCNE5                    KCNE5
#> LA16c-380H5.5    LA16c-380H5.5
#> RP11-160O5.1      RP11-160O5.1
#> MCAM                      MCAM
#> EBP                        EBP
#> IRX5                      IRX5
#> MALL1                     MALL
#> GAPDH                    GAPDH
#> NAT8L                    NAT8L
#> LCN21                     LCN2
#> DOK5                      DOK5
#> PDZD2                    PDZD2
#> CDH3                      CDH3
#> SLC7A5                  SLC7A5
#> ELF31                     ELF3
#> OAZ3                      OAZ3
#> NDRG11                   NDRG1
#> NFASC                    NFASC
#> TMPRSS13              TMPRSS13
#> EIF3B                    EIF3B
#> MAP3K8                  MAP3K8
#> ERBB31                   ERBB3
#> SRGAP3                  SRGAP3
#> GLRX1                     GLRX
#> KRTCAP3                KRTCAP3
#> STC2                      STC2
#> VGLL1                    VGLL1
#> SLC25A371             SLC25A37
#> ADAMTS9                ADAMTS9
#> PPFIA11                 PPFIA1
#> ABLIM1                  ABLIM1
#> IRX3                      IRX3
#> CCT6A                    CCT6A
#> SUSD2                    SUSD2
#> ST14                      ST14
#> GPRC5A1                 GPRC5A
#> SPON21                   SPON2
#> FSTL1                    FSTL1
#> TM4SF1-AS11         TM4SF1-AS1
#> EGFL7                    EGFL7
#> BARX1                    BARX1
#> SLC9A7                  SLC9A7
#> RASGRP1                RASGRP1
#> PRSS23                  PRSS23
#> HLA-C1                   HLA-C
#> ACTG2                    ACTG2
#> CTSV                      CTSV
#> C3orf38                C3orf38
#> KLF13                    KLF13
#> LGALSL                  LGALSL
#> PROX11                   PROX1
#> LINC00707            LINC00707
#> RAI14                    RAI14
#> PROM11                   PROM1
#> MTURN                    MTURN
#> CLUL1                    CLUL1
#> YES1                      YES1
#> SCD1                       SCD
#> CENPV                    CENPV
#> KLHDC31                 KLHDC3
#> CRNDE                    CRNDE
#> LOXL1                    LOXL1
#> ALCAM                    ALCAM
#> CAPN21                   CAPN2
#> TMEM154                TMEM154
#> MSLN                      MSLN
#> S100A8                  S100A8
#> TCF7L1                  TCF7L1
#> FKBP10                  FKBP10
#> CDC25B1                 CDC25B
#> TMPRSS3                TMPRSS3
#> MIR4435-2HG        MIR4435-2HG
#> ELF51                     ELF5
#> HOPX                      HOPX
#> CELF4                    CELF4
#> HSPA6                    HSPA6
#> SLC6A8                  SLC6A8
#> RP11-554I8.2      RP11-554I8.2
#> SOX41                     SOX4
#> SDF2L1                  SDF2L1
#> VTCN11                   VTCN1
#> CXADR                    CXADR
#> MIR210HG              MIR210HG
#> PCAT19                  PCAT19
#> THBS2                    THBS2
#> PHYH                      PHYH
#> SLPI1                     SLPI
#> COL6A1                  COL6A1
#> RAB25                    RAB25
#> DBNDD1                  DBNDD1
#> VASN1                     VASN
#> SQRDL1                   SQRDL
#> CASC15                  CASC15
#> UFD1L                    UFD1L
#> PDP1                      PDP1
#> PDGFA                    PDGFA
#> KLHL35                  KLHL35
#> INHBB                    INHBB
#> PDLIM41                 PDLIM4
#> MESP2                    MESP2
#> IFIT1                    IFIT1
#> NET11                     NET1
#> RP11-400K9.4      RP11-400K9.4
#> WEE1                      WEE1
#> DAPP11                   DAPP1
#> CLDN1                    CLDN1
#> CAMK2N1                CAMK2N1
#> MMP71                     MMP7
#> RARRES11               RARRES1
#> MB1                         MB
#> MEX3A                    MEX3A
#> S100A101               S100A10
#> PODXL2                  PODXL2
#> AC005152.31         AC005152.3
#> PLXNA2                  PLXNA2
#> NEDD9                    NEDD9
#> IFITM11                 IFITM1
#> ART3                      ART3
#> AARD                      AARD
#> OCLN                      OCLN
#> CLDN41                   CLDN4
#> HSP90AB11             HSP90AB1
#> GLYATL2                GLYATL2
#> SYCE1L                  SYCE1L
#> TEKT31                   TEKT3
#> TP53BP2                TP53BP2
#> ORM2                      ORM2
#> TINCR1                   TINCR
#> CELF21                   CELF2
#> FASN                      FASN
#> IFT1721                 IFT172
#> CSF3R1                   CSF3R
#> PSCA                      PSCA
#> S100P1                   S100P
#> SMTN1                     SMTN
#> SLC2A4RG              SLC2A4RG
#> CAPG1                     CAPG
#> DNTTIP11               DNTTIP1
#> RANBP1                  RANBP1
#> HLA-B1                   HLA-B
#> SELM1                     SELM
#> ISG201                   ISG20
#> NCCRP11                 NCCRP1
#> CSTB1                     CSTB
#> SMOC1                    SMOC1
#> PSTPIP2                PSTPIP2
#> NPR31                     NPR3
#> ANKRD371               ANKRD37
#> S100A92                 S100A9
#> LINC00152            LINC00152
#> MBD2                      MBD2
#> CHAF1B                  CHAF1B
#> CHI3L21                 CHI3L2
#> CAMK11                   CAMK1
#> PLOD21                   PLOD2
#> MTL5                      MTL5
#> EHD1                      EHD1
#> RDH101                   RDH10
#> PAM1                       PAM
#> HES11                     HES1
#> NR2F21                   NR2F2
#> LMO41                     LMO4
#> MLPH                      MLPH
#> PRSS81                   PRSS8
#> RCAN1                    RCAN1
#> AIF1L                    AIF1L
#> TTYH1                    TTYH1
#> NANOS1                  NANOS1
#> CTA-293F17.1      CTA-293F17.1
#> SDC41                     SDC4
#> SAC3D1                  SAC3D1
#> FADS2                    FADS2
#> KRT181                   KRT18
#> ELN                        ELN
#> GCNT1                    GCNT1
#> BACE21                   BACE2
#> THEM6                    THEM6
#> TPM11                     TPM1
#> ZBTB101                 ZBTB10
#> PMAIP1                  PMAIP1
#> GRB10                    GRB10
#> FH                          FH
#> BNIP3                    BNIP3
#> MAP2K3                  MAP2K3
#> SNX8                      SNX8
#> SCARB1                  SCARB1
#> TMEM25                  TMEM25
#> ZG16B                    ZG16B
#> IRAK1                    IRAK1
#> ADM                        ADM
#> POLD2                    POLD2
#> RPS3                      RPS3
#> HMGB3                    HMGB3
#> CKB                        CKB
#> ALDH1B11               ALDH1B1
#> MESP11                   MESP1
#> HLA-A1                   HLA-A
#> NFIB                      NFIB
#> BCL9L                    BCL9L
#> HES4                      HES4
#> MX11                       MX1
#> COL2A1                  COL2A1
#> TUBA4A                  TUBA4A
#> RNASET21               RNASET2
#> CREB5                    CREB5
#> KCTD11                   KCTD1
#> SMYD2                    SMYD2
#> GAL                        GAL
#> CENPQ                    CENPQ
#> KLRG2                    KLRG2
#> LY6E                      LY6E
#> PHLDA3                  PHLDA3
#> KLF101                   KLF10
#> ZNF83                    ZNF83
#> PRR15L                  PRR15L
#> TNFAIP21               TNFAIP2
#> PTPRS                    PTPRS
#> TMEM139                TMEM139
#> TJP3                      TJP3
#> RBBP7                    RBBP7
#> MBP                        MBP
#> SPEG                      SPEG
#> HLA-F1                   HLA-F
#> TM4SF11                 TM4SF1
#> MGP1                       MGP
#> PHKG1                    PHKG1
#> SELENBP11             SELENBP1
#> PPA1                      PPA1
#> TRIB11                   TRIB1
#> ADAM15                  ADAM15
#> TRIB3                    TRIB3
#> CP1                         CP
#> RAB6B                    RAB6B
#> CMSS1                    CMSS1
#> PAK1IP1                PAK1IP1
#> MAP7D3                  MAP7D3
#> SPHK1                    SPHK1
#> SCCPDH1                 SCCPDH
#> PITX1                    PITX1
#> MFAP21                   MFAP2
#> TNFRSF12A            TNFRSF12A
#> NFATC11                 NFATC1
#> RRAGD                    RRAGD
#> MAFF1                     MAFF
#> SOX81                     SOX8
#> ADHFE1                  ADHFE1
#> PDIA4                    PDIA4
#> USP18                    USP18
#> RAB11FIP1            RAB11FIP1
#> KRT82                     KRT8
#> PTP4A3                  PTP4A3
#> PMP221                   PMP22
#> PRSS221                 PRSS22
#> AQP5                      AQP5
#> QDPR                      QDPR
#> C1orf116              C1orf116
#> PLEKHB1                PLEKHB1
#> RAB30-AS1            RAB30-AS1
#> NUCKS1                  NUCKS1
#> NNMT1                     NNMT
#> RGS101                   RGS10
#> B4GALT1                B4GALT1
#> C10orf101             C10orf10
#> TTF2                      TTF2
#> COL4A2                  COL4A2
#> RMI2                      RMI2
#> HES6                      HES6
#> C1orf561               C1orf56
#> SERPINH1              SERPINH1
#> RAB3IP1                 RAB3IP
#> TMEM791                 TMEM79
#> MGLL1                     MGLL
#> PFKP1                     PFKP
#> PHLDA2                  PHLDA2
#> GRB141                   GRB14
#> TMEM106C1             TMEM106C
#> MARCKSL11             MARCKSL1
#> AZGP11                   AZGP1
#> IDH11                     IDH1
#> CLPSL1                  CLPSL1
#> IGFBP2                  IGFBP2
#> FRMD3                    FRMD3
#> PDGFRA1                 PDGFRA
#> SDC21                     SDC2
#> PCOLCE21               PCOLCE2
#> AQP51                     AQP5
#> PRSS331                 PRSS33
#> COL2A11                 COL2A1
#> LEFTY2                  LEFTY2
#> SSRP11                   SSRP1
#> HIBCH1                   HIBCH
#> MGST11                   MGST1
#> LNX11                     LNX1
#> MIA1                       MIA
#> PYCR11                   PYCR1
#> SYT81                     SYT8
#> SCRG11                   SCRG1
#> SLC9A3R21             SLC9A3R2
#> DBI1                       DBI
#> SOHLH11                 SOHLH1
#> S100B1                   S100B
#> TNFRSF211             TNFRSF21
#> FKBP41                   FKBP4
#> ROPN1B1                 ROPN1B
#> ARL6IP11               ARL6IP1
#> PEG101                   PEG10
#> HIST1H2AE1           HIST1H2AE
#> TBC1D11                 TBC1D1
#> CYP39A11               CYP39A1
#> EFHD1                    EFHD1
#> SERTAD4-AS11       SERTAD4-AS1
#> S100A11                 S100A1
#> ACTA21                   ACTA2
#> CTNND21                 CTNND2
#> HSP90AB12             HSP90AB1
#> RAC3                      RAC3
#> SERTAD41               SERTAD4
#> COL11A11               COL11A1
#> MYL91                     MYL9
#> LAPTM4B1               LAPTM4B
#> PBX11                     PBX1
#> TPM21                     TPM2
#> HAPLN1                  HAPLN1
#> IDH1-AS1              IDH1-AS1
#> SNHG251                 SNHG25
#> TIMM101                 TIMM10
#> EXTL1                    EXTL1
#> FBXL22                  FBXL22
#> TMX21                     TMX2
#> GDPD2                    GDPD2
#> C1QL41                   C1QL4
#> SYNGR11                 SYNGR1
#> GLS1                       GLS
#> KRT83                     KRT8
#> CSRP11                   CSRP1
#> ACTA1                    ACTA1
#> CA6                        CA6
#> PLEKHB11               PLEKHB1
#> SFRP11                   SFRP1
#> CKS1B1                   CKS1B
#> STMN11                   STMN1
#> FAM3C1                   FAM3C
#> NT5DC2                  NT5DC2
#> WFDC1                    WFDC1
#> ANGPT1                  ANGPT1
#> FABP7                    FABP7
#> DKK1                      DKK1
#> PPP1R12A              PPP1R12A
#> CA81                       CA8
#> ITGA101                 ITGA10
#> LIMCH11                 LIMCH1
#> PLOD31                   PLOD3
#> HRCT11                   HRCT1
#> CPED11                   CPED1
#> ROPN11                   ROPN1
#> NFIB1                     NFIB
#> C16orf74              C16orf74
#> HIST1H2BG1           HIST1H2BG
#> TAGLN1                   TAGLN
#> RP11-357H14.171 RP11-357H14.17
#> SPDEF                    SPDEF
#> PHGDH                    PHGDH
#> HACD11                   HACD1
#> PDIA41                   PDIA4
#> SNHG19                  SNHG19
#> NEGR11                   NEGR1
#> CRABP1                  CRABP1
#> SH3BGR1                 SH3BGR
#> MFAP22                   MFAP2
#> PPP1R1B1               PPP1R1B
#> MYOZ1                    MYOZ1
#> ENPP51                   ENPP5
#> LDHB                      LDHB
#> DTNB1                     DTNB
#> TUBB2B1                 TUBB2B
#> ITGA6                    ITGA6
#> IGFBP5                  IGFBP5
#> AIF1L1                   AIF1L
#> DSC31                     DSC3
#> RGCC                      RGCC
#> CACYBP1                 CACYBP
#> NDRG21                   NDRG2
#> SERPINH11             SERPINH1
#> AP1M21                   AP1M2
#> MATN31                   MATN3
#> NME1                      NME1
#> TTYH11                   TTYH1
#> MT1E                      MT1E
#> MYO1B1                   MYO1B
#> SUN31                     SUN3
#> DNM3OS1                 DNM3OS
#> PMP222                   PMP22
#> ANO1                      ANO1
#> MYLK1                     MYLK
#> HSPA1A1                 HSPA1A
#> CTGF1                     CTGF
#> VDR                        VDR
#> TUBB2A1                 TUBB2A
#> CKB1                       CKB
#> WIF1                      WIF1
#> SLC29A11               SLC29A1
#> HSPD1                    HSPD1
#> DNPH1                    DNPH1
#> ACTL8                    ACTL8
#> C1orf115              C1orf115
#> CMBL                      CMBL
#> FREM2                    FREM2
#> PODXL21                 PODXL2
#> CLU1                       CLU
#> NANOS11                 NANOS1
#> HIST2H2BE1           HIST2H2BE
#> PTGIS                    PTGIS
#> CPM                        CPM
#> MSRB31                   MSRB3
#> SERPINE21             SERPINE2
#> CCT51                     CCT5
#> EPHX11                   EPHX1
#> HMGB2                    HMGB2
#> KRT182                   KRT18
#> TSPAN21                 TSPAN2
#> NES                        NES
#> SYCP2                    SYCP2
#> NQO11                     NQO1
#> PFN2                      PFN2
#> VANGL11                 VANGL1
#> BAMBI1                   BAMBI
#> LAMB31                   LAMB3
#> SERPINA51             SERPINA5
#> CTHRC11                 CTHRC1
#> GOLM11                   GOLM1
#> IDI11                     IDI1
#> SOD31                     SOD3
#> IGFBP71                 IGFBP7
#> COL11A21               COL11A2
#> AARD1                     AARD
#> PRPS21                   PRPS2
#> HSPA1B1                 HSPA1B
#> TSPAN121               TSPAN12
#> PITX11                   PITX1
#> CRISPLD11             CRISPLD1
#> GTF3A                    GTF3A
#> HMGN2                    HMGN2
#> LRRC73                  LRRC73
#> METTL7A1               METTL7A
#> PRR7                      PRR7
#> HIST1H1C1             HIST1H1C
#> HMGB11                   HMGB1
#> XAGE2                    XAGE2
#> GOLT1A                  GOLT1A
#> RP11-89K21.1      RP11-89K21.1
#> NCAN                      NCAN
#> TUBB1                     TUBB
#> RP11-25K19.11     RP11-25K19.1
#> CITED41                 CITED4
#> TUBA1B                  TUBA1B
#> ST141                     ST14
#> GPM6B1                   GPM6B
#> MDFI1                     MDFI
#> DCXR                      DCXR
#> RAMP2                    RAMP2
#> PCBD1                    PCBD1
#> KCNQ1OT12             KCNQ1OT1
#> TEKT32                   TEKT3
#> CKS21                     CKS2
#> ERVMER34-1          ERVMER34-1
#> SMOC11                   SMOC1
#> NTHL1                    NTHL1
#> PRSS82                   PRSS8
#> KLK11                    KLK11
#> FXYD61                   FXYD6
#> TUBA1A1                 TUBA1A
#> TMEM204                TMEM204
#> VGF                        VGF
#> TNFSF13B1             TNFSF13B
#> TPM12                     TPM1
#> SQLE1                     SQLE
#> PHYH1                     PHYH
#> HEY21                     HEY2
#> COPZ2                    COPZ2
#> PLK2                      PLK2
#> TOX1                       TOX
#> SLBP1                     SLBP
#> DBNDD11                 DBNDD1
#> FBXO321                 FBXO32
#> IL17B1                   IL17B
#> TOB11                     TOB1
#> ITIH6                    ITIH6
#> KCNMB11                 KCNMB1
#> HN1                        HN1
#> RASL12                  RASL12
#> RUVBL1                  RUVBL1
#> SUSD5                    SUSD5
#> PTS1                       PTS
#> COL9A32                 COL9A3
#> MAPK13                  MAPK13
#> PPIL1                    PPIL1
#> P3H4                      P3H4
#> ELOVL5                  ELOVL5
#> LTBP1                    LTBP1
#> SNX221                   SNX22
#> GAS1                      GAS1
#> ACTL6A1                 ACTL6A
#> CTNNB11                 CTNNB1
#> LRRCC1                  LRRCC1
#> PAQR4                    PAQR4
#> VSNL1                    VSNL1
#> GAMT                      GAMT
#> PCP4                      PCP4
#> TMSB15A1               TMSB15A
#> FHL11                     FHL1
#> IMPA21                   IMPA2
#> SNAI21                   SNAI2
#> PRSS222                 PRSS22
#> C6orf141              C6orf141
#> NET12                     NET1
#> IFRD11                   IFRD1
#> DNM31                     DNM3
#> TPD52L11               TPD52L1
#> SCGB3A11               SCGB3A1
#> MLLT111                 MLLT11
#> NUDT41                   NUDT4
#> FRMD4A1                 FRMD4A
#> CRELD2                  CRELD2
#> ADAM151                 ADAM15
#> SNRNP25                SNRNP25
#> HSPH11                   HSPH1
#> BOP1                      BOP1
#> VASN2                     VASN
#> DCTPP1                  DCTPP1
#> B3GNT71                 B3GNT7
#> PALLD                    PALLD
#> SENCR                    SENCR
#> LIPH                      LIPH
#> SLC43A31               SLC43A3
#> HIST1H4E1             HIST1H4E
#> ELN1                       ELN
#> KLHDC32                 KLHDC3
#> SLC2A4RG1             SLC2A4RG
#> NPPC                      NPPC
#> PTRF1                     PTRF
#> CP2                         CP
#> NKD21                     NKD2
#> PPP1R14A              PPP1R14A
#> NREP                      NREP
#> NSG11                     NSG1
#> PLOD22                   PLOD2
#> CYR61                    CYR61
#> RP11-798M19.6    RP11-798M19.6
#> ODC1                      ODC1
#> SERPINB52             SERPINB5
#> TMEM1391               TMEM139
#> SOX82                     SOX8
#> CRNDE1                   CRNDE
#> DLX51                     DLX5
#> RAD211                   RAD21
#> POSTN1                   POSTN
#> FBLN1                    FBLN1
#> ATP6V0E2              ATP6V0E2
#> UCHL1                    UCHL1
#> LRP2                      LRP2
#> C9orf40                C9orf40
#> BGN1                       BGN
#> CA2                        CA2
#> SFN1                       SFN
#> MEST                      MEST
#> MT1G1                     MT1G
#> ZG16B1                   ZG16B
#> GJA1                      GJA1
#> CALD11                   CALD1
#> PTP4A31                 PTP4A3
#> PRKDC                    PRKDC
#> HSPB11                   HSPB1
#> SORBS22                 SORBS2
#> EPHX21                   EPHX2
#> NEXN                      NEXN
#> MT2A2                     MT2A
#> MAPK8IP21             MAPK8IP2
#> PRNP                      PRNP
#> KRT72                     KRT7
#> FAM89A1                 FAM89A
#> SLC43A1                SLC43A1
#> CD320                    CD320
#> CDH31                     CDH3
#> TMEM67                  TMEM67
#> TNC                        TNC
#> QPCT1                     QPCT
#> LY6E1                     LY6E
#> QPRT                      QPRT
#> RTN4RL21               RTN4RL2
#> UNG                        UNG
#> FAM46B1                 FAM46B
#> SLC26A7                SLC26A7
#> SMC1B                    SMC1B
#> HYLS1                    HYLS1
#> RP11-19E11.11     RP11-19E11.1
#> HIST1H2BJ            HIST1H2BJ
#> CDC42EP11             CDC42EP1
#> PRELP1                   PRELP
#> TUBA4A1                 TUBA4A
#> GCNT11                   GCNT1
#> SEPT4                    SEPT4
#> TMEM1581               TMEM158
#> C12orf75              C12orf75
#> GGH                        GGH
#> SCIN1                     SCIN
#> ASNS                      ASNS
#> PPDPF1                   PPDPF
#> MT1F1                     MT1F
#> UACA                      UACA
#> SBSPON1                 SBSPON
#> HIST1H2BN            HIST1H2BN
#> GAS61                     GAS6
#> HSPA51                   HSPA5
#> SCARB11                 SCARB1
#> C1orf1161             C1orf116
#> CELF41                   CELF4
#> KRT151                   KRT15
#> ELF52                     ELF5
#> RP1-313I6.12      RP1-313I6.12
#> TMEM611                 TMEM61
#> BYSL                      BYSL
#> YBX2                      YBX2
#> SKA2                      SKA2
#> CDCA7L                  CDCA7L
#> TSPAN51                 TSPAN5
#> TMEM97                  TMEM97
#> C2orf822               C2orf82
#> MBNL1-AS11           MBNL1-AS1
#> DGAT2                    DGAT2
#> SIGMAR1                SIGMAR1
#> ETV51                     ETV5
#> GMPR1                     GMPR
#> TUBB4B1                 TUBB4B
#> TUBA1C1                 TUBA1C
#> BOC1                       BOC
#> BARX11                   BARX1
#> PAICS                    PAICS
#> MESP12                   MESP1
#> CCND1                    CCND1
#> SAP30                    SAP30
#> MCM7                      MCM7
#> LOXL2                    LOXL2
#> GADD45G                GADD45G
#> TMEM792                 TMEM79
#> TSEN15                  TSEN15
#> S100A102               S100A10
#> MINCR                    MINCR
#> CTNNAL1                CTNNAL1
#> NUCKS11                 NUCKS1
#> NRM                        NRM
#> DHTKD1                  DHTKD1
#> MTL51                     MTL5
#> ANXA2R1                 ANXA2R
#> CHI3L11                 CHI3L1
#> CNTNAP3B1             CNTNAP3B
#> RRS1                      RRS1
#> FBLN22                   FBLN2
#> TRAF3IP31             TRAF3IP3
#> FKBP101                 FKBP10
#> TMEM45A                TMEM45A
#> SMOC2                    SMOC2
#> MT1X1                     MT1X
#> KCTD13                   KCTD1
#> DKC1                      DKC1
#> ID4                        ID4
#> CEP70                    CEP70
#> PEG32                     PEG3
#> BACE22                   BACE2
#> H2AFV                    H2AFV
#> PDZK1IP11             PDZK1IP1
#> RAB251                   RAB25
#> ALYREF                  ALYREF
#> VTCN12                   VTCN1
#> GCSH1                     GCSH
#> THBS11                   THBS1
#> KRT861                   KRT86
#> ACTG21                   ACTG2
#> TINAGL1                TINAGL1
#> PLA2G4A1               PLA2G4A
#> SDC42                     SDC4
#> TRIB31                   TRIB3
#> GPSM2                    GPSM2
#> FRZB                      FRZB
#> KLRG21                   KLRG2
#> C8orf461               C8orf46
#> KRT811                   KRT81
#> KRTCAP31               KRTCAP3
#> BARD11                   BARD1
#> TEX30                    TEX30
#> MEX3A1                   MEX3A
#> TUBB6                    TUBB6
#> CDK41                     CDK4
#> CDH11                     CDH1
#> MYBL11                   MYBL1
#> COL6A11                 COL6A1
#> CERCAM                  CERCAM
#> RUVBL2                  RUVBL2
#> A2M1                       A2M
#> HMGA11                   HMGA1
#> AC005152.32         AC005152.3
#> FGF1                      FGF1
#> HS3ST11                 HS3ST1
#> DNAJC9                  DNAJC9
#> PGP                        PGP
#> NEDD91                   NEDD9
#> PDLIM1                  PDLIM1
#> SELENBP12             SELENBP1
#> SDC12                     SDC1
#> IRX31                     IRX3
#> SLC12A2                SLC12A2
#> ANXA11                   ANXA1
#> BNIP31                   BNIP3
#> LBR                        LBR
#> PIR1                       PIR
#> WDR341                   WDR34
#> SMTN2                     SMTN
#> COL4A21                 COL4A2
#> KIF22                    KIF22
#> ANP32E                  ANP32E
#> PPA11                     PPA1
#> THEM61                   THEM6
#> ATAD2                    ATAD2
#> TOMM40                  TOMM40
#> FH1                         FH
#> CTD-3065J16.9    CTD-3065J16.9
#> GAPDH1                   GAPDH
#> MFSD31                   MFSD3
#> TCF7L11                 TCF7L1
#> ZBTB102                 ZBTB10
#> CTNNBIP1              CTNNBIP1
#> LYAR                      LYAR
#> MAP1B1                   MAP1B
#> HILPDA1                 HILPDA
#> FXYD32                   FXYD3
#> CLDN61                   CLDN6
#> SMIM221                 SMIM22
#> PCAT191                 PCAT19
#> MAGIX1                   MAGIX
#> RAB6B1                   RAB6B
#> PQBP11                   PQBP1
#> ILF21                     ILF2
#> KIF1A                    KIF1A
#> CT831                     CT83
#> PLP21                     PLP2
#> ORM1                      ORM1
#> CRABP21                 CRABP2
#> TIMM17B1               TIMM17B
#> CDK161                   CDK16
#> UXT1                       UXT
#> GAL1                       GAL
#> TRPM81                   TRPM8
#> NDUFB111               NDUFB11
#> MESP21                   MESP2
#> SLC34A2                SLC34A2
#> PRRX21                   PRRX2
#> EBP1                       EBP
#> GPKOW1                   GPKOW
#> KCNE51                   KCNE5
#> SLC35A21               SLC35A2
#> CENPV1                   CENPV
#> WDR451                   WDR45
#> WDR131                   WDR13
#> PGBD51                   PGBD5
#> RBM31                     RBM3
#> FTSJ11                   FTSJ1
#> UBA11                     UBA1
#> MAGEA4                  MAGEA4
#> C5orf66-AS1        C5orf66-AS1
#> RPS31                     RPS3
#> SOX181                   SOX18
#> VTCN13                   VTCN1
#> HDAC61                   HDAC6
#> SLC2A111               SLC2A11
#> IRX41                     IRX4
#> SLC6A141               SLC6A14
#> TIMM8B                  TIMM8B
#> NFASC1                   NFASC
#> PDXK1                     PDXK
#> OSR2                      OSR2
#> MISP1                     MISP
#> DPM11                     DPM1
#> CALB21                   CALB2
#> CALML51                 CALML5
#> RHOV1                     RHOV
#> PORCN1                   PORCN
#> ALCAM1                   ALCAM
#> SLC9A71                 SLC9A7
#> MORN31                   MORN3
#> RGN1                       RGN
#> SYCE1L1                 SYCE1L
#> PTGFRN1                 PTGFRN
#> RP11-817J15.3    RP11-817J15.3
#> PVRL11                   PVRL1
#> VGLL11                   VGLL1
#> MTL52                     MTL5
#> OTUD51                   OTUD5
#> ARMC11                   ARMC1
#> CLDN11                   CLDN1
#> HES41                     HES4
#> ORM21                     ORM2
#> ATP9A1                   ATP9A
#> LRP51                     LRP5
#> GS1-124K5.41       GS1-124K5.4
#> ADCYAP1                ADCYAP1
#> ELK11                     ELK1
#> HMGB31                   HMGB3
#> TTC9                      TTC9
#> RGS201                   RGS20
#> NPDC11                   NPDC1
#> CRYAB2                   CRYAB
#> ART31                     ART3
#> S100A141               S100A14
#> ACTG22                   ACTG2
#> CSF1                      CSF1
#> UBASH3B                UBASH3B
#> CFB1                       CFB
#> RP11-554I8.21     RP11-554I8.2
#> RP11-817J15.2    RP11-817J15.2
#> STC21                     STC2
#> EFNA5                    EFNA5
#> ARAF1                     ARAF
#> UBE2T                    UBE2T
#> TM4SF1-AS12         TM4SF1-AS1
#> PTX3                      PTX3
#> CLDN71                   CLDN7
#> S100A161               S100A16
#> PPP1R14B1             PPP1R14B
#> LA16c-380H5.51   LA16c-380H5.5
#> RBBP71                   RBBP7
#> CHST9                    CHST9
#> PRRT3-AS11           PRRT3-AS1
#> GPR271                   GPR27
#> CAGE1                    CAGE1
#> C1orf61                C1orf61
#> RAB30-AS11           RAB30-AS1
#> ARHGDIG                ARHGDIG
#> HYOU11                   HYOU1
#> IL322                     IL32
#> TBC1D251               TBC1D25
#> ATP1B11                 ATP1B1
#> TMEM2081               TMEM208
#> IFITM32                 IFITM3
#> LINC01133            LINC01133
#> MANEAL1                 MANEAL
#> ZNF674-AS11         ZNF674-AS1
#> KRT191                   KRT19
#> MFSD32                   MFSD3
#> GINS4                    GINS4
#> TMPRSS31               TMPRSS3
#> CDH12                     CDH1
#> PRAF21                   PRAF2
#> FDX1                      FDX1
#> HEBP22                   HEBP2
#> AARD2                     AARD
#> SEC131                   SEC13
#> RRP361                   RRP36
#> DDR11                     DDR1
#> BRK11                     BRK1
#> NUCKS12                 NUCKS1
#> TLDC11                   TLDC1
#> SPINT12                 SPINT1
#> RP11-490M8.1      RP11-490M8.1
#> PPP1R14C1             PPP1R14C
#> CEP57                    CEP57
#> METRN1                   METRN
#> CD242                     CD24
#> RBP11                     RBP1
#> TCF7L12                 TCF7L1
#> PPDPF2                   PPDPF
#> KRTCAP32               KRTCAP3
#> FSTL11                   FSTL1
#> PSAT1                    PSAT1
#> TM4SF12                 TM4SF1
#> LYNX1                    LYNX1
#> LSR1                       LSR
#> FBXO91                   FBXO9
#> MLLT41                   MLLT4
#> STARD101               STARD10
#> ANXA31                   ANXA3
#> PHGDH1                   PHGDH
#> HSPB21                   HSPB2
#> KRT231                   KRT23
#> POLD21                   POLD2
#> GRB101                   GRB10
#> KRBOX41                 KRBOX4
#> ECEL11                   ECEL1
#> TMEM251                 TMEM25
#> CCDC1671               CCDC167
#> KCNG11                   KCNG1
#> EEF1A2                  EEF1A2
#> SEMA6A1                 SEMA6A
#> LSM51                     LSM5
#> APOD1                     APOD
#> RANBP11                 RANBP1
#> KLK5                      KLK5
#> SERPINH12             SERPINH1
#> RPUSD31                 RPUSD3
#> CENPW                    CENPW
#> LCN22                     LCN2
#> LY6E2                     LY6E
#> DNAAF11                 DNAAF1
#> MESP13                   MESP1
#> MAD2L1BP1             MAD2L1BP
#> EGLN11                   EGLN1
#> IRX32                     IRX3
#> FADS21                   FADS2
#> KLK61                     KLK6
#> PDIA31                   PDIA3
#> SLC7A51                 SLC7A5
#> DSP1                       DSP
#> H2AFX                    H2AFX
#> MYBL2                    MYBL2
#> CAPN22                   CAPN2
#> PITX12                   PITX1
#> SCD5                      SCD5
#> CEBPB1                   CEBPB
#> VHL1                       VHL
#> CLDN9                    CLDN9
#> PPM1H                    PPM1H
#> RAI141                   RAI14
#> SECTM11                 SECTM1
#> DCUN1D5                DCUN1D5
#> LEMD12                   LEMD1
#> CENPF                    CENPF
#> PVRL41                   PVRL4
#> FAM92A11               FAM92A1
#> MMP151                   MMP15
#> SYNCRIP1               SYNCRIP
#> FASN1                     FASN
#> PPFIA12                 PPFIA1
#> DNAJC152               DNAJC15
#> ARHGAP291             ARHGAP29
#> SIPA1L21               SIPA1L2
#> TYMS                      TYMS
#> ZPR1                      ZPR1
#> TSEN151                 TSEN15
#> KNSTRN                  KNSTRN
#> ATP2A1-AS1          ATP2A1-AS1
#> PRICKLE11             PRICKLE1
#> CCT6A1                   CCT6A
#> CX3CL1                  CX3CL1
#> TTF21                     TTF2
#> QDPR1                     QDPR
#> TMEM200A1             TMEM200A
#> TMEM63A1               TMEM63A
#> XAGE21                   XAGE2
#> BARX12                   BARX1
#> CDH32                     CDH3
#> BCL9L1                   BCL9L
#> SNCG                      SNCG
#> GALNT21                 GALNT2
#> TM4SF18                TM4SF18
#> HMGA12                   HMGA1
#> TMEM132A1             TMEM132A
#> NLRP2                    NLRP2
#> PCSK1N1                 PCSK1N
#> SUSD21                   SUSD2
#> RGS102                   RGS10
#> GCNT2                    GCNT2
#> TNFRSF12A1           TNFRSF12A
#> OAZ31                     OAZ3
#> NRBP21                   NRBP2
#> CACNA1A1               CACNA1A
#> RASGRP11               RASGRP1
#> OGG11                     OGG1
#> NOA11                     NOA1
#> SUN11                     SUN1
#> WNT5B                    WNT5B
#> NUPR21                   NUPR2
#> DUT                        DUT
#> UFD1L1                   UFD1L
#> NT5E1                     NT5E
#> PRRT31                   PRRT3
#> PSORS1C11             PSORS1C1
#> RP11-273G15.2    RP11-273G15.2
#> AGT2                       AGT
#> C11orf73              C11orf73
#> CMSS11                   CMSS1
#> FERMT1                  FERMT1
#> PDIA61                   PDIA6
#> CUL71                     CUL7
#> MIS18A                  MIS18A
#> GABRP1                   GABRP
#> ARHGEF12              ARHGEF12
#> CXCL10                  CXCL10
#> SUSD4                    SUSD4
#> SELM2                     SELM
#> IRAK11                   IRAK1
#> LINC009601           LINC00960
#> NMB1                       NMB
#> NAT8L1                   NAT8L
#> NTN11                     NTN1
#> DEK                        DEK
#> APP2                       APP
#> H2AFY21                 H2AFY2
#> ELF32                     ELF3
#> PADI21                   PADI2
#> HES61                     HES6
#> PLK1                      PLK1
#> ZNF2171                 ZNF217
#> LSM4                      LSM4
#> GUCY1A3                GUCY1A3
#> EPPK11                   EPPK1
#> RRAGD1                   RRAGD
#> ECT2                      ECT2
#> LDOC11                   LDOC1
#> WEE11                     WEE1
#> SMYD21                   SMYD2
#> PSCA1                     PSCA
#> UBE2C                    UBE2C
#> CDC20                    CDC20
#> ITGB81                   ITGB8
#> CLDN32                   CLDN3
#> QPRT1                     QPRT
#> KLHL351                 KLHL35
#> PGP1                       PGP
#> PPA12                     PPA1
#> TJP31                     TJP3
#> IRX51                     IRX5
#> SORL11                   SORL1
#> PTP4A32                 PTP4A3
#> SETD51                   SETD5
#> CD551                     CD55
#> THUMPD3-AS11       THUMPD3-AS1
#> KIFC1                    KIFC1
#> PMAIP11                 PMAIP1
#> CENPA                    CENPA
#> THUMPD31               THUMPD3
#> FAM64A                  FAM64A
#> PXDN1                     PXDN
#> H2AFV1                   H2AFV
#> NOTCH3                  NOTCH3
#> CDKN3                    CDKN3
#> PLA2G16                PLA2G16
#> RIPK2                    RIPK2
#> DAPK2                    DAPK2
#> FANCD2                  FANCD2
#> PTPRF2                   PTPRF
#> DNPH11                   DNPH1
#> C1QBP                    C1QBP
#> NUF2                      NUF2
#> KLRG22                   KLRG2
#> SNX81                     SNX8
#> FH2                         FH
#> EIF3B1                   EIF3B
#> MUC11                     MUC1
#> CRNDE2                   CRNDE
#> CCND11                   CCND1
#> MUC5B1                   MUC5B
#> KRT812                   KRT81
#> GATA6-AS11           GATA6-AS1
#> LPIN11                   LPIN1
#> PTS2                       PTS
#> BCAS41                   BCAS4
#> GBA1                       GBA
#> TMPRSS131             TMPRSS13
#> BIRC5                    BIRC5
#> SLC1A31                 SLC1A3
#> FZD81                     FZD8
#> MND1                      MND1
#> LAD11                     LAD1
#> SEZ6L21                 SEZ6L2
#> DIAPH3                  DIAPH3
#> CHAF1B1                 CHAF1B
#> LINC007071           LINC00707
#> TACSTD22               TACSTD2
#> CDCA3                    CDCA3
#> FGF132                   FGF13
#> APMAP1                   APMAP
#> CHORDC1                CHORDC1
#> KRT141                   KRT14
#> YES11                     YES1
#> BYSL1                     BYSL
#> GLYATL21               GLYATL2
#> CARHSP11               CARHSP1
#> TARBP11                 TARBP1
#> EMC31                     EMC3
#> YDJC                      YDJC
#> MIR4435-2HG1       MIR4435-2HG
#> MCAM1                     MCAM
#> LMNB1                    LMNB1
#> TUNAR1                   TUNAR
#> ZNF695                  ZNF695
#> GINS2                    GINS2
#> BRPF11                   BRPF1
#> VEGFA1                   VEGFA
#> RCCD1                    RCCD1
#> TJP11                     TJP1
#> RNF81                     RNF8
#> TRIP13                  TRIP13
#> PRC1                      PRC1
#> TIMP11                   TIMP1
#> SDF2L11                 SDF2L1
#> IL17RC1                 IL17RC
#> GPC11                     GPC1
#> RAB252                   RAB25
#> SEC11C1                 SEC11C
#> ASPM                      ASPM
#> SDR16C51               SDR16C5
#> GGT5                      GGT5
#> MTMR141                 MTMR14
#> LTF                        LTF
#> CHST71                   CHST7
#> S100A21                 S100A2
#> SAPCD2                  SAPCD2
#> ZNRF21                   ZNRF2
#> AURKB                    AURKB
#> PTPN141                 PTPN14
#> OCLN1                     OCLN
#> NEK2                      NEK2
#> CP3                         CP
#> CEP55                    CEP55
#> CTA-293F17.11     CTA-293F17.1
#> ANP32E1                 ANP32E
#> FZD7                      FZD7
#> RP11-410L14.2    RP11-410L14.2
#> RP11-783K16.51   RP11-783K16.5
#> FKBP102                 FKBP10
#> SLC6A81                 SLC6A8
#> RAD51AP1              RAD51AP1
#> RIC31                     RIC3
#> ORC6                      ORC6
#> CRACR2B                CRACR2B
#> DHCR24                  DHCR24
#> GGCT1                     GGCT
#> NETO2                    NETO2
#> CENPN                    CENPN
#> IGSF23                  IGSF23
#> DGAT21                   DGAT2
#> KIAA0101              KIAA0101
#> SLC12A21               SLC12A2
#> MLPH1                     MLPH
#> TACC3                    TACC3
#> C16orf59              C16orf59
#> FAM129A                FAM129A
#> PBK                        PBK
#> CENPQ1                   CENPQ
#> FAM208B1               FAM208B
#> POC1A                    POC1A
#> LGALSL1                 LGALSL
#> MAD2L1                  MAD2L1
#> SLC38A1                SLC38A1
#> MTURN1                   MTURN
#> NPTXR1                   NPTXR
#> TTLL31                   TTLL3
#> RP11-19E11.12     RP11-19E11.1
#> APOBEC3B              APOBEC3B
#> LDLR                      LDLR
#> KDELC2                  KDELC2
#> CCNA2                    CCNA2
#> FLNB1                     FLNB
#> CENPU                    CENPU
#> MBD21                     MBD2
#> PROM12                   PROM1
#> TOP2A                    TOP2A
#> FLNA1                     FLNA
#> CCDC34                  CCDC34
#> NFIB2                     NFIB
#> BSPRY                    BSPRY
#> SAC3D11                 SAC3D1
#> GRXCR11                 GRXCR1
#> C6orf151               C6orf15
#> PAK1IP11               PAK1IP1
#> TBC1D71                 TBC1D7
#> NR2F22                   NR2F2
#> PIR2                       PIR
#> KLHDC33                 KLHDC3
#> SMC4                      SMC4
#> FOXM1                    FOXM1
#> RECQL4                  RECQL4
#> CD3201                   CD320
#> THEM62                   THEM6
#> E2F1                      E2F1
#> SMC2                      SMC2
#> MAP3K81                 MAP3K8
#> TPBG1                     TPBG
#> FKBP91                   FKBP9
#> ACTR3B1                 ACTR3B
#> CHPF1                     CHPF
#> SLC5A6                  SLC5A6
#> EGLN31                   EGLN3
#> EFNA12                   EFNA1
#> NUSAP1                  NUSAP1
#> COTL11                   COTL1
#> ARPC42                   ARPC4
#> SLC2A121               SLC2A12
#> KRT6B1                   KRT6B
#> EDN11                     EDN1
#> JAGN11                   JAGN1
#> C21orf58              C21orf58
#> ADHFE11                 ADHFE1
#> NBL11                     NBL1
#> DHFR                      DHFR
#> FBXL161                 FBXL16
#> LINC001521           LINC00152
#> NDC80                    NDC80
#> TROAP                    TROAP
#> TPX2                      TPX2
#> RDH102                   RDH10
#> CTD-3065J16.91   CTD-3065J16.9
#> MYO5B1                   MYO5B
#> TK1                        TK1
#> PSIP1                    PSIP1
#> RUVBL11                 RUVBL1
#> RMI21                     RMI2
#> PRR15L1                 PRR15L
#> CENPH                    CENPH
#> MEX3A2                   MEX3A
#> PCNA                      PCNA
#> ACOT7                    ACOT7
#> DSCC1                    DSCC1
#> PRSS83                   PRSS8
#> RARRES12               RARRES1
#> IGSF31                   IGSF3
#> CHMP2B1                 CHMP2B
#> MAFK1                     MAFK
#> MCM3                      MCM3
#> COL9A2                  COL9A2
#> CCNB1                    CCNB1
#> SLC25A51               SLC25A5
#> IRF2BP21               IRF2BP2
#> TRIB21                   TRIB2
#> SCCPDH2                 SCCPDH
#> ZG16B2                   ZG16B
#> C1S2                       C1S
#> KLF131                   KLF13
#> GMNN                      GMNN
#> MKI67                    MKI67
#> NUDT1                    NUDT1
#> RPA3                      RPA3
#> MELK                      MELK
#> PRR71                     PRR7
#> IFT1722                 IFT172
#> CYR611                   CYR61
#> RNASEH2A              RNASEH2A
#> CHAF1A                  CHAF1A
#> CYP27A11               CYP27A1
#> PKP11                     PKP1
#> EGFL71                   EGFL7
#> PARPBP                  PARPBP
#> PRKDC1                   PRKDC
#> PFN21                     PFN2
#> OBP2B                    OBP2B
#> CGGBP11                 CGGBP1
#> ATP6V0E21             ATP6V0E2
#> TADA31                   TADA3
#> TINAGL11               TINAGL1
#> DCTPP11                 DCTPP1
#> C6orf1322             C6orf132
#> LRP21                     LRP2
#> TRIB12                   TRIB1
#> GATA61                   GATA6
#> PKMYT1                  PKMYT1
#> CCNB2                    CCNB2
#> KRT171                   KRT17
#> MAP21                     MAP2
#> GAN1                       GAN
#> ASNS1                     ASNS
#> PDZD21                   PDZD2
#> MFI21                     MFI2
#> RCAN2                    RCAN2
#> HELLS                    HELLS
#> GGH1                       GGH
#> HSPA52                   HSPA5
#> SHISA91                 SHISA9
#> ANLN                      ANLN
#> MYC                        MYC
#> PRSS223                 PRSS22
#> TUBB61                   TUBB6
#> PSTPIP21               PSTPIP2
#> PCDHB91                 PCDHB9
#> PTTG1                    PTTG1
#> H2AFZ                    H2AFZ
#> MTHFD1                  MTHFD1
#> SPON22                   SPON2
#> PODXL22                 PODXL2
#> KRT73                     KRT7
#> CLMN1                     CLMN
#> HOMER21                 HOMER2
#> DTYMK                    DTYMK
#> MCM2                      MCM2
#> RCAN11                   RCAN1
#> ADAM152                 ADAM15
#> TNFRSF18              TNFRSF18
#> LINC016151           LINC01615
#> FANCI                    FANCI
#> SAA12                     SAA1
#> KCNN4                    KCNN4
#> FEN1                      FEN1
#> SPC25                    SPC25
#> FADS1                    FADS1
#> NDUFAF61               NDUFAF6
#> RAD212                   RAD21
#> SCD2                       SCD
#> DBNDD12                 DBNDD1
#> NANOS12                 NANOS1
#> RPP25                    RPP25
#> LAPTM4B2               LAPTM4B
#> CLCN42                   CLCN4
#> PAICS1                   PAICS
#> ABLIM11                 ABLIM1
#> MAL22                     MAL2
#> RASD11                   RASD1
#> MINCR1                   MINCR
#> PHLDA21                 PHLDA2
#> TMPO                      TMPO
#> CDT1                      CDT1
#> RAB11FIP11           RAB11FIP1
#> NPR32                     NPR3
#> LRR1                      LRR1
#> PHF19                    PHF19
#> SCGB3A12               SCGB3A1
#> YBX21                     YBX2
#> CRABP11                 CRABP1
#> ZWINT                    ZWINT
#> MCM6                      MCM6
#> MDC1                      MDC1
#> VANGL12                 VANGL1
#> ATAD21                   ATAD2
#> ERBB32                   ERBB3
#> LMNB2                    LMNB2
#> CDCA7L1                 CDCA7L
#> C31                         C3
#> CELF42                   CELF4
#> CCDC64B1               CCDC64B
#> TRIB32                   TRIB3
#> SAA21                     SAA2
#> SORBS23                 SORBS2
#> SGOL1                    SGOL1
#> CXADR1                   CXADR
#> LOXL11                   LOXL1
#> C1GALT11               C1GALT1
#> CSTB2                     CSTB
#> TCF19                    TCF19
#> RP11-400K9.41     RP11-400K9.4
#> BASP11                   BASP1
#> CD821                     CD82
#> BID1                       BID
#> TFAP2A2                 TFAP2A
#> ELN2                       ELN
#> KRT183                   KRT18
#> CEP701                   CEP70
#> EHD11                     EHD1
#> CDK1                      CDK1
#> SELENBP13             SELENBP1
#> PLXNA21                 PLXNA2
#> C19orf48              C19orf48
#> PHYH2                     PHYH
#> SOX42                     SOX4
#> TNF                        TNF
#> CAMK1D1                 CAMK1D
#> BACE23                   BACE2
#> SOX92                     SOX9
#> PLEKHB12               PLEKHB1
#> RRS11                     RRS1
#> TP53BP21               TP53BP2
#> TINCR2                   TINCR
#> COL4A22                 COL4A2
#> PDGFA1                   PDGFA
#> ARID5A1                 ARID5A
#> IRF2BPL1               IRF2BPL
#> RTP41                     RTP4
#> KRT161                   KRT16
#> MYO62                     MYO6
#> FADS31                   FADS3
#> HIST1H2BJ1           HIST1H2BJ
#> ACSL12                   ACSL1
#> MSLN1                     MSLN
#> KIAA05131             KIAA0513
#> CDKN2D                  CDKN2D
#> CENPM                    CENPM
#> CACYBP2                 CACYBP
#> MAP2K31                 MAP2K3
#> TMEM971                 TMEM97
#> THBS12                   THBS1
#> ADAMTS91               ADAMTS9
#> PEG102                   PEG10
#> RUNX31                   RUNX3
#> KIF20B                  KIF20B
#> SIGMAR11               SIGMAR1
#> RRM1                      RRM1
#> ARFGEF31               ARFGEF3
#> PTPRS1                   PTPRS
#> CTSF1                     CTSF
#> MTSS1L1                 MTSS1L
#> HIST1H2BN1           HIST1H2BN
#> CAMK12                   CAMK1
#> MBP1                       MBP
#> MCM4                      MCM4
#> CDCA4                    CDCA4
#> UBE2S                    UBE2S
#> ST6GAL11               ST6GAL1
#> SLPI2                     SLPI
#> HES12                     HES1
#> WDR342                   WDR34
#> GRB142                   GRB14
#> TUBB2                     TUBB
#> C21                         C2
#> SULF21                   SULF2
#> SLC2A4RG2             SLC2A4RG
#> RFC3                      RFC3
#> CKAP5                    CKAP5
#> SPDL1                    SPDL1
#> TPT1-AS11             TPT1-AS1
#> HLA-B2                   HLA-B
#> SMOC21                   SMOC2
#> TTC39A1                 TTC39A
#> LIMD21                   LIMD2
#> RACGAP1                RACGAP1
#> NEBL1                     NEBL
#> CKB2                       CKB
#> JHDM1D-AS11         JHDM1D-AS1
#> HIST1H4C1             HIST1H4C
#> LYAR1                     LYAR
#> MCM71                     MCM7
#> OAF                        OAF
#> NASP                      NASP
#> CCT52                     CCT5
#> VASN3                     VASN
#> EFHD11                   EFHD1
#> THOP1                    THOP1
#> MTHFD21                 MTHFD2
#> PHACTR11               PHACTR1
#> GSDMC2                   GSDMC
#> BOP11                     BOP1
#> CRLF11                   CRLF1
#> MAPK131                 MAPK13
#> ZNF6541                 ZNF654
#> S100A81                 S100A8
#> PSRC1                    PSRC1
#> MB2                         MB
#> FRZB1                     FRZB
#> CORO1A1                 CORO1A
#> HOPX1                     HOPX
#> CDK42                     CDK4
#> INHBB1                   INHBB
#> HLA-F2                   HLA-F
#> MARC1                    MARC1
#> SPHK11                   SPHK1
#> KIF221                   KIF22
#> NTHL11                   NTHL1
#> TOMM401                 TOMM40
#> IL17RE1                 IL17RE
#> C10orf102             C10orf10
#> NEAT12                   NEAT1
#> NDRG12                   NDRG1
#> DHTKD11                 DHTKD1
#> DONSON                  DONSON
#> MGLL2                     MGLL
#> C1R2                       C1R
#> SLC39A141             SLC39A14
#> NCAPD2                  NCAPD2
#> PDIA42                   PDIA4
#> RFC4                      RFC4
#> PRNP1                     PRNP
#> TPD52L12               TPD52L1
#> TPM13                     TPM1
#> DKC11                     DKC1
#> FBLN23                   FBLN2
#> GAPDH2                   GAPDH
#> PDP11                     PDP1
#> PDZK1IP12             PDZK1IP1
#> PHLDA31                 PHLDA3
#> HSP90AB13             HSP90AB1
#> ZNF528                  ZNF528
#> CRELD21                 CRELD2
#> SPEG1                     SPEG
#> FSCN11                   FSCN1
#> RPL39L1                 RPL39L
#> ALYREF1                 ALYREF
#> PLSCR12                 PLSCR1
#> NAALADL21             NAALADL2
#> AIF1L2                   AIF1L
#> LBR1                       LBR
#> RAB3IP2                 RAB3IP
#> GAS11                     GAS1
#> CXCL161                 CXCL16
#> C3orf381               C3orf38
#> S100A62                 S100A6
#> GPSM21                   GPSM2
#> CASP41                   CASP4
#> CREB51                   CREB5
#> EPAS1                    EPAS1
#> MBNL1-AS12           MBNL1-AS1
#> PPIL11                   PPIL1
#> B4GALT11               B4GALT1
#> TUBB4B2                 TUBB4B
#> ITGA61                   ITGA6
#> MARS                      MARS
#> CTSV1                     CTSV
#> CTNNBIP11             CTNNBIP1
#> SAP301                   SAP30
#> CHEK1                    CHEK1
#> CKAP2                    CKAP2
#> FGF11                     FGF1
#> BRCA2                    BRCA2
#> CLDN42                   CLDN4
#> C8orf462               C8orf46
#> UNG1                       UNG
#> TOX2                       TOX
#> SLC26A71               SLC26A7
#> SBSPON2                 SBSPON
#> COL11A22               COL11A2
#> S100P2                   S100P
#> PLA2G4A2               PLA2G4A
#> CNGA1                    CNGA1
#> DNM32                     DNM3
#> MATN32                   MATN3
#> SFN2                       SFN
#> SOD32                     SOD3
#> FXYD62                   FXYD6
#> RP1-27K12.2        RP1-27K12.2
#> EPHX12                   EPHX1
#> C2orf80                C2orf80
#> CDKN2A1                 CDKN2A
#> MDFI2                     MDFI
#> RAB253                   RAB25
#> RP11-25K19.12     RP11-25K19.1
#> CA82                       CA8
#> SLC29A12               SLC29A1
#> BOC2                       BOC
#> S100A12                 S100A1
#> HIST1H2AE2           HIST1H2AE
#> ISLR                      ISLR
#> ROPN12                   ROPN1
#> GJA11                     GJA1
#> CRISPLD12             CRISPLD1
#> SCRG12                   SCRG1
#> MYOZ11                   MYOZ1
#> GLS2                       GLS
#> TSPAN122               TSPAN12
#> TMEM793                 TMEM79
#> C1orf1162             C1orf116
#> DBI2                       DBI
#> LIPH1                     LIPH
#> CCDC64B2               CCDC64B
#> GPM6B2                   GPM6B
#> COL11A12               COL11A1
#> TSPAN52                 TSPAN5
#> HIBCH2                   HIBCH
#> MGST12                   MGST1
#> DMKN                      DMKN
#> CLDN33                   CLDN3
#> SH3BGR2                 SH3BGR
#> MARCKSL12             MARCKSL1
#> FBXO22                   FBXO2
#> CLU2                       CLU
#> GCNT12                   GCNT1
#> DIO2                      DIO2
#> HSPB12                   HSPB1
#> TTYH12                   TTYH1
#> AZGP12                   AZGP1
#> SLC43A32               SLC43A3
#> IGFBP21                 IGFBP2
#> HIST1H1C2             HIST1H1C
#> RP11-89K21.11     RP11-89K21.1
#> C1QL42                   C1QL4
#> KLK12                     KLK1
#> NET13                     NET1
#> CLDN43                   CLDN4
#> ROPN1B2                 ROPN1B
#> PAX11                     PAX1
#> TUBB2B2                 TUBB2B
#> IDI12                     IDI1
#> HIST1H2BG2           HIST1H2BG
#> GYLTL1B                GYLTL1B
#> DNM3OS2                 DNM3OS
#> SLC9A3R22             SLC9A3R2
#> RP11-268P4.5      RP11-268P4.5
#> C1orf1861             C1orf186
#> METTL7A2               METTL7A
#> MOG                        MOG
#> EXTL11                   EXTL1
#> PEG33                     PEG3
#> TAGLN2                   TAGLN
#> TIFA1                     TIFA
#> FXYD33                   FXYD3
#> SNX222                   SNX22
#> PDZK1IP13             PDZK1IP1
#> TUBB2A2                 TUBB2A
#> TUBA1A2                 TUBA1A
#> TOB12                     TOB1
#> GOLT1A1                 GOLT1A
#> C8orf463               C8orf46
#> BSPRY1                   BSPRY
#> LIMCH12                 LIMCH1
#> ENHO                      ENHO
#> PLEKHB13               PLEKHB1
#> TUBB4B3                 TUBB4B
#> FHL12                     FHL1
#> NAAA                      NAAA
#> TSPAN22                 TSPAN2
#> MYL92                     MYL9
#> MFAP23                   MFAP2
#> NGF                        NGF
#> IDH12                     IDH1
#> FBN2                      FBN2
#> MIA2                       MIA
#> PAM2                       PAM
#> IGFBP51                 IGFBP5
#> SERTAD42               SERTAD4
#> NCCRP12                 NCCRP1
#> VSNL11                   VSNL1
#> PFN22                     PFN2
#> PRSS332                 PRSS33
#> CYP26B11               CYP26B1
#> HIST2H2BE2           HIST2H2BE
#> NDRG22                   NDRG2
#> MUC15                    MUC15
#> SMTN3                     SMTN
#> CTSV2                     CTSV
#> NPPC1                     NPPC
#> GOLM12                   GOLM1
#> TACSTD23               TACSTD2
#> PRR15L2                 PRR15L
#> RIC32                     RIC3
#> BAMBI2                   BAMBI
#> LINC01048            LINC01048
#> SDC13                     SDC1
#> MYLK2                     MYLK
#> S100A103               S100A10
#> RTN4RL22               RTN4RL2
#> DTNB2                     DTNB
#> CAPN23                   CAPN2
#> ERBB33                   ERBB3
#> FAM84A1                 FAM84A
#> LRRC731                 LRRC73
#> PRSS84                   PRSS8
#> ENPP52                   ENPP5
#> SOHLH12                 SOHLH1
#> UCHL11                   UCHL1
#> SYT82                     SYT8
#> KCNQ1OT13             KCNQ1OT1
#> IFRD12                   IFRD1
#> HSPA1B2                 HSPA1B
#> WDR343                   WDR34
#> CITED42                 CITED4
#> PTGS21                   PTGS2
#> ALDH1B12               ALDH1B1
#> MBNL1-AS13           MBNL1-AS1
#> SORBS24                 SORBS2
#> PTGIS1                   PTGIS
#> RBP4                      RBP4
#> PMP223                   PMP22
#> ELF53                     ELF5
#> CYP39A12               CYP39A1
#> GCSH2                     GCSH
#> FZD51                     FZD5
#> CRABP12                 CRABP1
#> HSPA1A2                 HSPA1A
#> SERTAD4-AS12       SERTAD4-AS1
#> GATA31                   GATA3
#> S100B2                   S100B
#> PXDNL                    PXDNL
#> SCIN2                     SCIN
#> ARL6IP12               ARL6IP1
#> JHDM1D-AS12         JHDM1D-AS1
#> TIMM102                 TIMM10
#> TMEM1392               TMEM139
#> MFSD61                   MFSD6
#> LNX12                     LNX1
#> SMIM5                    SMIM5
#> FBXL162                 FBXL16
#> FOXC12                   FOXC1
#> SERPINB53             SERPINB5
#> ANXA2R2                 ANXA2R
#> LMTK3                    LMTK3
#> IL17B2                   IL17B
#> ADAM153                 ADAM15
#> SUN32                     SUN3
#> SENCR1                   SENCR
#> WFDC21                   WFDC2
#> C1orf1151             C1orf115
#> TMEM1582               TMEM158
#> PYCR12                   PYCR1
#> AP1M22                   AP1M2
#> SYNGR12                 SYNGR1
#> TMEM106C2             TMEM106C
#> KRT74                     KRT7
#> DEFB11                   DEFB1
#> PALLD1                   PALLD
#> TP53BP22               TP53BP2
#> MYBL12                   MYBL1
#> TEKT33                   TEKT3
#> MB3                         MB
#> TPM22                     TPM2
#> RHOBTB31               RHOBTB3
#> CLIC3                    CLIC3
#> TMX22                     TMX2
#> MAPK132                 MAPK13
#> POSTN2                   POSTN
#> OVOL12                   OVOL1
#> GRB143                   GRB14
#> KRT84                     KRT8
#> NUDT42                   NUDT4
#> RDH103                   RDH10
#> CA11                      CA11
#> CKS1B2                   CKS1B
#> SFRP12                   SFRP1
#> PDGFRA2                 PDGFRA
#> FBXO322                 FBXO32
#> MAP1B2                   MAP1B
#> ERVMER34-11         ERVMER34-1
#> C2orf823               C2orf82
#> CLDN72                   CLDN7
#> NUPR22                   NUPR2
#> ID41                       ID4
#> CHI3L12                 CHI3L1
#> RAB3IP3                 RAB3IP
#> HRCT12                   HRCT1
#> CRYAB3                   CRYAB
#> LEFTY21                 LEFTY2
#> SSRP12                   SSRP1
#> OPRK11                   OPRK1
#> AQP52                     AQP5
#> SOX83                     SOX8
#> MMP152                   MMP15
#> HEY22                     HEY2
#> RP11-554I8.22     RP11-554I8.2
#> PHLDA11                 PHLDA1
#> WNT7B                    WNT7B
#> IFI62                     IFI6
#> FAM3C2                   FAM3C
#> TUBB62                   TUBB6
#> PTRF2                     PTRF
#> AC005152.33         AC005152.3
#> TMEM612                 TMEM61
#> C6orf1323             C6orf132
#> PIM11                     PIM1
#> FAM46B2                 FAM46B
#> CEBPD1                   CEBPD
#> MLLT112                 MLLT11
#> SNHG252                 SNHG25
#> MGP2                       MGP
#> GADD45G1               GADD45G
#> RPL39L2                 RPL39L
#> TNFSF13B2             TNFSF13B
#> B3GNT72                 B3GNT7
#> BGN2                       BGN
#> CEP702                   CEP70
#> FRZB2                     FRZB
#> TINCR3                   TINCR
#> HN11                       HN1
#> DSC32                     DSC3
#> VGF1                       VGF
#> CSRP12                   CSRP1
#> PVRL42                   PVRL4
#> ZG16B3                   ZG16B
#> MARC11                   MARC1
#> NEXN1                     NEXN
#> CD243                     CD24
#> RP11-357H14.172 RP11-357H14.17
#> ARHGAP292             ARHGAP29
#> QPCT2                     QPCT
#> FBXL221                 FBXL22
#> VANGL13                 VANGL1
#> AIF1L3                   AIF1L
#> SAT12                     SAT1
#> CDC42EP12             CDC42EP1
#> IMPA22                   IMPA2
#> MMP72                     MMP7
#> PCBD11                   PCBD1
#> HIST1H4E2             HIST1H4E
#> DKK11                     DKK1
#> PCOLCE22               PCOLCE2
#> OCLN2                     OCLN
#> LA16c-380H5.52   LA16c-380H5.5
#> XAGE22                   XAGE2
#> GMPR2                     GMPR
#> STMN12                   STMN1
#> MUC12                     MUC1
#> PBX12                     PBX1
#> RP11-798M19.61   RP11-798M19.6
#> NEGR12                   NEGR1
#> MSRB32                   MSRB3
#> CTNND22                 CTNND2
#> KRT232                   KRT23
#> RAD213                   RAD21
#> LGALS32                 LGALS3
#> LSR2                       LSR
#> HAPLN11                 HAPLN1
#> BARD12                   BARD1
#> DSEL1                     DSEL
#> SNHG191                 SNHG19
#> ITGB42                   ITGB4
#> TNFRSF12A2           TNFRSF12A
#> NFIB3                     NFIB
#> STAT11                   STAT1
#> KCNMB12                 KCNMB1
#> ACTG23                   ACTG2
#> INAFM1                  INAFM1
#> KRT813                   KRT81
#> DLX52                     DLX5
#> AARD3                     AARD
#> SYCP21                   SYCP2
#> TPD52L13               TPD52L1
#> KRT192                   KRT19
#> EFHD12                   EFHD1
#> NSG12                     NSG1
#> ACTA22                   ACTA2
#> VASN4                     VASN
#> MTSS1L2                 MTSS1L
#> LAYN                      LAYN
#> IL17RE2                 IL17RE
#> CKS22                     CKS2
#> NQO12                     NQO1
#> TBC1D12                 TBC1D1
#> MESP14                   MESP1
#> MAL23                     MAL2
#> LIMA1                    LIMA1
#> LINC007072           LINC00707
#> PHYH3                     PHYH
#> NAB11                     NAB1
#> PPP1R14A1             PPP1R14A
#> SDC43                     SDC4
#> LEMD13                   LEMD1
#> FKBP42                   FKBP4
#> HACD12                   HACD1
#> CPM1                       CPM
#> OVOS21                   OVOS2
#> RGS161                   RGS16
#> SLC52A11               SLC52A1
#> TTC39A2                 TTC39A
#> SPEG2                     SPEG
#> GABRP2                   GABRP
#> HILPDA2                 HILPDA
#> MEST1                     MEST
#> MEX3A3                   MEX3A
#> ADHFE12                 ADHFE1
#> GPX11                     GPX1
#> NES1                       NES
#> DOK51                     DOK5
#> DSN1                      DSN1
#> CNTNAP3B2             CNTNAP3B
#> VGLL12                   VGLL1
#> SQLE2                     SQLE
#> CYR612                   CYR61
#> PMEPA1                  PMEPA1
#> ITGA102                 ITGA10
#> IRX33                     IRX3
#> CD822                     CD82
#> SLBP2                     SLBP
#> SOX43                     SOX4
#> NR2F23                   NR2F2
#> LRP22                     LRP2
#> GAS62                     GAS6
#> NME11                     NME1
#> EFNA51                   EFNA5
#> RAMP21                   RAMP2
#> RAC31                     RAC3
#> ELF33                     ELF3
#> RP11-19E11.13     RP11-19E11.1
#> DGAT22                   DGAT2
#> FAM89A2                 FAM89A
#> SLPI3                     SLPI
#> TUBA4A2                 TUBA4A
#> SHISA92                 SHISA9
#> PAQR41                   PAQR4
#> RP11-161M6.21     RP11-161M6.2
#> EPHX22                   EPHX2
#> GGCT2                     GGCT
#> AC022007.51         AC022007.5
#> HLA-C2                   HLA-C
#> P3H41                     P3H4
#> LTBP11                   LTBP1
#> LAPTM4B3               LAPTM4B
#> NANOS13                 NANOS1
#> HSP90AB14             HSP90AB1
#> ABCA51                   ABCA5
#> FKBP103                 FKBP10
#> HES13                     HES1
#> CACYBP3                 CACYBP
#> ACTL6A2                 ACTL6A
#> NEBL2                     NEBL
#> DSP2                       DSP
#> CHI3L22                 CHI3L2
#> SLC2A122               SLC2A12
#> CLUL11                   CLUL1
#> ZNF5281                 ZNF528
#> CLMN2                     CLMN
#> CIART1                   CIART
#> TFAP2A3                 TFAP2A
#> ST142                     ST14
#> NDRG13                   NDRG1
#> CARHSP12               CARHSP1
#> CALD12                   CALD1
#> TUBA1C2                 TUBA1C
#> KLRG23                   KLRG2
#> LINC014361           LINC01436
#> TUBB3                     TUBB
#> CKB3                       CKB
#> C1orf562               C1orf56
#> KRT184                   KRT18
#> CRNDE3                   CRNDE
#> MAP22                     MAP2
#> NKD22                     NKD2
#> PODXL23                 PODXL2
#> QPRT2                     QPRT
#> RGCC1                     RGCC
#> PRELP2                   PRELP
#> PROM13                   PROM1
#> CDH111                   CDH11
#> CLPSL11                 CLPSL1
#> IGFBP72                 IGFBP7
#> IGKC1                     IGKC
#> ATP1B12                 ATP1B1
#> C6orf1411             C6orf141
#> LAMB11                   LAMB1
#> APOBEC3B1             APOBEC3B
#> SOX93                     SOX9
#> MITF1                     MITF
#> CFH                        CFH
#> NTM                        NTM
#> PDGFRB                  PDGFRB
#> FBN1                      FBN1
#> ECM2                      ECM2
#> CHN1                      CHN1
#> FAP                        FAP
#> COL5A1                  COL5A1
#> ITGA1                    ITGA1
#> COL5A2                  COL5A2
#> LAMA4                    LAMA4
#> EMILIN1                EMILIN1
#> SULF1                    SULF1
#> SGIP1                    SGIP1
#> WNT2                      WNT2
#> COL10A1                COL10A1
#> VCAN                      VCAN
#> VCAM1                    VCAM1
#> NDN                        NDN
#> COL6A3                  COL6A3
#> PLAC9                    PLAC9
#> DCN                        DCN
#> DKK3                      DKK3
#> PODN                      PODN
#> LRRC32                  LRRC32
#> CAV1                      CAV1
#> NID2                      NID2
#> ITGBL1                  ITGBL1
#> ITGA11                  ITGA11
#> COX7A1                  COX7A1
#> THY1                      THY1
#> GNG11                    GNG11
#> EDNRA                    EDNRA
#> GPX8                      GPX8
#> MFAP5                    MFAP5
#> TFPI                      TFPI
#> SPARCL1                SPARCL1
#> PLAT                      PLAT
#> HEG1                      HEG1
#> ANGPTL2                ANGPTL2
#> RCN3                      RCN3
#> ADAMTS12              ADAMTS12
#> PRRX1                    PRRX1
#> CPXM1                    CPXM1
#> TIMP3                    TIMP3
#> ASPN                      ASPN
#> AEBP1                    AEBP1
#> C11orf96              C11orf96
#> ENPEP                    ENPEP
#> EBF1                      EBF1
#> OLFML3                  OLFML3
#> DPYSL3                  DPYSL3
#> WISP2                    WISP2
#> CD248                    CD248
#> CREB3L1                CREB3L1
#> CYGB                      CYGB
#> KCNE4                    KCNE4
#> RARRES2                RARRES2
#> SFRP2                    SFRP2
#> SFRP4                    SFRP4
#> COL15A1                COL15A1
#> COL3A1                  COL3A1
#> PLXDC1                  PLXDC1
#> SPRY1                    SPRY1
#> COL5A3                  COL5A3
#> ADAMTS2                ADAMTS2
#> FIBIN                    FIBIN
#> EDIL3                    EDIL3
#> SEMA5A                  SEMA5A
#> COL12A1                COL12A1
#> LHFP                      LHFP
#> ZFHX4                    ZFHX4
#> MEG3                      MEG3
#> TNFAIP6                TNFAIP6
#> FOXS1                    FOXS1
#> HEYL                      HEYL
#> IFI27                    IFI27
#> ADGRF5                  ADGRF5
#> SSPN                      SSPN
#> HIC1                      HIC1
#> NOX4                      NOX4
#> LUM                        LUM
#> NID1                      NID1
#> BICC1                    BICC1
#> OLFML2B1               OLFML2B
#> FILIP1                  FILIP1
#> COL1A1                  COL1A1
#> PRSS231                 PRSS23
#> FRMD6                    FRMD6
#> IGFBP4                  IGFBP4
#> FILIP1L1               FILIP1L
#> ECM11                     ECM1
#> COL8A1                  COL8A1
#> COL1A2                  COL1A2
#> AKAP12                  AKAP12
#> PCOLCE                  PCOLCE
#> MXRA5                    MXRA5
#> MAP1A                    MAP1A
#> ANGPT2                  ANGPT2
#> GJA4                      GJA4
#> CCDC102B              CCDC102B
#> COL18A1                COL18A1
#> INHBA                    INHBA
#> GGT51                     GGT5
#> ABCC9                    ABCC9
#> ANTXR1                  ANTXR1
#> VGLL3                    VGLL3
#> LRRC15                  LRRC15
#> EGFLAM                  EGFLAM
#> P3H3                      P3H3
#> EFEMP1                  EFEMP1
#> LGALS3BP1             LGALS3BP
#> NOTCH31                 NOTCH3
#> GUCY1A2                GUCY1A2
#> CTSK1                     CTSK
#> GUCY1A31               GUCY1A3
#> AXL                        AXL
#> C3orf80                C3orf80
#> SERPING1              SERPING1
#> CEMIP                    CEMIP
#> SYTL2                    SYTL2
#> COL4A1                  COL4A1
#> ID31                       ID3
#> EID11                     EID1
#> C1QTNF6                C1QTNF6
#> SPARC1                   SPARC
#> CRISPLD2              CRISPLD2
#> POSTN3                   POSTN
#> SERPINF11             SERPINF1
#> PDE1A                    PDE1A
#> CPE                        CPE
#> CDH112                   CDH11
#> UNC5B                    UNC5B
#> FMO1                      FMO1
#> ACKR4                    ACKR4
#> CERCAM1                 CERCAM
#> FN11                       FN1
#> PDGFRL                  PDGFRL
#> ZEB1                      ZEB1
#> FSTL12                   FSTL1
#> KCNJ8                    KCNJ8
#> GUCY1B3                GUCY1B3
#> PODNL1                  PODNL1
#> RUNX1T1                RUNX1T1
#> HTRA11                   HTRA1
#> PLPP4                    PLPP4
#> IFI161                   IFI16
#> CLEC11A                CLEC11A
#> SRPX2                    SRPX2
#> CCDC801                 CCDC80
#> RBMS3                    RBMS3
#> PTRF3                     PTRF
#> UACA1                     UACA
#> CAV2                      CAV2
#> COMP                      COMP
#> ADAM12                  ADAM12
#> GJC1                      GJC1
#> C5orf46                C5orf46
#> LIMA11                   LIMA1
#> SMOC22                   SMOC2
#> MMP111                   MMP11
#> TGFB3                    TGFB3
#> EHD2                      EHD2
#> XIST1                     XIST
#> CARMN                    CARMN
#> TSHZ2                    TSHZ2
#> S1PR3                    S1PR3
#> ADAMTS4                ADAMTS4
#> CALD13                   CALD1
#> SYDE1                    SYDE1
#> COL6A21                 COL6A2
#> RARRES3                RARRES3
#> IGFBP73                 IGFBP7
#> FGF7                      FGF7
#> INAFM11                 INAFM1
#> BMP1                      BMP1
#> COL6A12                 COL6A1
#> VAMP51                   VAMP5
#> TCF41                     TCF4
#> SPON1                    SPON1
#> LZTS1                    LZTS1
#> MTUS1                    MTUS1
#> LAMB12                   LAMB1
#> SALRNA2                SALRNA2
#> GLT8D2                  GLT8D2
#> ABCA6                    ABCA6
#> TIMP12                   TIMP1
#> THBS21                   THBS2
#> CFI                        CFI
#> SRPX                      SRPX
#> FOXP11                   FOXP1
#> TMEM119                TMEM119
#> TPM23                     TPM2
#> TBX2                      TBX2
#> ANKRD28                ANKRD28
#> COL4A23                 COL4A2
#> BGN3                       BGN
#> HTRA3                    HTRA3
#> C1R3                       C1R
#> F2R                        F2R
#> SCG5                      SCG5
#> IFITM21                 IFITM2
#> PSMB91                   PSMB9
#> UBB1                       UBB
#> HLA-A2                   HLA-A
#> IFITM12                 IFITM1
#> GRP                        GRP
#> SMIM3                    SMIM3
#> HEY1                      HEY1
#> IFITM33                 IFITM3
#> MIR4435-2HG2       MIR4435-2HG
#> RAMP1                    RAMP1
#> APOL1                    APOL1
#> HSPG2                    HSPG2
#> NRN1                      NRN1
#> IGFBP6                  IGFBP6
#> EPAS11                   EPAS1
#> BST21                     BST2
#> PCDH7                    PCDH7
#> B2M1                       B2M
#> HOXB2                    HOXB2
#> PLAU1                     PLAU
#> HIGD1B                  HIGD1B
#> STEAP1                  STEAP1
#> FNDC1                    FNDC1
#> TNFSF101               TNFSF10
#> NBL12                     NBL1
#> LMCD11                   LMCD1
#> PALLD2                   PALLD
#> TMEM45A1               TMEM45A
#> JAG1                      JAG1
#> NREP1                     NREP
#> OAF1                       OAF
#> COPZ21                   COPZ2
#> EVA1A                    EVA1A
#> HOPX2                     HOPX
#> IFI44L                  IFI44L
#> LXN1                       LXN
#> NNMT2                     NNMT
#> TMEM2041               TMEM204
#> NTRK2                    NTRK2
#> FMOD                      FMOD
#> CYP1B11                 CYP1B1
#> PTN                        PTN
#> NEXN2                     NEXN
#> LOXL12                   LOXL1
#> ENPP2                    ENPP2
#> PLEKHA4                PLEKHA4
#> FKBP11                  FKBP11
#> RAB311                   RAB31
#> PDPN1                     PDPN
#> MARCKS1                 MARCKS
#> SEPT41                   SEPT4
#> GEM                        GEM
#> LINC001522           LINC00152
#> SPOCK1                  SPOCK1
#> SLC2A31                 SLC2A3
#> IFI27L21               IFI27L2
#> OLFML2A                OLFML2A
#> CLEC2B1                 CLEC2B
#> TRIM221                 TRIM22
#> XAF11                     XAF1
#> PTGIR                    PTGIR
#> SPON23                   SPON2
#> C1S3                       C1S
#> ISG151                   ISG15
#> MMP21                     MMP2
#> RFTN1                    RFTN1
#> TPPP3                    TPPP3
#> CXCL121                 CXCL12
#> TMEM47                  TMEM47
#> CD200                    CD200
#> GBP11                     GBP1
#> HLA-B3                   HLA-B
#> HSD17B111             HSD17B11
#> PSMB81                   PSMB8
#> NR2F1-AS1            NR2F1-AS1
#> RP11-394O4.5      RP11-394O4.5
#> FRZB3                     FRZB
#> FAM198B1               FAM198B
#> MSRB33                   MSRB3
#> CST31                     CST3
#> LOXL21                   LOXL2
#> DAB21                     DAB2
#> MYL93                     MYL9
#> MARVELD1              MARVELD1
#> PLA2G161               PLA2G16
#> SELM3                     SELM
#> TPM14                     TPM1
#> PCDH18                  PCDH18
#> MX2                        MX2
#> EPSTI11                 EPSTI1
#> SAMD91                   SAMD9
#> COLEC121               COLEC12
#> DIO21                     DIO2
#> EPB41L21               EPB41L2
#> HES42                     HES4
#> ZEB21                     ZEB2
#> ISLR1                     ISLR
#> FBLN5                    FBLN5
#> C1orf541               C1orf54
#> CAMK2N11               CAMK2N1
#> ADA                        ADA
#> SLC40A11               SLC40A1
#> VIM1                       VIM
#> ACTA23                   ACTA2
#> PPP1R14A2             PPP1R14A
#> IFIT3                    IFIT3
#> HNMT1                     HNMT
#> DUSP11                   DUSP1
#> FBLN11                   FBLN1
#> GAS63                     GAS6
#> SERPINH13             SERPINH1
#> MVP1                       MVP
#> CTHRC12                 CTHRC1
#> LPAR61                   LPAR6
#> FKBP104                 FKBP10
#> EPS81                     EPS8
#> CD401                     CD40
#> FHL13                     FHL1
#> EGR1                      EGR1
#> BASP12                   BASP1
#> GJA12                     GJA1
#> PMEPA11                 PMEPA1
#> PLPP31                   PLPP3
#> LBH                        LBH
#> IFIT11                   IFIT1
#> IFI63                     IFI6
#> ARID5B1                 ARID5B
#> COL11A13               COL11A1
#> GREM1                    GREM1
#> IFI441                   IFI44
#> EMP31                     EMP3
#> TAGLN3                   TAGLN
#> DUT1                       DUT
#> IFI351                   IFI35
#> HSD17B41               HSD17B4
#> TGFB11                   TGFB1
#> PTPRE1                   PTPRE
#> PDLIM11                 PDLIM1
#> OAS2                      OAS2
#> HLA-C3                   HLA-C
#> CCDC85B1               CCDC85B
#> GJB2                      GJB2
#> TGM2                      TGM2
#> NEURL1B                NEURL1B
#> S100A162               S100A16
#> BHLHE411               BHLHE41
#> GADD45A                GADD45A
#> ANO11                     ANO1
#> SAMD9L1                 SAMD9L
#> CCL21                     CCL2
#> SULF22                   SULF2
#> HLA-F3                   HLA-F
#> CHPF2                     CHPF
#> NDUFA4L2              NDUFA4L2
#> PARP141                 PARP14
#> CRIP11                   CRIP1
#> BST1                      BST1
#> SPHK12                   SPHK1
#> TSC22D31               TSC22D3
#> RGS5                      RGS5
#> CTD-3193K9.4      CTD-3193K9.4
#> CTGF2                     CTGF
#> MYO1B2                   MYO1B
#> BICD1                    BICD1
#> FLNA2                     FLNA
#> TGFBI1                   TGFBI
#> IRF71                     IRF7
#> C12orf751             C12orf75
#> CTSB1                     CTSB
#> NR2F24                   NR2F2
#> NRP11                     NRP1
#> TRIB22                   TRIB2
#> MX12                       MX1
#> EMP11                     EMP1
#> ELOVL51                 ELOVL5
#> TAP11                     TAP1
#> SKA21                     SKA2
#> SDC22                     SDC2
#> PDIA32                   PDIA3
#> CYR613                   CYR61
#> TJP12                     TJP1
#> THBS13                   THBS1
#> LY961                     LY96
#> OAS11                     OAS1
#> SLC39A142             SLC39A14
#> RHOBTB32               RHOBTB3
#> PHLDA32                 PHLDA3
#> SNAI22                   SNAI2
#> ID42                       ID4
#> LY6E3                     LY6E
#> CTSL1                     CTSL
#> FAM46A                  FAM46A
#> RGS162                   RGS16
#> KDELC21                 KDELC2
#> STK17A                  STK17A
#> CD552                     CD55
#> CD1091                   CD109
#> MYOF1                     MYOF
#> CYBA1                     CYBA
#> AKR1B11                 AKR1B1
#> WLS                        WLS
#> SUSD22                   SUSD2
#> PGAM1                    PGAM1
#> SEPT61                   SEPT6
#> MEST2                     MEST
#> FCGRT1                   FCGRT
#> CASP42                   CASP4
#> PRDM11                   PRDM1
#> MT2A3                     MT2A
#> TPBG2                     TPBG
#> GBP21                     GBP2
#> PXDN2                     PXDN
#> MEF2C1                   MEF2C
#> LIMD22                   LIMD2
#> LDLR1                     LDLR
#> PSIP11                   PSIP1
#> PRRX22                   PRRX2
#> P3H42                     P3H4
#> SCCPDH3                 SCCPDH
#> DSE1                       DSE
#> PAPSS2                  PAPSS2
#> PDGFA2                   PDGFA
#> PRAF22                   PRAF2
#> MAF1                       MAF
#> GALM1                     GALM
#> NAAA1                     NAAA
#> BCAT11                   BCAT1
#> ARHGEF121             ARHGEF12
#> RARRES13               RARRES1
#> ARID5A2                 ARID5A
#> CCDC28B                CCDC28B
#> ARHGAP181             ARHGAP18
#> LMO42                     LMO4
#> ARL4C1                   ARL4C
#> UBE2L62                 UBE2L6
#> MMP91                     MMP9
#> MAP7D31                 MAP7D3
#> MAP1B3                   MAP1B
#> PHLDA12                 PHLDA1
#> PPP1R12A1             PPP1R12A
#> LDHB1                     LDHB
#> ZFP361                   ZFP36
#> FCER1A                  FCER1A
#> CD1E                      CD1E
#> CD1A                      CD1A
#> CLEC10A                CLEC10A
#> CD1C                      CD1C
#> ASGR2                    ASGR2
#> GAPT                      GAPT
#> ABI31                     ABI3
#> JAML                      JAML
#> IL1RN1                   IL1RN
#> CXorf211               CXorf21
#> MNDA1                     MNDA
#> NCKAP1L1               NCKAP1L
#> ALOX51                   ALOX5
#> HLA-DQA21             HLA-DQA2
#> CRIP12                   CRIP1
#> IL1B1                     IL1B
#> ADAM281                 ADAM28
#> HLA-DRB51             HLA-DRB5
#> PKIB1                     PKIB
#> CLEC4A1                 CLEC4A
#> TNFRSF9                TNFRSF9
#> HLA-DQA11             HLA-DQA1
#> SAMHD11                 SAMHD1
#> HCK1                       HCK
#> CD531                     CD53
#> LST11                     LST1
#> PYCARD1                 PYCARD
#> LCP11                     LCP1
#> NRG1                      NRG1
#> HCLS11                   HCLS1
#> CD481                     CD48
#> FGL21                     FGL2
#> RNASE61                 RNASE6
#> CSF1R1                   CSF1R
#> FCGBP                    FCGBP
#> PTPRC1                   PTPRC
#> CYTIP1                   CYTIP
#> PRKCB                    PRKCB
#> ALOX5AP1               ALOX5AP
#> WAS1                       WAS
#> SPI11                     SPI1
#> MS4A6A1                 MS4A6A
#> CD861                     CD86
#> DOK21                     DOK2
#> RGS11                     RGS1
#> FMNL11                   FMNL1
#> SNX20                    SNX20
#> SRGN1                     SRGN
#> BTK1                       BTK
#> HLA-DRB11             HLA-DRB1
#> HLA-DQB11             HLA-DQB1
#> FCN1                      FCN1
#> GMFG1                     GMFG
#> FCGR2B1                 FCGR2B
#> LILRB21                 LILRB2
#> CSTA1                     CSTA
#> LSP11                     LSP1
#> CYBB1                     CYBB
#> GPR1831                 GPR183
#> MYO1F1                   MYO1F
#> SLA1                       SLA
#> ITGB21                   ITGB2
#> TNFAIP81               TNFAIP8
#> HLA-DPA11             HLA-DPA1
#> PLEK21                   PLEK2
#> AIF11                     AIF1
#> CD831                     CD83
#> HLA-DPB11             HLA-DPB1
#> NAIP1                     NAIP
#> GNA151                   GNA15
#> NFKBID1                 NFKBID
#> PLEK1                     PLEK
#> IL7R1                     IL7R
#> EVI2A1                   EVI2A
#> CLECL11                 CLECL1
#> IFI162                   IFI16
#> PARVG1                   PARVG
#> LAPTM51                 LAPTM5
#> IL2RG1                   IL2RG
#> EID12                     EID1
#> EVI2B1                   EVI2B
#> IL10RA1                 IL10RA
#> RRM2                      RRM2
#> DUSP21                   DUSP2
#> RUNX32                   RUNX3
#> HLA-DMA1               HLA-DMA
#> LPXN1                     LPXN
#> RTN1                      RTN1
#> KYNU1                     KYNU
#> CST7                      CST7
#> TGFBI2                   TGFBI
#> CHST2                    CHST2
#> HLA-DMB1               HLA-DMB
#> EMP32                     EMP3
#> ESCO2                    ESCO2
#> HSD17B112             HSD17B11
#> C10orf541             C10orf54
#> FERMT31                 FERMT3
#> IFI301                   IFI30
#> CORO1A2                 CORO1A
#> XIST2                     XIST
#> IGSF61                   IGSF6
#> UBB2                       UBB
#> CASP11                   CASP1
#> ADAM81                   ADAM8
#> HIST1H1D              HIST1H1D
#> CD300C                  CD300C
#> HLA-DRA1               HLA-DRA
#> FAM111A                FAM111A
#> CD331                     CD33
#> MFNG                      MFNG
#> POU2F21                 POU2F2
#> SERPINB91             SERPINB9
#> LYZ1                       LYZ
#> CD741                     CD74
#> PALD1                    PALD1
#> CCL3L31                 CCL3L3
#> CD371                     CD37
#> CD41                       CD4
#> ARHGDIB1               ARHGDIB
#> NCF41                     NCF4
#> SLC8A11                 SLC8A1
#> TNFAIP31               TNFAIP3
#> AOAH1                     AOAH
#> RAC21                     RAC2
#> BIN2                      BIN2
#> HSD17B42               HSD17B4
#> HLA-DQB2              HLA-DQB2
#> CARD161                 CARD16
#> CTSS1                     CTSS
#> CD300A1                 CD300A
#> PSMB92                   PSMB9
#> MARCH11                 MARCH1
#> EFHD21                   EFHD2
#> UCP21                     UCP2
#> RGS21                     RGS2
#> SLC2A32                 SLC2A3
#> HCST1                     HCST
#> GPR841                   GPR84
#> HLA-DOA1               HLA-DOA
#> MEF2C2                   MEF2C
#> LY962                     LY96
#> CSF2RA1                 CSF2RA
#> SEPT62                   SEPT6
#> LY861                     LY86
#> CD69                      CD69
#> PSMB82                   PSMB8
#> DUSP12                   DUSP1
#> CASC5                    CASC5
#> PTGS11                   PTGS1
#> DSE2                       DSE
#> FCER1G1                 FCER1G
#> PTAFR1                   PTAFR
#> SLC7A71                 SLC7A7
#> CYBA2                     CYBA
#> CLEC5A1                 CLEC5A
#> CENPE                    CENPE
#> BST22                     BST2
#> CD300LF                CD300LF
#> GPR651                   GPR65
#> TBXAS11                 TBXAS1
#> PARM1                    PARM1
#> ADAP21                   ADAP2
#> CD841                     CD84
#> CASP43                   CASP4
#> RASAL3                  RASAL3
#> TYROBP1                 TYROBP
#> RNASET22               RNASET2
#> SAMSN11                 SAMSN1
#> FPR31                     FPR3
#> LGALS91                 LGALS9
#> CENPK                    CENPK
#> CD141                     CD14
#> IL4I11                   IL4I1
#> MS4A71                   MS4A7
#> CKLF1                     CKLF
#> GPR341                   GPR34
#> IFI352                   IFI35
#> CLSPN                    CLSPN
#> PRDM12                   PRDM1
#> FGR1                       FGR
#> GALM2                     GALM
#> PPM1N1                   PPM1N
#> CST6                      CST6
#> SCARF1                  SCARF1
#> MTHFD22                 MTHFD2
#> KIF11                    KIF11
#> SLC16A101             SLC16A10
#> RP11-1143G9.41   RP11-1143G9.4
#> ZFP362                   ZFP36
#> SH2B31                   SH2B3
#> KCTD121                 KCTD12
#> GTSE1                    GTSE1
#> HMMR                      HMMR
#> HLA-A3                   HLA-A
#> FYB1                       FYB
#> CXCR41                   CXCR4
#> SAMD9L2                 SAMD9L
#> C1orf1621             C1orf162
#> SGK11                     SGK1
#> CLEC2B2                 CLEC2B
#> TFEC1                     TFEC
#> TGFB12                   TGFB1
#> CCR11                     CCR1
#> ZEB22                     ZEB2
#> FAM26F1                 FAM26F
#> PLAUR1                   PLAUR
#> CENPM1                   CENPM
#> NR4A21                   NR4A2
#> SERPINF12             SERPINF1
#> DOCK101                 DOCK10
#> F13A1                    F13A1
#> PLAU2                     PLAU
#> OSM1                       OSM
#> CTSC1                     CTSC
#> CD681                     CD68
#> FOXP12                   FOXP1
#> HNMT2                     HNMT
#> CYTH41                   CYTH4
#> LTB                        LTB
#> FBP11                     FBP1
#> CDKN2C                  CDKN2C
#> RPPH1                    RPPH1
#> PDE4B1                   PDE4B
#> B2M2                       B2M
#> LAIR11                   LAIR1
#> LIMS11                   LIMS1
#> C3AR11                   C3AR1
#> ARHGAP182             ARHGAP18
#> SNX101                   SNX10
#> ALDH21                   ALDH2
#> GPAT3                    GPAT3
#> PTPRE2                   PTPRE
#> PTPN221                 PTPN22
#> HJURP                    HJURP
#> TCF42                     TCF4
#> ADAM19                  ADAM19
#> NCAPG                    NCAPG
#> LRRC251                 LRRC25
#> HAVCR21                 HAVCR2
#> LILRB11                 LILRB1
#> HIST1H1E              HIST1H1E
#> ARRB21                   ARRB2
#> MKI671                   MKI67
#> IFI27L22               IFI27L2
#> KIF15                    KIF15
#> ARHGAP11B            ARHGAP11B
#> CPVL1                     CPVL
#> PTPN61                   PTPN6
#> MPEG11                   MPEG1
#> BCL2A11                 BCL2A1
#> AC092484.11         AC092484.1
#> CCL4L21                 CCL4L2
#> SPC251                   SPC25
#> CDK11                     CDK1
#> HLX                        HLX
#> PID1                      PID1
#> C1QC1                     C1QC
#> CCNA21                   CCNA2
#> CST32                     CST3
#> PRR11                    PRR11
#> CKAP2L                  CKAP2L
#> COTL12                   COTL1
#> CD402                     CD40
#> FCGR2A1                 FCGR2A
#> CAPG2                     CAPG
#> SLC25A52               SLC25A5
#> ANPEP2                   ANPEP
#> TOP2A1                   TOP2A
#> CXCL81                   CXCL8
#> NCAPH                    NCAPH
#> OLFML2B2               OLFML2B
#> IRF72                     IRF7
#> LMNB11                   LMNB1
#> BMP2K1                   BMP2K
#> C15orf481             C15orf48
#> ASF1B                    ASF1B
#> STX111                   STX11
#> RHOF                      RHOF
#> FAM105A1               FAM105A
#> RASGRP3                RASGRP3
#> AURKB1                   AURKB
#> NFKBIA2                 NFKBIA
#> ID21                       ID2
#> MS4A4A1                 MS4A4A
#> CDCA8                    CDCA8
#> CELF22                   CELF2
#> EBI31                     EBI3
#> NUSAP11                 NUSAP1
#> RAB312                   RAB31
#> CD931                     CD93
#> GPX12                     GPX1
#> PLXNC11                 PLXNC1
#> FCGRT2                   FCGRT
#> LIMD23                   LIMD2
#> IL181                     IL18
#> PTGER41                 PTGER4
#> SLC16A31               SLC16A3
#> SERPINA11             SERPINA1
#> SELPLG1                 SELPLG
#> HBEGF1                   HBEGF
#> CACNA2D41             CACNA2D4
#> STEAP11                 STEAP1
#> SLAMF81                 SLAMF8
#> GIMAP41                 GIMAP4
#> SHCBP1                  SHCBP1
#> KIAA01011             KIAA0101
#> GK1                         GK
#> FILIP1L2               FILIP1L
#> MMP92                     MMP9
#> MYO1G1                   MYO1G
#> MAFB1                     MAFB
#> CCL31                     CCL3
#> TYMS1                     TYMS
#> VSIG41                   VSIG4
#> C1QA1                     C1QA
#> DLGAP5                  DLGAP5
#> MCOLN21                 MCOLN2
#> CDKN31                   CDKN3
#> NCF11                     NCF1
#> CEP551                   CEP55
#> FTL1                       FTL
#> FBXO5                    FBXO5
#> CD300E1                 CD300E
#> LILRB41                 LILRB4
#> APOBEC3C1             APOBEC3C
#> MRC1                      MRC1
#> RGS103                   RGS10
#> RAB7B                    RAB7B
#> ASPM1                     ASPM
#> IER51                     IER5
#> CENPA1                   CENPA
#> TSC22D32               TSC22D3
#> GNG2                      GNG2
#> VIM2                       VIM
#> RASSF41                 RASSF4
#> ARL4C2                   ARL4C
#> SPC24                    SPC24
#> OSCAR1                   OSCAR
#> CD361                     CD36
#> EPSTI12                 EPSTI1
#> OAS12                     OAS1
#> ANLN1                     ANLN
#> EGR2                      EGR2
#> NR4A3                    NR4A3
#> DEPDC1                  DEPDC1
#> HIST1H4C2             HIST1H4C
#> MVP2                       MVP
#> C1orf211               C1orf21
#> DHFR1                     DHFR
#> OLR11                     OLR1
#> PLTP1                     PLTP
#> GATM1                     GATM
#> KIF20B1                 KIF20B
#> IL27RA                  IL27RA
#> ADAMDEC11             ADAMDEC1
#> BRCA1                    BRCA1
#> BIRC3                    BIRC3
#> SLC39A81               SLC39A8
#> KIF2C                    KIF2C
#> AKR1B12                 AKR1B1
#> NCF21                     NCF2
#> APOBEC3G              APOBEC3G
#> ARHGAP11A            ARHGAP11A
#> GHRL                      GHRL
#> BIRC51                   BIRC5
#> RFTN11                   RFTN1
#> AXL1                       AXL
#> HIVEP3                  HIVEP3
#> ACP51                     ACP5
#> AP1S21                   AP1S2
#> IGSF211                 IGSF21
#> CDCA5                    CDCA5
#> VASH11                   VASH1
#> LPAR62                   LPAR6
#> FADS11                   FADS1
#> NDC801                   NDC80
#> TRPV21                   TRPV2
#> MND11                     MND1
#> RPS32                     RPS3
#> C1QB1                     C1QB
#> SGOL2                    SGOL2
#> TUBA1B1                 TUBA1B
#> PILRA1                   PILRA
#> NUF21                     NUF2
#> UBE2C1                   UBE2C
#> SMC21                     SMC2
#> CAMK13                   CAMK1
#> FOXM11                   FOXM1
#> SLC25A191             SLC25A19
#> AURKA                    AURKA
#> CDCA2                    CDCA2
#> FAM198B2               FAM198B
#> BUB1                      BUB1
#> CFD1                       CFD
#> MXD3                      MXD3
#> HIST2H2AC            HIST2H2AC
#> PRC11                     PRC1
#> ADGRE21                 ADGRE2
#> PSTPIP11               PSTPIP1
#> HLA-B4                   HLA-B
#> CCL41                     CCL4
#> CCDC85B2               CCDC85B
#> EPB41L31               EPB41L3
#> COQ21                     COQ2
#> TAP12                     TAP1
#> C19orf481             C19orf48
#> KIF4A                    KIF4A
#> TACC31                   TACC3
#> PTPN7                    PTPN7
#> STAB11                   STAB1
#> PLXND11                 PLXND1
#> CKAP21                   CKAP2
#> ZWINT1                   ZWINT
#> PLK11                     PLK1
#> NUDT11                   NUDT1
#> CCNB11                   CCNB1
#> CCL51                     CCL5
#> BCAT12                   BCAT1
#> LAT22                     LAT2
#> SKA1                      SKA1
#> CENPF1                   CENPF
#> ANKRD281               ANKRD28
#> SQRDL2                   SQRDL
#> ATF51                     ATF5
#> GAS2L3                  GAS2L3
#> DNMT1                    DNMT1
#> KIAA1524              KIAA1524
#> SLCO2B11               SLCO2B1
#> PTTG11                   PTTG1
#> GBP22                     GBP2
#> MAD2L11                 MAD2L1
#> LDHB2                     LDHB
#> SGOL11                   SGOL1
#> KIF23                    KIF23
#> GIMAP71                 GIMAP7
#> UPP11                     UPP1
#> IQGAP21                 IQGAP2
#> CTSB2                     CTSB
#> MCM5                      MCM5
#> LINC009361           LINC00936
#> CD1631                   CD163
#> TK11                       TK1
#> UBE2S1                   UBE2S
#> PGAM11                   PGAM1
#> S100A42                 S100A4
#> FEN11                     FEN1
#> CD521                     CD52
#> DTYMK1                   DTYMK
#> SKA22                     SKA2
#> FCGR1A1                 FCGR1A
#> SOCS32                   SOCS3
#> TROAP1                   TROAP
#> GBP12                     GBP1
#> CCNB21                   CCNB2
#> LXN2                       LXN
#> LINC010941           LINC01094
#> CDCA31                   CDCA3
#> NASP1                     NASP
#> ETV52                     ETV5
#> SAMD92                   SAMD9
#> H2AFZ1                   H2AFZ
#> PRIM1                    PRIM1
#> OAS3                      OAS3
#> TPX21                     TPX2
#> ID32                       ID3
#> TRIM222                 TRIM22
#> IFI442                   IFI44
#> CXCL21                   CXCL2
#> IRF81                     IRF8
#> HMGB21                   HMGB2
#> HMGN21                   HMGN2
#> CDC201                   CDC20
#> HS3ST12                 HS3ST1
#> LMO43                     LMO4
#> XAF12                     XAF1
#> FANCI1                   FANCI
#> ITGAX1                   ITGAX
#> RAD51AP11             RAD51AP1
#> LRR11                     LRR1
#> IER31                     IER3
#> PARPBP1                 PARPBP
#> KIFC11                   KIFC1
#> THBD1                     THBD
#> G0S21                     G0S2
#> APOC11                   APOC1
#> TNFSF102               TNFSF10
#> ZNF3311                 ZNF331
#> PHACTR12               PHACTR1
#> FABP51                   FABP5
#> DAB22                     DAB2
#> SMC41                     SMC4
#> PPIF1                     PPIF
#> TCIRG11                 TCIRG1
#> TRAC                      TRAC
#> HLA-F4                   HLA-F
#> WARS                      WARS
#> TMPO1                     TMPO
#> FAM129A1               FAM129A
#> DUSP61                   DUSP6
#> ICAM12                   ICAM1
#> DNAJC153               DNAJC15
#> DMXL21                   DMXL2
#> DUT2                       DUT
#> LINC001523           LINC00152
#> MPP11                     MPP1
#> LIMA12                   LIMA1
#> CENPH1                   CENPH
#> GPX31                     GPX3
#> PARVB1                   PARVB
#> ERO1A1                   ERO1A
#> SAC3D12                 SAC3D1
#> CDT11                     CDT1
#> C1orf542               C1orf54
#> BTG22                     BTG2
#> CENPU1                   CENPU
#> EZH2                      EZH2
#> HMOX11                   HMOX1
#> SPP11                     SPP1
#> KCNN41                   KCNN4
#> DAPP12                   DAPP1
#> THOP11                   THOP1
#> IRF1                      IRF1
#> ORC61                     ORC6
#> NCAPD21                 NCAPD2
#> SORL12                   SORL1
#> BRCA21                   BRCA2
#> MASTL                    MASTL
#> PARP142                 PARP14
#> CXCL31                   CXCL3
#> RACGAP11               RACGAP1
#> ARPC43                   ARPC4
#> PKMYT11                 PKMYT1
#> MSR11                     MSR1
#> CLEC4E1                 CLEC4E
#> KNSTRN1                 KNSTRN
#> MELK1                     MELK
#> MARCKS2                 MARCKS
#> NEK21                     NEK2
#> NRP12                     NRP1
#> EPB41L22               EPB41L2
#> OTOA1                     OTOA
#> LSM41                     LSM4
#> NETO21                   NETO2
#> TUBA1C3                 TUBA1C
#> H2AFX1                   H2AFX
#> DCXR1                     DCXR
#> C5AR11                   C5AR1
#> LGALS3BP2             LGALS3BP
#> HTRA12                   HTRA1
#> MYBL21                   MYBL2
#> C1QBP1                   C1QBP
#> CD1092                   CD109
#> GPSM22                   GPSM2
#> GMNN1                     GMNN
#> RRM11                     RRM1
#> ACOT71                   ACOT7
#> RPA31                     RPA3
#> ANP32E2                 ANP32E
#> BID2                       BID
#> CAMK1D2                 CAMK1D
#> PPA13                     PPA1
#> CKAP51                   CKAP5
#> RFC41                     RFC4
#> PAPSS21                 PAPSS2
#> MX13                       MX1
#> CMSS12                   CMSS1
#> CLDN12                   CLDN1
#> TREM21                   TREM2
#> VAMP52                   VAMP5
#> CBR3                      CBR3
#> CDKN1A1                 CDKN1A
#> BASP13                   BASP1
#> NKG71                     NKG7
#> AHCY                      AHCY
#> C32                         C3
#> SDS1                       SDS
#> NEK22                     NEK2
#> PBK1                       PBK
#> CENPA2                   CENPA
#> CDC25C                  CDC25C
#> NMU                        NMU
#> KIF20A                  KIF20A
#> KIF231                   KIF23
#> HMMR1                     HMMR
#> LGR6                      LGR6
#> BUB11                     BUB1
#> CCNA22                   CCNA2
#> PRSS3                    PRSS3
#> DLGAP51                 DLGAP5
#> CENPE1                   CENPE
#> NCAPG1                   NCAPG
#> PLK12                     PLK1
#> NUF22                     NUF2
#> FAM72D                  FAM72D
#> FAM83D                  FAM83D
#> FAM72C                  FAM72C
#> KIF151                   KIF15
#> FAM64A1                 FAM64A
#> KIF2C1                   KIF2C
#> DEPDC11                 DEPDC1
#> SGOL12                   SGOL1
#> GTSE11                   GTSE1
#> KIF4A1                   KIF4A
#> CDCA32                   CDCA3
#> SHCBP11                 SHCBP1
#> SPC252                   SPC25
#> CDKN32                   CDKN3
#> BIRC52                   BIRC5
#> TPX22                     TPX2
#> RACGAP12               RACGAP1
#> AURKB2                   AURKB
#> ANLN2                     ANLN
#> TROAP2                   TROAP
#> KIF14                    KIF14
#> CDCA81                   CDCA8
#> CDC202                   CDC20
#> NUSAP12                 NUSAP1
#> CCNB12                   CCNB1
#> KIFC12                   KIFC1
#> AURKA1                   AURKA
#> TOP2A2                   TOP2A
#> KIAA15241             KIAA1524
#> MKI672                   MKI67
#> CEP552                   CEP55
#> MAD2L12                 MAD2L1
#> ECT21                     ECT2
#> CENPF2                   CENPF
#> CCNB22                   CCNB2
#> SPAG5                    SPAG5
#> ASPM2                     ASPM
#> KIF20B2                 KIF20B
#> PKMYT12                 PKMYT1
#> PARPBP2                 PARPBP
#> FOXM12                   FOXM1
#> ASF1B1                   ASF1B
#> PRC12                     PRC1
#> TACC32                   TACC3
#> SGOL21                   SGOL2
#> FAM72B                  FAM72B
#> TK12                       TK1
#> KIF111                   KIF11
#> CKAP22                   CKAP2
#> UBE2C2                   UBE2C
#> SKA11                     SKA1
#> PTTG12                   PTTG1
#> CDCA21                   CDCA2
#> LMNB12                   LMNB1
#> NDC802                   NDC80
#> TTK                        TTK
#> SKA23                     SKA2
#> MELK2                     MELK
#> DBF4                      DBF4
#> KIF18A                  KIF18A
#> MYBL22                   MYBL2
#> GPSM23                   GPSM2
#> POC1A1                   POC1A
#> TRIP131                 TRIP13
#> SKA3                      SKA3
#> ZMYND10                ZMYND10
#> IQGAP3                  IQGAP3
#> LBR2                       LBR
#> ATAD22                   ATAD2
#> TRAF3IP32             TRAF3IP3
#> DEPDC1B                DEPDC1B
#> UBE2S2                   UBE2S
#> ACTL81                   ACTL8
#> CDKN2C1                 CDKN2C
#> UBE2T1                   UBE2T
#> DIAPH31                 DIAPH3
#> CKS23                     CKS2
#> SMC42                     SMC4
#> NCAPD22                 NCAPD2
#> HMGB22                   HMGB2
#> PSRC11                   PSRC1
#> TEX301                   TEX30
#> CDKN2D1                 CDKN2D
#> RIBC2                    RIBC2
#> REEP41                   REEP4
#> HN12                       HN1
#> CDCA51                   CDCA5
#> RAD214                   RAD21
#> MXD31                     MXD3
#> ARL6IP13               ARL6IP1
#> ANP32E3                 ANP32E
#> CENPW1                   CENPW
#> KPNA2                    KPNA2
#> PRPH                      PRPH
#> KIF222                   KIF22
#> MCM72                     MCM7
#> AC004381.6          AC004381.6
#> WFDC11                   WFDC1
#> GGH2                       GGH
#> HMGB12                   HMGB1
#> ZWINT2                   ZWINT
#> CKAP2L1                 CKAP2L
#> GAS2L31                 GAS2L3
#> CA83                       CA8
#> CKS1B3                   CKS1B
#> HSPD11                   HSPD1
#> SAPCD21                 SAPCD2
#> CENPM2                   CENPM
#> BRCA22                   BRCA2
#> SSRP13                   SSRP1
#> CENPU2                   CENPU
#> CDT12                     CDT1
#> KIAA01012             KIAA0101
#> RAC32                     RAC3
#> STMN13                   STMN1
#> CDK12                     CDK1
#> EME1                      EME1
#> YBX22                     YBX2
#> MT1E1                     MT1E
#> CCNF                      CCNF
#> TMX23                     TMX2
#> DTYMK2                   DTYMK
#> HSP90AB15             HSP90AB1
#> TUBA1B2                 TUBA1B
#> RAD51AP12             RAD51AP1
#> CHEK11                   CHEK1
#> CENPN1                   CENPN
#> CKAP52                   CKAP5
#> ACOT72                   ACOT7
#> FANCI2                   FANCI
#> CCT53                     CCT5
#> CACYBP4                 CACYBP
#> RTKN2                    RTKN2
#> ODC11                     ODC1
#> LSM52                     LSM5
#> SNRNP251               SNRNP25
#> TYMS2                     TYMS
#> TMPO2                     TMPO
#> NME12                     NME1
#> DKC12                     DKC1
#> HYLS11                   HYLS1
#> KRT85                     KRT8
#> SMC22                     SMC2
#> KNSTRN2                 KNSTRN
#> MT2A4                     MT2A
#> H2AFV2                   H2AFV
#> ACTL6A3                 ACTL6A
#> NUCKS13                 NUCKS1
#> RUVBL21                 RUVBL2
#> HSPH12                   HSPH1
#> CKB4                       CKB
#> RUVBL12                 RUVBL1
#> PIF1                      PIF1
#> TIMM103                 TIMM10
#> H2AFX2                   H2AFX
#> CEP703                   CEP70
#> DNMT11                   DNMT1
#> RAD51                    RAD51
#> DSC33                     DSC3
#> MDC11                     MDC1
#> GMNN2                     GMNN
#> PCBD12                   PCBD1
#> LYAR2                     LYAR
#> AZGP13                   AZGP1
#> H2AFZ2                   H2AFZ
#> SLBP3                     SLBP
#> CASC51                   CASC5
#> OBP2B1                   OBP2B
#> FKBP43                   FKBP4
#> RRS12                     RRS1
#> HMGN22                   HMGN2
#> TOMM402                 TOMM40
#> DCXR2                     DCXR
#> CLSPN1                   CLSPN
#> CA61                       CA6
#> BSPRY2                   BSPRY
#> HIST1H2BN2           HIST1H2BN
#> EFHD13                   EFHD1
#> HSPA1A3                 HSPA1A
#> ZNF6951                 ZNF695
#> NCAPH1                   NCAPH
#> MTHFD23                 MTHFD2
#> SPDL11                   SPDL1
#> MT1G2                     MT1G
#> NTHL12                   NTHL1
#> RFC42                     RFC4
#> TMEM106C3             TMEM106C
#> PCNA1                     PCNA
#> LRRCC11                 LRRCC1
#> HJURP1                   HJURP
#> TUBA1C4                 TUBA1C
#> PPA14                     PPA1
#> PHGDH2                   PHGDH
#> SYCP22                   SYCP2
#> NES2                       NES
#> SLC2A4RG3             SLC2A4RG
#> DNAJC91                 DNAJC9
#> MCM21                     MCM2
#> RPL39L3                 RPL39L
#> PBX13                     PBX1
#> ALYREF2                 ALYREF
#> TTF22                     TTF2
#> ASNS2                     ASNS
#> EZH21                     EZH2
#> CENPH2                   CENPH
#> CDC45                    CDC45
#> ORC62                     ORC6
#> IGFBP22                 IGFBP2
#> PDGFRA3                 PDGFRA
#> PRNP2                     PRNP
#> SMC1B1                   SMC1B
#> ZNF367                  ZNF367
#> SLC9A3R23             SLC9A3R2
#> TUBB4B4                 TUBB4B
#> CHAF1A1                 CHAF1A
#> CDCA7                    CDCA7
#> CCT6A2                   CCT6A
#> DSN11                     DSN1
#> VANGL14                 VANGL1
#> C9orf401               C9orf40
#> UHRF1                    UHRF1
#> BOP12                     BOP1
#> MSX1                      MSX1
#> GAL2                       GAL
#> RAD54B                  RAD54B
#> PLOD32                   PLOD3
#> PAICS2                   PAICS
#> GINS21                   GINS2
#> CHORDC11               CHORDC1
#> RASL11B                RASL11B
#> HIST1H4C3             HIST1H4C
#> TUBB2A3                 TUBB2A
#> VGF2                       VGF
#> C12orf752             C12orf75
#> CTHRC13                 CTHRC1
#> MIS18A1                 MIS18A
#> ARHGAP11A1           ARHGAP11A
#> PRR111                   PRR11
#> FBLN12                   FBLN1
#> LAPTM4B4               LAPTM4B
#> CD3202                   CD320
#> TMEM972                 TMEM97
#> HSPA1B3                 HSPA1B
#> RFC31                     RFC3
#> HIST1H2BJ2           HIST1H2BJ
#> RMI22                     RMI2
#> ST143                     ST14
#> FRMD4A2                 FRMD4A
#> GTF3A1                   GTF3A
#> NT5DC21                 NT5DC2
#> PPIL12                   PPIL1
#> RNASEH2A1             RNASEH2A
#> PPP1R12A2             PPP1R12A
#> UCHL12                   UCHL1
#> CDC42EP13             CDC42EP1
#> TNFRSF181             TNFRSF18
#> SH3BGR3                 SH3BGR
#> DCTPP12                 DCTPP1
#> SLC43A33               SLC43A3
#> PYCR13                   PYCR1
#> SOX84                     SOX8
#> CCDC341                 CCDC34
#> SIGMAR12               SIGMAR1
#> CCDC85B3               CCDC85B
#> EXOSC8                  EXOSC8
#> BYSL2                     BYSL
#> RCCD11                   RCCD1
#> SERTAD43               SERTAD4
#> MCM10                    MCM10
#> MGST13                   MGST1
#> MND12                     MND1
#> SQLE3                     SQLE
#> KRT185                   KRT18
#> NAV2                      NAV2
#> TIMELESS              TIMELESS
#> ADAMTS41               ADAMTS4
#> PRKDC2                   PRKDC
#> DUSP9                    DUSP9
#> CENPQ2                   CENPQ
#> C6orf1412             C6orf141
#> SUSD51                   SUSD5
#> LMNB21                   LMNB2
#> PARVB2                   PARVB
#> DSCC11                   DSCC1
#> IL17B3                   IL17B
#> CTNNAL11               CTNNAL1
#> CMSS13                   CMSS1
#> HIBCH3                   HIBCH
#> GAMT1                     GAMT
#> MYC1                       MYC
#> GOLM13                   GOLM1
#> HACD13                   HACD1
#> CENPK1                   CENPK
#> PITX13                   PITX1
#> PSIP12                   PSIP1
#> FRMD31                   FRMD3
#> APOLD1                  APOLD1
#> BARD13                   BARD1
#> SDC23                     SDC2
#> APOBEC3B2             APOBEC3B
#> MCM31                     MCM3
#> GINS1                    GINS1
#> YDJC1                     YDJC
#> PSAT11                   PSAT1
#> MCM41                     MCM4
#> LSM42                     LSM4
#> CDC6                      CDC6
#> COL2A12                 COL2A1
#> NUDT12                   NUDT1
#> LY6E4                     LY6E
#> GINS41                   GINS4
#> RANBP12                 RANBP1
#> DKK12                     DKK1
#> DEK1                       DEK
#> SAP302                   SAP30
#> THOP12                   THOP1
#> ADAM154                 ADAM15
#> RBBP72                   RBBP7
#> CTSV3                     CTSV
#> ELOVL52                 ELOVL5
#> LRR12                     LRR1
#> ACTA24                   ACTA2
#> C21orf581             C21orf58
#> SERPINE22             SERPINE2
#> TCF191                   TCF19
#> TM4SF13                 TM4SF1
#> AHCY1                     AHCY
#> DCUN1D51               DCUN1D5
#> NEMP2                    NEMP2
#> RECQL41                 RECQL4
#> MYBL13                   MYBL1
#> TUBB2B3                 TUBB2B
#> DACT3                    DACT3
#> ERVMER34-12         ERVMER34-1
#> MT1X2                     MT1X
#> ITGA62                   ITGA6
#> TBC1D13                 TBC1D1
#> NQO13                     NQO1
#> FEN12                     FEN1
#> CDCA41                   CDCA4
#> NTF3                      NTF3
#> CPNE7                    CPNE7
#> FABP71                   FABP7
#> POLD22                   POLD2
#> RRM12                     RRM1
#> SEPT42                   SEPT4
#> EXO1                      EXO1
#> HAPLN12                 HAPLN1
#> LEFTY22                 LEFTY2
#> TOB13                     TOB1
#> PCOLCE23               PCOLCE2
#> PALLD3                   PALLD
#> LDHB3                     LDHB
#> MAPK133                 MAPK13
#> YES12                     YES1
#> PTS3                       PTS
#> DSP3                       DSP
#> TPM24                     TPM2
#> CRABP13                 CRABP1
#> ART32                     ART3
#> ACHE                      ACHE
#> SMYD22                   SMYD2
#> TNFRSF212             TNFRSF21
#> KLHL352                 KLHL35
#> KCTD14                   KCTD1
#> SYNCRIP2               SYNCRIP
#> SLC38A11               SLC38A1
#> CDCA7L2                 CDCA7L
#> ANO12                     ANO1
#> E2F11                     E2F1
#> NRM1                       NRM
#> BARX13                   BARX1
#> MARVELD11             MARVELD1
#> KLK111                   KLK11
#> RP11-357H14.173 RP11-357H14.17
#> PRSS333                 PRSS33
#> NEURL1B1               NEURL1B
#> DHTKD12                 DHTKD1
#> IDH13                     IDH1
#> HMGB32                   HMGB3
#> RPA32                     RPA3
#> RRM21                     RRM2
#> CYP39A13               CYP39A1
#> ACTA11                   ACTA1
#> PRIM11                   PRIM1
#> AP1M23                   AP1M2
#> UACA2                     UACA
#> PHF191                   PHF19
#> AIF1L4                   AIF1L
#> PDIA43                   PDIA4
#> EPHX13                   EPHX1
#> LNX13                     LNX1
#> CPED12                   CPED1
#> YEATS4                  YEATS4
#> MEX3A4                   MEX3A
#> DTL                        DTL
#> TUBB63                   TUBB6
#> PPIF2                     PPIF
#> DNPH12                   DNPH1
#> ANXA2R3                 ANXA2R
#> SCARB12                 SCARB1
#> CRISPLD13             CRISPLD1
#> COQ22                     COQ2
#> SLC25A192             SLC25A19
#> P3H43                     P3H4
#> IDI13                     IDI1
#> NASP2                     NASP
#> GLS3                       GLS
#> MTHFD11                 MTHFD1
#> HAUS1                    HAUS1
#> PIM12                     PIM1
#> TMPRSS32               TMPRSS3
#> INSIG11                 INSIG1
#> SERPINH14             SERPINH1
#> LRP23                     LRP2
#> EBP2                       EBP
#> MARCKSL13             MARCKSL1
#> EIF3B2                   EIF3B
#> FAM3C3                   FAM3C
#> KDELC22                 KDELC2
#> RASL121                 RASL12
#> WEE12                     WEE1
#> CITED43                 CITED4
#> EPB41L23               EPB41L2
#> SLC7A52                 SLC7A5
#> HSPB13                   HSPB1
#> IMPA23                   IMPA2
#> FAM208B2               FAM208B
#> MTL53                     MTL5
#> AARD4                     AARD
#> PAQR42                   PAQR4
#> DONSON1                 DONSON
#> MCM51                     MCM5
#> CSRP13                   CSRP1
#> SNHG192                 SNHG19
#> RTN4RL23               RTN4RL2
#> NKD23                     NKD2
#> CTNNBIP12             CTNNBIP1
#> SFN3                       SFN
#> TNFRSF12A3           TNFRSF12A
#> RAMP22                   RAMP2
#> TTC39A3                 TTC39A
#> CDK43                     CDK4
#> ENPP53                   ENPP5
#> STAT12                   STAT1
#> TUBA1A3                 TUBA1A
#> NDUFAF62               NDUFAF6
#> COL9A21                 COL9A2
#> FABP52                   FABP5
#> PPP1R1B2               PPP1R1B
#> C1QBP2                   C1QBP
#> PHYH4                     PHYH
#> UCP22                     UCP2
#> COL11A14               COL11A1
#> C1QL43                   C1QL4
#> HMGA13                   HMGA1
#> MARS1                     MARS
#> SYT83                     SYT8
#> BNIP32                   BNIP3
#> MASTL1                   MASTL
#> DDIT41                   DDIT4
#> LDHA1                     LDHA
#> PGP2                       PGP
#> CLPSL12                 CLPSL1
#> ATP6V0E22             ATP6V0E2
#> HELLS1                   HELLS
#> IGSF32                   IGSF3
#> GGCT3                     GGCT
#> RIBC21                   RIBC2
#> MCM101                   MCM10
#> DTL1                       DTL
#> EXO11                     EXO1
#> CDC451                   CDC45
#> CDC61                     CDC6
#> PKMYT13                 PKMYT1
#> DSCC12                   DSCC1
#> MCM22                     MCM2
#> CLSPN2                   CLSPN
#> CCNE2                    CCNE2
#> UHRF11                   UHRF1
#> ASF1B2                   ASF1B
#> EME11                     EME1
#> RRM22                     RRM2
#> ESCO21                   ESCO2
#> SMC1B2                   SMC1B
#> ZNF3671                 ZNF367
#> SGOL13                   SGOL1
#> GINS11                   GINS1
#> SPC253                   SPC25
#> NCAPG2                   NCAPG
#> CDCA52                   CDCA5
#> HELLS2                   HELLS
#> E2F12                     E2F1
#> CDCA71                   CDCA7
#> C16orf591             C16orf59
#> FAM111B                FAM111B
#> CDT13                     CDT1
#> CDK13                     CDK1
#> DHFR2                     DHFR
#> KIFC13                   KIFC1
#> PBK2                       PBK
#> ZWINT3                   ZWINT
#> MYBL23                   MYBL2
#> RAD511                   RAD51
#> RAD54L                  RAD54L
#> CHEK12                   CHEK1
#> RAD51AP13             RAD51AP1
#> NUF23                     NUF2
#> TK13                       TK1
#> RMI23                     RMI2
#> TCF192                   TCF19
#> MELK3                     MELK
#> RECQL42                 RECQL4
#> CENPM3                   CENPM
#> TRIP132                 TRIP13
#> RAD54B1                 RAD54B
#> TPX23                     TPX2
#> CENPU3                   CENPU
#> ELFN1-AS1            ELFN1-AS1
#> GINS22                   GINS2
#> KIAA01013             KIAA0101
#> AC004381.61         AC004381.6
#> BRCA11                   BRCA1
#> BRCA23                   BRCA2
#> NCAPH2                   NCAPH
#> PDZK1                    PDZK1
#> GMNN3                     GMNN
#> DONSON2                 DONSON
#> CNIH2                    CNIH2
#> LMNB13                   LMNB1
#> SPC241                   SPC24
#> KIF152                   KIF15
#> ATAD23                   ATAD2
#> FBXO51                   FBXO5
#> TYMS3                     TYMS
#> MAD2L13                 MAD2L1
#> BIRC53                   BIRC5
#> UBE2T2                   UBE2T
#> TOP2A3                   TOP2A
#> RFC32                     RFC3
#> PRC13                     PRC1
#> NUSAP13                 NUSAP1
#> CCNA23                   CCNA2
#> RTKN21                   RTKN2
#> FEN13                     FEN1
#> RFC43                     RFC4
#> CDCA22                   CDCA2
#> MASTL2                   MASTL
#> TIMELESS1             TIMELESS
#> MCM42                     MCM4
#> RNASEH2A2             RNASEH2A
#> C9orf135              C9orf135
#> CENPK2                   CENPK
#> C21orf582             C21orf58
#> CENPH3                   CENPH
#> RASL11B1               RASL11B
#> ORC63                     ORC6
#> UBE2C3                   UBE2C
#> MCM73                     MCM7
#> ERVMER34-13         ERVMER34-1
#> RRM13                     RRM1
#> ACTL82                   ACTL8
#> DSN12                     DSN1
#> DIAPH32                 DIAPH3
#> PCNA2                     PCNA
#> TMEM973                 TMEM97
#> DNAJC92                 DNAJC9
#> GINS42                   GINS4
#> MND13                     MND1
#> FANCI3                   FANCI
#> CENPF3                   CENPF
#> UNG2                       UNG
#> WFDC12                   WFDC1
#> RAC33                     RAC3
#> SLBP4                     SLBP
#> TEX302                   TEX30
#> CHAF1A2                 CHAF1A
#> RACGAP13               RACGAP1
#> MCM32                     MCM3
#> CENPW2                   CENPW
#> TTF23                     TTF2
#> DHTKD13                 DHTKD1
#> PRIM12                   PRIM1
#> POC1A2                   POC1A
#> SHCBP12                 SHCBP1
#> FAM83D1                 FAM83D
#> WDR76                    WDR76
#> GTSE12                   GTSE1
#> MKI673                   MKI67
#> ATP2A1-AS11         ATP2A1-AS1
#> CHAF1B2                 CHAF1B
#> KIF20B3                 KIF20B
#> DBF41                     DBF4
#> CCNB13                   CCNB1
#> NRM2                       NRM
#> SSRP14                   SSRP1
#> SYCP23                   SYCP2
#> MTHFD12                 MTHFD1
#> HMMR2                     HMMR
#> TACC33                   TACC3
#> BARD14                   BARD1
#> PAICS3                   PAICS
#> AURKB3                   AURKB
#> ANLN3                     ANLN
#> CDCA82                   CDCA8
#> CDC203                   CDC20
#> KIF18A1                 KIF18A
#> CKS1B4                   CKS1B
#> GGH3                       GGH
#> CENPN2                   CENPN
#> SPSB4                    SPSB4
#> C6orf1413             C6orf141
#> LRRCC12                 LRRCC1
#> CDCA42                   CDCA4
#> KIF4A2                   KIF4A
#> DNMT12                   DNMT1
#> NEMP21                   NEMP2
#> SKA12                     SKA1
#> SLC29A2                SLC29A2
#> PRKDC3                   PRKDC
#> TMEM106C4             TMEM106C
#> MCM52                     MCM5
#> KIF232                   KIF23
#> EZH22                     EZH2
#> YEATS41                 YEATS4
#> DTYMK3                   DTYMK
#> CDCA7L3                 CDCA7L
#> CA62                       CA6
#> TUBA1B3                 TUBA1B
#> CKS24                     CKS2
#> C16orf741             C16orf74
#> THEM63                   THEM6
#> SGOL22                   SGOL2
#> H2AFZ3                   H2AFZ
#> CCNB23                   CCNB2
#> HAPLN13                 HAPLN1
#> DKC13                     DKC1
#> KIAA15242             KIAA1524
#> DEK2                       DEK
#> PRPS22                   PRPS2
#> SKA31                     SKA3
#> RPA33                     RPA3
#> TMSB15A2               TMSB15A
#> HMGB23                   HMGB2
#> PARPBP3                 PARPBP
#> HACD14                   HACD1
#> KIF2C2                   KIF2C
#> FOXM13                   FOXM1
#> GTF3A2                   GTF3A
#> DCTPP13                 DCTPP1
#> LSM43                     LSM4
#> ALYREF3                 ALYREF
#> CDK44                     CDK4
#> TROAP3                   TROAP
#> HMGB13                   HMGB1
#> HIST1H1A              HIST1H1A
#> C9orf402               C9orf40
#> RTN4RL24               RTN4RL2
#> MIS18A2                 MIS18A
#> PRR72                     PRR7
#> VGF3                       VGF
#> TOMM403                 TOMM40
#> SNRNP252               SNRNP25
#> NME13                     NME1
#> CTNNAL12               CTNNAL1
#> AURKA2                   AURKA
#> GPR12                    GPR12
#> SMC23                     SMC2
#> HSPD12                   HSPD1
#> HYLS12                   HYLS1
#> CACYBP5                 CACYBP
#> CCNE1                    CCNE1
#> WDR344                   WDR34
#> CDCA33                   CDCA3
#> TIMM104                 TIMM10
#> APOBEC3B3             APOBEC3B
#> CENPA3                   CENPA
#> KIF141                   KIF14
#> CCT54                     CCT5
#> HIST1H1D1             HIST1H1D
#> MCM61                     MCM6
#> LAPTM4B5               LAPTM4B
#> EXOSC81                 EXOSC8
#> HMGN23                   HMGN2
#> NES3                       NES
#> BUB12                     BUB1
#> PTTG13                   PTTG1
#> STMN14                   STMN1
#> YBX23                     YBX2
#> NETO22                   NETO2
#> HIST1H4C4             HIST1H4C
#> GYLTL1B1               GYLTL1B
#> SLC43A34               SLC43A3
#> SMC43                     SMC4
#> KLK13                     KLK1
#> KIF223                   KIF22
#> SFN4                       SFN
#> LDHB4                     LDHB
#> ODC12                     ODC1
#> DUT3                       DUT
#> PDGFRA4                 PDGFRA
#> GGCT4                     GGCT
#> HAUS11                   HAUS1
#> LYAR3                     LYAR
#> PLEKHB14               PLEKHB1
#> RAD215                   RAD21
#> NASP3                     NASP
#> SOSTDC1                SOSTDC1
#> TBC1D14                 TBC1D1
#> NTHL13                   NTHL1
#> GDPD21                   GDPD2
#> NDC803                   NDC80
#> SPDEF1                   SPDEF
#> CRABP14                 CRABP1
#> LRR13                     LRR1
#> MXD32                     MXD3
#> AZGP14                   AZGP1
#> FANCD21                 FANCD2
#> IGFBP23                 IGFBP2
#> FKBP44                   FKBP4
#> CA84                       CA8
#> CKB5                       CKB
#> RAMP23                   RAMP2
#> PAQR43                   PAQR4
#> DNPH13                   DNPH1
#> P3H44                     P3H4
#> RANBP13                 RANBP1
#> BYSL3                     BYSL
#> PHGDH3                   PHGDH
#> LEFTY23                 LEFTY2
#> GAMT2                     GAMT
#> CCNF1                     CCNF
#> NUCKS14                 NUCKS1
#> NXPH4                    NXPH4
#> CLPSL13                 CLPSL1
#> PYCR14                   PYCR1
#> TUBB4                     TUBB
#> LNX14                     LNX1
#> RP11-357H14.174 RP11-357H14.17
#> DUSP91                   DUSP9
#> DEPDC12                 DEPDC1
#> ANP32E4                 ANP32E
#> FAM111A1               FAM111A
#> RCCD12                   RCCD1
#> VSNL12                   VSNL1
#> MYBL14                   MYBL1
#> GOLT1A2                 GOLT1A
#> PLK13                     PLK1
#> PPIL13                   PPIL1
#> HIST2H2AC1           HIST2H2AC
#> SUN33                     SUN3
#> NCAPD23                 NCAPD2
#> AIF1L5                   AIF1L
#> MYEOV                    MYEOV
#> TMPO3                     TMPO
#> SYT84                     SYT8
#> SUSD52                   SUSD5
#> B4GALNT1              B4GALNT1
#> PSRC12                   PSRC1
#> SLC2A4RG4             SLC2A4RG
#> AQP53                     AQP5
#> SLC25A193             SLC25A19
#> RPL39L4                 RPL39L
#> FAM19A3                FAM19A3
#> CACNB4                  CACNB4
#> C12orf753             C12orf75
#> CDKN2C2                 CDKN2C
#> CMBL1                     CMBL
#> AHCY2                     AHCY
#> HIST1H2BN3           HIST1H2BN
#> FBLN13                   FBLN1
#> SYNGR13                 SYNGR1
#> SHISA2                  SHISA2
#> C19orf482             C19orf48
#> HJURP2                   HJURP
#> FH3                         FH
#> PITX14                   PITX1
#> AP1M24                   AP1M2
#> COL2A13                 COL2A1
#> SCRG13                   SCRG1
#> TMX24                     TMX2
#> CKAP23                   CKAP2
#> VDR1                       VDR
#> CBR31                     CBR3
#> HMGA14                   HMGA1
#> SKA24                     SKA2
#> HIST1H1E1             HIST1H1E
#> ADAM155                 ADAM15
#> KPNA21                   KPNA2
#> PRSS334                 PRSS33
#> NT5DC22                 NT5DC2
#> HSD11B2                HSD11B2
#> CENPQ3                   CENPQ
#> MARC12                   MARC1
#> POLD23                   POLD2
#> EFHD14                   EFHD1
#> PSAT12                   PSAT1
#> BOP13                     BOP1
#> RUVBL13                 RUVBL1
#> H2AFX3                   H2AFX
#> PHLDA22                 PHLDA2
#> SLC43A11               SLC43A1
#> FAM64A2                 FAM64A
#> DHCR241                 DHCR24
#> ACOT73                   ACOT7
#> MGST14                   MGST1
#> ABCA7                    ABCA7
#> MDC12                     MDC1
#> HSP90AB16             HSP90AB1
#> HN13                       HN1
#> NUDT13                   NUDT1
#> MAPK134                 MAPK13
#> LMNB22                   LMNB2
#> SIGMAR13               SIGMAR1
#> MARCKSL14             MARCKSL1
#> PODXL24                 PODXL2
#> SERTAD44               SERTAD4
#> RUVBL22                 RUVBL2
#> KRT186                   KRT18
#> TTYH13                   TTYH1
#> SAPCD22                 SAPCD2
#> SEPT43                   SEPT4
#> CMSS14                   CMSS1
#> CEP553                   CEP55
#> CTHRC14                 CTHRC1
#> ITGA103                 ITGA10
#> CDKN33                   CDKN3
#> PRSS31                   PRSS3
#> PGP3                       PGP
#> MSX11                     MSX1
#> PSIP13                   PSIP1
#> NKD24                     NKD2
#> SERPINA52             SERPINA5
#> H2AFV3                   H2AFV
#> KLHDC34                 KLHDC3
#> TRAF3IP33             TRAF3IP3
#> SPDL12                   SPDL1
#> DEPDC1B1               DEPDC1B
#> CD3203                   CD320
#> HRCT13                   HRCT1
#> PCBD13                   PCBD1
#> RBBP73                   RBBP7
#> SERPINE1              SERPINE1
#> THOP13                   THOP1
#> IDI14                     IDI1
#> CYP39A14               CYP39A1
#> YDJC2                     YDJC
#> TTK1                       TTK
#> WIF11                     WIF1
#> CCDC342                 CCDC34
#> MYOZ12                   MYOZ1
#> CPNE71                   CPNE7
#> HIBCH4                   HIBCH
#> NANOS14                 NANOS1
#> CKAP2L2                 CKAP2L
#> RRS13                     RRS1
#> DGAT23                   DGAT2
#> EBP3                       EBP
#> PPA15                     PPA1
#> SLC9A3R24             SLC9A3R2
#> SAC3D13                 SAC3D1
#> MT1G3                     MT1G
#> UCHL13                   UCHL1
#> DCXR3                     DCXR
#> GOLM14                   GOLM1
#> ACTL6A4                 ACTL6A
#> ROPN1B3                 ROPN1B
#> HSPH13                   HSPH1
#> TUBB2B4                 TUBB2B
#> UBE2S3                   UBE2S
#> TUBB2A4                 TUBB2A
#> MT1E2                     MT1E
#> TUBA4A3                 TUBA4A
#> PPP1R12A3             PPP1R12A
#> ZNF6952                 ZNF695
#> CLDN111                 CLDN11
#> VANGL15                 VANGL1
#> GAPDH3                   GAPDH
#> ANO13                     ANO1
#> CA21                       CA2
#> PARVB3                   PARVB
#> ROPN13                   ROPN1
#> ENPP54                   ENPP5
#> LGR61                     LGR6
#> CCDC1672               CCDC167
#> LINC010481           LINC01048
#> GAL3                       GAL
#> LY6E5                     LY6E
#> MAPK8IP22             MAPK8IP2
#> S100B3                   S100B
#> DBI3                       DBI
#> PDLIM12                 PDLIM1
#> ASPM3                     ASPM
#> CSRP14                   CSRP1
#> PFN23                     PFN2
#> HIST1H2BJ3           HIST1H2BJ
#> NSG13                     NSG1
#> TPM25                     TPM2
#> PBX14                     PBX1
#> LSM53                     LSM5
#> HILPDA3                 HILPDA
#> NOTCH4                  NOTCH4
#> CKLF2                     CKLF
#> XAGE23                   XAGE2
#> ACTA25                   ACTA2
#> CTNND23                 CTNND2
#> FXYD63                   FXYD6
#> KLK112                   KLK11
#> SOHLH13                 SOHLH1
#> ELOVL53                 ELOVL5
#> CCND12                   CCND1
#> ATP6V0E23             ATP6V0E2
#> S100A13                 S100A1
#> TUBB4B5                 TUBB4B
#> GPSM24                   GPSM2
#> PCOLCE24               PCOLCE2
#> KIF112                   KIF11
#> KLHL353                 KLHL35
#> ARL6IP14               ARL6IP1
#> PDIA44                   PDIA4
#> FBXL222                 FBXL22
#> TUBA1C5                 TUBA1C
#> PTS4                       PTS
#> PRR112                   PRR11
#> MIA3                       MIA
#> KDELC23                 KDELC2
#> FKBP111                 FKBP11
#> NEK23                     NEK2
#> PPIF3                     PPIF
#> CCT6A3                   CCT6A
#> MTHFD24                 MTHFD2
#> COL11A23               COL11A2
#> CENPE2                   CENPE
#> COPZ22                   COPZ2
#> KRT87                     KRT8
#> ST144                     ST14
#> EXTL12                   EXTL1
#> DKK13                     DKK1
#> DMKN1                     DMKN
#> IDH14                     IDH1
#> UCP23                     UCP2
#> COL4A24                 COL4A2
#> PGAM12                   PGAM1
#> TEKT34                   TEKT3
#> MYC2                       MYC
#> MYL94                     MYL9
#> AARD5                     AARD
#> KCNN42                   KCNN4
#> FABP72                   FABP7
#> SNHG253                 SNHG25
#> CITED44                 CITED4
#> RIPPLY3                RIPPLY3
#> PTPRS2                   PTPRS
#> LBR3                       LBR
#> TRIB33                   TRIB3
#> PALLD4                   PALLD
#> HES62                     HES6
#> SERPINE23             SERPINE2
#> TSPAN23                 TSPAN2
#> PPP1R14A3             PPP1R14A
#> HOXB9                    HOXB9
#> PLOD33                   PLOD3
#> IGFBP74                 IGFBP7
#> SFRP13                   SFRP1
#> PPP1R1B3               PPP1R1B
#> NDUFAF63               NDUFAF6
#> GAS12                     GAS1
#> LAMB32                   LAMB3
#> LRP24                     LRP2
#> SH3BGR4                 SH3BGR
#> DLGAP52                 DLGAP5
#> LIMCH13                 LIMCH1
#> CTNNB12                 CTNNB1
#> CDKN2D2                 CDKN2D
#> NQO14                     NQO1
#> TNFRSF213             TNFRSF21
#> UACA3                     UACA
#> COL11A15               COL11A1
#> TTC39A4                 TTC39A
#> C1QL44                   C1QL4
#> LTBP12                   LTBP1
#> RPP251                   RPP25
#> QDPR2                     QDPR
#> ACTA12                   ACTA1
#> BARX14                   BARX1
#> DNM33                     DNM3
#> MT2A5                     MT2A
#> RP11-294J22.6    RP11-294J22.6
#> IMPA24                   IMPA2
#> FREM21                   FREM2
#> SNHG193                 SNHG19
#> TINAGL12               TINAGL1
#> DSC34                     DSC3
#> PEG103                   PEG10
#> NTF31                     NTF3
#> DACT31                   DACT3
#> FAM3C4                   FAM3C
#> RAB254                   RAB25
#> CPED13                   CPED1
#> BAMBI3                   BAMBI
#> NUDT43                   NUDT4
#> SELENBP14             SELENBP1
#> CENPV2                   CENPV
#> FHL14                     FHL1
#> IL17B4                   IL17B
#> MTL54                     MTL5
#> KIF20A1                 KIF20A
#> SCARB13                 SCARB1
#> FABP53                   FABP5
#> KRTCAP33               KRTCAP3
#> RGCC2                     RGCC
#> B3GNT73                 B3GNT7
#> CCDC28B1               CCDC28B
#> N4BP2                    N4BP2
#> SERTAD4-AS13       SERTAD4-AS1
#> STRADB                  STRADB
#> MYO1B3                   MYO1B
#> C2orf801               C2orf80
#> OBP2B2                   OBP2B
#> MTSS11                   MTSS1
#> LOXL13                   LOXL1
#> PHF192                   PHF19
#> RP11-273G15.21   RP11-273G15.2
#> NEGR13                   NEGR1
#> MYLK3                     MYLK
#> SLC29A13               SLC29A1
#> GCHFR1                   GCHFR
#> ITGA63                   ITGA6
#> SMOC12                   SMOC1
#> LIPH2                     LIPH
#> CDC42EP14             CDC42EP1
#> SERPINH15             SERPINH1
#> IDH1-AS11             IDH1-AS1
#> ETV53                     ETV5
#> NREP2                     NREP
#> SAP303                   SAP30
#> FKBP105                 FKBP10
#> GMPR3                     GMPR
#> CCDC85B4               CCDC85B
#> TNFRSF182             TNFRSF18
#> ZG16B4                   ZG16B
#> ANGPT11                 ANGPT1
#> TSEN152                 TSEN15
#> HSPA61                   HSPA6
#> COL4A11                 COL4A1
#> USP181                   USP18
#> GLS4                       GLS
#> WEE13                     WEE1
#> TPD52L14               TPD52L1
#> SDC24                     SDC2
#> FASN2                     FASN
#> HSPA1A4                 HSPA1A
#> FRMD4A3                 FRMD4A
#> TNFSF13B3             TNFSF13B
#> ARRB22                   ARRB2
#> HES43                     HES4
#> NAB12                     NAB1
#> QPCT3                     QPCT
#> MESP22                   MESP2
#> TMEM794                 TMEM79
#> TMPRSS132             TMPRSS13
#> NPPC2                     NPPC
#> IGF2                      IGF2
#> GJA13                     GJA1
#> FMOD1                     FMOD
#> EPHX14                   EPHX1
#> PIR3                       PIR
#> CEP704                   CEP70
#> SNX223                   SNX22
#> BSPRY3                   BSPRY
#> NFIB4                     NFIB
#> PIM13                     PIM1
#> TAGLN4                   TAGLN
#> DBNDD13                 DBNDD1
#> APOLD11                 APOLD1
#> MDFI3                     MDFI
#> C11orf731             C11orf73
#> SMYD23                   SMYD2
#> ASNS3                     ASNS
#> SOX85                     SOX8
#> KLRG24                   KLRG2
#> SMOC23                   SMOC2
#> PRSS85                   PRSS8
#> INSR                      INSR
#> CHORDC12               CHORDC1
#> CLU3                       CLU
#> INHBB2                   INHBB
#> DUSP22                   DUSP2
#> UPP12                     UPP1
#> CELF43                   CELF4
#> QPRT3                     QPRT
#> BARX2                    BARX2
#> PPDPF3                   PPDPF
#> FRMD32                   FRMD3
#> REEP42                   REEP4
#> TMEM2042               TMEM204
#> CRNDE4                   CRNDE
#> LAD12                     LAD1
#> CRELD22                 CRELD2
#> SLC25A53               SLC25A5
#> CTC-425F1.4        CTC-425F1.4
#> C1QBP3                   C1QBP
#> SCGB3A13               SCGB3A1
#> HSPA1B4                 HSPA1B
#> PMAIP12                 PMAIP1
#> PRNP3                     PRNP
#> TPM15                     TPM1
#> GRB102                   GRB10
#> DTNB3                     DTNB
#> FADS22                   FADS2
#> NAAA2                     NAAA
#> SQLE4                     SQLE
#> FAM105A2               FAM105A
#> MATN33                   MATN3
#> ELN3                       ELN
#> ATP1B13                 ATP1B1
#> MARVELD12             MARVELD1
#> LMTK31                   LMTK3
#> PTP4A33                 PTP4A3
#> TUBA1A4                 TUBA1A
#> EPPK12                   EPPK1
#> ERBB34                   ERBB3
#> ITIH61                   ITIH6
#> COL9A33                 COL9A3
#> MSRB34                   MSRB3
#> BACE24                   BACE2
#> TUBB64                   TUBB6
#> CTGF3                     CTGF
#> MAL24                     MAL2
#> ATP6V0D2              ATP6V0D2
#> RUFY4                    RUFY4
#> SIGLEC15              SIGLEC15
#> STRIP2                  STRIP2
#> DPP4                      DPP4
#> SLC9B2                  SLC9B2
#> PAGE2                    PAGE2
#> AKAP6                    AKAP6
#> AK5                        AK5
#> SUCNR1                  SUCNR1
#> SH3GL2                  SH3GL2
#> F5                          F5
#> SLC37A2                SLC37A2
#> ITGB3                    ITGB3
#> DGKI                      DGKI
#> TM4SF19                TM4SF19
#> MATK                      MATK
#> HTRA31                   HTRA3
#> RP11-702B10.2    RP11-702B10.2
#> ANGPTL6                ANGPTL6
#> CD1093                   CD109
#> SCD51                     SCD5
#> ECSCR.1                ECSCR.1
#> PTPN222                 PTPN22
#> LAT                        LAT
#> MAP1A1                   MAP1A
#> RAC22                     RAC2
#> GNG21                     GNG2
#> RFTN2                    RFTN2
#> NRIP3                    NRIP3
#> RAB38                    RAB38
#> CTSK2                     CTSK
#> HPGD                      HPGD
#> ITGA2                    ITGA2
#> CD842                     CD84
#> STRADB1                 STRADB
#> HSD3B7                  HSD3B7
#> SCARF11                 SCARF1
#> PARM11                   PARM1
#> THOP14                   THOP1
#> TNFRSF11A            TNFRSF11A
#> DUSP4                    DUSP4
#> CCR12                     CCR1
#> GAPLINC1               GAPLINC
#> CFI1                       CFI
#> ACP52                     ACP5
#> FJX1                      FJX1
#> SELPLG2                 SELPLG
#> OSCAR2                   OSCAR
#> SPINK2                  SPINK2
#> IGFBP61                 IGFBP6
#> SNX102                   SNX10
#> SPN                        SPN
#> ANPEP3                   ANPEP
#> F2R1                       F2R
#> EVI2A2                   EVI2A
#> CYTH42                   CYTH4
#> CLEC4A2                 CLEC4A
#> GNPTAB1                 GNPTAB
#> EID13                     EID1
#> MMP93                     MMP9
#> CD332                     CD33
#> FAM43A                  FAM43A
#> TFEC2                     TFEC
#> CSF1R2                   CSF1R
#> JDP2                      JDP2
#> EVI2B2                   EVI2B
#> TCIRG12                 TCIRG1
#> COL5A11                 COL5A1
#> VASH12                   VASH1
#> CAV11                     CAV1
#> ID33                       ID3
#> MITF2                     MITF
#> CKLF3                     CKLF
#> SPP12                     SPP1
#> FAM134B                FAM134B
#> CD682                     CD68
#> KIAA00401             KIAA0040
#> RGS104                   RGS10
#> EDIL31                   EDIL3
#> F13A11                   F13A1
#> RNASE62                 RNASE6
#> CST33                     CST3
#> CLEC11A1               CLEC11A
#> PRDM13                   PRDM1
#> CACNA2D42             CACNA2D4
#> IFI302                   IFI30
#> SPI12                     SPI1
#> SEPT63                   SEPT6
#> ERO1A2                   ERO1A
#> GPX13                     GPX1
#> WAS2                       WAS
#> SPARC2                   SPARC
#> LAPTM52                 LAPTM5
#> EHD21                     EHD2
#> UBB3                       UBB
#> PSTPIP12               PSTPIP1
#> FERMT32                 FERMT3
#> TGFB13                   TGFB1
#> MVP3                       MVP
#> LY963                     LY96
#> DOCK102                 DOCK10
#> APOBEC3C2             APOBEC3C
#> LSP12                     LSP1
#> LIMS12                   LIMS1
#> TSC22D33               TSC22D3
#> RCAN12                   RCAN1
#> GPR1832                 GPR183
#> XIST3                     XIST
#> EMP33                     EMP3
#> PTPRC2                   PTPRC
#> FILIP1L3               FILIP1L
#> CLECL12                 CLECL1
#> LCP12                     LCP1
#> CD300C1                 CD300C
#> CYTIP2                   CYTIP
#> PGAM13                   PGAM1
#> ITGB22                   ITGB2
#> CSTB3                     CSTB
#> N4BP21                   N4BP2
#> CD300A2                 CD300A
#> IGSF212                 IGSF21
#> BST23                     BST2
#> OAF2                       OAF
#> FAM162A1               FAM162A
#> SYCE1L2                 SYCE1L
#> TBXAS12                 TBXAS1
#> NFATC12                 NFATC1
#> FCER1G2                 FCER1G
#> TYROBP2                 TYROBP
#> MAP2K32                 MAP2K3
#> HLX1                       HLX
#> NCKAP1L2               NCKAP1L
#> CTNNBIP13             CTNNBIP1
#> B2M3                       B2M
#> DUSP62                   DUSP6
#> ADM1                       ADM
#> CITED21                 CITED2
#> SLC25A54               SLC25A5
#> CYBA3                     CYBA
#> STEAP12                 STEAP1
#> VIM3                       VIM
#> FCGRT3                   FCGRT
#> ADAM282                 ADAM28
#> CLEC2B3                 CLEC2B
#> S100A43                 S100A4
#> APOE1                     APOE
#> RHOF1                     RHOF
#> LILRB42                 LILRB4
#> LAT23                     LAT2
#> TRPV22                   TRPV2
#> BTK2                       BTK
#> ZEB11                     ZEB1
#> ABI32                     ABI3
#> RPS33                     RPS3
#> PDLIM42                 PDLIM4
#> LDHB5                     LDHB
#> TGM21                     TGM2
#> GNA152                   GNA15
#> APOC12                   APOC1
#> SQRDL3                   SQRDL
#> HNMT3                     HNMT
#> HLA-A4                   HLA-A
#> SLC16A32               SLC16A3
#> HSD17B113             HSD17B11
#> SLC39A82               SLC39A8
#> PALLD5                   PALLD
#> GAPDH4                   GAPDH
#> C1orf212               C1orf21
#> CD532                     CD53
#> IL27RA1                 IL27RA
#> PTTG14                   PTTG1
#> FTL2                       FTL
#> CKB6                       CKB
#> FAM46A1                 FAM46A
#> AIF12                     AIF1
#> ANKRD282               ANKRD28
#> PSMB83                   PSMB8
#> CD42                       CD4
#> HAVCR22                 HAVCR2
#> ID22                       ID2
#> ZEB23                     ZEB2
#> ADAM82                   ADAM8
#> C10orf542             C10orf54
#> COL18A11               COL18A1
#> CRIP13                   CRIP1
#> CA22                       CA2
#> BCAT13                   BCAT1
#> MPP12                     MPP1
#> SLC2A33                 SLC2A3
#> HTRA13                   HTRA1
#> COX7A11                 COX7A1
#> IFI163                   IFI16
#> RAI142                   RAI14
#> HIVEP31                 HIVEP3
#> FAM26F2                 FAM26F
#> VWA5A1                   VWA5A
#> RNASET23               RNASET2
#> EFHD22                   EFHD2
#> FAM198B3               FAM198B
#> MYO1F2                   MYO1F
#> TRAC1                     TRAC
#> GALM3                     GALM
#> ENO2                      ENO2
#> LST12                     LST1
#> OTOA2                     OTOA
#> PTPN71                   PTPN7
#> CASP44                   CASP4
#> MYO1B4                   MYO1B
#> OAS13                     OAS1
#> PLXND12                 PLXND1
#> FKBP112                 FKBP11
#> GPNMB2                   GPNMB
#> SAMHD12                 SAMHD1
#> ABCC32                   ABCC3
#> MTSS12                   MTSS1
#> GBP23                     GBP2
#> LYZ2                       LYZ
#> GCHFR2                   GCHFR
#> TCF43                     TCF4
#> LGALS92                 LGALS9
#> BMP2K2                   BMP2K
#> ABL21                     ABL2
#> CTSC2                     CTSC
#> NCF12                     NCF1
#> PLEK3                     PLEK
#> FUCA11                   FUCA1
#> SRGAP31                 SRGAP3
#> ARFGEF32               ARFGEF3
#> RAB313                   RAB31
#> RNASE11                 RNASE1
#> PFKP2                     PFKP
#> IFI27L23               IFI27L2
#> DHTKD14                 DHTKD1
#> DOK22                     DOK2
#> ARL4C3                   ARL4C
#> FOXP13                   FOXP1
#> NCF42                     NCF4
#> PYCARD2                 PYCARD
#> RRAGD2                   RRAGD
#> GMFG2                     GMFG
#> DCXR4                     DCXR
#> CTSB3                     CTSB
#> SAMD9L3                 SAMD9L
#> LGALS33                 LGALS3
#> ECM12                     ECM1
#> HSD17B43               HSD17B4
#> COL1A21                 COL1A2
#> MBP2                       MBP
#> NRP22                     NRP2
#> DUT4                       DUT
#> LRRC252                 LRRC25
#> ARHGAP183             ARHGAP18
#> DNMT13                   DNMT1
#> SLC9A72                 SLC9A7
#> PSMB93                   PSMB9
#> [1] ">>>>> Features that will be displayed..."
#>  [1] "IFRD1"    "CEBPD"    "SAT1"     "CLU"      "MMP7"     "CD163"   
#>  [7] "MS4A7"    "CYBB"     "FPR3"     "SLCO2B1"  "PCSK1N"   "CALML5"  
#> [13] "CFB"      "SEZ6L2"   "S100A6"   "TMEM106C" "MARCKSL1" "AZGP1"   
#> [19] "IDH1"     "CLPSL1"   "CLDN6"    "SMIM22"   "PCAT19"   "MAGIX"   
#> [25] "RAB6B"    "TOX"      "SLC26A7"  "SBSPON"   "COL11A2"  "S100P"   
#> [31] "CFH"      "NTM"      "PDGFRB"   "FBN1"     "ECM2"     "FCER1A"  
#> [37] "CD1E"     "CD1A"     "CLEC10A"  "CD1C"     "NEK2"     "PBK"     
#> [43] "CENPA"    "CDC25C"   "NMU"      "RIBC2"    "MCM10"    "DTL"     
#> [49] "EXO1"     "CDC45"    "ATP6V0D2" "RUFY4"    "SIGLEC15" "STRIP2"  
#> [55] "DPP4"    
#> [1] ">>>>> Head of feature data..."
#>       AAACCTGCACCTTGTC AAACGGGAGTCCTCCT AAACGGGTCCAGAGGA AAAGATGCAGTTTACG
#> IFRD1      -0.53124386       -0.2362754        0.9024895       -1.2434526
#> CEBPD      -0.04417121       -1.0730323        0.8303619       -1.8354638
#> SAT1       -0.15677948        0.4396964        0.9859063       -0.7343847
#> CLU         0.67627503       -1.2361823        1.2489661        0.4217206
#> MMP7       -0.40480716       -0.9602354        1.9093987        0.5072913
#> CD163      -0.40891468        1.1342452       -0.4089147       -0.4089147
#>       AAAGCAACAGGAATGC AAAGCAATCGGAATCT AAAGTAGAGTGTACTC AAAGTAGCAGCCTATA
#> IFRD1        1.6231150      -1.24345261        0.2350531        0.1851641
#> CEBPD        1.4661604      -1.83546381        0.1362996       -0.8819674
#> SAT1        -0.7006112       0.07341751        0.7776911       -1.0054831
#> CLU          1.0069610      -1.23618235       -0.7546259        0.1834973
#> MMP7         0.9266660      -0.20044974       -0.2596847       -0.6248048
#> CD163       -0.4089147      -0.40891468        2.4434918       -0.4089147
#>       AAAGTAGGTACAGTTC AAAGTAGTCGCATGAT AAAGTAGTCTATCCCG AAATGCCAGTACGTTC
#> IFRD1       -1.2434526       0.66497243        0.2966520        1.1256303
#> CEBPD        0.3816619       1.72797779        0.1955085       -0.4030097
#> SAT1        -0.3702289      -1.63165506        1.0652489        0.4745282
#> CLU         -0.7098273      -0.07356937        0.9272701        1.4411393
#> MMP7        -0.4753143      -0.28420995        0.7258350        1.8308129
#> CD163       -0.4089147      -0.40891468       -0.4089147       -0.4089147
#>       AAATGCCCATTGGTAC AAATGCCGTTGGTTTG AACACGTTCTTAGAGC AACCATGAGATGTCGG
#> IFRD1     0.4662465999       1.84250368       -0.5812834        0.5358595
#> CEBPD     2.0898613191      -0.05740542        1.2749966        0.9232562
#> SAT1     -0.0006470345       1.71916726        1.6553902        0.8579519
#> CLU      -1.2361823459       1.28299588        1.6287690        1.2235082
#> MMP7     -0.9602354113       2.23102725        1.5961554        1.0218828
#> CD163    -0.4089146752      -0.40891468       -0.4089147       -0.4089147
#>       AACCATGGTGTGGTTT AACCATGTCACCACCT AACCATGTCACCAGGC AACCATGTCGTATCAG
#> IFRD1        1.3174301      -0.06922206        0.4251841       -1.2434526
#> CEBPD        0.8346139       0.42477227       -1.8354638       -1.1940711
#> SAT1         0.2018444       1.07886645       -1.9623118        0.6317864
#> CLU          1.3523846       0.83287283       -1.2361823       -1.2361823
#> MMP7         1.4922774      -0.96023541        0.5483970       -0.9602354
#> CD163       -0.4089147      -0.40891468       -0.4089147        3.6110860
#>       AACCATGTCGTTTATC AACCGCGCACCTATCC AACCGCGGTGTCAATC AACCGCGTCAAAGACA
#> IFRD1       0.78704020        2.6525917       -1.0746040        0.8395387
#> CEBPD       0.03248883        1.8297231        0.3950349        1.3197322
#> SAT1        0.22067414        0.9225546       -1.8742395        1.4827249
#> CLU         1.08332259        1.2372284       -0.7352289        1.8425355
#> MMP7       -0.72476204        1.7893906       -0.4220233        0.9533575
#> CD163      -0.40891468       -0.4089147       -0.4089147       -0.4089147
#>       AACGTTGTCGACCAGC AACTCAGAGCTAGGCA AACTCAGCAAGAAAGG AACTCAGTCAGGATCT
#> IFRD1       -0.9341160      -0.07039225        2.0958270        1.5296400
#> CEBPD       -0.8867366      -0.94745907        0.6331833        2.1831390
#> SAT1        -1.3823762      -0.78804131        0.5294523        0.6747323
#> CLU         -0.9482958      -1.23618235        1.7106821        1.3252061
#> MMP7        -0.5953226      -0.96023541        0.2998086        1.3693671
#> CD163       -0.4089147       0.68321121       -0.4089147       -0.4089147
#>       AACTCAGTCCCAAGAT AACTCCCCAAGCCCAC AACTCCCCAGCAGTTT AACTCCCGTAGGAGTC
#> IFRD1       -1.2434526        0.8538223        2.1590559        1.4886823
#> CEBPD        0.4531251        1.1059551        1.7639695        1.6320972
#> SAT1        -1.3816437        1.1480840        1.6949564        1.0965638
#> CLU         -0.6195335        0.9803124        1.7425097        1.5413377
#> MMP7        -0.6137642        1.4002255        1.2348435        1.6605702
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       AACTCCCGTGCACCAC AACTCCCTCTCTAGGA AACTCTTCAATAACGA AACTCTTTCTTGTCAT
#> IFRD1        1.1925154       -0.2489341       -1.2434526       -0.1649797
#> CEBPD        1.4376869       -0.6582198       -1.8354638       -0.2270046
#> SAT1        -0.1493864       -0.3302763        0.6327619       -0.1355941
#> CLU          1.4827477       -0.7246923       -0.7320943       -0.3804238
#> MMP7         0.9651699       -0.9602354       -0.9602354        0.1209269
#> CD163       -0.4089147        3.4137069       -0.4089147       -0.4089147
#>       AACTGGTGTTTCCACC AAGACCTCAGCTGCTG AAGCCGCAGAAACCGC AAGGAGCAGACAAGCC
#> IFRD1       -0.6115918       0.09682116       -0.5223900      -0.02100873
#> CEBPD       -1.3571462       0.98593783       -0.4738173       0.55210760
#> SAT1        -0.5767323       1.23286873        0.4579239      -2.82713424
#> CLU         -0.9112105       1.28284903       -0.8653332      -0.60746819
#> MMP7        -0.7927893       0.34344694       -0.3988944       0.38711081
#> CD163       -0.4089147      -0.40891468        0.6958710      -0.40891468
#>       AAGGAGCTCGAACTGT AAGGCAGGTCACCTAA AAGGCAGTCGCATGGC AAGGTTCGTAGCGCTC
#> IFRD1        1.1711579        1.2103741       -1.2434526       -1.2434526
#> CEBPD        2.0226538        1.2154958       -0.1382978       -1.2765256
#> SAT1        -0.1143882        1.1132039        0.8613470       -0.3798359
#> CLU          1.4142692        1.5819640       -1.2361823       -0.8564364
#> MMP7        -0.2542190        1.8781335       -0.9602354       -0.6103825
#> CD163       -0.4089147       -0.4089147        1.6636445        3.2959635
#>       AAGGTTCGTTTAAGCC AAGGTTCTCGAGCCCA AAGTCTGCAGCTATTG AAGTCTGCATCCGTGG
#> IFRD1        0.4068617       -1.2434526       -1.2434526       -0.4255869
#> CEBPD       -0.8024120       -0.3342581        0.1476519       -0.1442722
#> SAT1        -0.7364219       -0.8867419       -1.6540685       -1.5654929
#> CLU          0.4810302        0.8907676       -1.2361823       -0.2254018
#> MMP7        -0.9602354       -0.7070569       -0.9602354       -0.3031706
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       AATCCAGGTGCAACGA AATCGGTTCGCATGGC AATCGGTTCTCTAGGA ACACCAAAGCACGCCT
#> IFRD1       -1.2434526        1.2230900        1.3452966       -0.4920625
#> CEBPD       -1.8354638        0.6886944        0.6561230       -0.3966093
#> SAT1        -2.2991059        0.8322489        0.9409753       -0.5518887
#> CLU         -1.2361823        1.7867904        1.1527366       -1.1392991
#> MMP7        -0.6066657        1.2427583        0.7513276       -0.7234733
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       ACACCCTAGAGCTTCT ACACCCTGTAGCTAAA ACACCGGCATGGTCTA ACAGCTACACAGACTT
#> IFRD1        0.4522827       -1.2434526        1.1753992      -1.24345261
#> CEBPD       -1.0845526       -1.4756727        1.7051311      -0.48343087
#> SAT1        -0.2863854       -1.2635849        1.2288799      -0.02332579
#> CLU          0.5762614       -0.9917381        2.0386254      -1.23618235
#> MMP7        -0.6113741       -0.7350335        1.6897631      -0.54044644
#> CD163       -0.4089147       -0.4089147       -0.4089147      -0.40891468
#>       ACAGCTATCCAAGCCG ACATCAGCAGTGACAG ACCAGTAAGTTCGATC ACCAGTACAGTTAACC
#> IFRD1       -0.1924045       0.83766963        0.8479075       -1.2434526
#> CEBPD       -0.6026184       1.22006522        1.2575685       -0.6585050
#> SAT1         1.8715030       0.26202308        0.5815784       -2.3048367
#> CLU          0.9898242       1.18243778        0.5921743        1.1102185
#> MMP7         1.3457704      -0.80865371       -0.5715798        0.5665004
#> CD163       -0.4089147       0.08124232       -0.4089147       -0.4089147
#>       ACCAGTATCGTAGGAG ACCCACTCATCGATGT ACCCACTCATTGTGCA ACCGTAAAGTATCTCG
#> IFRD1       -1.2434526       0.25813408        1.8589813       -1.2434526
#> CEBPD        0.6365599      -0.69876519        1.5577574       -1.1188499
#> SAT1        -2.0872438      -0.02596362        1.5939228        0.1636551
#> CLU         -0.4021947      -1.23618235        1.1349504       -0.7493106
#> MMP7        -0.9602354      -0.50832200        2.1090095        0.0701476
#> CD163       -0.4089147      -0.40891468       -0.4089147       -0.4089147
#>       ACCGTAACATGGAATA ACCTTTACATCCGTGG ACCTTTATCACCTTAT ACGAGCCAGGCAAAGA
#> IFRD1       -0.8290510        0.3219412       -1.2434526      -0.30412965
#> CEBPD       -1.5217625        0.4745444        0.3512586      -0.07957978
#> SAT1        -0.8842545       -1.7194337        0.7539403      -1.86038799
#> CLU          0.6190907       -0.7206589       -0.4331711      -0.84550185
#> MMP7        -0.6146775        1.6577997       -0.9602354      -0.31289855
#> CD163       -0.4089147       -0.4089147        1.1221087      -0.40891468
#>       ACGAGCCTCCCTCAGT ACGAGGACATGGTCTA ACGATACGTCTAGTCA ACGATACGTGACAAAT
#> IFRD1       0.40484361        1.7168581        0.2048805       -1.2434526
#> CEBPD      -0.66162245        1.3167337        0.5980048       -1.8354638
#> SAT1        0.03315462        1.4530965        0.1171678       -0.5485277
#> CLU         0.96512229        1.2237313       -0.1447451       -1.2361823
#> MMP7       -0.82135637        1.8142218        0.2751482       -0.9602354
#> CD163      -0.40891468       -0.4089147       -0.4089147        2.1914980
#>       ACGATACGTTGTTTGG ACGATACTCCGTAGTA ACGATGTTCTGTCAAG ACGCAGCAGGCAATTA
#> IFRD1        0.6053256       -0.3041471       -1.2434526       0.24594382
#> CEBPD        0.2271069       -0.7132385       -0.7781608      -0.04233179
#> SAT1         0.3020321       -0.3307350       -0.2133313       0.51177220
#> CLU          2.3746261       -0.2764145       -0.8553208       1.15803338
#> MMP7        -0.3799449        1.5566235       -0.9602354      -0.58012519
#> CD163       -0.4089147       -0.4089147        1.7310646      -0.40891468
#>       ACGCAGCCAACACCCG ACGCAGCCAGATAATG ACGCAGCGTTCAACCA ACGCCAGTCACAGTAC
#> IFRD1        2.1564145       -1.2434526    -1.939462e-01       -0.2557361
#> CEBPD        1.4145072       -1.1001392     1.052514e+00       -0.8133421
#> SAT1         1.4450455       -0.2896648     1.143524e+00       -0.4141173
#> CLU          1.6569108       -1.2361823     1.115674e+00       -0.6601910
#> MMP7         0.9998798       -0.2377535     9.336994e-05        0.4451005
#> CD163       -0.4089147       -0.4089147    -4.089147e-01       -0.4089147
#>       ACGGAGAAGCAGGTCA ACGGAGACATCCGTGG ACGGCCAAGGGATGGG ACGGCCACAAAGGAAG
#> IFRD1        0.2790804        0.9627711       -1.2434526       -0.4398987
#> CEBPD       -0.1619583        2.1004929       -1.8354638       -0.5796384
#> SAT1        -0.5184712        1.4768921       -2.0917235       -0.8706893
#> CLU         -1.2361823        1.9112625       -1.2361823       -0.6834946
#> MMP7        -0.9602354        2.1681784       -0.9602354       -0.7413202
#> CD163        2.5190250       -0.4089147       -0.4089147       -0.4089147
#>       ACGGCCACACTTAACG ACGGCCAGTCAATACC ACGGCCAGTTTGCATG ACGGCCATCAGATAAG
#> IFRD1       -0.3504882       -0.4849304        0.2893728        0.3858285
#> CEBPD        0.9924762        0.9670956       -1.2869677        0.1925570
#> SAT1         0.5515978        0.2512257        0.2948578        0.3881888
#> CLU          1.1833675        0.6183016        0.4060776        0.6873283
#> MMP7        -0.7135553        0.9934394       -0.9602354       -0.9602354
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       ACGGCCATCCCTAACC ACGGGCTAGCTGGAAC ACGGGCTGTATAGGTA ACGGGTCCACCAGATT
#> IFRD1       -0.3132586       -0.3212495       -0.4426141      0.377765069
#> CEBPD       -0.2090079       -0.7304311       -1.4870475     -0.779019340
#> SAT1        -0.8662719       -0.2281965       -0.1224175      0.009861347
#> CLU          0.4259765       -1.2361823       -0.9994662      0.419672805
#> MMP7        -0.9602354       -0.9602354       -0.5807793     -0.811537006
#> CD163       -0.4089147        2.4112316       -0.4089147     -0.408914675
#>       ACGGGTCCATACGCCG ACGGGTCGTGTGAAAT ACGTCAAAGACCTAGG ACGTCAAGTGCAACGA
#> IFRD1        2.0615050       -0.2337882      -0.64342230        0.1232426
#> CEBPD        0.6398247       -0.1676897       0.64167311        1.1147719
#> SAT1         1.2105011       -0.6590789      -0.23999836        1.1098689
#> CLU          0.2807782       -1.0647740       0.06050867        1.3274575
#> MMP7         0.8801465       -0.4818327       1.36626893        2.2787024
#> CD163       -0.4089147       -0.4089147      -0.40891468       -0.4089147
#>       ACTATCTCAGTTTACG ACTATCTGTAAGTGGC ACTGAACAGATGCCTT ACTGAACGTTCGTGAT
#> IFRD1       0.02755483        1.1332686       -1.2434526        1.1859383
#> CEBPD       0.28812830        0.8346943        0.5552537        2.0021585
#> SAT1        0.34257983        1.2606365       -0.1724025        1.8615954
#> CLU        -0.73913636        1.3938481        1.3947881        1.9935252
#> MMP7       -0.50231619        1.3513535       -0.9602354        1.8539277
#> CD163      -0.07427060       -0.4089147        1.2249349       -0.4089147
#>       ACTGAACTCTTTAGTC ACTGAGTAGCGGCTTC ACTGAGTCATCGGAAG ACTGAGTTCCCGACTT
#> IFRD1       -1.2434526        1.0531090       -1.2434526       -0.1780087
#> CEBPD       -0.3859846        0.7522276       -1.8354638       -1.1256967
#> SAT1        -0.3054983        0.5551851       -0.2947288       -0.9052957
#> CLU         -0.3955994        1.6780271       -1.2361823       -0.9083353
#> MMP7        -0.9602354       -0.9602354       -0.9602354       -0.8425951
#> CD163       -0.4089147       -0.4089147        2.1562215       -0.4089147
#>       ACTGATGGTCTCCACT ACTGATGGTGAGTGAC ACTGCTCAGGGCACTA ACTGCTCCATACTACG
#> IFRD1     -0.042169313       -0.4942667        2.5698546       -0.2593467
#> CEBPD     -0.926094355        0.9038293        1.4895778       -0.3728087
#> SAT1       0.006930587       -0.9949854        1.2896630        0.4098250
#> CLU       -1.236182346        0.4760373        0.4143938       -1.2361823
#> MMP7      -0.529939863        1.2831082        1.4441037       -0.4939429
#> CD163      1.798932544       -0.4089147       -0.4089147        1.0988962
#>       ACTGTCCGTAGAAAGG ACTGTCCGTCATGCAT ACTGTCCTCAAACGGG ACTTACTAGTTTGCGT
#> IFRD1       -0.9218249        2.8714132       -1.2434526       -1.2434526
#> CEBPD        0.3900382        1.7562610       -0.0415927       -0.6841794
#> SAT1         0.5204261        1.1758007       -1.0062314       -0.2467483
#> CLU          0.4679311        1.2872836       -0.7641174       -1.0532063
#> MMP7        -0.8078406        1.9270005        1.1034235       -0.7916630
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       ACTTACTCACGAAGCA ACTTGTTAGTTCGCAT ACTTGTTTCAGTACGT ACTTTCACACAACTGT
#> IFRD1        0.8241887       -0.6615296       0.03023389        1.0209519
#> CEBPD       -0.6840207        0.5059696      -1.40299034       -1.4258472
#> SAT1        -0.6095178       -0.9270734      -1.92585098       -0.2261561
#> CLU         -0.2541897       -0.1604171       0.22298976        1.1954477
#> MMP7        -0.2395194       -0.9602354      -0.68953984       -0.3821624
#> CD163       -0.4089147       -0.4089147      -0.40891468        0.4201487
#>       ACTTTCACATGTCTCC ACTTTCAGTGCTGTAT ACTTTCATCCAGTAGT AGAATAGAGAAGGCCT
#> IFRD1       -1.2434526        1.0014109        0.2829842        1.2506467
#> CEBPD       -1.1337309        1.2284750        0.4035940        1.3865423
#> SAT1         0.6799563        1.4281192        1.9620849        0.2645878
#> CLU         -1.2361823        1.2155217        1.7438283        1.6318686
#> MMP7        -0.9602354        1.4704830        2.1386610        1.0279268
#> CD163        1.8374322       -0.4089147       -0.4089147       -0.4089147
#>       AGAATAGAGAAGGGTA AGAATAGTCAATCTCT AGACGTTAGCACGCCT AGACGTTAGCGTTGCC
#> IFRD1       -0.1862938       -0.9041952      -0.09400614        0.1250625
#> CEBPD        0.1281860       -1.5786465      -0.01118264       -0.3019798
#> SAT1        -1.1360290       -0.8812194       1.39528868       -1.8587486
#> CLU         -0.9112792       -1.0616991      -0.64501143       -0.6517867
#> MMP7        -0.8438078       -0.9602354      -0.51382248       -0.5554769
#> CD163       -0.4089147       -0.4089147      -0.40891468       -0.4089147
#>       AGACGTTGTAAAGTCA AGAGCGACATTTGCCC AGAGCTTTCGAACTGT AGAGTGGTCACTCTTA
#> IFRD1       0.44025706      -0.19852843       -1.2434526       -0.5455159
#> CEBPD      -0.34579315      -0.10942998       -0.9894046       -0.9633729
#> SAT1        0.35696549      -0.14368921        0.1780473       -2.3332618
#> CLU         0.08686456      -0.06350459       -1.2361823       -1.2361823
#> MMP7       -0.33427530      -0.05997619       -0.6409856       -0.4143727
#> CD163       1.61519986      -0.40891468        1.3035076       -0.4089147
#>       AGATCTGAGATGTGTA AGATCTGCACAGCCCA AGATCTGGTCCGAAGA AGATCTGGTCTAGCCG
#> IFRD1       -0.6853269        1.9399251       -0.8785299        1.0165460
#> CEBPD        0.1763666        1.4965692        0.6053762        1.0038007
#> SAT1         1.4596240        1.1531784        0.2324527        0.6085462
#> CLU         -0.9491332        1.6935393       -0.5491118        0.6190587
#> MMP7        -0.5105109        1.1965319       -0.4041685        0.2616281
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       AGATTGCCACCACCAG AGATTGCTCGTATCAG AGATTGCTCTTGAGAC AGCAGCCGTCCGACGT
#> IFRD1        0.4577900       -1.2434526        0.1415229       -0.4456838
#> CEBPD        0.3471140       -1.5069666        1.1813453       -0.1977870
#> SAT1         1.6696864       -0.9799214        2.4601768        0.6561266
#> CLU          0.7123785       -1.0129994        1.4339017       -1.2361823
#> MMP7         1.3968631       -0.6001489        2.2227717       -0.9602354
#> CD163       -0.4089147       -0.4089147       -0.4089147        2.5503776
#>       AGCAGCCTCAGCTCTC AGCAGCCTCGACAGCC AGCAGCCTCTCCGGTT AGCATACGTACCGAGA
#> IFRD1       -0.3869458      -1.24345261        1.7044670       -1.2434526
#> CEBPD       -1.8354638       0.56733278        1.4061225       -0.9773235
#> SAT1        -2.2210547      -0.08568643        1.5888173        0.1190208
#> CLU         -1.2361823      -0.84653928        1.9436080       -1.2361823
#> MMP7        -0.9602354      -0.96023541        1.1017203       -0.6356832
#> CD163        1.6925788       2.44499287       -0.4089147       -0.4089147
#>       AGCATACTCGTCACGG AGCCTAAAGCTGCAAG AGCCTAACAGCTGCTG AGCCTAAGTCGCTTTC
#> IFRD1        1.0267308      -0.07601338       -1.2434526        0.5069408
#> CEBPD        0.9520242       1.51904186        1.0542933       -0.7457408
#> SAT1         1.2967416      -1.41120051        0.5577727       -1.4953417
#> CLU          0.7470517      -0.87163025       -1.2361823       -0.5925334
#> MMP7         1.6518706       1.51697185       -0.5886905       -0.6437330
#> CD163       -0.4089147      -0.40891468        3.8709723       -0.4089147
#>       AGCGGTCAGCTCCTCT AGCGGTCTCATCTGTT AGCGGTCTCGATGAGG AGCGTATAGGACAGCT
#> IFRD1       -1.2434526      -0.19095281        0.3523592       -0.6668948
#> CEBPD       -0.4161606      -1.53669809        1.6046335       -0.2241152
#> SAT1        -1.6967442      -0.54988826        1.2261004       -0.3290730
#> CLU         -0.8428950       0.02404683        1.6184375       -0.4485291
#> MMP7        -0.5979071      -0.69722743        1.8087927        1.2321742
#> CD163       -0.4089147      -0.40891468       -0.4089147       -0.4089147
#>       AGCGTATGTTCAGACT AGCGTCGAGGGATGGG AGCGTCGAGTGATCGG AGCGTCGTCGAGCCCA
#> IFRD1       -0.4952048        0.8173580      -1.24345261        0.8796383
#> CEBPD       -0.1184417        0.3599381      -0.98075690        1.4185141
#> SAT1         0.3486573       -0.1065461      -0.09822076        0.9797430
#> CLU         -1.2361823        1.0098440      -0.48398951        1.9283642
#> MMP7        -0.6056980       -0.9602354      -0.26725415        0.8108534
#> CD163        0.7375231       -0.4089147      -0.40891468       -0.4089147
#>       AGCTCCTCACATCCGG AGCTCTCCAGCCACCA AGCTCTCTCAGCGATT AGCTTGAGTTACGGAG
#> IFRD1      -0.40562845      -0.68594918       -1.2434526      -0.77319910
#> CEBPD       0.68990083      -1.83546381        0.4702561      -0.02297766
#> SAT1       -0.40264250       0.02642327        1.2675001      -0.26442485
#> CLU         0.01523773       0.20500519       -0.8545733      -0.81637915
#> MMP7        1.47175804      -0.51094944       -0.9602354       0.12903208
#> CD163      -0.40891468      -0.40891468        3.3068283      -0.40891468
#>       AGGCCACCAGTAAGCG AGGCCGTAGTTCCACA AGGCCGTTCCGTACAA AGGGAGTTCACGAAGG
#> IFRD1        0.7127888       -1.2434526        1.4012431      -0.17109126
#> CEBPD        1.4263530       -1.8354638        1.5527717       0.05598846
#> SAT1        -0.5152408       -0.9871422       -0.8743079      -0.13196836
#> CLU          0.3162650       -1.2361823        0.4198944       0.74604097
#> MMP7         1.6091602       -0.5178300       -0.6374308       0.63925778
#> CD163       -0.4089147        1.0216545       -0.4089147      -0.40891468
#>       AGGGATGAGAAGATTC AGGGATGCAATCCGAT AGGGATGTCGCTTAGA AGGTCATAGACACTAA
#> IFRD1        0.5096497        0.1448749        1.3564008        2.3051082
#> CEBPD       -1.0327431       -1.1777119        1.0213300        1.9189021
#> SAT1         0.3484375       -0.2998530        0.2039190        1.5183802
#> CLU          0.1565545       -1.2361823        1.4239305        0.7639740
#> MMP7        -0.6598479       -0.9602354        0.8975538        1.2992402
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       AGGTCATAGCAGACTG AGGTCATAGCGTGAGT AGGTCATGTATAAACG AGGTCATGTGCTAGCC
#> IFRD1       -0.1059298       0.59985263        1.4806184       -1.2434526
#> CEBPD        0.9371342      -0.06782915        0.4357401       -0.2191491
#> SAT1        -0.5853625       0.77292002        0.9012116       -0.9164226
#> CLU          1.2632935       1.61462069        0.1648325       -0.4370420
#> MMP7        -0.1361549      -0.31078853       -0.9602354       -0.9602354
#> CD163       -0.4089147      -0.40891468       -0.4089147       -0.4089147
#>       AGGTCATTCTCTAAGG AGGTCCGTCTGTCTCG AGTCTTTAGCGCTTAT AGTGAGGAGGTGCACA
#> IFRD1       -1.2434526       -1.2434526       -1.2434526       -1.2434526
#> CEBPD       -0.6006314        0.3350811       -1.8354638       -1.2577647
#> SAT1        -2.8271342       -1.0825824       -1.8866327       -0.7931930
#> CLU         -0.3972298       -0.6304226       -0.3655326       -0.8436901
#> MMP7        -0.9602354        0.1794569       -0.9602354       -0.9602354
#> CD163       -0.4089147       -0.4089147       -0.4089147        2.8105211
#>       AGTGGGAGTAGAGCTG AGTGGGAGTCAATGTC AGTGGGAGTCGAAAGC AGTGGGAGTGGTCCGT
#> IFRD1       -1.2434526       -1.2434526       -1.2434526        1.4286646
#> CEBPD       -1.8354638       -1.0130565       -0.8176508        1.4130621
#> SAT1         0.6665133        0.4363842       -1.1986001        1.5032761
#> CLU         -0.6978326       -1.2361823       -1.2361823        0.2006005
#> MMP7        -0.9602354       -0.1670211       -0.9602354        1.6596352
#> CD163       -0.4089147        4.2355721       -0.4089147       -0.4089147
#>       AGTGTCAAGACTAGGC AGTGTCATCTGCTGTC AGTTGGTAGACCACGA AGTTGGTCAAGGTTTC
#> IFRD1        0.4519322        0.1754407       -0.3592790        0.3005036
#> CEBPD        0.5213809       -1.8354638       -1.1661459       -0.4233331
#> SAT1        -0.3021361       -0.7066169       -0.1693735       -0.9809724
#> CLU         -0.1547255        0.5794217       -1.2361823        0.5273756
#> MMP7        -0.7341362       -0.5377663       -0.2926653       -0.5626020
#> CD163       -0.4089147       -0.4089147        0.9457836        0.3345170
#>       AGTTGGTCAGGACCCT AGTTGGTGTCTGATTG AGTTGGTTCCCTAATT AGTTGGTTCGCCATAA
#> IFRD1       -0.4490498        1.3315690       -0.1158754        0.6239739
#> CEBPD       -1.2341022        1.8924746       -0.9818897        1.0307763
#> SAT1         0.1126506        0.9940910       -1.4781533        1.5846710
#> CLU         -1.0016066        1.9185384       -0.9431112        0.8384729
#> MMP7        -0.9602354        2.3195198       -0.6902344        1.9428984
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       ATAACGCCACGTCAGC ATAACGCCAGTGAGTG ATAACGCGTGTGGTTT ATAACGCTCCACTCCA
#> IFRD1        1.0386732       -0.5707756       -0.3793345      -0.74224459
#> CEBPD        1.3115040       -0.2987510       -1.8354638      -0.03201339
#> SAT1         0.4624685       -0.5221308       -0.9489413      -2.14414723
#> CLU          1.3123037       -0.7124069       -1.2361823      -0.79203469
#> MMP7         1.9479534        1.0549510       -0.9602354      -0.60321495
#> CD163       -0.4089147       -0.4089147       -0.4089147      -0.40891468
#>       ATAACGCTCCGCATAA ATAAGAGAGATATGGT ATAAGAGCATGTAGTC ATAAGAGTCCCTAATT
#> IFRD1       -1.2434526        0.8030148       -0.1308513       -0.5883645
#> CEBPD       -1.5241379        1.5327429       -1.5499058       -0.1325424
#> SAT1        -0.5988325        0.3041871       -0.4147060       -0.9546186
#> CLU         -1.2361823        1.4308484        0.2991002        1.4510283
#> MMP7        -0.6170216       -0.6940222       -0.7814977        1.4999172
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       ATAAGAGTCTTGAGAC ATAGACCAGCAGGTCA ATAGACCCATCCTTGC ATCATCTAGCTTATCG
#> IFRD1       0.04571267        1.2818296       -0.3948151      -0.21548399
#> CEBPD       0.14947150        0.6920332        0.4670445      -0.54118158
#> SAT1        0.37380686        1.0644252        1.2798826       0.03088762
#> CLU         0.99657802        1.8588355       -0.7997203       0.97660308
#> MMP7       -0.96023541        1.0482296       -0.5581310      -0.47315971
#> CD163      -0.40891468       -0.4089147        2.6807054      -0.40891468
#>       ATCATCTGTTCCCGAG ATCATCTTCGAGAACG ATCATCTTCTGTGCAA ATCATGGCACCACCAG
#> IFRD1      -0.48409141       1.16805889       -0.1839512        0.4583546
#> CEBPD       0.46706686       1.08198932        0.8394187        0.8119562
#> SAT1        0.09817722       0.01235415        0.8536372        1.4063205
#> CLU        -0.84563588       1.55652299       -1.2361823        1.3737822
#> MMP7        0.01190532       0.01989405       -0.9602354       -0.9602354
#> CD163      -0.40891468      -0.40891468        1.2144144       -0.4089147
#>       ATCATGGCACGACGAA ATCCACCAGCGTAATA ATCCACCAGGTGATAT ATCCACCTCGAGAGCA
#> IFRD1       -0.3495614        0.1640446        0.7173141      0.007082315
#> CEBPD        0.9044735       -0.7699908       -0.8397903      0.405509944
#> SAT1        -1.7469649        0.8704684       -0.9045450     -0.732343680
#> CLU          0.5586789        0.5244580       -1.2361823      0.599634827
#> MMP7        -0.3991547        1.0941037       -0.5739964     -0.814106623
#> CD163       -0.4089147       -0.4089147       -0.4089147     -0.408914675
#>       ATCGAGTAGGAACTGC ATCGAGTAGTGTCCCG ATCGAGTCAAGAGGCT ATCTACTAGCACCGTC
#> IFRD1       0.34329586       -0.7689090       0.19921636        1.8426296
#> CEBPD      -0.10450174       -0.0139357      -0.01958258        0.9165340
#> SAT1        0.03217096       -1.5258067      -0.04630022        0.8228854
#> CLU         0.65117509       -1.2361823       0.15316632        1.3293298
#> MMP7       -0.96023541       -0.5703534       2.56374968        1.3218661
#> CD163      -0.40891468       -0.4089147      -0.40891468       -0.4089147
#>       ATCTACTAGCGCCTTG ATCTACTGTGCTGTAT ATCTGCCGTGGCGAAT ATGAGGGAGTGTTTGC
#> IFRD1       -0.7723168       -0.1772667       -1.2434526        0.9554925
#> CEBPD       -1.8354638       -1.8354638       -1.8354638        0.6443693
#> SAT1        -1.0681510        0.7173524       -1.0151587        0.3518199
#> CLU         -0.9938728       -1.2361823       -1.2361823        0.5490019
#> MMP7        -0.9602354       -0.9602354       -0.5011091        2.9054279
#> CD163       -0.4089147        2.1162573       -0.4089147       -0.4089147
#>       ATGAGGGGTATGAAAC ATGAGGGTCACTTATC ATGAGGGTCTTAACCT ATGCGATGTAGAGGAA
#> IFRD1       -0.3934573       -1.2434526        0.7967466     -0.614064322
#> CEBPD       -1.8354638       -0.6797333       -0.2910364     -0.233238122
#> SAT1         0.1250913       -1.5546137        0.6060087     -0.004135506
#> CLU         -0.9829960       -0.6270531        1.4746518     -0.606779571
#> MMP7        -0.7269795        0.4617863        0.1332648      0.042636283
#> CD163       -0.4089147       -0.4089147       -0.4089147     -0.408914675
#>       ATGCGATGTTCCGTCT ATGGGAGGTAAGTGTA ATGTGTGAGCGTTGCC ATGTGTGAGGCGTACA
#> IFRD1       -1.2434526        1.0338250       -0.8015868       -0.7370202
#> CEBPD       -1.2859881        1.3211636       -1.2508272       -0.6310463
#> SAT1        -0.3495690        1.3170398       -0.3033833        0.8193227
#> CLU         -0.6234067        1.4559433       -1.2361823       -1.2361823
#> MMP7        -0.9602354        2.1630809        0.3210222       -0.9602354
#> CD163        3.2582503       -0.4089147       -0.4089147        0.9263849
#>       ATTATCCCAGCCAGAA ATTATCCGTCTGCAAT ATTCTACAGTGTCTCA ATTCTACGTTTACTCT
#> IFRD1      -0.53238403      -0.51185410        0.4669589        1.5102268
#> CEBPD      -0.08634301       0.09795073        0.7824349        1.8576320
#> SAT1       -1.08299298       0.34556881        0.6839554        1.5697887
#> CLU        -0.51262206       1.10262778        0.8690807       -1.2361823
#> MMP7       -0.47007695      -0.96023541        0.3730334       -0.5242257
#> CD163      -0.40891468      -0.40891468       -0.4089147       -0.4089147
#>       ATTGGACGTCGGCTCA ATTGGACGTTACGCGC ATTGGTGCATGGGACA ATTGGTGTCAGTCCCT
#> IFRD1        1.6771139      -1.24345261        0.8051698        0.7287109
#> CEBPD        0.3754002      -1.01658553        0.7104668        1.1941170
#> SAT1         1.7177965      -0.62930029        1.2178050        0.9769523
#> CLU          1.9922808       0.05623615        0.5563908        0.7542499
#> MMP7         0.7967995       1.49235912        0.2182086        1.1589569
#> CD163       -0.4089147      -0.40891468       -0.4089147       -0.4089147
#>       ATTTCTGAGGACCACA ATTTCTGCAAAGCAAT ATTTCTGGTCCATCCT CAACCAACAAACTGTC
#> IFRD1       -0.4478157        0.4310795        0.6095200       -0.5681990
#> CEBPD       -1.2331680        0.5338752        1.1283796       -0.9878044
#> SAT1        -1.6620333        0.1988566       -0.4115217       -2.3493127
#> CLU          0.2754528       -0.2837840        1.4069599       -0.6602775
#> MMP7        -0.9602354       -0.7019669        2.0103596       -0.4296650
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       CAACCAACAATTCCTT CAACCTCCATCACGTA CAACCTCGTTGCTCCT CAACCTCTCCTAGGGC
#> IFRD1       -0.5770486       -0.9669877      -0.33423325       0.05574882
#> CEBPD        0.8107234       -0.4432017       0.70624097       0.29346082
#> SAT1         0.2139689       -0.7592254      -0.48831407      -0.23296318
#> CLU          0.9389013       -0.8763700      -0.49431778      -1.09045751
#> MMP7         0.5675122       -0.5483780       0.08222572      -0.40325757
#> CD163       -0.4089147       -0.4089147      -0.40891468       0.02520984
#>       CAACTAGAGGTACTCT CAAGAAAAGACAAAGG CAAGAAACAAACTGTC CAAGAAATCATAGCAC
#> IFRD1      -0.71907792        2.0948302       0.36235987       -1.2434526
#> CEBPD      -0.09396100        0.1503357       0.78888069        0.4599085
#> SAT1       -0.04374574        2.1153218      -0.07460918       -2.1101289
#> CLU        -1.04505196        0.9161327       1.79132752       -0.1780913
#> MMP7       -0.78415055        1.3272555       0.89826543       -0.5700812
#> CD163      -0.40891468       -0.4089147      -0.40891468       -0.4089147
#>       CAAGATCAGCTCTCGG CAAGATCAGGATCGCA CAAGATCAGTCGTACT CAAGATCCACCAGATT
#> IFRD1       -0.5361707        1.7014834       -1.2434526       -1.2434526
#> CEBPD       -1.3000526        0.4655228       -0.3519588       -1.0644430
#> SAT1        -1.7624499        0.6604459       -2.8271342       -0.1656179
#> CLU         -0.8724207        0.7114134       -1.2361823       -1.2361823
#> MMP7        -0.9602354        0.8340483       -0.4849903       -0.4776349
#> CD163       -0.4089147       -0.4089147       -0.4089147        1.6331653
#>       CAAGATCCATGTCCTC CAAGATCGTCGGCATC CAAGATCTCGATCCCT CAAGATCTCTAACTGG
#> IFRD1       -1.2434526        1.1180423        0.5411139      -1.24345261
#> CEBPD       -1.8354638        1.1628613        0.2266234      -0.02544352
#> SAT1        -0.7795320        1.2445814        0.1599157      -0.46678353
#> CLU         -0.8061286        0.7574761       -1.2361823      -1.23618235
#> MMP7         1.5743554        1.6085000        2.0804447      -0.96023541
#> CD163       -0.4089147       -0.4089147       -0.4089147      -0.40891468
#>       CAAGATCTCTGTACGA CAAGGCCGTGCGATAG CAAGGCCTCGGTTAAC CAAGTTGAGCTCCTCT
#> IFRD1        0.4970403       0.06624319        0.9001627      -1.24345261
#> CEBPD       -1.1873945       0.31874053       -0.5239298      -0.78383540
#> SAT1         0.3270390       0.67625727       -0.7598681      -0.04603116
#> CLU         -1.2361823      -0.81904217        1.7005095      -1.23618235
#> MMP7        -0.5545932      -0.96023541       -0.9602354      -0.54822324
#> CD163        2.6995564       3.95119085       -0.4089147       2.28933347
#>       CACAAACCAAAGCAAT CACACAACAAGTTCTG CACACAACATTTGCTT CACACAATCAGTCCCT
#> IFRD1       -1.2434526        0.4243313       -0.5168946       -0.2113614
#> CEBPD       -1.2670883       -1.0689050        0.1916444       -1.8354638
#> SAT1        -1.5173968       -0.6347138        0.3258087        0.3689242
#> CLU         -0.6054290        0.2071061        0.2877345       -0.7053683
#> MMP7        -0.6044755       -0.7229124        0.7754072       -0.9602354
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       CACACAATCGTGGACC CACACCTCAAGAAGAG CACACCTTCTTACCGC CACACTCAGGGCTTCC
#> IFRD1       2.21677832        0.3491136      -0.05965665       -1.2434526
#> CEBPD       0.91526015       -0.6298939       1.17180543       -1.2245850
#> SAT1        2.14274644       -0.4245176       0.06615284        0.1639321
#> CLU         0.02532127        1.1127485      -0.42951280       -1.2361823
#> MMP7        1.42663858       -0.5993584      -0.25789506       -0.9602354
#> CD163      -0.40891468       -0.4089147      -0.40891468        0.8275028
#>       CACACTCAGTCCTCCT CACACTCCAACTGCTA CACAGGCAGTACACCT CACAGTAAGGTAAACT
#> IFRD1      -0.60858606       -1.2434526        0.5024443      -1.24345261
#> CEBPD      -0.42858465       -0.6291552        0.2521907      -0.47597520
#> SAT1       -0.08172575        0.2964366        0.8653409       0.02896492
#> CLU         0.59727747       -1.2361823        0.6920670      -0.86616220
#> MMP7       -0.55015588       -0.3716061       -0.1666453      -0.96023541
#> CD163      -0.40891468        1.4944864       -0.4089147       0.69340127
#>       CACAGTACAAACCCAT CACAGTAGTAGAGTGC CACAGTATCGGTTCGG CACATAGGTGGTTTCA
#> IFRD1       -0.2447041      0.895663746        0.6256469       -0.5270044
#> CEBPD       -0.3064773      1.281132836        1.6335333       -1.2931137
#> SAT1         0.2434876      0.557805456        0.5501591       -0.3642669
#> CLU         -0.7225168      0.920631960        1.2297871       -1.2361823
#> MMP7        -0.4870049     -0.005053363        1.4245801       -0.6207654
#> CD163       -0.4089147     -0.408914675       -0.4089147        2.9729186
#>       CACCACTAGACTGTAA CACCACTAGCTAAACA CACCACTTCACAACGT CACCACTTCCTCATTA
#> IFRD1       -0.4831560       -1.2434526        0.8961413      -0.02279485
#> CEBPD       -1.0621223       -0.5874936        0.6241014       0.78841567
#> SAT1        -0.8395898        0.7178064        0.5509738       0.43669909
#> CLU         -1.0128882       -1.2361823        1.5345331       0.04995461
#> MMP7         0.5241067       -0.9602354        2.0095872      -0.83722780
#> CD163       -0.4089147        2.7336725       -0.4089147      -0.40891468
#>       CACCAGGGTCGCATCG CACCTTGCATTACCTT CACCTTGGTACAGCAG CACCTTGGTCTTCTCG
#> IFRD1        0.3175759       0.73236309        0.3552432       0.70587451
#> CEBPD        0.5768155       0.15142281        0.8797238       0.76866292
#> SAT1         0.9579527       0.28857032        0.2484900      -0.06052908
#> CLU          1.4297321       0.53605220        1.6640207       1.08249610
#> MMP7         1.4958220       0.09602187        0.8126516      -0.96023541
#> CD163       -0.4089147      -0.40891468       -0.4089147      -0.40891468
#>       CACTCCACATGATCCA CAGAATCAGGACAGAA CAGAATCAGGATGCGT CAGAATCTCCATGAAC
#> IFRD1       -0.3505572       -1.2434526      -0.00845389      -0.57336703
#> CEBPD        0.2846834       -0.6697974      -0.10592597      -0.99339591
#> SAT1        -1.1199288        0.4399325      -2.13132848      -0.56662553
#> CLU         -0.7133001       -1.2361823      -0.89055592      -0.06317796
#> MMP7         0.2206283       -0.6850000      -0.27495018      -0.52920965
#> CD163       -0.4089147        1.5683734      -0.40891468      -0.40891468
#>       CAGAGAGCAGTATAAG CAGCAGCAGGACGAAA CAGCAGCCAATCTACG CAGCAGCTCCCGACTT
#> IFRD1       -0.3510966       -1.2434526       -1.2434526       -1.2434526
#> CEBPD       -1.5089707       -1.0672769       -0.3227961       -0.1335399
#> SAT1        -1.0373344        0.5120057       -1.0221451       -0.7895067
#> CLU         -0.8474567       -0.5529833       -1.2361823       -0.4368310
#> MMP7        -0.7558755       -0.9602354       -0.9602354       -0.9602354
#> CD163       -0.4089147        2.0143276       -0.4089147       -0.4089147
#>       CAGCCGAGTCTAGTCA CAGCGACAGCATCATC CAGCGACCACTCGACG CAGCGACCATGCGCAC
#> IFRD1       -0.3626070      -0.50096899       2.22000663       -0.2981100
#> CEBPD        0.3671342       0.77101838       0.87497561       -0.4160218
#> SAT1        -0.1154473       0.09683671      -0.02054718       -1.2898007
#> CLU         -0.1285870       1.20993847      -0.16598478       -0.7499839
#> MMP7         0.3116265       0.79675580      -0.44493229       -0.5123099
#> CD163       -0.4089147      -0.40891468      -0.40891468       -0.4089147
#>       CAGCTAAAGATGGGTC CAGCTAAAGTACATGA CAGCTAACAGGTGGAT CAGCTAACATGCAATC
#> IFRD1        0.3810552        1.3834100       -0.3657961       -1.2434526
#> CEBPD       -0.1424629        1.4112630        0.9176179       -0.6740033
#> SAT1        -1.3165778        0.8094536       -1.1981125       -1.4655898
#> CLU         -0.2568798        1.1679188        0.5747004       -0.4470791
#> MMP7        -0.5898463       -0.7125058        0.6657149       -0.9602354
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       CAGCTAAGTAGTAGTA CAGCTAAGTGCTAGCC CAGCTAATCCACGACG CAGCTGGAGAGGGATA
#> IFRD1      -0.03469055        0.4530140        0.6348822      -0.39813415
#> CEBPD      -1.57348884        0.3426833        0.8007147       0.01660679
#> SAT1       -0.11571329        1.0341112        0.2366407      -0.45614993
#> CLU        -0.15785047       -1.2361823        0.7080575      -0.68382046
#> MMP7       -0.66647368       -0.4371732       -0.5969654      -0.55970365
#> CD163      -0.40891468        1.2824677       -0.4089147      -0.40891468
#>       CAGCTGGAGGAACTGC CAGCTGGGTACCGTAT CAGGTGCAGCAGCCTC CAGGTGCAGCGCCTCA
#> IFRD1       -1.2434526       -0.3444518       -1.2434526        0.5704299
#> CEBPD       -0.9042852       -1.8354638       -1.8354638        1.2373661
#> SAT1         0.6414335       -0.9049029       -1.7380397        0.9691562
#> CLU         -0.4238283       -1.2361823       -0.7309944        2.1420386
#> MMP7        -0.3773883       -0.5342677       -0.4948152        1.2011581
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       CAGGTGCAGTGACTCT CATATGGAGACCACGA CATATGGAGCCTTGAT CATATGGAGGTGATAT
#> IFRD1       0.78431938       -1.2434526        1.0623586        1.5698316
#> CEBPD       0.64238476       -1.8354638        1.2125786        1.5177810
#> SAT1        0.01885519       -2.0425036        1.3483131        1.1387534
#> CLU         0.58726967       -0.6658988        1.9576485       -0.2340856
#> MMP7       -0.55489273       -0.9602354        1.2369466       -0.4257154
#> CD163      -0.40891468       -0.4089147       -0.4089147       -0.4089147
#>       CATATGGCAAAGTCAA CATATGGCACGAAACG CATATGGGTCAAAGAT CATATGGTCTCTTGAT
#> IFRD1        1.8575160       -0.7298758       -0.2571126       -1.2434526
#> CEBPD        1.0098817        0.7186376       -0.1422981       -0.5842044
#> SAT1         1.8234318       -1.0496778        0.3103943       -0.3143250
#> CLU          0.9216887       -0.6343229       -0.7288986       -1.2361823
#> MMP7         1.4172822        1.0100287       -0.4928843       -0.4532784
#> CD163       -0.4089147       -0.4089147        1.1023192       -0.4089147
#>       CATATTCGTACCGAGA CATCAAGAGAGTGAGA CATCAAGAGTCTCCTC CATCAAGAGTCTTGCA
#> IFRD1        0.8536910       -1.2434526       0.02740798       -1.2434526
#> CEBPD        1.3557161       -0.3326666      -0.87342446       -1.2553152
#> SAT1         0.2943380       -1.8822846      -0.19303239       -0.7091920
#> CLU          1.3790100       -0.3619353       0.83800676       -1.2361823
#> MMP7         1.1848709       -0.9602354      -0.96023541       -0.5971064
#> CD163       -0.4089147       -0.4089147      -0.40891468       -0.4089147
#>       CATCAAGCAATCGGTT CATCAAGGTGTGTGCC CATCAGATCCTCGCAT CATCCACCAAGCGCTC
#> IFRD1      -0.09397948       -0.1448895       -1.2434526       0.07731218
#> CEBPD      -0.96531457        0.6502903       -1.2564909       0.12890020
#> SAT1       -1.29334231        0.7309134        0.3258547      -0.10498452
#> CLU         0.65748560       -1.2361823       -1.2361823      -0.96480150
#> MMP7       -0.96023541       -0.9602354       -0.5978423      -0.39363468
#> CD163      -0.40891468        2.8031015        2.0389867      -0.40891468
#>       CATCCACCAGCCTATA CATCCACCAGGCTCAC CATCGGGAGAATAGGG CATCGGGAGCGTCAAG
#> IFRD1        1.6776290        0.5410952       -1.2434526       -1.2434526
#> CEBPD        0.9533245        1.0362349       -0.8286685        0.2235060
#> SAT1         1.0302695        0.5403733        0.6003652       -0.1686905
#> CLU          0.1886185        1.5344493       -0.8114139       -0.4542872
#> MMP7         0.1791160        2.0127297       -0.5689041       -0.3155727
#> CD163       -0.4089147       -0.4089147        3.2863170       -0.4089147
#>       CATCGGGAGTACTTGC CATCGGGCAATGCCAT CATCGGGTCACGGTTA CATCGGGTCATACGGT
#> IFRD1      1.774818704      -0.53165711       1.72433741       -0.1541215
#> CEBPD      1.089724443      -0.58333756      -0.01589387        0.7248315
#> SAT1       0.849279900      -0.08387271      -0.51984948       -0.7624565
#> CLU        0.005199599      -1.23618235       1.43757377       -0.1073208
#> MMP7      -0.748882932      -0.76913307      -0.29796339        1.4025551
#> CD163     -0.408914675      -0.40891468      -0.40891468       -0.4089147
#>       CATCGGGTCGGTGTTA CATCGGGTCTCCTATA CATGACAGTTGTACAC CATGCCTGTCGAGTTT
#> IFRD1     -0.026935088      -1.24345261        1.0764081       -1.2434526
#> CEBPD     -0.195583280       0.38401685        1.6329848       -0.5780629
#> SAT1       0.001531876      -0.52606363        0.8297055        0.1766810
#> CLU       -0.515332833      -0.08808249        1.0398371       -1.2361823
#> MMP7       0.793301727       0.93171505        0.5160075       -0.9602354
#> CD163     -0.408914675      -0.40891468       -0.4089147       -0.4089147
#>       CATGGCGTCTTTCCTC CATTCGCTCCAGGGCT CCAATCCAGGTGATAT CCAATCCAGTTGTAGA
#> IFRD1       -0.3132450       0.32298458        0.3120031       -1.2434526
#> CEBPD        1.0627573      -0.12262039        0.5015509       -0.1596640
#> SAT1        -1.3849973       0.03759641        0.9338291        0.4558840
#> CLU         -0.5778754      -0.22816061        1.5018723       -0.7364829
#> MMP7        -0.2635265       2.33906395       -0.2232239       -0.9602354
#> CD163       -0.4089147      -0.40891468       -0.4089147        1.9277663
#>       CCAATCCCATAACCTG CCAATCCGTATATCCG CCAATCCGTCTGCCAG CCAATCCGTGTAATGA
#> IFRD1       -1.2434526        1.0740517     -0.240260351        0.8217883
#> CEBPD       -1.8354638       -0.3125835     -0.006704607        1.2996066
#> SAT1        -0.1027505        0.5428505     -0.586252561        1.0001925
#> CLU          1.2447276        0.9789430      0.318511133        1.3772226
#> MMP7        -0.9602354        0.9539957      1.145489490        2.0993928
#> CD163       -0.4089147       -0.4089147     -0.408914675       -0.4089147
#>       CCACCTAAGGGTTCCC CCACCTACAGTCGTGC CCACGGAAGGAGTAGA CCACGGACACTTAACG
#> IFRD1       -0.8414331        0.4758079      -0.39237918       -0.4328636
#> CEBPD       -0.4006241       -1.4389654       0.02436131        0.1394162
#> SAT1        -1.7596039       -1.6105548      -1.60099964       -0.5268299
#> CLU          0.5102895        0.1382987      -0.34500591        0.4107570
#> MMP7         1.7157968       -0.8239131      -0.31383329        0.3562013
#> CD163       -0.4089147       -0.4089147      -0.40891468       -0.4089147
#>       CCACTACCAGGATTGG CCACTACTCTATCCCG CCAGCGAAGGCATGGT CCAGCGACAGTCACTA
#> IFRD1        1.2412937       -0.3892804       0.20479272       -1.2434526
#> CEBPD        1.0800047       -0.6974992      -0.04701346       -1.1208026
#> SAT1         1.4554130       -1.0009158       0.59227800       -1.2912197
#> CLU          0.4584780        0.5117896      -0.76632829       -1.2361823
#> MMP7         1.5469399       -0.9602354      -0.96023541       -0.9602354
#> CD163       -0.4089147       -0.4089147       2.39139928       -0.4089147
#>       CCAGCGATCTGCTGCT CCATGTCCACCACGTG CCATGTCGTTCGTTGA CCATTCGAGCTCTCGG
#> IFRD1      -1.24345261        0.6118002       -1.2434526        0.1047531
#> CEBPD      -0.59258146        0.4541997       -0.2995372        0.1396329
#> SAT1        0.41748148       -0.4338484       -0.1256670        1.2865297
#> CLU        -1.23618235        0.5713573       -0.3985798       -1.2361823
#> MMP7       -0.04915944       -0.5186333       -0.4622241        1.0062999
#> CD163      -0.40891468       -0.4089147        2.0863621       -0.4089147
#>       CCATTCGAGGTGCTTT CCATTCGCATAGGATA CCATTCGGTTAAAGAC CCCAATCAGTGTGGCA
#> IFRD1       -1.2434526       -0.2606131       -0.0615830        1.3911152
#> CEBPD       -0.5536543        0.6965922       -0.3877454        0.6391345
#> SAT1        -0.4110794        0.2110782       -1.0377807        1.4653224
#> CLU         -0.6169423        1.3754107       -0.7369699        1.3732939
#> MMP7         0.4758388       -0.3487564        0.7674998        1.3523531
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       CCCAATCTCCGTACAA CCCAGTTAGCTCCCAG CCCAGTTCACCACGTG CCCAGTTGTATAGGGC
#> IFRD1      -0.19841221       -0.5023887        0.1118468       -1.2434526
#> CEBPD       0.93073656       -0.4472079       -0.2598620       -0.6124467
#> SAT1       -0.84018085       -0.4280956        0.1813742        0.8302549
#> CLU        -0.37321765        0.4460283        0.4931004       -1.2361823
#> MMP7        0.09162471       -0.6091019       -0.7016364       -0.9602354
#> CD163      -0.40891468       -0.4089147       -0.4089147        1.1860720
#>       CCCATACAGTCGTACT CCCATACCATTTGCCC CCCATACGTCTGGTCG CCCTCCTCACTCTGTC
#> IFRD1        0.7907053       -0.4785342       0.76815918        1.2899177
#> CEBPD        0.6731462       -0.5153782       0.84075121        1.4638109
#> SAT1        -1.0782388       -0.7501363       0.04880243        1.2547098
#> CLU          0.6706857       -1.0113655       1.15969920        0.5888295
#> MMP7         1.1549256       -0.5977991       0.29006893        2.6590058
#> CD163       -0.4089147       -0.4089147      -0.40891468       -0.4089147
#>       CCGGGATGTACTCGCG CCGGGATTCCAGGGCT CCGGTAGAGCGATAGC CCGGTAGGTACCAGTT
#> IFRD1       -0.7609983       -0.2607281       -1.2434526        0.1113697
#> CEBPD       -1.3279337       -1.0915430       -1.2033072       -0.8098656
#> SAT1        -1.5637353       -1.4613210        0.2820445        0.3932204
#> CLU         -0.9880516       -1.2361823       -0.9880149        1.5235254
#> MMP7        -0.7316372       -0.4945974       -0.9602354        0.6690954
#> CD163       -0.4089147       -0.4089147        2.8342559       -0.4089147
#>       CCGGTAGGTATGAATG CCGGTAGGTCATCGGC CCGGTAGTCAATCTCT CCGGTAGTCCGATATG
#> IFRD1      -0.02113029        0.7161805       -0.6492214      -1.24345261
#> CEBPD       0.40159767       -0.8404901        0.1952534       0.10416758
#> SAT1        0.73215774       -0.9054495       -0.8404645      -2.25121049
#> CLU        -1.23618235        0.1605323       -1.2361823      -0.22883262
#> MMP7       -0.96023541       -0.9602354       -0.4852365      -0.03218286
#> CD163       4.55980395       -0.4089147       -0.4089147      -0.40891468
#>       CCGGTAGTCCTGCAGG CCGTACTGTTCCGTCT CCGTACTGTTCCTCCA CCGTGGATCTGCCCTA
#> IFRD1       -1.2434526       -1.2434526       -0.2269893        0.2712445
#> CEBPD       -0.7726674       -0.6833539       -1.8354638        0.8542029
#> SAT1         0.1681212        0.7589349        0.8634587        0.9611357
#> CLU         -1.2361823       -0.7377113       -1.2361823        0.7647207
#> MMP7        -0.5430338       -0.5010033       -0.9602354        0.7587027
#> CD163        2.7606752       -0.4089147        3.0891828       -0.4089147
#>       CCGTTCAAGATGTTAG CCGTTCATCCCGGATG CCGTTCATCCGGCACA CCTAAAGAGGTCATCT
#> IFRD1       -1.2434526       -0.3433365        1.2385719       -0.1897208
#> CEBPD       -0.7889640        0.2223992        0.6566324        0.4609956
#> SAT1        -1.8488997       -1.5481467       -0.6869018       -0.2577916
#> CLU         -1.2361823       -0.3966766       -1.2361823        0.2028247
#> MMP7        -0.9602354        1.2581711       -0.4928111       -0.7465492
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       CCTAAAGCAATGCCAT CCTAAAGCAATTCCTT CCTAAAGCATGGAATA CCTAAAGGTGATGTGG
#> IFRD1        1.2719676        1.0017130        0.8800898    -1.2434526134
#> CEBPD        1.8764102        0.3418226        1.0223109    -0.0008522143
#> SAT1         1.6677795        0.2024976        1.2111420    -0.1912208411
#> CLU          1.9984886        1.2804152        1.5849435    -1.2361823459
#> MMP7         0.8532990        1.5617689        2.2205218    -0.4751956573
#> CD163       -0.4089147       -0.4089147       -0.4089147    -0.4089146752
#>       CCTAAAGTCCGAACGC CCTACACCAGATAATG CCTACACGTAAGGATT CCTACCATCTTCCTTC
#> IFRD1       -1.2434526      0.442372224       -1.2434526       -1.2434526
#> CEBPD       -1.8354638      0.867969331       -1.8354638       -1.1886979
#> SAT1         0.4657272     -0.047141837        0.3535529       -0.0896509
#> CLU          1.2588905      1.441804290       -1.2361823       -0.9815158
#> MMP7        -0.9602354     -0.002957117       -0.9602354       -0.7256158
#> CD163       -0.4089147     -0.408914675        0.8429875        1.6881709
#>       CCTAGCTAGCCACGCT CCTAGCTGTCACTGGC CCTAGCTGTGTGACGA CCTATTAAGTGGAGAA
#> IFRD1      -0.01401190        0.2663105      -1.24345261      -0.03872325
#> CEBPD      -1.18214054       -1.1087790      -1.83546381      -0.33731357
#> SAT1       -0.30787920       -0.2394763       0.05937639      -0.95017735
#> CLU         0.03829627        0.6669232      -0.44445415       0.28835076
#> MMP7       -0.66033261       -0.6923314      -0.96023541      -0.38940632
#> CD163      -0.40891468       -0.4089147       3.00890527      -0.40891468
#>       CCTCAGTTCGGAATCT CCTCTGATCATGCAAC CCTTACGCACGAAACG CCTTACGCACGGCTAC
#> IFRD1       -1.2434526        0.9687905      -1.24345261        0.1110359
#> CEBPD       -1.8354638        0.4955581       0.07195947       -1.1064421
#> SAT1         0.9342097        0.5450104      -1.82764991       -0.6163078
#> CLU         -1.2361823       -0.5223447      -0.50973936        0.4952157
#> MMP7        -0.9602354       -0.9602354      -0.54003641       -0.6224379
#> CD163       -0.4089147        2.1145554      -0.40891468       -0.4089147
#>       CCTTACGCATAGAAAC CCTTACGCATTCTCAT CCTTCCCAGTTACCCA CCTTCCCGTTGCGCAC
#> IFRD1       -1.2434526       -1.2434526       0.74163304       -1.2434526
#> CEBPD        0.5055738       -0.4684970      -0.57592638        0.4210315
#> SAT1        -0.7188853       -2.1909853       0.07583167       -0.1505664
#> CLU         -0.8164060       -0.9663649       0.66888482       -0.5584445
#> MMP7        -0.9602354        0.1476494      -0.84910928        1.5615012
#> CD163       -0.4089147       -0.4089147      -0.40891468       -0.4089147
#>       CCTTCCCGTTTGTTTC CCTTCGACAACACGCC CCTTCGACACCGGAAA CCTTTCTCATCGGGTC
#> IFRD1       -0.8586860       -0.2295595        1.6231937        0.8729482
#> CEBPD        0.7663810        1.1892758        1.2897661        0.7202808
#> SAT1        -1.0913476       -1.2143520        1.8479274        0.9561709
#> CLU         -0.7604767       -1.0149446        1.2290181        0.7167628
#> MMP7         1.0004006       -0.9602354        1.2485078       -0.4329644
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       CCTTTCTGTAACGTTC CGAATGTGTCTTCTCG CGAATGTGTGACGCCT CGAATGTGTGCTGTAT
#> IFRD1        0.3606059        0.9916754        0.4484063       -1.2434526
#> CEBPD        1.0445393        0.8613633       -0.5547294       -1.8354638
#> SAT1         0.4043634        1.1074461       -0.9443637       -1.0255457
#> CLU          0.3757775        0.8812814       -0.1567234       -0.8947127
#> MMP7        -0.9602354       -0.1296367        0.3004747       -0.9602354
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       CGACCTTCACCCTATC CGACCTTCACCTATCC CGACTTCTCTACCAGA CGAGAAGAGAACTCGG
#> IFRD1       1.36541671     -0.009402868       0.53454456        0.7020580
#> CEBPD       0.03893169     -1.835463810       0.29929096        2.8826585
#> SAT1        1.29417971     -1.510664927      -0.37541061        1.4781387
#> CLU         1.79791361     -0.847204853       0.08935779       -0.1282508
#> MMP7        1.03981832     -0.960235411      -0.40670950        0.5235275
#> CD163      -0.40891468     -0.408914675      -0.40891468       -0.4089147
#>       CGAGAAGAGCGCTCCA CGAGCACCACCACCAG CGAGCACTCAACTCTT CGAGCCACACATGGGA
#> IFRD1       0.33151453       -0.3218531        0.5013402       -0.6394262
#> CEBPD       1.33866166       -0.7310393        1.2447498       -0.3367796
#> SAT1       -0.61784805        0.5541654       -0.3387846        0.6351124
#> CLU         0.08219622       -1.2361823        1.3536560       -0.9255260
#> MMP7       -0.55253855       -0.5235599       -0.3694660       -0.8009000
#> CD163      -0.40891468        1.0031262       -0.4089147        2.6244193
#>       CGAGCCATCCCAAGTA CGAGCCATCTTGCATT CGATCGGCACGGTGTC CGATCGGCATCTACGA
#> IFRD1        1.1566886      -0.93214755       -0.5783180       0.70871008
#> CEBPD        1.1752803       0.54563506       -1.8354638      -0.91746684
#> SAT1         0.6411100       0.28012069       -1.2204943      -0.07373254
#> CLU          0.8510509      -0.01393669       -0.8940975      -1.01922498
#> MMP7         2.2072712       0.21846246       -0.6450790      -0.76035663
#> CD163       -0.4089147      -0.40891468       -0.4089147      -0.40891468
#>       CGATCGGTCTGCGTAA CGATGGCCACCAGTTA CGATGTAAGAACTGTA CGATGTAGTAAATGTG
#> IFRD1       -1.2434526        1.0896973       -0.5437507       -0.4345128
#> CEBPD       -1.0652549        1.3257170       -1.2267798        1.2201817
#> SAT1        -1.4255802        1.1802595       -0.3530571       -1.2822902
#> CLU         -1.2361823        0.9343698        0.3232482        0.2888685
#> MMP7        -0.9602354        1.0394539       -0.8060141       -0.5769407
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       CGATGTATCAGTGTTG CGCCAAGAGAGTCGGT CGCCAAGCAATCAGAA CGCGTTTAGGTGATAT
#> IFRD1        0.6971757       -0.7152209        1.3209604      -0.98419417
#> CEBPD        1.3853287       -1.1508085        1.4197607       0.05631385
#> SAT1         0.6302670       -1.7607248        1.8101730      -1.47462417
#> CLU          1.5575199       -0.7171066        1.2513822      -0.07195296
#> MMP7         2.0349025       -0.1549194        2.1612955       1.37400724
#> CD163       -0.4089147       -0.4089147       -0.4089147      -0.40891468
#>       CGCGTTTCACCGATAT CGCGTTTTCATCGGAT CGCTATCAGCTCTCGG CGCTATCAGTACATGA
#> IFRD1      -0.66087682        1.1094938        0.9612284       -0.3924360
#> CEBPD       0.39487426        1.0498910        2.1667049        0.1558453
#> SAT1        0.16305003        1.0181775        0.7758921        0.5948638
#> CLU        -1.06996825        1.7569925        0.8793187       -0.5345882
#> MMP7        0.03153646        1.9428152       -0.9602354       -0.5570037
#> CD163       0.48368651        0.2029871       -0.4089147        0.8949817
#>       CGCTATCCAGCTTCGG CGCTGGAGTAAATGAC CGCTTCAGTGACGCCT CGCTTCATCTCTTATG
#> IFRD1       1.28244729        2.2222261        0.8749569       -0.3353180
#> CEBPD      -0.01237523        1.4619789        0.4088566       -1.8354638
#> SAT1        0.65110031        2.1283205        0.6471600        0.3398529
#> CLU         1.22216966        1.3229775        1.2202188       -1.2361823
#> MMP7        0.33503412        1.9169093        0.9694105       -0.9602354
#> CD163      -0.40891468       -0.4089147       -0.4089147       -0.4089147
#>       CGGACGTTCCACGCAG CGGACTGGTTTAGCTG CGGACTGTCAGTTTGG CGGAGCTCAATCTGCA
#> IFRD1       -0.4572578        0.9357561        1.8986401       0.67060438
#> CEBPD        0.3753398        0.8844800        0.7925067      -0.45237697
#> SAT1        -0.5842347        1.3749747        0.3045028      -0.41465771
#> CLU         -0.4527639        1.0833872        1.3120683       0.08313015
#> MMP7        -0.6890569        1.3460032        0.4920504      -0.58708429
#> CD163       -0.4089147       -0.4089147       -0.4089147      -0.40891468
#>       CGGAGTCAGCAACGGT CGGAGTCCATGGTCTA CGGAGTCTCAGTTGAC CGGCTAGGTCCGAATT
#> IFRD1       -1.2434526        0.5654237        1.0159046       -0.1178370
#> CEBPD       -0.1175654        1.7101469        0.6287055       -0.1677698
#> SAT1         1.1366046        1.7688583       -0.2193578        0.9033436
#> CLU         -1.2361823        1.2536809        1.4612071       -0.8342516
#> MMP7        -0.9602354        1.9813728       -0.9602354       -0.5899441
#> CD163        1.9997965       -0.4089147       -0.4089147        0.7884651
#>       CGGGTCAGTACCGGCT CGGGTCATCATTGCGA CGGTTAAAGGCAATTA CGGTTAAGTTACCGAT
#> IFRD1       0.01848972       -0.7152138       -1.2434526        0.2558881
#> CEBPD       0.31515299       -1.8354638       -1.4535198        1.1046366
#> SAT1        0.96881279       -0.6012442       -0.7840330        0.4596767
#> CLU         0.83092434       -0.9645043       -1.2361823       -0.8740704
#> MMP7        0.46757074       -0.9602354       -0.9602354        0.6917374
#> CD163      -0.40891468       -0.4089147       -0.4089147       -0.4089147
#>       CGGTTAATCGCCATAA CGGTTAATCTAGAGTC CGTAGCGAGATGCCAG CGTAGCGAGCGTGAAC
#> IFRD1       -0.2383684       -1.2434526       -0.1149944       -1.2434526
#> CEBPD       -1.0746167        0.6480918        0.8824945       -1.1048072
#> SAT1         1.1783658       -0.4807410       -1.7927292       -1.4795198
#> CLU         -1.2361823       -0.6721012        0.4396749        1.0984687
#> MMP7        -0.4840029        1.7851215       -0.9602354        0.5219315
#> CD163        3.0624152       -0.4089147       -0.4089147       -0.4089147
#>       CGTAGCGAGGATTCGG CGTAGCGCAAGCCGTC CGTAGGCCATGTCCTC CGTAGGCTCTGTCTCG
#> IFRD1        2.3847092        0.3106681        0.5940720      -1.24345261
#> CEBPD        0.7489771       -0.5761213        1.0596204      -0.03901386
#> SAT1        -0.2597823        0.1947098        1.7539980       0.34815440
#> CLU          0.5769892        0.2003019        0.9537934      -1.23618235
#> MMP7        -0.9602354       -0.9602354       -0.9602354      -0.96023541
#> CD163       -0.4089147       -0.4089147       -0.4089147       2.13542387
#>       CGTCACTGTACAGTGG CGTCACTGTTAAGGGC CGTCACTTCTCCAACC CGTCAGGAGCCGGTAA
#> IFRD1       -1.2434526        1.5890216      -0.61380418      -0.86358254
#> CEBPD       -1.8354638        1.9961941      -1.35882100       0.87040575
#> SAT1        -0.2158712        1.4911534      -0.25729288      -1.20683570
#> CLU         -0.4471817        1.7376577       0.63162407       0.05179216
#> MMP7        -0.2333438        2.3245777      -0.08453627       0.30000443
#> CD163       -0.4089147       -0.4089147      -0.40891468      -0.40891468
#>       CGTCAGGAGGGTGTTG CGTCAGGTCAGAGCTT CGTCCATCAGCTGCTG CGTCCATCATGTCTCC
#> IFRD1        0.9310491        1.0325396       -1.2434526        0.7796833
#> CEBPD        1.2748965       -1.1826743       -1.0868718       -1.8354638
#> SAT1         0.4811204       -0.9746032       -1.4549387        0.3730393
#> CLU          1.4902299        0.3291106       -1.2361823       -1.2361823
#> MMP7        -0.5117231       -0.7639039       -0.9602354       -0.9602354
#> CD163       -0.4089147       -0.4089147       -0.4089147        1.4053759
#>       CGTCCATTCGTTGACA CGTGAGCAGACTTGAA CGTGAGCGTCAATGTC CGTGAGCGTCCGAACC
#> IFRD1       -0.2431639        0.2854154      -0.01912117        0.5179812
#> CEBPD       -1.8354638       -1.2888091      -0.64480338       -1.2608190
#> SAT1        -0.1276124       -0.3553727       0.61718082        0.5427516
#> CLU          0.3686157        0.5056423      -1.23618235        0.2498474
#> MMP7        -0.4862751       -0.9602354      -0.96023541       -0.8698556
#> CD163       -0.4089147       -0.4089147       2.42301615       -0.4089147
#>       CGTGTAACAAGCTGAG CGTGTAACAGCGTAAG CGTGTAACAGGCTGAA CGTGTAAGTTAAGATG
#> IFRD1       -1.2434526       -0.7884182       -1.2434526       -0.1332173
#> CEBPD       -0.9203221        0.1332136       -1.0017949       -0.9045316
#> SAT1        -0.9775578        0.4864291        0.5892546       -0.7840856
#> CLU          0.1311304        1.2831625       -1.2361823        0.3424005
#> MMP7        -0.9602354        0.3399314       -0.9602354       -0.7323155
#> CD163       -0.4089147       -0.4089147        1.2784295       -0.4089147
#>       CGTGTAATCAGTGTTG CGTGTAATCGCCCTTA CGTGTCTCATCGGGTC CGTTAGAAGACACTAA
#> IFRD1       -0.4656865       -0.7021211        0.3800364        2.1023392
#> CEBPD       -1.2466961        0.1205591        0.8235194        1.5395258
#> SAT1         0.4853131        0.5575698       -0.8547884        0.3946216
#> CLU         -0.5861815       -0.9577706       -0.9148481        0.7678400
#> MMP7        -0.9602354        0.6016261       -0.3127486       -0.9602354
#> CD163        1.5274833        1.0069364       -0.4089147       -0.4089147
#>       CGTTAGAAGCTGTTCA CGTTAGACACCGAAAG CGTTAGATCAAACAAG CGTTAGATCAGCGACC
#> IFRD1       -0.3122325       -1.2434526        2.1081001        0.5321268
#> CEBPD        0.7105065       -1.8354638        0.2114147        1.4042986
#> SAT1        -1.6865481        0.1903180        1.5961222        1.0060636
#> CLU         -0.3116715       -0.3993406        0.8192446        1.8177368
#> MMP7        -0.7764288        0.6510861        2.1507251        0.5080492
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       CGTTAGATCGTCCAGG CGTTCTGAGTTGTCGT CGTTCTGCACGAAGCA CGTTCTGGTAAGAGAG
#> IFRD1      -0.05469683       -1.2434526      1.650099326        1.3912105
#> CEBPD      -0.29932702       -0.2961091     -0.083641658       -0.6455651
#> SAT1       -1.74319885        0.2394671      1.096327906       -0.7098625
#> CLU        -0.86385224       -0.8626756      1.386851163        0.3375644
#> MMP7       -0.96023541       -0.6161306      0.009324389       -0.8184436
#> CD163      -0.40891468        0.7037881     -0.408914675       -0.4089147
#>       CGTTGGGAGACCCACC CGTTGGGCAGTGGAGT CGTTGGGCATGAAGTA CTAACTTCATCATCCC
#> IFRD1       -0.5186951        0.4567961      -0.02569043        0.9081882
#> CEBPD       -1.6676083       -0.1386838       1.23968902        1.0194752
#> SAT1        -1.0295983        1.9477276      -0.17199233        0.2740879
#> CLU          0.7329659        1.4719511      -0.20964520        0.9655564
#> MMP7        -0.4969690        2.1648410       1.84598598        0.8016870
#> CD163       -0.4089147       -0.4089147      -0.40891468       -0.4089147
#>       CTAACTTTCTCCAACC CTAAGACAGAAACGAG CTAAGACTCCCTAACC CTAATGGCATACAGCT
#> IFRD1       -0.3296479        1.3148210        1.6787614       -0.2524565
#> CEBPD       -0.2287931        0.2948819        0.9840290        0.5663779
#> SAT1         0.4460097        0.2253912        1.9661789        1.3081680
#> CLU          1.4826560        0.9824474        1.3692861       -1.2361823
#> MMP7        -0.4309282        0.3731993        1.3953418       -0.4906781
#> CD163       -0.4089147       -0.4089147       -0.4089147        3.4050439
#>       CTAATGGGTATTAGCC CTAATGGGTTCCTCCA CTAATGGTCATTGCCC CTACACCTCCTATGTT
#> IFRD1       -0.8048037      -0.38690738       1.12341095       2.14853576
#> CEBPD       -0.6230432      -1.03952144      -0.04375204       1.20415813
#> SAT1        -2.2841813      -0.02096361       1.25303906       1.67423312
#> CLU         -1.0105812      -0.34068255       1.16519750       0.08794857
#> MMP7        -0.9602354       1.52325697      -0.04232296       2.22662139
#> CD163       -0.4089147      -0.40891468      -0.40891468      -0.40891468
#>       CTACATTAGTAGATGT CTACCCAGTCTCACCT CTACGTCTCAATACCG CTAGAGTTCTGCCAGG
#> IFRD1      -0.85126845       -0.4630349       -1.2434526      -0.27767783
#> CEBPD       0.59403348       -0.3957167       -0.5369737      -0.49343287
#> SAT1        0.04172818       -0.5477632       -1.3236767      -0.15692713
#> CLU         0.10637192       -1.2361823       -1.2361823       0.03727288
#> MMP7        0.24596281       -0.3596639       -0.9602354      -0.81091619
#> CD163      -0.40891468        0.7868126       -0.4089147      -0.40891468
#>       CTAGCCTAGCAATCTC CTAGCCTCAAGCCCAC CTAGCCTCACCTCGGA CTAGTGAAGAGGACGG
#> IFRD1       -0.2221125        0.4111992      -0.16637512      0.762778691
#> CEBPD        0.3790577        0.2006401       0.06143463      1.011081930
#> SAT1        -0.1012989        0.3808835      -0.03311952      0.009257211
#> CLU         -0.7108977        1.6083424      -0.32315177      1.360343717
#> MMP7         0.3407543        0.8060369       0.09613444     -0.594748408
#> CD163       -0.4089147       -0.4089147      -0.40891468     -0.408914675
#>       CTAGTGACAACCGCCA CTAGTGATCCTCCTAG CTCACACGTGTGGTTT CTCACACTCAAGATCC
#> IFRD1      -0.60293046        1.7398660        1.4665733       -0.2877459
#> CEBPD       0.05312138        1.2079374        1.0974369        0.2955058
#> SAT1       -2.37388935        1.3107167       -0.1286369        0.2702549
#> CLU        -1.23618235        1.5776240        1.8862087        0.6845751
#> MMP7       -0.96023541        2.0825532        0.3238386        0.1921521
#> CD163      -0.40891468       -0.4089147       -0.4089147       -0.4089147
#>       CTCACACTCACATGCA CTCACACTCCTAGTGA CTCAGAAAGCGGCTTC CTCAGAAAGGGCTTGA
#> IFRD1        0.1126862      -0.24289870      -1.24345261        0.9525271
#> CEBPD        0.9681365       0.06050273      -0.16949275        1.1007114
#> SAT1        -0.3135727      -0.45215799      -0.01355092        0.3052460
#> CLU          0.9130826       1.01239842      -1.23618235        1.2526575
#> MMP7        -0.6670664      -0.60807660      -0.96023541        1.0747018
#> CD163       -0.4089147      -0.40891468       1.06699912       -0.4089147
#>       CTCAGAAGTGTGACGA CTCATTACAGGCGATA CTCATTAGTCTACCTC CTCATTAGTCTAGCCG
#> IFRD1        2.1386630      -0.05599331        1.1055190       -0.4436400
#> CEBPD        1.6400144      -0.75627125        1.1520387        0.3956914
#> SAT1         0.9128611      -0.35399163        1.0299470        0.3539312
#> CLU          0.9050383       1.03693941        1.2885106       -0.2410415
#> MMP7         0.1045803      -0.61765139        0.4818936        1.6352008
#> CD163       -0.4089147      -0.40891468       -0.4089147       -0.4089147
#>       CTCATTATCAGCTTAG CTCCTAGAGGCATGGT CTCCTAGAGTAGGCCA CTCGAAAAGTTAGCGG
#> IFRD1        0.5255253        0.3811685      -0.08953655        0.9997207
#> CEBPD        0.4095716       -1.0421527       0.77020844       -0.4550958
#> SAT1         0.8711165        0.9835924      -0.44077181       -1.2398248
#> CLU          0.9851425       -1.2361823      -0.46912689        0.2126872
#> MMP7         0.1534709       -0.9602354       1.30318151        0.3168662
#> CD163       -0.4089147        2.0802699      -0.40891468       -0.4089147
#>       CTCGGAGGTCCGAATT CTCGGGATCGTGGTCG CTCGGGATCTGAGGGA CTCGTACAGAGCAATT
#> IFRD1        0.4063550       -1.2434526       -0.7455111        1.5115720
#> CEBPD       -0.4225868       -0.2529235       -0.5149993        0.5924726
#> SAT1         0.5785402        1.0417498       -1.1013794        0.9869921
#> CLU          0.8072959       -0.3699190       -0.9800864        1.3860564
#> MMP7        -0.6575201       -0.4416760       -0.1868700       -0.9602354
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       CTCGTACCAGACGCTC CTCGTCAGTCTAAAGA CTCTAATAGGATCGCA CTCTAATAGTGAATTG
#> IFRD1       -0.4653088       -0.2954877      -1.24345261        1.6110677
#> CEBPD       -0.0759514        0.9162479      -0.01215722       -1.0408005
#> SAT1         0.5685033        0.3293083       0.48816861        0.7877442
#> CLU         -0.8359758       -0.3616254      -1.23618235        0.8253658
#> MMP7        -0.5915326       -0.5110674      -0.96023541       -0.7475401
#> CD163       -0.4089147       -0.4089147       2.18214852       -0.4089147
#>       CTCTAATGTTTAGCTG CTCTACGAGTGTTTGC CTCTACGCACTTGGAT CTCTACGGTGACTACT
#> IFRD1       0.07148321      -1.24345261        0.8728086        0.6197807
#> CEBPD      -1.21859049      -0.15101804        1.3285713       -0.4249993
#> SAT1       -2.25050080      -0.01330791        1.4344990        0.8936707
#> CLU        -1.23618235      -0.54724816        1.4299396       -1.2361823
#> MMP7       -0.96023541      -0.96023541        2.0686599       -0.9602354
#> CD163      -0.40891468      -0.40891468       -0.4089147        3.2713393
#>       CTCTGGTAGTGCCATT CTCTGGTCATTAACCG CTCTGGTGTCGGCATC CTCTGGTGTTAAGATG
#> IFRD1       -1.2434526        0.4389484       -0.2054667       -0.6534898
#> CEBPD       -0.1651725       -0.8362756        0.6268856       -1.3888630
#> SAT1        -0.9869107       -1.6366373        0.9898087       -0.9132651
#> CLU         -1.2361823       -0.5573279       -1.2361823        0.7004427
#> MMP7        -0.9602354       -0.7367098       -0.9602354       -0.9602354
#> CD163       -0.4089147       -0.4089147        1.1814492       -0.4089147
#>       CTCTGGTTCGTGGACC CTGAAACAGAAACCTA CTGAAACAGGCATTGG CTGAAACGTCCAGTGC
#> IFRD1       -1.2434526       -1.2434526        0.8710163         1.390338
#> CEBPD       -1.1311016       -0.1470862        1.4377274         1.205103
#> SAT1        -2.1687189       -1.2488928        1.8407895         1.506535
#> CLU         -0.2831814       -1.2361823        1.5849992         1.681909
#> MMP7        -0.5193581       -0.4948152        1.2624810         1.754566
#> CD163       -0.4089147       -0.4089147       -0.4089147         1.329388
#>       CTGAAACTCATGTCCC CTGAAACTCCTGCCAT CTGAAGTGTGAAGGCT CTGAAGTTCAACACTG
#> IFRD1        1.9004652       1.46478316       -0.4915699       -0.8862969
#> CEBPD        1.0753581       0.15486900       -1.2662898       -0.5878702
#> SAT1         2.0856750      -0.05113497        0.5469342       -0.5410269
#> CLU          1.5167124       0.44217161       -1.2361823       -1.0524938
#> MMP7         2.6604421      -0.96023541       -0.6039757       -0.9602354
#> CD163       -0.4089147      -0.40891468       -0.4089147        1.2246214
#>       CTGAAGTTCCCATTTA CTGATAGGTCCGTCAG CTGATCCAGCCGATTT CTGATCCGTAGTACCT
#> IFRD1       -0.3148615        0.7394200     -0.001910622       -0.3727747
#> CEBPD       -0.3166482        0.9060821     -0.329697166       -0.7826815
#> SAT1        -1.9597240        0.9368593     -0.118695489       -0.2974134
#> CLU         -0.4810478        0.6005505      0.030780815       -1.2361823
#> MMP7         0.2677227       -0.6733386     -0.599326315       -0.3012736
#> CD163       -0.4089147       -0.4089147     -0.408914675        0.9251060
#>       CTGCGGAAGAGTCGGT CTGCGGATCGGCTACG CTGCTGTGTAATCGTC CTGCTGTGTAGCACGA
#> IFRD1        0.9489207       -0.8437523        1.3418370      -0.14286533
#> CEBPD        0.5025640       -0.6979971        2.3880471       0.07479507
#> SAT1        -0.8563014       -1.4148689        2.1918548      -1.45443048
#> CLU          1.1541658       -0.7452601        0.7109079      -1.02145396
#> MMP7        -0.7380584       -0.6256871        1.3160647      -0.04106789
#> CD163       -0.4089147       -0.4089147       -0.4089147      -0.40891468
#>       CTGGTCTTCGCTGATA CTGTGCTGTGGGTCAA CTGTGCTGTGTGACCC CTTAACTAGAGGTTAT
#> IFRD1       -0.8044945       -1.2434526       0.49506067       -0.3001548
#> CEBPD       -0.8892214        0.2434406      -0.51941215       -0.1282938
#> SAT1        -1.4966492       -1.3228536      -0.07067892       -0.9756373
#> CLU         -1.0104221       -1.2361823       0.62025683       -0.8776928
#> MMP7        -0.5964437       -0.9602354      -0.70793715       -0.4149811
#> CD163       -0.4089147       -0.4089147      -0.40891468       -0.4089147
#>       CTTAACTAGGGTGTTG CTTAACTTCAACGAAA CTTACCGGTAGCGTAG CTTACCGTCCATTCTA
#> IFRD1       1.20280578        0.2469903       -1.2434526        0.3692358
#> CEBPD       1.12177695        0.7378580       -0.7444977       -0.9458124
#> SAT1       -0.01211792        0.3714192        0.3083855       -0.8373595
#> CLU         0.77298481        0.1492322       -1.2361823        0.4904803
#> MMP7        0.08335487       -0.7415136       -0.5298783       -0.6217634
#> CD163      -0.40891468       -0.4089147        0.9826950       -0.4089147
#>       CTTCTCTCAATCAGAA CTTCTCTCAATGAAAC CTTCTCTCAATGGAAT CTTCTCTGTACGCACC
#> IFRD1       -1.2434526        0.7264392        1.2055459       0.16013991
#> CEBPD       -1.2991399        0.3054019        0.5502914      -0.37385321
#> SAT1        -0.1627417        0.6854468        1.2872797       0.02262656
#> CLU         -0.8718007        0.9086919        0.5783936      -0.46090742
#> MMP7        -0.4072971        0.2238867        1.1985717      -0.71755060
#> CD163        1.8994113       -0.4089147       -0.4089147      -0.40891468
#>       CTTCTCTGTGTATGGG CTTGGCTAGCGCCTCA CTTGGCTCAGTAGAGC CTTGGCTGTTGAGTTC
#> IFRD1       -0.7393701       0.04597621        0.7171202        0.5220290
#> CEBPD        1.4058715      -0.67070968        1.4217366        0.5451079
#> SAT1         1.0669751      -0.74639214        0.3205457        0.1745867
#> CLU         -0.6432558      -0.57301722        1.3978270       -0.8749584
#> MMP7         1.3100988      -0.96023541       -0.3825350        0.8664472
#> CD163       -0.4089147      -0.40891468       -0.4089147       -0.4089147
#>       CTTTGCGTCCGAAGAG GAAACTCCACCCAGTG GAAACTCGTCTGCCAG GAAACTCTCATGCAAC
#> IFRD1       -1.2434526        0.2978596        0.1137009       -0.6036588
#> CEBPD       -0.7866177        0.8215242        0.1068312       -0.4212197
#> SAT1        -0.6505482        0.6121403       -0.5727504       -0.7639854
#> CLU         -0.7903667        2.4011429        0.7821008       -0.5202265
#> MMP7        -0.9602354        1.8583813       -0.6667861        1.4410370
#> CD163        1.7139479       -0.4089147       -0.4089147       -0.4089147
#>       GAAATGAAGTACGCCC GAACATCCAATGGTCT GAACATCGTAATAGCA GAACATCGTACCGTTA
#> IFRD1       -1.2434526        0.9563607        2.7475831        1.4703918
#> CEBPD       -0.4377816        0.6451120        1.9044553        0.6767400
#> SAT1        -1.7154337        0.6954942        2.0111687        1.0407320
#> CLU         -0.6073118        1.1747181        0.6900702        1.0454996
#> MMP7        -0.3808687        1.6721013        1.8372759       -0.3380963
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       GAACCTAGTCGCATCG GAACCTATCAGAGACG GAACCTATCGTTTATC GAACGGAAGCCTATGT
#> IFRD1      -0.74995768       -1.2434526        1.4010763       0.64669853
#> CEBPD       0.27323213       -1.8354638        0.5076501      -0.12256442
#> SAT1        0.03371221       -1.3344108        2.3531002       0.04749791
#> CLU        -0.43260553       -0.8413420        1.4122967       0.75371796
#> MMP7       -0.55661302       -0.9602354        1.9163712      -0.96023541
#> CD163      -0.40891468       -0.4089147       -0.4089147      -0.40891468
#>       GAACGGAAGGAATTAC GAACGGAAGTTACCCA GAATAAGAGATCCGAG GAATAAGGTGAAGGCT
#> IFRD1       -1.2434526       -1.2434526       0.01124565       -1.2434526
#> CEBPD       -0.1937734       -0.8159907      -0.23383776       -1.8354638
#> SAT1         0.1030214        0.2914687       0.22963791       -1.6501948
#> CLU         -0.9188425       -0.9207796       0.41074485       -1.2361823
#> MMP7        -0.4694494       -0.6696607      -0.59483630       -0.9602354
#> CD163        1.6644173        2.0452465      -0.40891468       -0.4089147
#>       GAATAAGGTTATGTGC GACACGCTCTCGCATC GACAGAGAGCAGGCTA GACAGAGGTGTCTGAT
#> IFRD1       -0.2328750       -0.8519717       -1.2434526        0.8042021
#> CEBPD       -1.0704582        1.1036974        0.6785506       -0.1057403
#> SAT1         1.1708090        0.5160079        0.9735459       -0.1418590
#> CLU         -1.2361823        0.9805459       -0.6812262        0.4327436
#> MMP7        -0.9602354        1.2282449        0.9317150       -0.9602354
#> CD163        1.1394552       -0.4089147       -0.4089147       -0.4089147
#>       GACAGAGTCAACACGT GACCAATAGACTAAGT GACCTGGAGGAGTACC GACCTGGCACCTCGTT
#> IFRD1      -0.55253836       -0.5905231        0.6188647        0.3054923
#> CEBPD       0.07643427       -1.0120256        0.9652749        0.3951495
#> SAT1       -1.20541104       -0.5946419        1.5810108       -1.0248417
#> CLU        -0.47620764       -1.1354291        0.3960767       -0.2874483
#> MMP7        0.04492058       -0.6508621        1.8269401        1.4431134
#> CD163      -0.40891468        0.9519206       -0.4089147       -0.4089147
#>       GACGCGTGTGACCAAG GACGGCTGTTGTTTGG GACGGCTTCTACTATC GACGTGCAGACGCACA
#> IFRD1       -0.6333725       -0.3277121       -1.2434526        2.9211506
#> CEBPD       -1.0590263       -0.7369470       -1.1758808        1.3923058
#> SAT1        -0.5887283       -0.9443637       -1.3738736        2.0260818
#> CLU         -1.2361823       -0.3851716       -1.0519973        0.9927140
#> MMP7        -0.6711651       -0.5263360       -0.7905492        1.7225264
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       GACGTGCCACTACAGT GACGTTAGTTACGCGC GACGTTATCTGGAGCC GACTAACTCCCATTAT
#> IFRD1      0.002948534       -0.3969559      -0.21212513       0.07738604
#> CEBPD     -0.625584509        1.0566301      -0.08996857       0.73027817
#> SAT1       0.398708011       -1.2358869       0.33111152      -0.11742910
#> CLU       -0.842637019        0.6190907      -0.70576112       0.78274538
#> MMP7      -0.960235411        1.7141761      -0.96023541      -0.63276225
#> CD163     -0.408914675       -0.4089147       1.17124734      -0.40891468
#>       GACTAACTCGGAATCT GACTACAGTCGTCTTC GACTACATCATCGATG GACTACATCCTTTCTC
#> IFRD1       -0.6294712        0.6765949       -1.2434526       -1.2434526
#> CEBPD       -1.3706809        1.2360879       -1.8354638       -0.6658643
#> SAT1        -0.3804970        1.2061459       -1.5452443        0.8786464
#> CLU          0.2981796        1.7863636       -0.4987042       -0.8596575
#> MMP7        -0.7980068        1.2837190        0.3301038       -0.6133501
#> CD163       -0.4089147       -0.4089147       -0.4089147        4.3720748
#>       GACTGCGCAATGCCAT GAGCAGAAGCCCAGCT GAGCAGAAGTGGCACA GAGCAGACAATGACCT
#> IFRD1       1.18912272        0.9176806        0.5016946       -0.5814191
#> CEBPD       0.78376399        1.4530315       -0.7167588       -1.3343056
#> SAT1        0.19571602        1.9733765       -1.0072195       -1.1215926
#> CLU         1.77918574        1.3856558       -0.9537396       -1.2361823
#> MMP7       -0.06430115        1.9655349       -0.2258585       -0.9602354
#> CD163      -0.40891468       -0.4089147       -0.4089147       -0.4089147
#>       GAGTCCGGTAATCACC GAGTCCGGTCTTCAAG GAGTCCGGTGGTAACG GAGTCCGTCTGACCTC
#> IFRD1       -0.1728601        1.4560407       -0.2098305       -1.2434526
#> CEBPD       -0.5835576        1.4073109       -1.8354638       -1.8354638
#> SAT1        -0.1179094        0.6109689       -1.4090805        0.3836183
#> CLU          0.5147469        1.3439265       -1.2361823       -1.2361823
#> MMP7        -0.9602354        0.8606503       -0.1992486       -0.9602354
#> CD163       -0.4089147       -0.4089147       -0.4089147        1.1903235
#>       GATCAGTGTAGGAGTC GATCGATCATTTCACT GATCGATGTAGTACCT GATCGATTCACCGGGT
#> IFRD1        1.1479401      -0.43571321       -1.2434526        0.4494499
#> CEBPD        1.0008246      -0.03469502        0.5878024        0.3967418
#> SAT1         0.9481301       0.59367782       -0.4718508        0.1972949
#> CLU          1.4900870      -1.23618235        0.5353944        1.5114852
#> MMP7         1.6977261      -0.96023541        1.7890013        1.1021493
#> CD163       -0.4089147       3.73239471       -0.4089147       -0.4089147
#>       GATCGTAAGTACGACG GATCGTAAGTAGTGCG GATCGTACACGAAAGC GATCGTACAGATCGGA
#> IFRD1       -1.2434526        1.9229341        1.7859037       -1.2434526
#> CEBPD       -0.8639892        1.4014698        1.0985626       -1.8354638
#> SAT1        -0.5534060        0.3658498        1.6938360       -1.8587486
#> CLU         -1.2361823        0.2928046        1.4290643       -1.2361823
#> MMP7        -0.4974116        0.8122960        2.0670230       -0.9602354
#> CD163        1.0876796       -0.4089147       -0.4089147       -0.4089147
#>       GATCGTAGTTCCCGAG GATCTAGAGCTGCGAA GATCTAGCAGTCAGCC GATCTAGGTAGCAAAT
#> IFRD1        1.3288262       -0.8491790        0.4459274       -0.4957708
#> CEBPD        0.7303893       -0.1926756       -1.8354638       -0.6468203
#> SAT1         0.5589842        0.2177687       -0.6685777       -0.7363188
#> CLU          1.8018475       -0.4706647       -1.2361823       -1.2361823
#> MMP7         0.7921899       -0.6297691       -0.9602354       -0.6059662
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       GATCTAGTCAAACCGT GATCTAGTCTATCCCG GATGAAATCAACACAC GATGAGGCACAACGCC
#> IFRD1       -0.4909973       -0.2096601       -0.6559998        0.1543880
#> CEBPD       -1.2658564       -0.3182669        0.6645669        0.1898899
#> SAT1        -0.0936037        0.6229191       -0.1374062        0.9413918
#> CLU         -1.2361823       -0.4100597       -0.4377007       -0.7856720
#> MMP7        -0.9602354       -0.9602354       -0.6818865       -0.5451886
#> CD163        2.4313139        3.5093713       -0.4089147        2.7493360
#>       GATGAGGTCTCTGTCG GATGCTATCATCGGAT GATTCAGGTCCATGAT GATTCAGGTCTGATTG
#> IFRD1        1.7850305       -1.2434526      -0.15737788      -0.13235818
#> CEBPD        1.1681034        0.1412923      -0.66254653       0.06823689
#> SAT1         0.6672162        0.4437872      -0.42238051      -1.21455491
#> CLU          1.4779145        0.8865188      -0.04429539      -0.54683860
#> MMP7         0.8526091        0.8696600       0.17139550      -0.56530904
#> CD163       -0.4089147       -0.4089147      -0.40891468      -0.40891468
#>       GATTCAGGTTAAGATG GATTCAGTCCTATTCA GCAAACTAGCGTGAAC GCAAACTTCAACACGT
#> IFRD1       -1.2434526       -0.1091182      -0.03470242       -0.3109294
#> CEBPD       -1.1933885       -0.3541717      -0.65840129       -0.7200478
#> SAT1        -0.2721781       -1.1908909       0.46525478        1.7501136
#> CLU         -1.2361823       -0.6527837      -1.23618235        1.9140912
#> MMP7        -0.5583450       -0.6354418      -0.38750115        2.2058650
#> CD163       -0.4089147       -0.4089147       0.72211539       -0.4089147
#>       GCAATCACAACTGCTA GCACATAGTAACGCGA GCACATAGTTCCGGCA GCACTCTTCCATGAGT
#> IFRD1       0.54122603        0.8842705        0.0309602       1.38708045
#> CEBPD      -0.01926895        0.4877897       -0.1828851      -0.02372273
#> SAT1       -1.22063876        1.2843113        0.5875957      -0.70009129
#> CLU         0.21254656        0.8702007        1.4916278       0.81788738
#> MMP7       -0.48257113        0.2584447        0.7846485      -0.96023541
#> CD163      -0.40891468       -0.4089147       -0.4089147      -0.40891468
#>       GCACTCTTCCCTAATT GCACTCTTCTCGCTTG GCAGCCAAGGAGTTTA GCAGCCAAGGGCTTGA
#> IFRD1       -1.2434526       -0.3707088        0.3487164        0.1335534
#> CEBPD       -1.2445810       -1.2764694        0.7922431        0.8712730
#> SAT1         0.4677704       -0.5019882        1.0931609       -0.1040119
#> CLU         -1.2361823        0.5698517        1.3805065       -0.4751322
#> MMP7        -0.9602354       -0.8527908        0.5020607        1.2549505
#> CD163        3.8601836       -0.4089147       -0.4089147       -0.4089147
#>       GCAGCCACATCCGGGT GCAGCCACATTCTCAT GCAGCCATCATCTGTT GCAGTTAAGCTAACTC
#> IFRD1        0.2727385        0.4047734       -0.8239837       -0.2915615
#> CEBPD       -1.8354638        2.2133887       -0.6035331        0.8524938
#> SAT1        -0.4427172        1.4947034       -0.9167725        0.4297558
#> CLU         -1.2361823        1.2269053       -1.1197777       -0.7466160
#> MMP7        -0.9602354        0.0886407       -0.6108991       -0.5092070
#> CD163       -0.4089147       -0.4089147       -0.4089147        2.4788947
#>       GCATACAAGTACGTTC GCATACACACATGACT GCATACATCCACTCCA GCATGATAGAGGACGG
#> IFRD1       -1.2434526       -1.2434526       -0.6121618       -1.2434526
#> CEBPD       -0.3287093       -1.1053096       -0.9067684       -0.7304987
#> SAT1         0.4285445       -0.2336027       -1.0762290       -1.1452539
#> CLU         -1.2361823        0.1278626       -0.9115036       -0.4854624
#> MMP7        -0.5655860       -0.9602354       -0.6611150       -0.5233057
#> CD163        2.6407535       -0.4089147        0.5583259       -0.4089147
#>       GCATGATCACACCGAC GCATGATGTATGCTTG GCATGATGTCGCGAAA GCATGATTCCGAAGAG
#> IFRD1       -0.6396512       -1.2434526       -0.3888413        2.8417534
#> CEBPD        0.3073548       -0.5473986        0.2206365        2.4147927
#> SAT1        -1.5560257       -0.1575727       -1.0841852        1.7909977
#> CLU         -1.2361823       -0.6663810       -0.4305278        0.4819651
#> MMP7        -0.9602354        0.7250875        1.4709640        1.9898806
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       GCATGCGAGAGAGCTC GCATGCGAGATCCTGT GCATGCGAGTTCCACA GCATGCGGTGCCTGTG
#> IFRD1       -0.1451800      -0.03807476       -0.5280989      -1.24345261
#> CEBPD        0.7024103      -0.92299478       -1.1041548      -0.62583247
#> SAT1         0.8322383       0.02258984       -0.9946516       0.06804202
#> CLU          1.1721535      -1.23618235       -0.1020765      -0.70780241
#> MMP7         1.7466364      -0.38909905       -0.6212840      -0.47344881
#> CD163       -0.4089147       0.71842635       -0.4089147       1.16516618
#>       GCATGCGGTTTAGGAA GCATGTACAATCGAAA GCATGTAGTAGCTAAA GCATGTAGTGTAATGA
#> IFRD1        0.4681986       -1.2434526       -1.2434526       -0.7190568
#> CEBPD       -0.5397467       -1.4315731       -1.8354638        0.6105220
#> SAT1         1.3614749       -1.6596263       -0.6810396       -0.5276720
#> CLU         -0.8055079       -0.6157089       -0.6769468        0.2147844
#> MMP7        -0.9602354       -0.9602354       -0.4450221        1.3593941
#> CD163        2.2136144       -0.4089147       -0.4089147       -0.4089147
#>       GCCAAATAGAATGTTG GCCAAATCATCCTTGC GCCAAATTCCGTCAAA GCGAGAACAAGACACG
#> IFRD1       -0.2584428       -0.1626922        0.8830083        0.9031331
#> CEBPD       -0.2513463       -1.8354638        0.9706281        1.1262280
#> SAT1         0.2238105       -0.1551645        1.2786092        0.4906886
#> CLU         -0.5711298       -0.3789143        0.7557877        1.2917606
#> MMP7        -0.1298534       -0.1704503        1.8092394        1.6151115
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       GCGAGAACACGTTGGC GCGAGAATCATGCAAC GCGAGAATCGGTCTAA GCGAGAATCGTCACGG
#> IFRD1      -0.08779588        0.9959654      -0.49886845       -0.6220513
#> CEBPD      -0.38578822        0.2581691       1.15475097        0.2685349
#> SAT1       -0.16297218        0.4150326       0.06881915       -1.2858536
#> CLU        -0.18081073        0.7231849       1.77716912       -0.5363526
#> MMP7        1.10390342        1.9774781       1.49862708       -0.5582024
#> CD163      -0.40891468       -0.4089147      -0.40891468       -0.4089147
#>       GCGAGAATCTGCGTAA GCGAGAATCTGGTTCC GCGCAACAGGCGCTCT GCGCAACCACAACGTT
#> IFRD1      -0.55350659      -1.24345261      -1.24345261      -0.75725928
#> CEBPD      -0.51583824      -0.25818511      -0.05572497       0.14159903
#> SAT1        0.07149599      -0.06368788       0.50802955       0.01783524
#> CLU        -0.12766461      -1.23618235      -1.23618235      -0.08479953
#> MMP7       -0.96023541      -0.96023541      -0.96023541       0.41489176
#> CD163       1.33881681      -0.40891468       1.21703109      -0.40891468
#>       GCGCAACTCCAAAGTC GCGCAGTAGCAGCGTA GCGCAGTCATGCAACT GCGCAGTGTATCACCA
#> IFRD1       0.08930698       -1.2434526       -1.2434526       -0.3351229
#> CEBPD      -1.11969835       -0.2703591       -1.4172538        0.1001282
#> SAT1        0.19476919        0.9152258       -0.8534096        0.4431301
#> CLU         0.80375510       -1.2361823       -1.0793533       -1.2361823
#> MMP7       -0.51222055       -0.9602354       -0.6984677       -0.7087190
#> CD163      -0.40891468        1.7406768       -0.4089147       -0.4089147
#>       GCGCAGTGTGGTAACG GCGCAGTTCACAACGT GCGCCAACAACTTGAC GCGCCAAGTTAAGGGC
#> IFRD1       -0.1660034       -1.2434526       -1.2434526        0.5799248
#> CEBPD        1.6820473       -1.8354638       -1.8354638        1.1637583
#> SAT1         1.5172528       -1.0956305        0.0197250        0.8746932
#> CLU          1.1785062       -0.6572259       -1.2361823        1.1430715
#> MMP7         2.1118860       -0.9602354       -0.9602354        0.8605541
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       GCGCCAATCAAGAAGT GCGCGATAGGGTATCG GCGCGATAGTCCGTAT GCGCGATCACATAACC
#> IFRD1       -0.4121216       -1.2434526       -1.2434526      -0.33115979
#> CEBPD       -0.6754355       -0.6649435        0.8238651       0.67085558
#> SAT1        -0.3234950       -1.3795171        1.2536513       0.08385603
#> CLU          0.3695070       -0.8592714       -1.2361823       0.56859285
#> MMP7        -0.7327640       -0.9602354       -0.5420493       0.84486274
#> CD163       -0.4089147       -0.4089147        1.7464644      -0.40891468
#>       GCGCGATTCATTGCGA GCGCGATTCTGTCTCG GCGGGTTAGGATATAC GCGGGTTAGGGAAACA
#> IFRD1        0.4644370        0.6315353        0.6810054       -0.9440319
#> CEBPD       -1.1465114       -0.0073004        0.6465030        0.1243760
#> SAT1        -0.5390108        0.1417037        0.1440447       -0.0611739
#> CLU          0.0473660       -0.7500203       -1.2361823        0.3711270
#> MMP7        -0.7814411       -0.9602354        0.5932875       -0.2151720
#> CD163       -0.4089147        3.8767356       -0.4089147       -0.4089147
#>       GCGGGTTGTCAAAGAT GCTCCTAAGAGTAATC GCTCCTAAGCTCCTCT GCTCCTAAGGAGCGAG
#> IFRD1        1.3330227       -1.2434526       -1.2434526       -0.1214652
#> CEBPD        0.7789076        0.6600440       -1.0642554       -0.2236134
#> SAT1         0.9714244        0.5637084        0.7080184        0.1831990
#> CLU          1.0588242       -0.8936566       -0.7122187       -0.8880922
#> MMP7         0.5873407       -0.9602354       -0.9602354       -0.6395464
#> CD163       -0.4089147        2.1954323        2.0222841        2.2256229
#>       GCTCCTAGTATTCTCT GCTCCTATCAGTTAGC GCTCTGTAGATCACGG GCTCTGTCAAGAGTCG
#> IFRD1       -0.3271937       -1.2434526        1.1415365       -0.7647673
#> CEBPD       -0.6310954       -0.1089692        1.3806349       -0.9983154
#> SAT1        -0.8834451       -0.9099430        0.2319805       -0.1717520
#> CLU         -0.9604581       -1.0535679        1.4817328        0.6583359
#> MMP7        -0.6402383       -0.1290211       -0.9602354       -0.5673416
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       GCTCTGTCATTCCTGC GCTCTGTTCGATAGAA GCTGCAGAGAAGGGTA GCTGCAGCAAAGGTGC
#> IFRD1        1.6295097       -1.2434526       -0.3688251        0.2618515
#> CEBPD        1.7244200        1.9350100        0.4094246        0.8444798
#> SAT1         1.5352412        0.8711265        0.3237891        1.8084921
#> CLU          1.3182755       -1.2361823       -1.2361823        1.1780329
#> MMP7         1.8979385        1.2914817       -0.9602354        1.7147492
#> CD163       -0.4089147       -0.4089147        3.6867627       -0.4089147
#>       GCTGCAGTCAGTACGT GCTGCGAGTCTCTCTG GCTGCGATCAGCAACT GCTGGGTCACCACGTG
#> IFRD1        1.3217062       -1.2434526      0.776793817        1.9416974
#> CEBPD        1.7066610       -1.0660028      0.958600311        1.4945436
#> SAT1         0.6828326       -1.0369033      0.979760222        1.0274783
#> CLU          0.9468824       -0.4216335      0.960591383       -0.2170941
#> MMP7        -0.2242937       -0.9602354     -0.002995155       -0.2088546
#> CD163       -0.4089147       -0.4089147     -0.408914675       -0.4089147
#>       GCTGGGTCATATACGC GCTGGGTGTGCCTTGG GCTTCCACAGGTCCAC GCTTCCATCGTCCAGG
#> IFRD1        0.5659784       -0.5557502        1.2804404        0.6785434
#> CEBPD        0.2241488       -0.2072249        0.7881498        0.1241893
#> SAT1         0.3612084        0.6789136        1.8816218       -0.5989078
#> CLU          0.3142302       -1.2361823        0.8261061        0.9850657
#> MMP7        -0.9602354       -0.9602354        0.8354680       -0.9602354
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       GCTTCCATCTTTCCTC GCTTGAAGTAGGCATG GGAAAGCAGCTAACTC GGAACTTAGATCTGAA
#> IFRD1        0.8752309       -0.8598294        0.4970403        0.9762064
#> CEBPD       -0.3245207       -0.1852499       -0.5179136       -0.3824759
#> SAT1        -0.7344723       -0.9569142       -1.8569614       -1.0448325
#> CLU         -0.6630214       -0.9582600        0.2713210       -0.9400259
#> MMP7        -0.1500594       -0.5778041        1.1896079       -0.6873919
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       GGAACTTCAACTGCTA GGAACTTCAATACGCT GGAATAACAGTAAGAT GGAATAATCCGTAGTA
#> IFRD1        1.0404698       -0.3532847       -1.2434526        1.5559470
#> CEBPD        0.9321277        0.8738774       -1.8354638        0.8161547
#> SAT1         1.4691060       -1.1232259       -0.2930229        1.7353899
#> CLU          1.6922538       -0.9266313       -1.2361823        1.0023305
#> MMP7         1.7504882       -0.6030295       -0.9602354        2.1017429
#> CD163       -0.4089147       -0.4089147        0.8114202       -0.4089147
#>       GGACAAGGTAGCGTGA GGACAAGGTGGCTCCA GGACAGAAGTGAAGAG GGACAGATCAAAGTAG
#> IFRD1       -0.3763923       -1.2434526        0.3776326        2.2422684
#> CEBPD        0.9257745       -1.8354638       -0.8991405        2.1442607
#> SAT1        -0.3623356       -1.7212012       -0.1509249        1.5262896
#> CLU         -0.7902453       -0.4323705       -0.2227583        1.5842739
#> MMP7         1.4945873       -0.4861997       -0.6320412        1.9334776
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       GGACAGATCACCTCGT GGACAGATCAGTGTTG GGACAGATCCAAACAC GGACATTAGAACTCGG
#> IFRD1        0.7769811        1.1842031      -0.75754897        0.2673215
#> CEBPD        1.8925482        0.6853074      -0.36944718        0.9112743
#> SAT1         2.0190417       -0.7364995       0.02703163        0.7403346
#> CLU          2.2252123       -0.5459701      -0.66051081        0.6142937
#> MMP7         1.3164703       -0.2317815       0.41449567       -0.5050226
#> CD163       -0.4089147       -0.4089147      -0.40891468       -0.4089147
#>       GGACATTAGTCATGCT GGACATTCAGGGCATA GGACATTCATTACCTT GGACGTCAGCACAGGT
#> IFRD1        2.3187449        1.5280141       -0.6527917       -0.5092649
#> CEBPD        1.5401076        0.1707801       -1.3883345       -0.9246883
#> SAT1         1.3456454        0.9813401       -1.3333011       -1.3777665
#> CLU          0.6738109        1.4761142       -0.7232942        0.3709619
#> MMP7         1.8675297        1.2172008       -0.9602354        0.9462632
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       GGACGTCCACAGCGTC GGACGTCCACCGTTGG GGACGTCGTCAACATC GGACGTCTCCGTTGCT
#> IFRD1       -0.1766297       -1.2434526     -0.006055867       -0.7338024
#> CEBPD       -0.9228310       -0.7128886      1.090854586       -0.2263543
#> SAT1        -1.5357841       -1.2961613      0.128319817       -0.4649896
#> CLU         -0.2177530       -1.2361823      1.149809954       -0.1117237
#> MMP7        -0.7076479       -0.9602354      0.419229180       -0.3509056
#> CD163       -0.4089147       -0.4089147     -0.408914675       -0.4089147
#>       GGACGTCTCCTAGGGC GGAGCAAGTCCTCTTG GGAGCAAGTCTAGTCA GGAGCAATCCGAGCCA
#> IFRD1       -0.5167031        0.2694928        1.2828376       1.31736715
#> CEBPD        0.4945932       -0.3969660        1.6860853      -0.38314958
#> SAT1        -1.4908873        0.7606290        1.7136262      -0.93832588
#> CLU         -0.9656838       -0.9445499        0.7799975       0.08087086
#> MMP7        -0.9602354       -0.9602354        2.1338104      -0.96023541
#> CD163       -0.4089147        1.5310583        0.5137591      -0.40891468
#>       GGATGTTGTCATCCCT GGATTACAGCAATATG GGCAATTAGCTGTCTA GGCAATTAGTCCAGGA
#> IFRD1      -0.44545999       0.07497798        0.5140489       -0.6570878
#> CEBPD       0.57579900       0.58634304        0.6800042       -0.6617856
#> SAT1        0.06698579      -2.16620933       -1.4126794        0.4811870
#> CLU        -0.82576738      -0.88152644        0.5755863       -0.9346095
#> MMP7        0.03929566      -0.96023541        0.3667450       -0.3442114
#> CD163      -0.40891468      -0.40891468       -0.4089147        1.1093336
#>       GGCAATTAGTCGTACT GGCCGATCAGCGTAAG GGCCGATGTGCTTCTC GGCGACTAGCTGTTCA
#> IFRD1        2.3564903       -0.4972038       -1.2434526       -0.6789525
#> CEBPD        1.1109701       -0.4403592        0.5890972       -0.6939212
#> SAT1         0.6665094        0.5571090        0.4669209        0.3762944
#> CLU          0.7116867       -0.6087161        1.0463689       -1.2361823
#> MMP7        -0.3975333       -0.9602354        0.8301791       -0.6927620
#> CD163       -0.4089147        2.4147743       -0.4089147        2.9517673
#>       GGCGACTAGTACGCCC GGCGACTAGTCCGTAT GGCGACTGTGTAAGTA GGCGACTTCTGGGCCA
#> IFRD1        1.3343929       -0.4691152       -0.3020503       -1.2434526
#> CEBPD        0.3016947        0.5208688       -1.8354638        0.2547441
#> SAT1         0.7355190       -0.6245099       -0.4724052       -0.2382407
#> CLU          1.3445422       -0.8379335       -1.2361823       -0.5932504
#> MMP7         0.8669450        1.4615094       -0.2564913       -0.5964110
#> CD163       -0.4089147       -0.4089147       -0.4089147        0.7675536
#>       GGCGTGTCACCGTTGG GGCTCGAAGTAGCGGT GGCTCGACAAGTAGTA GGCTCGACACTGTGTA
#> IFRD1        1.2448095     -0.330661876       0.05707672      -0.55953199
#> CEBPD        0.7906023     -0.007229316       0.25872290       0.05003512
#> SAT1         0.2974131     -0.804091038       0.37995723       0.27285905
#> CLU          1.0610073     -0.961647171      -0.74914648      -0.70463101
#> MMP7        -0.3445181     -0.960235411      -0.39502584      -0.47052706
#> CD163       -0.4089147     -0.408914675      -0.40891468      -0.40891468
#>       GGCTCGATCATTATCC GGCTCGATCGGTTAAC GGCTGGTCACGGTAAG GGGAATGTCGCGGATC
#> IFRD1        0.3786577      -0.11591634       -0.7270602      -1.24345261
#> CEBPD        0.3301033      -0.98192074        0.4260716       0.05074065
#> SAT1        -0.9723841       0.02687716       -1.3876654       0.94940691
#> CLU          0.8246657       1.56214599       -0.7270019      -0.52083307
#> MMP7        -0.9602354       0.71184014        0.3292978      -0.54763130
#> CD163       -0.4089147      -0.40891468       -0.4089147       1.72216094
#>       GGGACCTAGCAGCCTC GGGACCTCATCCCATC GGGACCTTCGACAGCC GGGAGATAGATATGGT
#> IFRD1        0.2874960        1.3244341       -1.2434526       1.18882205
#> CEBPD        0.7020519       -1.0478187       -1.8354638       0.07941684
#> SAT1         0.5390175        0.1786209       -2.1118014       0.57486191
#> CLU          1.3662208        2.0642389       -1.2361823       1.69357833
#> MMP7         1.6309235        0.4408176       -0.9602354       1.71840562
#> CD163       -0.4089147       -0.4089147       -0.4089147      -0.40891468
#>       GGGAGATAGCGTTTAC GGGAGATCATACGCCG GGGATGAAGAGGACGG GGGATGACAAGGTTCT
#> IFRD1       1.67478025       -0.6215442       -1.2434526        0.6888854
#> CEBPD       0.89095256       -0.3083582        0.1562275        0.8545193
#> SAT1       -0.04959208       -0.0456284       -0.7290769        1.5145225
#> CLU         0.85866025       -0.9163291       -0.6392917        0.8128409
#> MMP7       -0.47835127       -0.6655606        0.8681960        1.7747058
#> CD163      -0.40891468       -0.4089147       -0.4089147       -0.4089147
#>       GGGATGACACCAGGTC GGGATGAGTTCTGTTT GGGCACTAGCCTATGT GGGCATCCATGCTAGT
#> IFRD1        0.2629037        0.5501209        0.2030833       -0.2181444
#> CEBPD       -1.8354638       -0.1621149       -1.2779982        0.5037810
#> SAT1        -0.7223670       -0.7075056       -1.8035396        0.1467533
#> CLU          0.1216961       -0.5069069       -0.8574368       -0.7088569
#> MMP7        -0.5066098       -0.9602354       -0.9602354       -0.3062083
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       GGGCATCCATTCCTCG GGGCATCGTCGGCATC GGGCATCGTGACGGTA GGGTCTGAGGGCTCTC
#> IFRD1       -1.2434526       2.46686553        1.0903686       -1.2434526
#> CEBPD       -1.8354638       1.65908557        1.6559310       -1.8354638
#> SAT1        -1.3994687       0.04726138        0.5121003       -2.3674020
#> CLU         -1.2361823       0.87698175        1.1066590       -0.9020407
#> MMP7        -0.9602354      -0.07380911        0.2785593       -0.9602354
#> CD163       -0.4089147      -0.40891468       -0.4089147       -0.4089147
#>       GGGTCTGGTAGAGTGC GGGTCTGTCATCATTC GGGTTGCCAGCTCGCA GGGTTGCCAGTCACTA
#> IFRD1        0.7113667        0.2762170        1.0342677       -1.2434526
#> CEBPD       -0.4359994        0.4294139        0.5477151        0.1415070
#> SAT1        -0.2210775       -0.1534308        1.7396932        0.8422653
#> CLU          0.1274929        0.6176666        1.0714949       -0.6866613
#> MMP7        -0.7772900        0.7396219        1.0400827       -0.1779338
#> CD163       -0.4089147       -0.4089147       -0.4089147        1.2281471
#>       GGGTTGCGTACGCTGC GGGTTGCGTGGCAAAC GGTATTGCATCCTAGA GGTATTGTCTGTCTAT
#> IFRD1        0.9072030        0.6694233        0.5995153        1.4863974
#> CEBPD        0.3819749        0.7746020       -0.5050381        1.0373192
#> SAT1        -0.1592458        1.1601584       -0.2864517        0.7281455
#> CLU          0.6555822        1.4554957        0.5513213        1.1594294
#> MMP7        -0.9602354        1.4999326       -0.9602354        1.6253422
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       GGTGAAGGTCCAGTAT GGTGCGTGTGGTCCGT GGTGTTAGTAACGCGA GTAACGTAGTGCTGCC
#> IFRD1        1.7336347       -1.2434526        0.9334855      -1.24345261
#> CEBPD        1.1141599       -0.3519588       -0.5774265       0.73455317
#> SAT1        -0.1803190       -1.4404013        0.8729086      -0.15435029
#> CLU         -0.2614910       -0.2282802        0.2648022       0.03426455
#> MMP7         0.4503785       -0.9602354       -0.3433442       0.13535316
#> CD163       -0.4089147       -0.4089147       -0.4089147      -0.40891468
#>       GTAACTGCACCAGTTA GTAACTGGTACCGTAT GTAACTGGTTAAAGTG GTAACTGTCGCCTGAG
#> IFRD1       -0.2089779        1.2552215        1.7178523     -1.243452613
#> CEBPD       -0.3175244       -0.0561059        1.5238326     -0.009426542
#> SAT1         1.1243224        0.9889149        1.4583298      0.508406015
#> CLU         -0.5410196        0.9851885        1.5688250     -0.895220776
#> MMP7         1.6636717       -0.9602354        1.3938008     -0.960235411
#> CD163       -0.4089147       -0.4089147       -0.4089147      4.400285494
#>       GTACGTACAGTAAGCG GTACTCCAGGAATGGA GTACTCCAGTCGTACT GTACTCCCACAGACAG
#> IFRD1       0.09261077        0.8329146        0.1834626     -1.243452613
#> CEBPD       0.12723972       -0.2636573       -1.1559801      0.740269788
#> SAT1        0.31840579       -0.7625684        0.1171609      0.001000293
#> CLU        -0.54903263        0.8009353       -1.2361823      0.878236381
#> MMP7       -0.56679233       -0.7159857       -0.5349301      1.240675458
#> CD163       2.63425927       -0.4089147        0.9663592     -0.408914675
#>       GTACTTTCACGGTGTC GTACTTTGTGAAAGAG GTACTTTGTGTCGCTG GTAGGCCAGCGTTCCG
#> IFRD1      -1.24345261       -0.2089210       -1.2434526       -0.2414909
#> CEBPD       0.09668525        0.4371107        0.4885028       -1.0769805
#> SAT1        0.41119598       -1.4929430        0.7156299       -1.7198024
#> CLU        -0.86034199        0.3349471       -1.2361823       -1.2361823
#> MMP7       -0.39243210        0.4622238       -0.9602354       -0.9602354
#> CD163       1.42714309       -0.4089147       -0.4089147       -0.4089147
#>       GTAGGCCCAAGAGGCT GTAGGCCGTAATCACC GTAGGCCGTCGCCATG GTAGTCACACTTAAGC
#> IFRD1       -1.0456223        1.0579990       -1.2434526      -1.24345261
#> CEBPD       -1.1579354        1.0145267       -1.0309697      -0.65868630
#> SAT1        -0.9294670        1.6440998        0.7638003       0.02601663
#> CLU         -0.6264616        1.9330759       -1.2361823      -1.23618235
#> MMP7        -0.1902835        1.5340845       -0.9602354      -0.96023541
#> CD163       -0.4089147       -0.4089147       -0.4089147      -0.40891468
#>       GTAGTCACATACAGCT GTAGTCACATATGGTC GTAGTCAGTTCCACAA GTATCTTGTACTTGAC
#> IFRD1       -0.8636046      -1.24345261       -0.1025789       -1.2434526
#> CEBPD       -0.9886815       0.02016971        0.5245529        0.4530681
#> SAT1        -1.0551137      -2.02847198       -0.2315527        1.1507101
#> CLU         -1.0408230      -1.23618235       -0.5080931       -0.3785437
#> MMP7        -0.9602354      -0.42544821        0.4575873       -0.9602354
#> CD163       -0.4089147      -0.40891468       -0.4089147        3.9452066
#>       GTATCTTTCTTGTACT GTATTCTAGTGTACTC GTATTCTCAGCCAGAA GTCAAGTCACTTCGAA
#> IFRD1      -1.24345261       -1.2434526        1.2359640       -1.2434526
#> CEBPD      -1.37789711       -0.2536667        0.7178926       -0.1222845
#> SAT1       -1.03579331        2.2000374        0.4514707        0.9563760
#> CLU        -0.31164394        0.4842413        1.6069858       -0.7965299
#> MMP7       -0.02174962        1.6393490        1.2201579       -0.9602354
#> CD163      -0.40891468       -0.4089147       -0.4089147        4.0784262
#>       GTCAAGTCAGACAAGC GTCACGGAGCTGAACG GTCACGGAGTTGCAGG GTCACGGCACTCTGTC
#> IFRD1        0.1079494       -1.2434526        0.8830201       -1.2434526
#> CEBPD       -0.1074319        0.1459837       -0.4487599       -0.9230394
#> SAT1        -0.1399147        0.1163218        0.3966556        0.5211504
#> CLU         -1.2361823       -0.6582556        0.4795608       -1.2361823
#> MMP7         0.4032448       -0.7288828       -0.5260538       -0.6116236
#> CD163       -0.4089147       -0.4089147       -0.4089147        1.4378305
#>       GTCACGGGTCAGCTAT GTCATTTCACCGAAAG GTCCTCAAGACATAAC GTCGGGTAGCAGGCTA
#> IFRD1        1.3201676        0.1564308       0.03699997        2.2965158
#> CEBPD        1.7350745        0.7865499      -0.38360589        0.8442864
#> SAT1         1.0616259        1.4087202      -0.82969847        1.9690291
#> CLU          1.5170991        0.8280586       0.75472774        0.1928689
#> MMP7         1.9017031        0.5778146      -0.82870884        2.0637576
#> CD163       -0.4089147       -0.4089147      -0.40891468       -0.4089147
#>       GTCGGGTGTGGGTCAA GTCGGGTTCACAACGT GTCGGGTTCGACAGCC GTCGTAATCCGTAGTA
#> IFRD1        1.5954529       -0.6525228       -0.7007891     -0.796530963
#> CEBPD        0.4926861        0.2137313       -0.9097996     -0.607549906
#> SAT1         0.6893555       -0.1608203        0.3326856     -0.145915530
#> CLU          1.3327966        1.0548405       -0.9570855     -0.008273431
#> MMP7         1.4637266        1.3546857       -0.7031087     -0.748473426
#> CD163       -0.4089147       -0.4089147        1.8355686     -0.408914675
#>       GTCGTAATCTAACTTC GTCGTAATCTTTACGT GTCTCGTGTTCAGACT GTCTTCGAGCTCCTTC
#> IFRD1       -1.2434526       -0.9644000        1.5970525       -0.1279413
#> CEBPD       -0.1494229        0.1897061        1.2657764        1.0424225
#> SAT1        -0.0423420       -0.5973716        1.6571789        0.4955663
#> CLU         -1.2361823       -0.7080618        1.0904260       -0.3560719
#> MMP7        -0.9602354       -0.5452529        1.3824213        1.8731922
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       GTGAAGGCAATGAATG GTGAAGGGTCGCGAAA GTGAAGGTCGTTTGCC GTGCAGCAGAGTAATC
#> IFRD1       -0.6760937        0.6149418        0.8533952        0.3925106
#> CEBPD       -0.5323394        0.3476030        0.7031675        0.8359110
#> SAT1        -0.4258540        0.9416287        1.0279211        0.5549871
#> CLU         -0.5846892        1.1560936        1.2745233        1.0253813
#> MMP7        -0.8114837        1.8207526        1.5898114       -0.9602354
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       GTGCAGCAGGCATTGG GTGCAGCCACCAGCAC GTGCAGCGTTGTTTGG GTGCAGCTCATGCATG
#> IFRD1      -0.08447858       -1.2434526        0.0781744      -1.24345261
#> CEBPD       1.19118800       -0.9791609       -1.6157075      -0.55660596
#> SAT1        1.11369959        0.1161665       -0.4638213      -0.05089445
#> CLU         2.15214599       -0.6544050        0.4815156      -0.52412223
#> MMP7        1.02627221       -0.9602354       -0.9602354       2.25715643
#> CD163      -0.40891468        1.3242409       -0.4089147      -0.40891468
#>       GTGCATACACCACCAG GTGCATACACGCTTTC GTGCTTCCACATCCGG GTGCTTCTCAACACTG
#> IFRD1      -0.03119043       -0.3232583       -0.5341068        0.3276289
#> CEBPD      -0.12930361       -1.1388783       -0.8142791       -1.5126430
#> SAT1       -1.23227027        0.7555856       -1.1567445       -0.4172381
#> CLU         0.50708041       -1.2361823       -1.0295469        0.7995068
#> MMP7       -0.60927389       -0.5242257       -0.9602354       -0.6057083
#> CD163      -0.40891468        1.0009731        1.6579612       -0.4089147
#>       GTGGGTCGTATTCTCT GTGTGCGCAAAGTGCG GTGTGCGTCGGAATCT GTGTTAGCAATCTACG
#> IFRD1        1.0783325       -1.2434526       -0.6574848        2.0084619
#> CEBPD        0.8978065       -1.0319256       -0.3037694        1.3760263
#> SAT1         0.6662568       -1.3808630       -1.0276526        2.0656182
#> CLU          0.6208179       -1.2361823       -0.3829466        1.2334204
#> MMP7        -0.6925527       -0.9602354       -0.6825901        1.8409555
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       GTGTTAGCACCCTATC GTGTTAGGTTAGAACA GTTAAGCCAATGTAAG GTTACAGGTAATCGTC
#> IFRD1        2.0028241        0.2134131        1.7536573        0.1895517
#> CEBPD        1.6177986       -1.5826137        1.6879373       -0.7506820
#> SAT1         0.9380942       -0.4670641        1.3527667        0.3622120
#> CLU          0.3471951        0.3599994        1.3256905        0.4926697
#> MMP7         1.2948069       -0.5708341        1.9514539        2.1946746
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       GTTACAGGTACGACCC GTTACAGGTCAAGCGA GTTACAGTCTTGTTTG GTTCATTAGAGACTAT
#> IFRD1       -0.4617260       -0.3951257     -1.243452613        1.1730293
#> CEBPD       -0.3029344        0.8378347      0.003491901        1.3363635
#> SAT1        -0.3146609        0.2401779     -0.852081293        0.9101443
#> CLU          0.7396053        1.3226255      0.387630649        1.8888270
#> MMP7        -0.9602354        1.2514480     -0.275750840        2.2885238
#> CD163       -0.4089147       -0.4089147     -0.408914675       -0.4089147
#>       GTTCATTGTTATCGGT GTTCGGGAGATGTCGG GTTCGGGCACTTGGAT GTTCTCGAGAAACCGC
#> IFRD1       -0.6705903       -0.1508603        1.7915626       -0.1453736
#> CEBPD       -0.6815716       -1.8354638        0.7980698       -0.2490140
#> SAT1         0.6079701        0.3685191        1.1373110        0.9663626
#> CLU         -0.7367906        0.5493926        1.0111873       -0.8966837
#> MMP7        -0.5001551       -0.1632678       -0.1782454       -0.6474617
#> CD163        1.5467952       -0.4089147       -0.4089147        1.7764848
#>       GTTCTCGAGAGGTTGC GTTCTCGCAATCACAC GTTCTCGGTAATAGCA GTTCTCGGTCCCTTGT
#> IFRD1       -1.2434526        0.3454034      -0.57811248       -0.6368265
#> CEBPD       -1.8354638        0.0258587      -0.00558617        0.3124890
#> SAT1        -0.9493003        2.3635331       0.08584964       -0.4558990
#> CLU         -0.6726220        1.5078446      -0.36258017       -1.2361823
#> MMP7        -0.9602354        0.7848841       1.12061857       -0.9602354
#> CD163       -0.4089147       -0.4089147      -0.40891468        2.0259952
#>       GTTTCTATCCGTAGTA TAAACCGAGGAACTGC TAAACCGCACGGATAG TAAACCGCAGATTGCT
#> IFRD1       -1.2434526      -0.37297320        0.4739397       -0.9751230
#> CEBPD       -0.8717206       0.46765471       -1.1087652       -0.4679263
#> SAT1         0.3045613      -0.96032744       -0.6201106       -2.0281331
#> CLU         -1.2361823      -0.45110734        0.2959192       -0.7998655
#> MMP7        -0.5885312       0.05203737       -0.8120398        0.8116234
#> CD163        2.5157566      -0.40891468       -0.1551466       -0.4089147
#>       TAAACCGGTATATGAG TAAACCGTCAAACCAC TAAGAGATCTGAGTGT TAAGCGTAGCCCAGCT
#> IFRD1     -1.243452613      -0.35956177        0.5902507       -0.1057907
#> CEBPD      0.002868304       0.06808837        0.2925389       -1.3147909
#> SAT1       1.458857037      -2.11468498        0.6939885       -0.1588502
#> CLU       -1.236182346      -1.01682732        0.6079837        0.4849926
#> MMP7      -0.960235411      -0.85106238        0.4370691       -0.6343337
#> CD163      2.939404004      -0.40891468       -0.4089147       -0.4089147
#>       TAAGCGTCAGGCAGTA TAAGCGTGTCTAACGT TAAGTGCAGAATAGGG TAAGTGCCAATCACAC
#> IFRD1       -0.6470088       0.07660142        1.1929765        1.9008466
#> CEBPD        0.9947525       0.45288066        0.9953098        0.2876574
#> SAT1         0.3171765      -0.97042219        1.1759556        1.2316678
#> CLU         -0.1438772      -0.06884030        1.5973175       -0.6920774
#> MMP7        -0.6776263       1.20569289        0.2501336       -0.4589617
#> CD163       -0.4089147      -0.40891468       -0.4089147       -0.4089147
#>       TACACGAAGTGAATTG TACACGAAGTTAACGA TACACGACACAAGCCC TACACGAGTCAAAGAT
#> IFRD1       -0.1512295       -0.2336462       -1.2434526       -0.5440771
#> CEBPD       -0.1824084        0.4261088       -0.1093432        1.1851456
#> SAT1        -1.0005161       -0.2414584       -2.8271342       -1.1709294
#> CLU         -0.8618252        1.1288055        0.9869494       -0.8764871
#> MMP7        -0.6854016       -0.4817654       -0.9602354       -0.5118965
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       TACACGAGTTGAACTC TACAGTGCAGACGCTC TACAGTGGTAACGTTC TACAGTGGTTCTCATT
#> IFRD1        1.0785077        0.9665617       -0.1610051        2.0930400
#> CEBPD        1.1690164        1.3066298        0.4544368        1.5756977
#> SAT1         1.0134267        1.3567952        0.5661916        1.4801808
#> CLU          1.0689932        0.8143553       -0.3778014        1.3619744
#> MMP7         1.8188722        1.7279886        0.1688312       -0.3322829
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       TACCTATAGATGAGAG TACCTATAGATGGCGT TACCTATAGGACACCA TACCTATCAGATGGGT
#> IFRD1       -1.2434526       -1.0424891       -1.2434526       -1.2434526
#> CEBPD       -0.7707392        0.4361513        1.1341839       -1.8354638
#> SAT1        -1.5671275       -1.3340457        1.2067554       -0.4845583
#> CLU         -0.7823589       -0.7162000       -0.7137863       -1.2361823
#> MMP7        -0.9602354        1.0528232       -0.2102627       -0.9602354
#> CD163       -0.4089147       -0.4089147        3.7839908        1.0189166
#>       TACCTTAGTCAGAAGC TACGGATCATTGGGCC TACGGGCCACCACCAG TACGGGCTCAATCACG
#> IFRD1       -0.1430459       -0.5616727        2.2662272       -1.2434526
#> CEBPD       -0.2465328       -1.3193576        0.8213580        0.2815173
#> SAT1        -0.9442583       -0.8353471        1.3526818        0.9457318
#> CLU         -1.2361823       -1.2361823        1.0398400       -0.7950206
#> MMP7        -0.1585341       -0.9602354        1.3597837       -0.9602354
#> CD163       -0.4089147       -0.4089147       -0.4089147        3.0662841
#>       TACGGTAAGCGTAGTG TACGGTATCAGTTAGC TACTCGCAGGGTGTTG TACTCGCGTTCCAACA
#> IFRD1     -0.543441030       -0.9052225       -1.2434526       -1.2434526
#> CEBPD      0.524895473       -0.3986451       -1.8354638       -0.8387477
#> SAT1       0.005374336       -0.5570162       -0.1438303       -1.8954359
#> CLU       -0.809676311       -0.8640281       -0.8657067       -1.2361823
#> MMP7       0.379489122       -0.4762266       -0.9602354       -0.9602354
#> CD163     -0.408914675       -0.4089147       -0.4089147       -0.4089147
#>       TACTTACAGACACGAC TACTTACCAAGAGTCG TACTTACGTTAAGAAC TACTTACTCTGACCTC
#> IFRD1       0.63415917        0.3994097       0.64997687        1.1384198
#> CEBPD      -0.33323768       -0.5918198       1.12636940        0.9831702
#> SAT1        0.02872469        0.4184506       0.08405695        0.9274720
#> CLU         1.24036973        0.7955608       1.62932329        1.4995690
#> MMP7        0.63850846       -0.5843417       1.67970160        1.9386657
#> CD163      -0.40891468       -0.4089147      -0.40891468       -0.4089147
#>       TACTTGTAGAGGTAGA TACTTGTCATGCCCGA TACTTGTGTAGCGTCC TACTTGTTCACTGGGC
#> IFRD1      -0.02088714      -0.08937528        1.5561282       -0.7897927
#> CEBPD      -0.06548044      -0.63385433        1.4285329        0.9197177
#> SAT1       -0.41746855      -0.49268158        1.6740131        0.8070041
#> CLU         0.46644873      -0.78927121        0.5698092       -0.4790459
#> MMP7       -0.96023541      -0.96023541        1.8004344        1.7277000
#> CD163      -0.40891468       1.71835228       -0.4089147       -0.4089147
#>       TAGACCAAGACAGACC TAGACCATCTCTGTCG TAGAGCTTCCATGAGT TAGCCGGCATAACCTG
#> IFRD1       -1.2434526        0.8107387       -1.2434526       -0.1750796
#> CEBPD       -1.1584204       -0.5853623       -0.8321136       -0.2808378
#> SAT1        -0.1050919        0.3777588       -0.5561670       -0.2402251
#> CLU         -0.7761950       -1.2361823       -1.2361823       -0.3870976
#> MMP7        -0.2861758       -0.9602354       -0.9602354       -0.4540151
#> CD163       -0.4089147       -0.4089147        2.0093537       -0.4089147
#>       TAGTGGTCAATGTAAG TAGTTGGCACCGAAAG TATCTCAGTCCTCTTG TATGCCCCACCCATGG
#> IFRD1       -0.5105619       -0.2024495       -0.4043614       -0.7741400
#> CEBPD       -1.8354638       -1.8354638       -1.2002731       -1.0115895
#> SAT1        -1.2428334       -2.8271342       -0.3410967       -0.9221471
#> CLU         -1.2361823       -0.4052572       -1.2361823       -0.5599658
#> MMP7        -0.6129745       -0.9602354       -0.9602354       -0.8035750
#> CD163       -0.4089147       -0.4089147        2.2176578       -0.4089147
#>       TATTACCAGTACGATA TATTACCGTTTAAGCC TATTACCTCGATAGAA TCAACGAGTTTGGCGC
#> IFRD1       0.46096161       -0.5039872        1.1809031        0.3853725
#> CEBPD      -0.06233862        0.6058973        0.9451518       -0.6024459
#> SAT1        0.54936858        0.8707691        0.2862007       -1.6745485
#> CLU         0.85335147       -0.7239774        1.4848698       -0.6955278
#> MMP7        2.33717851        0.7443640       -0.6433988       -0.9602354
#> CD163      -0.40891468       -0.4089147       -0.4089147       -0.4089147
#>       TCAACGATCCGCATAA TCAATCTAGTGTGGCA TCAATCTCACGACTCG TCACAAGCAGCTTAAC
#> IFRD1       0.42062573     -0.008695274        1.0636875       -1.2434526
#> CEBPD       1.09015369      0.552014607        1.1561239       -0.6693280
#> SAT1        0.15567256     -0.209983759        1.0861517       -0.2524463
#> CLU        -0.08839762     -0.790143997        1.1747231       -0.8611090
#> MMP7        1.61725209      1.024927056        1.5324662       -0.2303229
#> CD163      -0.40891468     -0.408914675       -0.4089147        2.3692462
#>       TCACAAGTCGAATCCA TCACAAGTCGTGGACC TCACGAAAGTTAACGA TCACGAACACAGGCCT
#> IFRD1       -0.5444374       -1.2434526       -0.5951954        1.5514369
#> CEBPD       -1.5362518       -1.8354638        0.5733987        1.6405909
#> SAT1        -0.8872107       -1.1498308       -0.8165310        0.8259298
#> CLU         -1.0328959       -0.6837798        0.1170013        1.0779537
#> MMP7        -0.7729514       -0.9602354       -0.6530759        1.4858698
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       TCACGAAGTTTGTTGG TCAGATGAGGTAAACT TCAGATGGTTACGCGC TCAGATGGTTGATTCG
#> IFRD1        0.1930486       -0.1770848        1.4378508       -1.2434526
#> CEBPD       -0.7480348       -0.7838556        1.7594893       -0.4564662
#> SAT1        -1.2467855       -1.5075214        1.7600662        1.2096367
#> CLU         -0.0875568       -0.9080066        1.4080118        0.2351102
#> MMP7        -0.9602354       -0.4549653        2.0502956        0.5115075
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       TCAGCAAAGGCTATCT TCAGCAACAGACGCCT TCAGCAAGTACAGCAG TCAGCAATCATATCGG
#> IFRD1        1.1254406      -0.19043740       1.35966530      -0.03160824
#> CEBPD        0.5117962      -1.56934365       1.22589244       0.49728549
#> SAT1         0.9980801      -0.91918176      -0.05012698       0.37344780
#> CLU          1.8798442       0.03362965      -0.02562109       0.49090046
#> MMP7         1.1306593      -0.66226265       0.83692115      -0.14062993
#> CD163       -0.4089147      -0.40891468      -0.40891468      -0.40891468
#>       TCAGCTCCAAGAGGCT TCAGCTCTCGTAGGTT TCAGGATGTTCCATGA TCAGGATTCCCAACGG
#> IFRD1        2.0871999        1.6879825      1.568516367     -0.759001208
#> CEBPD        1.8676188        0.1053331      0.981316132      0.357682264
#> SAT1         2.0808020        0.8591868      0.249717738      0.004087429
#> CLU          0.3148638        1.5679196      0.597982029     -0.890027540
#> MMP7         1.2048466        1.4990556      0.008893283      0.993018300
#> CD163       -0.4089147       -0.4089147     -0.408914675     -0.408914675
#>       TCAGGTACACGAAGCA TCAGGTACACGGTAGA TCAGGTATCATACGGT TCATTACAGTGCCATT
#> IFRD1      -0.24898122        1.4605360        0.0459853       0.25424783
#> CEBPD       1.01913185       -0.0228533        0.2773583      -0.06720359
#> SAT1       -0.02423444       -0.9811993       -0.3245444       0.22129771
#> CLU         0.04693811        0.5539757        0.4466128       0.34371543
#> MMP7        1.30113954       -0.9602354       -0.7812814      -0.86614830
#> CD163      -0.40891468       -0.4089147       -0.4089147      -0.40891468
#>       TCATTACCACTAAGTC TCATTACGTACCAGTT TCATTACGTACCTACA TCATTACGTACGCACC
#> IFRD1        0.3779656      -1.24345261      -0.24101185        0.9638525
#> CEBPD       -0.6080530      -0.01602795       0.06278637       -0.4799466
#> SAT1         1.0539681      -1.77909369       0.03668446       -0.2438613
#> CLU         -1.2361823      -1.23618235      -0.72061782        0.9623166
#> MMP7        -0.5907701      -0.96023541      -0.03205911       -0.9602354
#> CD163        1.5316946      -0.40891468      -0.40891468       -0.4089147
#>       TCATTACTCGAATCCA TCATTTGCAATCCGAT TCCACACGTACTTCTT TCCCGATAGGCTAGCA
#> IFRD1       -0.1357781       -1.2434526       -0.4652616      -1.24345261
#> CEBPD       -0.5476128       -0.2669289       -0.8782995      -0.66724230
#> SAT1        -0.6431793        0.8582770       -0.8587870      -0.08944962
#> CLU         -0.3612086       -1.2361823       -1.2361823      -0.72937471
#> MMP7         1.3352168       -0.9602354       -0.9602354      -0.96023541
#> CD163       -0.1012265        2.7657970       -0.4089147       4.22237784
#>       TCCCGATGTCGCGGTT TCGAGGCCACATTCGA TCGAGGCCATTGGCGC TCGCGAGAGAAGGTTT
#> IFRD1        0.1829530      -1.24345261        1.8502187       -1.2434526
#> CEBPD        0.4069924      -0.54093555        0.9776492       -0.5227051
#> SAT1         0.2931678      -0.04960188        1.1333521        0.1630319
#> CLU         -0.2548639      -0.66293585        1.3279999       -1.2361823
#> MMP7         1.0844929      -0.43211408        2.2371653       -0.5565456
#> CD163        1.7765715       1.29882681       -0.4089147        0.8964631
#>       TCGCGAGAGACATAAC TCGCGAGAGATCGATA TCGCGAGCACCTATCC TCGCGAGCATGGATGG
#> IFRD1      -1.24345261       0.13357308       -0.4494438        1.5125459
#> CEBPD       0.05029181       0.57032083       -1.8354638        1.9910552
#> SAT1       -1.20678097       0.66127520       -1.4649200        1.1627882
#> CLU        -0.90762594       0.08227641       -1.2361823        1.3891284
#> MMP7       -0.45445828       0.46275700       -0.9602354        1.7341346
#> CD163      -0.40891468      -0.40891468       -0.4089147       -0.4089147
#>       TCGCGAGTCCGCTGTT TCGCGTTCACCAGGCT TCGCGTTGTGACCAAG TCGCGTTTCGAATGCT
#> IFRD1      -0.01575687       -1.2434526        1.5274517        0.4270515
#> CEBPD      -0.11341917       -0.5062527        1.2566298        1.2526297
#> SAT1       -0.41191950       -0.9756933        0.8896972        0.6428967
#> CLU         0.51923869       -0.8777122        1.5324751        0.2755107
#> MMP7       -0.21313659       -0.9602354        1.4387443       -0.9602354
#> CD163      -0.40891468       -0.4089147       -0.4089147       -0.4089147
#>       TCGGGACTCACTCCTG TCGGTAAAGCAATCTC TCGGTAAGTGACGGTA TCGGTAAGTTTACTCT
#> IFRD1       0.10945041      -0.11593097        1.2291156       -0.4844820
#> CEBPD       0.01311138      -1.06747331        0.2845409       -1.8354638
#> SAT1        1.10940668      -0.25834443        0.5858461        0.1822827
#> CLU         1.71433819       0.36885978        1.5455160        0.2672006
#> MMP7        0.30417809      -0.60270943        1.5425052       -0.7549208
#> CD163      -0.40891468      -0.05234569       -0.4089147       -0.4089147
#>       TCGGTAATCCGCGCAA TCGTACCCAAGTCTGT TCGTACCGTCGAATCT TCGTAGACAACAACCT
#> IFRD1       -0.9214542       -0.6384681        1.3511446        1.6532768
#> CEBPD       -0.3746142       -1.3774915        0.9417899        1.3342060
#> SAT1        -0.4813978        0.5011764        0.9898078        1.2236212
#> CLU         -0.5055495       -1.2361823        1.7371541        1.4737056
#> MMP7         1.2165434       -0.4777718        0.2691457        1.2143897
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       TCGTAGATCGCGTAGC TCTATTGCAAACGCGA TCTATTGGTGCGAAAC TCTATTGTCATGTAGC
#> IFRD1      -0.46620511       -1.2434526      -0.47195659       -1.2434526
#> CEBPD      -0.61106156       -1.8354638       0.04359946       -1.3741985
#> SAT1        0.10315791       -0.3054568       1.35490891       -0.4090545
#> CLU        -0.14239800        0.4164000      -0.66387935       -0.7092258
#> MMP7        0.07198406       -0.4795058       1.32826213       -0.6715184
#> CD163      -0.40891468       -0.4089147      -0.40891468       -0.4089147
#>       TCTCATATCGGTGTCG TCTCTAAAGCGTAATA TCTCTAAGTCAATACC TCTGAGACAGTACACT
#> IFRD1        0.3088750       0.01822946       -1.2434526       -0.6079627
#> CEBPD        1.1876998       0.57876819       -1.2616705       -1.8354638
#> SAT1         0.3876731      -0.77461441       -0.4025515       -0.9855693
#> CLU         -0.7258010      -0.70071070       -0.7122329       -0.9093440
#> MMP7         1.0688888      -0.06247399       -0.7552197       -0.4567473
#> CD163       -0.4089147      -0.40891468        1.1519675       -0.4089147
#>       TCTGAGACATATGCTG TCTGAGATCTAACTGG TCTGGAACAGACAAGC TCTGGAACAGTATGCT
#> IFRD1       -1.2434526       -1.2434526       -1.2434526       -0.3353180
#> CEBPD       -1.3297187       -1.0623110       -1.2122324       -0.7446292
#> SAT1        -1.1124970       -0.2971486       -0.4229982        1.0175571
#> CLU         -1.2361823       -1.2361823       -1.2361823       -1.2361823
#> MMP7        -0.2780584       -0.9602354       -0.9602354       -0.5299399
#> CD163       -0.4089147       -0.4089147        2.9746959        2.8290603
#>       TCTGGAACATGGTCTA TCTGGAAGTCGGATCC TCTTCGGGTGATGTGG TCTTCGGGTTTGTGTG
#> IFRD1       -1.2434526     -0.423961267       -1.2434526        0.5761767
#> CEBPD       -1.4392044     -1.835463810       -0.2868843        1.0965400
#> SAT1         0.3554308      0.006777542       -0.2041936        1.3861193
#> CLU         -0.9669613     -1.236182346       -0.8592928        1.2940006
#> MMP7        -0.5350954     -0.960235411       -0.9602354        1.7890345
#> CD163        3.1690604     -0.408914675        1.4315278       -0.4089147
#>       TCTTTCCAGTACGACG TCTTTCCAGTGGTCCC TCTTTCCTCACCGTAA TCTTTCCTCGGCGCTA
#> IFRD1        1.6273299        0.3228951        1.8461354       -0.2933265
#> CEBPD        0.8347383        0.7782188        1.3350330        1.3618322
#> SAT1         0.8513775        1.2914878        1.3650443       -0.4613032
#> CLU          1.4914544        0.1017960        1.6210213        0.8081396
#> MMP7        -0.2469869        1.3814482        2.3527536       -0.9602354
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       TGAAAGAAGCCACCTG TGAAAGATCGAATGCT TGACAACCACGCGAAA TGACAACCACGGCTAC
#> IFRD1       1.51704182      -0.44153524       -0.5548325       0.41583136
#> CEBPD      -1.05541329      -1.22841369       -0.2058501       1.00736877
#> SAT1        0.08827576       0.01882139        0.2030732      -0.01626901
#> CLU         0.40685303      -1.23618235        0.8689929       0.58801548
#> MMP7       -0.47198298      -0.96023541        0.0597791       1.68664401
#> CD163      -0.40891468       2.12893818       -0.4089147      -0.40891468
#>       TGACGGCTCCTAGTGA TGACGGCTCTTTACAC TGACTAGAGGCAGTCA TGACTAGCAACACCCG
#> IFRD1       -0.2552169       -1.2434526        1.4568414       -1.2434526
#> CEBPD       -0.5054032       -0.7220454        0.7487294       -0.5601542
#> SAT1         0.1772685       -1.3063911        1.5588289        0.7439139
#> CLU         -0.9355465       -0.7575992        1.1708548       -0.8144400
#> MMP7        -0.6832652       -0.9602354        1.0857792       -0.3338665
#> CD163        1.5781922       -0.4089147       -0.4089147        0.8474855
#>       TGACTAGCACGTCTCT TGACTTTAGATGTTAG TGACTTTTCCGCATAA TGAGAGGAGCTAACAA
#> IFRD1        1.7534663       -0.3547042       -0.1405124       -0.5653261
#> CEBPD        1.6656677       -1.1626827        0.6002450       -1.8354638
#> SAT1         1.1469030        0.3435984       -0.4818619        0.4817350
#> CLU          1.0619539       -1.2361823       -0.1548197       -1.2361823
#> MMP7         0.5351488       -0.9602354        0.4885920       -0.9602354
#> CD163       -0.4089147       -0.4089147       -0.4089147        0.6300858
#>       TGAGAGGAGTATCGAA TGAGCCGAGACAAGCC TGAGCCGAGCACGCCT TGAGGGAAGTTAAGTG
#> IFRD1      -0.53770734        0.5334678        0.3454753      1.003522428
#> CEBPD       0.30734683       -0.3382809        1.5284308      0.003971956
#> SAT1       -1.09065462        0.5611214       -0.2743672     -0.186176750
#> CLU        -0.36818547        0.8168128        1.2787419      1.150045377
#> MMP7       -0.96023541       -0.9602354        2.2658290     -0.587114875
#> CD163       0.01900282       -0.4089147       -0.4089147     -0.408914675
#>       TGAGGGATCGTCGTTC TGATTTCGTAACGTTC TGATTTCTCACAGGCC TGCACCTAGAGGTACC
#> IFRD1       -1.0268579      0.008632085        1.0666290       -1.2434526
#> CEBPD       -0.3993483     -1.148469612        1.9934706        1.0797855
#> SAT1        -1.1275070     -1.332393805        1.5298585       -0.5337357
#> CLU         -0.9452903      0.638617846        0.8179169        0.3848676
#> MMP7        -0.6922420     -0.595729264       -0.4702386       -0.6062957
#> CD163       -0.4089147     -0.408914675       -0.4089147       -0.4089147
#>       TGCACCTGTCGCGGTT TGCACCTTCAACACAC TGCACCTTCAACACCA TGCCAAAGTTAGAACA
#> IFRD1       -0.1424623        0.3475263      -0.61496547        0.9166478
#> CEBPD        0.4104202       -0.6310954       0.12263308        1.4955283
#> SAT1         0.0246595       -0.2263815       0.02161845        0.7587229
#> CLU          0.7233058       -1.2361823      -1.23618235        0.2757420
#> MMP7        -0.9602354       -0.9602354      -0.66244343        1.4980686
#> CD163       -0.4089147        1.1569640      -0.40891468       -0.4089147
#>       TGCCCATAGCCTTGAT TGCCCATGTACTTGAC TGCCCATTCAGAGCTT TGCCCTAAGACGACGT
#> IFRD1        2.0369013        1.7315168       -1.2434526       2.57368346
#> CEBPD        0.9751705        1.0551610       -1.2527044       1.22695899
#> SAT1         1.7550629        1.3026053        0.3583381      -0.04177941
#> CLU         -0.6796466        1.1005468       -1.2361823      -0.32300924
#> MMP7         1.6158764        1.7855474       -0.5954722      -0.48569027
#> CD163       -0.4089147       -0.4089147        0.7705893      -0.40891468
#>       TGCGCAGAGAGACTTA TGCGCAGGTACCATCA TGCGCAGGTATATGGA TGCGCAGGTCTTCGTC
#> IFRD1       0.04456899       -0.2716348       -1.2434526       -1.2434526
#> CEBPD       0.35265774       -1.0997993        0.4470053       -0.9350137
#> SAT1        0.19630678       -1.4726356       -1.1092472       -1.6739031
#> CLU        -0.76833891       -1.2361823       -1.2361823       -0.7741614
#> MMP7        0.08256122       -0.4997653       -0.1511424       -0.7118677
#> CD163      -0.40891468       -0.4089147       -0.4089147       -0.4089147
#>       TGCGGGTTCATTGCGA TGCGTGGAGCGAAGGG TGCGTGGCATATGAGA TGCGTGGGTACAGACG
#> IFRD1        1.4155278      -0.10398838       -1.2434526       -0.6274381
#> CEBPD        1.1386946       0.67262980       -1.8354638       -0.2752419
#> SAT1         1.0973368       0.07005769       -1.6852412       -0.4685459
#> CLU          0.7538383      -0.65014539       -1.2361823       -0.7686384
#> MMP7         0.5120192      -0.51797133        0.2557828        0.1125706
#> CD163       -0.4089147      -0.40891468       -0.4089147       -0.4089147
#>       TGCTACCAGATGTCGG TGCTACCCACGAGGTA TGCTACCTCTTACCTA TGCTGCTCAAGCGAGT
#> IFRD1       -0.2687059        0.2071229       -1.2434526        1.6797497
#> CEBPD        0.1851104        0.3321364       -1.8354638        0.7035031
#> SAT1         1.5554348       -0.7241307       -1.2440216       -0.3408651
#> CLU          0.3811828        0.5555114       -1.2361823        0.7703830
#> MMP7         1.5277633       -0.5265396       -0.6133896       -0.8097498
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       TGCTGCTCACCCATTC TGCTGCTGTCTAGAGG TGGACGCAGAAGGTGA TGGACGCCACCAACCG
#> IFRD1       1.03140957      0.717314136       -0.2783682      -0.04830402
#> CEBPD      -1.40289523      0.443607713        0.2013585       0.81936558
#> SAT1       -0.56319908     -0.005085086       -0.2329044      -0.06913308
#> CLU        -0.06619976      1.089861591        0.8246684      -0.77174597
#> MMP7       -0.81031543      1.086399522       -0.3335255       0.15085309
#> CD163      -0.40891468     -0.408914675       -0.4089147      -0.40891468
#>       TGGCCAGAGCAGGTCA TGGCCAGCAGTATGCT TGGCCAGGTCTCATCC TGGCCAGTCAAAGTAG
#> IFRD1       -0.7031774       -0.1069108       -0.3729831      -1.24345261
#> CEBPD       -1.8354638       -0.1434687        0.4036079       0.03574588
#> SAT1        -0.1440609       -1.7860678       -1.5243160       0.37185373
#> CLU         -0.9583138       -0.7574468       -0.6416880      -0.73205514
#> MMP7        -0.7042403        0.8496762       -0.3014067      -0.23225567
#> CD163        2.1410914       -0.4089147       -0.4089147       3.69459347
#>       TGGCCAGTCACAATGC TGGCGCAAGAATAGGG TGGCGCATCGAGAACG TGGCTGGAGCACAGGT
#> IFRD1        1.0626168        0.8411253      -0.05760967        1.7521393
#> CEBPD        0.3483610        0.7126333       0.27628785        1.2687952
#> SAT1        -0.7320799       -1.4028316      -0.19200006        1.6104027
#> CLU         -0.6659598       -0.7564534       0.48608914        1.1906054
#> MMP7         1.3500284       -0.8172714      -0.45369039        0.7339010
#> CD163       -0.4089147       -0.4089147      -0.40891468       -0.4089147
#>       TGGCTGGTCAGGCCCA TGGGCGTAGTATCTCG TGGGCGTGTTCCACGG TGGGCGTTCAAAGACA
#> IFRD1        1.5493612        0.6120473      -1.24345261        1.1883927
#> CEBPD        1.3110700       -0.3407766      -0.09542055        0.8039890
#> SAT1         1.0199960       -0.1515396      -1.31766589       -0.6171047
#> CLU          1.2390841        0.8551601      -0.35000738        0.4001819
#> MMP7         2.1127747       -0.7565873       0.12889962       -0.7137819
#> CD163       -0.4089147       -0.4089147      -0.40891468       -0.4089147
#>       TGGTTAGGTCGGGTCT TGGTTCCCATGCTAGT TGGTTCCTCATGGTCA TGTATTCCATGATCCA
#> IFRD1        2.1136691        0.1243416       -0.2828581        0.4176421
#> CEBPD        1.3137332        0.5438391        0.9167844       -0.9319393
#> SAT1         1.5510641        0.8415871       -0.3550405       -1.2147823
#> CLU          0.4904163        1.3799194       -0.4592640       -0.7488399
#> MMP7         1.6699108       -0.4988292        0.1959304       -0.7409012
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       TGTATTCTCTGAAAGA TGTCCCACATGGTCTA TGTGGTAAGAGCCTAG TGTGGTAGTGTAAGTA
#> IFRD1        1.4226948        1.0573455        1.6687706        1.5335073
#> CEBPD        1.3714543        0.1911440        1.5836654        1.2794552
#> SAT1         1.0654894        0.1388644        0.9690057        0.5265701
#> CLU          1.4473245        0.8342309        0.3554670        0.9615178
#> MMP7         1.3218027       -0.7892512        0.8648983        1.2337011
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       TGTGGTATCCATTCTA TGTTCCGCACCTCGTT TGTTCCGGTATGCTTG TTAACTCGTTAAGGGC
#> IFRD1      -0.36252565      1.611573809       -0.2822915      -1.24345261
#> CEBPD      -0.08855977      0.008210959       -1.1078665      -1.06079865
#> SAT1        0.52247367     -0.004407746        0.2407405       0.06301348
#> CLU        -1.23618235      0.556556112       -1.2361823      -1.23618235
#> MMP7       -0.96023541     -0.471548116       -0.9602354      -0.96023541
#> CD163       0.94080919     -0.408914675       -0.4089147      -0.40891468
#>       TTAGGACAGCCAGTTT TTAGGACAGTGAATTG TTAGGACCACAGGCCT TTAGGACGTGTATGGG
#> IFRD1      0.748345892       -0.6327114        2.3657963       -0.2506594
#> CEBPD     -0.135242912       -0.8192676        1.8071172        1.1607521
#> SAT1       0.226217441       -0.6389020        1.7994338       -0.1176760
#> CLU       -0.006061776        0.7280301        1.6191090        0.0455532
#> MMP7      -0.539125513       -0.6708518        1.1051870        0.4033579
#> CD163     -0.408914675       -0.4089147       -0.4089147       -0.4089147
#>       TTAGTTCAGATCCCGC TTAGTTCAGCGATCCC TTATGCTAGATGTCGG TTATGCTGTTGATTGC
#> IFRD1        1.0197535       -1.2434526       -1.2434526       1.14688997
#> CEBPD        1.8180035       -0.5115976        0.7704522      -0.02597846
#> SAT1         1.5315481        0.2431520        0.6765814      -2.05763700
#> CLU          1.5379431       -1.2361823       -1.2361823      -0.50899780
#> MMP7         1.1450551       -0.5520146       -0.5369940      -0.44497721
#> CD163       -0.4089147        0.9111146        4.7159877      -0.40891468
#>       TTCGAAGAGTGACTCT TTCGAAGCAATGGTCT TTCGAAGTCAAACAAG TTCGAAGTCACCTCGT
#> IFRD1       -1.2434526       -1.2434526       -1.2434526        1.6974448
#> CEBPD       -0.8083698       -1.0352527       -0.5652198        1.5220687
#> SAT1        -1.6069645       -1.6691988        0.1105442        0.3851043
#> CLU         -1.2361823       -0.6925141       -1.2361823        0.9231253
#> MMP7        -0.9602354       -0.9602354       -0.9602354        0.9975801
#> CD163       -0.4089147       -0.4089147        0.8409104       -0.4089147
#>       TTCGAAGTCAGGTAAA TTCGAAGTCTTGCCGT TTCGGTCAGTGGCACA TTCGGTCTCCTGCAGG
#> IFRD1       1.31655885       -0.2407721        1.5530798        1.2527069
#> CEBPD       0.29623637        0.5815657        1.0568702        1.1380375
#> SAT1        0.74164517       -0.7790077        0.6954166        1.4459954
#> CLU         0.08045515       -0.7204945        1.2258322        2.0631016
#> MMP7       -0.61218120       -0.9602354        0.8501462        1.2186421
#> CD163      -0.40891468        1.9901633       -0.4089147       -0.4089147
#>       TTCGGTCTCTGGCGAC TTCTACAAGATGGGTC TTCTACACAATACGCT TTCTCAAAGCGTTGCC
#> IFRD1       -1.2434526        0.2665190      -1.24345261       -1.2434526
#> CEBPD       -1.8354638        1.4450241      -0.09137456        0.1597633
#> SAT1        -0.1832067        0.2639991       1.24813888       -0.1512475
#> CLU         -0.7663960        1.1950587      -1.23618235       -1.2361823
#> MMP7        -0.9602354        1.6449231      -0.54386976       -0.4465889
#> CD163       -0.4089147       -0.4089147       2.75627915        3.2457085
#>       TTCTCAAGTCCGAACC TTCTCAATCGCCTGAG TTCTCAATCTGATTCT TTCTTAGGTTTGGCGC
#> IFRD1        0.3319868       -1.2434526        1.3755660        0.9730387
#> CEBPD        0.8031048        0.3713986        1.6222564        0.3794322
#> SAT1         0.5259014        0.1469381        0.5729475       -0.8748288
#> CLU          1.4888180       -1.2361823        2.1003377        0.8182531
#> MMP7        -0.7573539       -0.9602354        1.2533623       -0.6880436
#> CD163       -0.4089147        4.0577708       -0.4089147       -0.4089147
#>       TTCTTAGTCATTGCGA TTGACTTCATATACCG TTGACTTTCTTTAGGG TTGCGTCAGAGACTTA
#> IFRD1       -0.4849528        2.1315758       -1.2434526       -0.3736077
#> CEBPD        0.5687436        1.2715817       -1.3202828       -0.7835320
#> SAT1        -1.2747793        0.9856461       -1.0939312       -1.4724911
#> CLU          0.1743302        1.0054315       -0.6563853       -0.7888131
#> MMP7        -0.4772302        1.6052376       -0.9602354       -0.7207967
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       TTGCGTCAGGTAGCCA TTGGAACAGGGCATGT TTGGAACAGTAAGTAC TTGGAACCACAGTCGC
#> IFRD1       -1.2434526        1.1573454       -0.4926247       -1.2434526
#> CEBPD       -1.8354638        0.8314371        0.8029270        0.2720784
#> SAT1        -2.1138091        0.4956199        1.1294902        1.1696423
#> CLU         -1.2361823        1.6595215       -0.8500246       -0.6986860
#> MMP7        -0.9602354        1.2267494        1.4768230       -0.4650500
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       TTGGAACCACCGTTGG TTGGCAATCAGTTTGG TTGGCAATCGTCTGAA TTGGCAATCGTTTGCC
#> IFRD1        0.2642041      -0.12050215       -1.2434526        0.8967556
#> CEBPD       -0.6941702       0.94712201       -0.5319091        0.3075596
#> SAT1        -0.5572961       0.50451150       -0.5264554        0.4020245
#> CLU          0.1225457       0.80488476       -0.8020602        0.3456640
#> MMP7        -0.6625085       0.04929864       -0.1443091       -0.3297585
#> CD163       -0.4089147      -0.40891468       -0.4089147       -0.4089147
#>       TTGGCAATCTTGTTTG TTGTAGGAGAGCCCAA TTGTAGGGTCGACTGC TTGTAGGTCTTCAACT
#> IFRD1        1.9492157       -1.2434526       -1.2434526        0.1243151
#> CEBPD        1.6845652       -0.2648168       -1.3420133        1.6243674
#> SAT1         1.7918685        0.5353441        0.4540881        0.3277215
#> CLU          1.2248490       -1.2361823       -1.2361823        0.1765795
#> MMP7         2.3243213       -0.9602354       -0.6513728        1.8738632
#> CD163       -0.4089147        4.2279205        3.2506690       -0.4089147
#>       TTTACTGAGGACGAAA TTTACTGCATGCTAGT TTTACTGGTACTCGCG TTTATGCGTAAAGTCA
#> IFRD1       0.28597618       -1.2434526       0.97445692        0.7304483
#> CEBPD      -1.09746811        0.2361317       0.06169439        1.7657957
#> SAT1       -0.04924902        0.7521621      -0.66074440        0.7002965
#> CLU        -0.73478363       -1.2361823       1.35587506        1.1195020
#> MMP7       -0.96023541       -0.9602354      -0.56755833        1.9890352
#> CD163      -0.40891468        3.0874466      -0.40891468       -0.4089147
#>       TTTATGCTCCTCATTA TTTATGCTCTGTTGAG TTTCCTCTCGGAAACG TTTGCGCCAATCACAC
#> IFRD1        1.5560445        1.4711768        0.7873104        0.2075723
#> CEBPD        1.4842779        0.3228759        1.1944232        0.1377859
#> SAT1         2.1969442       -0.1251083        0.4838643       -0.2725630
#> CLU          1.6289273        0.2302066        1.8522125        0.6876170
#> MMP7         1.2259184        0.7870085        1.2434796       -0.8460843
#> CD163       -0.4089147       -0.4089147       -0.4089147       -0.4089147
#>       TTTGTCATCTTGTATC
#> IFRD1        1.5677657
#> CEBPD        0.4547755
#> SAT1         0.1227886
#> CLU         -0.3009281
#> MMP7        -0.6689371
#> CD163        0.1164242
#> [1] ">>>>> Features not exited in matrix data..."
#>    Mode    TRUE 
#> logical      55 
#> character(0)
#> 载入需要的程辑包：AnnotationDbi
#> 
#> 载入程辑包：'AnnotationDbi'
#> The following object is masked from 'package:dplyr':
#> 
#>     select
#> 
#> 'select()' returned 1:1 mapping between keys and columns
#> 'select()' returned 1:1 mapping between keys and columns
#> 'select()' returned 1:1 mapping between keys and columns
#> 'select()' returned 1:1 mapping between keys and columns
#> 'select()' returned 1:1 mapping between keys and columns
#> 'select()' returned 1:1 mapping between keys and columns
#> 'select()' returned 1:1 mapping between keys and columns
#> 'select()' returned 1:1 mapping between keys and columns
#> 'select()' returned 1:1 mapping between keys and columns
#> 'select()' returned 1:1 mapping between keys and columns
#> 'select()' returned 1:1 mapping between keys and columns
#> >>>>--- Processsing terms: enrichGO
#>   Cluster         ID                                 Description GeneRatio
#> 1       0 GO:0005201 extracellular matrix structural constituent    13/268
#> 2       0 GO:0044548                        S100 protein binding     4/268
#> 3       0 GO:0019887           protein kinase regulator activity    12/268
#> 4       0 GO:0004857                   enzyme inhibitor activity    17/268
#> 5       0 GO:0004860           protein kinase inhibitor activity     7/268
#> 6       0 GO:0019210                   kinase inhibitor activity     7/268
#>     BgRatio       pvalue     p.adjust       qvalue
#> 1 172/18410 1.466386e-06 0.0007170629 0.0006282307
#> 2  14/18410 3.918438e-05 0.0068689337 0.0060179873
#> 3 207/18410 5.362759e-05 0.0068689337 0.0060179873
#> 4 390/18410 6.187818e-05 0.0068689337 0.0060179873
#> 5  70/18410 7.023450e-05 0.0068689337 0.0060179873
#> 6  74/18410 1.003521e-04 0.0074791242 0.0065525854
#>                                                                                geneID
#> 1                  157869/4256/1116/633/4148/5549/3490/115908/1301/176/4237/1299/3914
#> 2                                                                6271/6285/27101/6277
#> 3                  2810/1026/3315/2621/1164/9021/84152/9052/1163/1029/100133941/10221
#> 4 2810/144203/1026/3315/6590/2621/3958/5270/9021/84152/183/2/5104/1029/301/5268/10221
#> 5                                                2810/1026/3315/9021/84152/1029/10221
#> 6                                                2810/1026/3315/9021/84152/1029/10221
#>   Count
#> 1    13
#> 2     4
#> 3    12
#> 4    17
#> 5     7
#> 6     7
#> Scale for 'colour' is already present. Adding another scale for 'colour',
#> which will replace the existing scale.
#> Scale for 'y' is already present. Adding another scale for 'y', which will
#> replace the existing scale.
#> Scale for 'size' is already present. Adding another scale for 'size', which
#> will replace the existing scale.
#> >>>>Options for `theme`: light, bw, classic and classic2
#> >>>>Options for 'legend.position' : none, left, right, bottom, top
#> >>>>Options for 'legend.direction' : horizontal, vertical
```

<img src="man/figuresunnamed-chunk-6-1.png" width="100%" />

    #> >>>>--- Processsing terms: enrichKEGG
    #>   Cluster       ID                         Description GeneRatio  BgRatio
    #> 1       0 hsa04310               Wnt signaling pathway    10/142 170/8217
    #> 2       1 hsa05323                Rheumatoid arthritis    31/274  93/8217
    #> 3       1 hsa04145                           Phagosome    36/274 152/8217
    #> 4       1 hsa05140                       Leishmaniasis    27/274  77/8217
    #> 5       1 hsa05150     Staphylococcus aureus infection    29/274  96/8217
    #> 6       1 hsa04612 Antigen processing and presentation    23/274  78/8217
    #>         pvalue     p.adjust       qvalue
    #> 1 6.980089e-04 1.786903e-01 1.741348e-01
    #> 2 2.163980e-23 5.388310e-21 4.054615e-21
    #> 3 2.397883e-21 2.985364e-19 2.246437e-19
    #> 4 3.673899e-21 3.049336e-19 2.294576e-19
    #> 5 1.425411e-20 8.873182e-19 6.676924e-19
    #> 6 2.855679e-16 1.422128e-14 1.070128e-14
    #>                                                                                                                                                                                  geneID
    #> 1                                                                                                                               4316/25805/6422/7855/27101/147111/1501/81839/85409/4772
    #> 2                            3117/3118/54/3127/3109/3606/942/3113/3689/3108/6348/2921/3123/3119/2920/3122/3115/3576/1514/7040/3111/414062/6352/3553/7099/6347/2919/3383/10312/1513/6387
    #> 3 1536/4481/2209/2212/3117/3118/2214/3127/3109/1520/929/3113/3689/3108/3123/3119/1535/3122/3115/4688/2213/4689/1514/3105/653361/3111/4973/7099/3106/948/11151/81035/10312/718/3134/6890
    #> 4                                               1536/2209/2212/3117/3118/2214/3127/3109/3113/3689/3108/3123/3119/1535/3122/3115/4688/4689/653361/7040/3111/3586/3553/7099/4792/5777/718
    #> 5                                             2359/714/719/2209/728/2212/3117/3118/712/2214/3127/3109/713/3113/3689/3108/3123/3119/3122/3115/5724/2213/1675/3111/717/3586/6404/3383/718
    #> 6                                                                      920/3117/3118/3127/3109/1520/3113/3108/10437/1508/3123/3119/972/3122/3115/567/1514/3105/3111/5641/3106/3134/6890
    #>   Count
    #> 1    10
    #> 2    31
    #> 3    36
    #> 4    27
    #> 5    29
    #> 6    23
    #> Scale for 'colour' is already present. Adding another scale for 'colour',
    #> which will replace the existing scale.
    #> Scale for 'y' is already present. Adding another scale for 'y', which will
    #> replace the existing scale.
    #> Scale for 'size' is already present. Adding another scale for 'size', which
    #> will replace the existing scale.
    #> >>>>Options for `theme`: light, bw, classic and classic2
    #> >>>>Options for 'legend.position' : none, left, right, bottom, top
    #> >>>>Options for 'legend.direction' : horizontal, vertical

<img src="man/figuresunnamed-chunk-6-2.png" width="100%" />

    #> >>>--- Processing cluster: 0
    #>               p_val avg_diff pct.1 pct.2    p_val_adj
    #> IFRD1  3.772950e-83 1.559992 0.963 0.315 1.131885e-79
    #> CEBPD  8.688609e-81 1.428352 0.968 0.387 2.606583e-77
    #> SAT1   8.133918e-77 1.386217 0.959 0.412 2.440175e-73
    #> CLU    1.793408e-72 1.417723 0.945 0.354 5.380225e-69
    #> MMP7   8.665662e-72 1.514401 0.886 0.246 2.599699e-68
    #> SBSPON 1.308976e-70 1.400558 0.836 0.190 3.926929e-67
    #> [1] "Default assay is RNA"
    #> [1] ">>>-- Colors could be change by parameter: 'cols'"
    #> Warning: 'ncol' is ignored with 'stack' is TRUE
    #> >>>--- Processing cluster: 1
    #>                 p_val avg_diff pct.1 pct.2     p_val_adj
    #> CD163   2.265804e-151 2.097133 0.839 0.040 6.797411e-148
    #> MS4A7   1.269907e-139 2.179377 0.911 0.076 3.809722e-136
    #> CYBB    2.875015e-139 2.079585 0.906 0.076 8.625044e-136
    #> FPR3    5.081822e-136 2.072075 0.856 0.068 1.524547e-132
    #> SLCO2B1 1.028761e-132 2.016369 0.756 0.043 3.086282e-129
    #> MSR1    1.424992e-132 2.080712 0.811 0.063 4.274976e-129
    #> [1] "Default assay is RNA"
    #> [1] ">>>-- Colors could be change by parameter: 'cols'"
    #> Warning: 'ncol' is ignored with 'stack' is TRUE
    #> >>>--- Processing cluster: 10
    #>                  p_val avg_diff pct.1 pct.2     p_val_adj
    #> ATP6V0D2 7.748189e-119 5.681203 0.810 0.010 2.324457e-115
    #> RUFY4    4.422049e-117 5.317588 0.714 0.007 1.326615e-113
    #> SIGLEC15 2.321180e-106 6.131332 0.952 0.022 6.963540e-103
    #> STRIP2   4.587842e-103 4.320153 0.571 0.004  1.376353e-99
    #> DPP4     3.257920e-101 5.366486 0.810 0.015  9.773761e-98
    #> SLC9B2    2.241746e-89 6.456339 1.000 0.035  6.725238e-86
    #> [1] "Default assay is RNA"
    #> [1] ">>>-- Colors could be change by parameter: 'cols'"
    #> Warning: 'ncol' is ignored with 'stack' is TRUE
    #> >>>--- Processing cluster: 2
    #>               p_val avg_diff pct.1 pct.2    p_val_adj
    #> PCSK1N 1.290257e-82 1.954099 0.982 0.304 3.870770e-79
    #> CALML5 6.381724e-81 1.959466 1.000 0.214 1.914517e-77
    #> CFB    4.616001e-79 1.600839 0.833 0.168 1.384800e-75
    #> SEZ6L2 2.384760e-74 1.551717 0.673 0.097 7.154279e-71
    #> S100A6 3.319170e-74 1.628335 0.970 0.367 9.957511e-71
    #> APOD   4.931076e-65 1.819956 0.851 0.159 1.479323e-61
    #> [1] "Default assay is RNA"
    #> [1] ">>>-- Colors could be change by parameter: 'cols'"
    #> Warning: 'ncol' is ignored with 'stack' is TRUE
    #> >>>--- Processing cluster: 3
    #>                 p_val avg_diff pct.1 pct.2    p_val_adj
    #> TMEM106C 9.634129e-58 1.550635 0.959 0.342 2.890239e-54
    #> MARCKSL1 1.079709e-54 1.231641 0.973 0.519 3.239128e-51
    #> AZGP1    1.356838e-52 1.530683 0.939 0.293 4.070513e-49
    #> IDH1     1.519149e-50 1.427358 0.959 0.331 4.557448e-47
    #> CLPSL1   2.063225e-50 1.544419 0.797 0.263 6.189676e-47
    #> IGFBP2   1.193092e-46 1.323517 0.878 0.342 3.579275e-43
    #> [1] "Default assay is RNA"
    #> [1] ">>>-- Colors could be change by parameter: 'cols'"
    #> Warning: 'ncol' is ignored with 'stack' is TRUE
    #> >>>--- Processing cluster: 4
    #>               p_val avg_diff pct.1 pct.2    p_val_adj
    #> CLDN6  7.172546e-86 1.898127 0.860 0.145 2.151764e-82
    #> SMIM22 1.385336e-78 1.850544 0.860 0.165 4.156007e-75
    #> PCAT19 2.292897e-73 1.443169 0.676 0.079 6.878691e-70
    #> MAGIX  1.984059e-67 1.835863 0.956 0.278 5.952177e-64
    #> RAB6B  3.165224e-63 1.653225 0.647 0.099 9.495672e-60
    #> PQBP1  1.131886e-60 1.681424 0.985 0.367 3.395659e-57
    #> [1] "Default assay is RNA"
    #> [1] ">>>-- Colors could be change by parameter: 'cols'"
    #> Warning: 'ncol' is ignored with 'stack' is TRUE
    #> >>>--- Processing cluster: 5
    #>                 p_val avg_diff pct.1 pct.2    p_val_adj
    #> TOX      3.680723e-46 2.080471 0.955 0.226 1.104217e-42
    #> SLC26A7  1.205001e-35 1.652019 0.612 0.100 3.615003e-32
    #> SBSPON   9.894931e-31 1.642054 0.940 0.279 2.968479e-27
    #> C10orf99 1.732970e-30 1.478415 0.239 0.012 5.198909e-27
    #> COL11A2  3.837653e-30 1.594112 0.761 0.195 1.151296e-26
    #> S100P    5.688115e-30 1.534660 1.000 0.435 1.706435e-26
    #> [1] "Default assay is RNA"
    #> [1] ">>>-- Colors could be change by parameter: 'cols'"
    #> Warning: 'ncol' is ignored with 'stack' is TRUE
    #> >>>--- Processing cluster: 6
    #>                p_val avg_diff pct.1 pct.2     p_val_adj
    #> CFH    1.549368e-162 3.855307 0.852 0.012 4.648103e-159
    #> NTM    1.774126e-160 3.644234 0.741 0.004 5.322378e-157
    #> PDGFRB 6.010230e-159 3.979637 0.907 0.018 1.803069e-155
    #> FBN1   1.553514e-150 3.901692 0.852 0.016 4.660541e-147
    #> ECM2   5.919022e-149 3.562662 0.722 0.006 1.775707e-145
    #> CHN1   7.750026e-148 3.859079 0.852 0.017 2.325008e-144
    #> [1] "Default assay is RNA"
    #> [1] ">>>-- Colors could be change by parameter: 'cols'"
    #> Warning: 'ncol' is ignored with 'stack' is TRUE
    #> >>>--- Processing cluster: 7
    #>                p_val avg_diff pct.1 pct.2    p_val_adj
    #> FCER1A  6.304196e-77 2.772433 0.408 0.005 1.891259e-73
    #> CD1E    1.332495e-69 2.541608 0.347 0.003 3.997484e-66
    #> CD1A    1.046284e-63 2.135411 0.306 0.002 3.138853e-60
    #> CLEC10A 2.176725e-54 1.811671 0.265 0.002 6.530175e-51
    #> CD1C    6.947607e-51 1.854845 0.265 0.003 2.084282e-47
    #> ASGR2   4.367738e-49 1.843508 0.286 0.005 1.310321e-45
    #> [1] "Default assay is RNA"
    #> [1] ">>>-- Colors could be change by parameter: 'cols'"
    #> Warning: 'ncol' is ignored with 'stack' is TRUE
    #> >>>--- Processing cluster: 8
    #>               p_val avg_diff pct.1 pct.2    p_val_adj
    #> NEK2   1.803522e-61 4.430332 1.000 0.084 5.410566e-58
    #> PBK    2.064247e-60 4.235198 1.000 0.086 6.192741e-57
    #> CENPA  4.710395e-58 3.448689 0.964 0.077 1.413118e-54
    #> CDC25C 3.077076e-57 2.845775 0.607 0.023 9.231229e-54
    #> NMU    1.576128e-54 3.614909 0.786 0.051 4.728384e-51
    #> KIF20A 6.266232e-54 3.619765 0.786 0.051 1.879870e-50
    #> [1] "Default assay is RNA"
    #> [1] ">>>-- Colors could be change by parameter: 'cols'"
    #> Warning: 'ncol' is ignored with 'stack' is TRUE
    #> >>>--- Processing cluster: 9
    #>              p_val avg_diff pct.1 pct.2    p_val_adj
    #> RIBC2 9.166234e-77 3.465913 0.704 0.021 2.749870e-73
    #> MCM10 7.017606e-59 3.548696 0.889 0.059 2.105282e-55
    #> DTL   2.417974e-58 3.473248 0.852 0.053 7.253921e-55
    #> EXO1  9.729034e-52 2.967511 0.741 0.043 2.918710e-48
    #> CDC45 5.391993e-48 3.232693 0.778 0.056 1.617598e-44
    #> CDC6  4.520882e-47 2.414208 0.741 0.048 1.356265e-43
    #> [1] "Default assay is RNA"
    #> [1] ">>>-- Colors could be change by parameter: 'cols'"
    #> Warning: 'ncol' is ignored with 'stack' is TRUE

<img src="man/figuresunnamed-chunk-6-3.png" width="100%" />

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
                         method             = "doubletfinder",
                         save_path          = "./test")
#> An object of class Seurat 
#> 33694 features across 1097 samples within 1 assay 
#> Active assay: RNA (33694 features, 3000 variable features)
#>  3 dimensional reductions calculated: pca, tsne, umap
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
    #>      25    1072

``` r
head(sce_df@meta.data)
#>                  orig.ident nCount_RNA nFeature_RNA percent.mt percent.rp
#> AAACCTGCACCTTGTC       TNBC      13976         3590  16.463938  10.997424
#> AAACGGGAGTCCTCCT       TNBC       8732         2629   5.817682   9.757215
#> AAACGGGTCCAGAGGA       TNBC      17138         4166  10.228731   9.738593
#> AAAGATGCAGTTTACG       TNBC       9519         1703   5.200126  16.482824
#> AAAGCAACAGGAATGC       TNBC      10285         2778  18.395722  13.349538
#> AAAGCAATCGGAATCT       TNBC       8436         2822  16.797060   5.263158
#>                      S.Score   G2M.Score Phase RNA_snn_res.0.4 seurat_clusters
#> AAACCTGCACCTTGTC  0.07723506  0.42077041   G2M               3               8
#> AAACGGGAGTCCTCCT -0.07738382 -0.13026069    G1               2               1
#> AAACGGGTCCAGAGGA -0.03538864 -0.09533497    G1               0               0
#> AAAGATGCAGTTTACG -0.04374442 -0.09006494    G1               1               2
#> AAAGCAACAGGAATGC -0.07575323 -0.07794991    G1               0               5
#> AAAGCAATCGGAATCT  0.12829573  0.35263693   G2M               1               4
#>                  RNA_snn_res.0.8 RNA_snn_res.1.2 RNA_snn_res.1.6 RNA_snn_res.2
#> AAACCTGCACCTTGTC               7               8              11            11
#> AAACGGGAGTCCTCCT               1               4               3             3
#> AAACGGGTCCAGAGGA               0               0               6             7
#> AAAGATGCAGTTTACG               2               1               0             0
#> AAAGCAACAGGAATGC               0               6               5             5
#> AAAGCAATCGGAATCT               4               3              10            10
#>                  RNA_snn_res.1 pANN_0.25_0.01_27
#> AAACCTGCACCTTGTC             8        0.33333333
#> AAACGGGAGTCCTCCT             1        0.00000000
#> AAACGGGTCCAGAGGA             0        0.06666667
#> AAAGATGCAGTTTACG             2        0.00000000
#> AAAGCAACAGGAATGC             5        0.00000000
#> AAAGCAATCGGAATCT             4        0.53333333
#>                  DF.classifications_0.25_0.01_27 DoubletFinder_final
#> AAACCTGCACCTTGTC                         Singlet             Singlet
#> AAACGGGAGTCCTCCT                         Singlet             Singlet
#> AAACGGGTCCAGAGGA                         Singlet             Singlet
#> AAAGATGCAGTTTACG                         Singlet             Singlet
#> AAAGCAACAGGAATGC                         Singlet             Singlet
#> AAAGCAATCGGAATCT                         Doublet             Doublet
```

``` r

p1<-dong_dimplot(sce              = sce_df,
                 reduction        = "umap",
                 groups           = "DoubletFinder_final",
                 split.by         = NULL,
                 label            = T,
                 label.size       = 5,
                 pt.size          = 0.5,
                 cols             = "normal", 
                 seed             = 54321, 
                 show_col         = F,
                 width            = 8, 
                 height           = 8, 
                 w_index          = 7, 
                 w_add            = 2,
                 max_category     = 14,
                 show_plot        = T,
                 path             = "./test",
                 index            = paste0("1-DoubletFinder"),
                 legend.position  = "right",
                 legend.direction = "vertical",
                 legend.size      = 10)
#> >>>>=== Palette option for random: 1: palette1; 2: palette2; 3: palette3;  4: palette4
#> >>> Processing group:: DoubletFinder_final
#> >>>>Options for `theme`: light, bw, classic and classic2
#> >>>>Options for 'legend.position' : none, left, right, bottom, top
#> >>>>Options for 'legend.direction' : horizontal, vertical
#> >>>>>> Figure name is:: 1-DoubletFinder-1-umap-DoubletFinder_final.pdf
```

<img src="man/figuresunnamed-chunk-9-1.png" width="100%" />

### Citation

If you use scsig in published research, please cite:

### Contact

E-mail any questions to <dongqiangzeng0808@gmail.com>
