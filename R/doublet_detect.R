








#' Integration of  doublet detection
#'
#' @param eset expression set, default is null
#' @param file_type if file is from 10X, data_path must be provide
#' @param data_path data_path of 10X data
#' @param project project name
#' @param check_data if true
#' @param index folder name of save_path
#' @param minFeature minimum feature of cell
#' @param minCount minimum count of cell
#' @param percent.mt max of percentage of mitochondrial gens
#' @param filter_data if TRUE, cells will be filtered by aforementioned criterion
#' @param method default is `doubletfinder`: Benchmarking Computational Doublet-Detection Methods for Single-Cell RNA Sequencing Data
#' @param sce default is null, a seurat object
#' @param cores default is 1
#' @param already_normalized
#' @param propotion default is 0.025
#' @param save_path default is NULL
#'
#' @return
#' @export
#'
#' @examples
doublet_detect<-function(sce                = NULL,
                         already_normalized = FALSE,
                         eset               = NULL,
                         file_type          = "10X",
                         data_path          = NULL,
                         project            = NULL,
                         index              = 1,
                         propotion          = 0.025,
                         check_data         = TRUE,
                         filter_data        = TRUE,
                         minFeature         = 2000,
                         minCount           = 1000,
                         percent.mt         = 20,
                         method             = "doubletfinder",
                         cores              = 1,
                         save_path          = NULL){

  if(is.null(eset) & file_type =="10X" & is.null(sce)){
    sce<-CreateSeuratObject(Read10X(paste0(data_path)),project)
  }

  if(!is.null(eset)){
    sce <- CreateSeuratObject(eset,project)
  }

  if(!is.null(sce)){
    sce<-sce
  }

  print(sce)


  if(!is.null(save_path)){
    file_name<-creat_folder(save_path)
  }else{
    file_name<-mydb:: creat_folder(paste0(index,"-",project,"-doublet-detection"))
  }

  abspath<-file_name$abspath
  message(" >>>---  Result will be deposite in ",abspath)
  #############################################
  #############################################

  if(already_normalized) message(">>>--- Input data has been normalized, QC and data filtering will not be proceed...")


  if(c(check_data|filter_data) & !already_normalized){

    print(">>>---- Quality control")
    if(length(unique(grepl('^mt-',rownames(sce))))<=1){
      message(paste0("Mitochondrial genes have been removed"))
      mito<-TRUE
    }else{
      print("Mitochondrial genes")
      print(rownames(sce)[grepl('^mt-',rownames(sce))])
    }

    if(length(unique(grepl('^Rp[sl]',rownames(sce))))<=1){
      message(paste0("Ribosomal genes have been removed"))
    }else{
      print("Ribosomal genes")
      print(rownames(sce)[grepl('^Rp[sl]',rownames(sce))])
    }

    #################################################
    sce[["percent.mt"]] <- Seurat:: PercentageFeatureSet(sce, pattern = "^mt-")
    rb.genes <- rownames(sce)[grep("^Rp[sl]",rownames(sce))]
    C<-Seurat:: GetAssayData(object = sce, slot = "counts")
    percent.ribo <- Matrix::colSums(C[rb.genes,])/Matrix::colSums(C)*100
    sce <-Seurat:: AddMetaData(sce, percent.ribo, col.name = "percent.ribo")
    print(head(sce@meta.data))
    ##################################################

    if(check_data){
      plot1 <- Seurat::FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "percent.mt")
      plot2 <- Seurat:: FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
      pp<-Seurat:: CombinePlots(plots = list(plot1, plot2))
      print(pp)
      # ggsave(pp,filename = paste0("1-Mitocondrial-genes-and-feature-count.pdf"),width = 8.5, height = 5.8,path = path1$folder_name)
      pp<- VlnPlot(sce, features = c("nFeature_RNA", "nCount_RNA","percent.ribo", "percent.mt"), ncol = 2)
      print(pp)
      # ggsave(pp,filename = paste0("2-Mitocondrial-and-ribosomal-gene-count.pdf"),width = 8.5, height = 5.8,path = path1$folder_name)
      # Filtering data
      # Idents(sce)<-c("nCount_RNA","nFeature_RNA","percent.mt")
    }
    if(filter_data){
      message("Default parameters are : minFeature = 2000, minCount = 1000, percent.mt =20.")
      if(mito){
        sce <-subset(sce, subset = nFeature_RNA > minFeature & nCount_RNA > minCount)
      }else{
        sce <-subset(sce, subset = nFeature_RNA > minFeature & nCount_RNA > minCount & percent.mt < percent.mt)
      }
      message(">> After data filtering...")

    }
    print(sce)
  }



  if(method == "scds"){

    # devtools::install_github('kostkalab/scds',ref="master")
    message(">>>--- Running scds doublet detection...")
    db_res <- Seurat:: as.SingleCellExperiment(sce) #convert seurat subject into SingleCellExperiment object
    logcounts(db_res) = log1p(counts(db_res))

    db_res = scds:: cxds(db_res, retRes=T)  #calculate cxds score
    db_res = scds:: bcds(db_res, retRes=T, estNdbl=T, verb=T)
    db_res = scds:: cxds_bcds_hybrid(db_res)


    print("Results of bcds: ")
    print(table(colData(db_res)$bcds_call))


    sce[["cxds_score"]] <- colData(db_res)$cxds_score
    # sce[["cxds_call"]] <- colData(db_res)$cxds_call
    sce[["bcds_score"]] <- colData(db_res)$bcds_score
    sce[["bcds_call"]] <- colData(db_res)$bcds_call
    sce[["hybrid_score"]] <- colData(db_res)$hybrid_score

    meta<-sce@meta.data
    meta<-rownames_to_column(meta, var = "rowname")
    writexl::write_xlsx(meta, paste0(file_name$abspath, "0-Doublet-detected-by-",method,".xlsx") )
    # distribution of bcds to define the best cutoff of doublet call

    pdf(file = paste0(file_name$abspath, "4-Doublet-detection-by-",method,".pdf"),width =9.71, height = 5.39 )
    par(mfrow=c(1,4))
    h1 <- hist(sce$cxds_score, main="cxds_score", xlab="")
    h2 <- hist(sce$bcds_score, main="bcds_score", xlab="")
    #define doublets as those with bcds > 0.6
    cutoff_bcds <- 0.4
    abline(v=cutoff_bcds, col="red")
    h3 <- hist(sce$hybrid_score, main="hybrid_score", xlab="")
    boxplot(sce$bcds_score ~ sce$bcds_call, main="bcds")
    invisible(dev.off())

  }


  if(method == "doubletfinder"){

    # remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
    message(">>>----- Running DoubletFinder doublet detection...")
    ## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
    DefaultAssay(sce) <- "RNA"
    if(!already_normalized){

      sce <- NormalizeData(sce)
      sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
      sce <- ScaleData(sce)
      sce <- RunPCA(sce)
      sce <- RunUMAP(sce, dims = 1:30)
    }

    ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
    sweep.res.list_sce <- paramSweep_v3(sce, PCs = 1:15, sct = FALSE,num.cores = cores)
    sweep.stats_sce <- summarizeSweep(sweep.res.list_sce, GT = FALSE)
    bcmvn_sce <- find.pK(sweep.stats_sce)

    pK <- bcmvn_sce$pK[which.max(bcmvn_sce$BCmetric)]; pK <- as.numeric(levels(pK))[pK]
    pK

    ##  Doublet Proportion Estimate -------------------------------------------------------------------------------------

    nExp_poi <- round(propotion*nrow(sce@meta.data))  ## Assuming 2.5% doublet formation rate - tailor for your dataset

    ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
    sce <- doubletFinder_v3(sce, PCs = 1:15, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

    attribute <- paste('pANN', 0.25, pK ,nExp_poi, sep = '_')
    attribute
    score <- sce@meta.data[[attribute]]; score

    t <- sort(score, decreasing = TRUE)[nExp_poi]
    # remove predicted doublet
    pred.index <- which(score > t)
    length(pred.index)

    sce@meta.data[["DoubletFinder_final"]]<-ifelse(sce@meta.data[[attribute]]>t, "Doublet","Singlet" )

    message(paste0(">>> Doublet detected by : DoubleFinder R package, with parametar: pN = 0.25, proportion = ", propotion, ", pK = ", pK))
    print(summary(as.factor( sce@meta.data[["DoubletFinder_final"]])))

    message(">>> The column name of result is : DoubletFinder_final.")
    meta<-sce@meta.data
    meta<-rownames_to_column(meta, var = "rowname")
    writexl::write_xlsx(meta, paste0(file_name$abspath, "0-Doublet-detected-by-",method,".xlsx") )


  }



  return(sce)

}
