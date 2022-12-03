




#' modified_copykat
#'
#' @param eset expression set, default is null
#' @param file_type if file is from 10X, data_path must be provide
#' @param data_path data_path of 10X data
#' @param project project name
#' @param id_type symbol or ensembl
#' @param check_data if true
#' @param index folder name of save_path
#' @param minFeature minimum feature of cell
#' @param minCount minimum count of cell
#' @param percent.mt max of percentage of mitochondrial gens
#' @param filter_data if TRUE, cells will be filtered by aforementioned criterion
#' @param sce seurat object
#' @param normal_cell_id
#'
#' @return
#' @export
#'
#' @examples
#' copykat_res<-copykat_plus(eset = sc_tnbc, file_type = "10X", data_path = NULL, index = 1, project = "TNBC", id_type = "symbol", check_data = TRUE, filter_data = FALSE, minFeature = 2000, minCount = 1000, percent.mt =20)
#'
copykat_plus<-function(sce,
                       eset        = NULL,
                       file_type   = "10X",
                       data_path   = NULL,
                       project     = NULL,
                       index       = 1,
                       id_type     = "symbol",
                       normal_cell_id = "",
                       check_data  = TRUE,
                       filter_data = TRUE,
                       minFeature  = 2000,
                       minCount    = 1000,
                       percent.mt  =20){

  if(is.null(eset) & file_type =="10X"){
    sce<-CreateSeuratObject(Read10X(paste0(data_path)),project)
  }

  if(!is.null(eset)){
    sce <- CreateSeuratObject(eset, project)
  }
  print(sce)

  file_name<-mydb:: creat_folder(paste0(index,"-",project))
  abspath<-file_name$abspath
  message(" >>>  Data will be deposite in ",abspath)
  ##############################################
  path1<-creat_folder(paste0(file_name$folder_name,"/2-CNV-estimation"))
  #############################################
  print(">>> Quality control")
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
    plot1 <-Seurat::  FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "percent.mt")
    plot2 <- Seurat:: FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    pp<-Seurat:: CombinePlots(plots = list(plot1, plot2))
    print(pp)
    # ggsave(pp,filename = paste0("1-Mitocondrial-genes-and-feature-count.pdf"),width = 8.5, height = 5.8,path = path1$folder_name)

    pp<- VlnPlot(sce, features = c("percent.ribo", "percent.mt"), ncol = 2)
    print(pp)
    # ggsave(pp,filename = paste0("2-Mitocondrial-and-ribosomal-gene-count.pdf"),width = 8.5, height = 5.8,path = path1$folder_name)

    pp<- VlnPlot(sce, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
    print(pp)
    # ggsave(pp,filename = paste0("3-nFeature_RNA-and-nCount_RNA.pdf"),width = 8.5, height = 5.8,path = path1$folder_name)
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
    message("After data filtering")
    print(sce)
  }

  eset<-as.matrix(GetAssayData(object = sce, assay = "RNA", slot = 'counts'))

  print(eset[1:10,1:3])
  if(id_type=="symbol") id.type<-"S"
  if(id_type=="emsembl") id.type<-"E"
  copykat<-mycopykat(rawmat            = eset,
                    id.type            = id.type,
                    cell.line          ="no",
                    ngene.chr          =5,
                    win.size           =25,
                    KS.cut             =0.15,
                    sam.name           = project,
                    norm.cell.names    = normal_cell_id,
                    distance           ="euclidean",
                    n.cores            =1,
                    save_path          = path1$abspath)
  return(copykat)
}
