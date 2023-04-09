





#' single cell annotation by scPred
#'
#'  Cell type annotation using a combination of unbiased feature selection from a reduced-dimension space, and machine-learning probability-based prediction method.
#'  Refer to [scPred](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1862-5)
#' @param sce A seurat object containing cells to be classified
#' @param sce_ref A Seurat object with trained model(s) using scPred or a scPred object
#' @param name A character used to name a new column of metadata containing cell type annotation by `scPred`
#' @param mini_cluster Minimal cell counts of the cluster
#' @param threshold Threshold used for probabilities to classify cells into classes. All cells below this threshold value will be labels as "unassigned". In the case of binary classification (two cell tyoes), a threshold of 0.5 will force all cells to be classified to any of the two cell types. For multi-class classification, if there's no probability higher than the threshold associated to a cell type, this will be labelled as "unassigned"
#' @param merge_seurat_cluster Whether to merge prediction results by`scPred` with clustering results
#' @param name_merge_seurat A character used to name a new column of metadata containing the results of `merge_seurat_cluster`
#' @param source Character string related to computer system name choose from win and linux
#' @param path_ref The path or the name of the file where the `sce_ref` is read from
#'
#' @return Seurat object with updated metadata containing cell type annotation by `scPred`
#' @export
#'
#' @examples
dong_cell_anno<-function(sce, sce_ref = NULL, name = "Model1",merge_seurat_cluster = T, name_merge_seurat = "scpred_seurat", mini_cluster = 50, threshold = 0.65, source = "win", path_ref = NULL){


  if(is.null(sce_ref)){

    message(">>>--- Reading reference data from ")

    if(!is.null(path_ref)){
      sce_ref<-readRDS(path_ref)
      message(paste0("  ", path_ref))
    }else{
      if(source=="win"){
        sce_ref<-base:: readRDS("E:/03-NSCLC/13-NSCLC-scRNA-Immunotherapy/4-Data-analyses/0-cell-annotation-model/4-Cell_type.refined-cell-annotation-model.rds")
        message("    E:/03-NSCLC/13-NSCLC-scRNA-Immunotherapy/4-Data-analyses/0-cell-annotation-model/4-Cell_type.refined-cell-annotation-model.rds ")
      }
      if(source =="linux"){
        sce_ref<-base:: readRDS(paste0("/hdd/nsclc/IO/1-data-analysis/0-Cell-annotation-modles/","4-Cell_type.refined-cell-annotation-model.rds"))
        message("  /hdd/nsclc/IO/1-data-analysis/0-Cell-annotation-modles/","4-Cell_type.refined-cell-annotation-model.rds ")
      }
    }

  }

  sce <-scPred:: scPredict(sce, sce_ref, recompute_alignment = FALSE, threshold = threshold)

  print(head(sce@meta.data))
  ################################################################
  colnames(sce@meta.data)<-gsub(colnames(sce@meta.data), pattern = "scpred", replacement = name)

  prefix<-colnames(sce@meta.data)[str_detect(colnames(sce@meta.data),name)]

  prefix_contain<-c(paste0(name,"_prediction"), paste0(name,"_no_rejection"))

  sce@meta.data<-sce@meta.data[!colnames(sce@meta.data)%in%setdiff(prefix, prefix_contain)]

  sce@meta.data[,paste0(name,"_prediction")] <-gsub(sce@meta.data[,paste0(name,"_prediction")], pattern = "\\/", replacement = "-")
  sce@meta.data[,paste0(name,"_no_rejection")] <-gsub(sce@meta.data[,paste0(name,"_no_rejection")], pattern = "\\/", replacement = "-")


  print(table(sce@meta.data[,paste0(name,"_no_rejection")]))
  ###################################################################


  if(merge_seurat_cluster){


    if(name_merge_seurat%in%colnames(sce@meta.data)) stop(paste0(name_merge_seurat, " already exist..."))

    if("seurat_clusters"%in%colnames(sce@meta.data)){

      # print(paste0("The Identity of seurat object is: ", Idents(sce)))

      sce@meta.data[,name_merge_seurat]<-paste0(sce@meta.data[,paste0(name,"_no_rejection")]," ", sce$seurat_clusters)

      clu_com<-as.data.frame(table(sce@meta.data[,name_merge_seurat]))
      clu_certein<-as.character(clu_com[clu_com$Freq>= mini_cluster, ]$Var1)

      sce@meta.data[,name_merge_seurat]<-ifelse(!sce@meta.data[,name_merge_seurat]%in%clu_certein, "unassigned", sce@meta.data[,name_merge_seurat])

      print(table(sce@meta.data[,name_merge_seurat]))

      message(">>>--- Final head of seurat meta.data: ")
      print(head(sce@meta.data))

    }else{
      stop("seurat_clusters was not exixt in meta.data.... Please performe TSNE process...")
    }

  }

  return(sce)
}
