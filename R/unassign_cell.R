








#' Title
#'
#' @param sce
#' @param cluster
#' @param ignore_cell_prefix
#' @param min_cell_count
#' @param new_col
#' @param delete_unassigned
#' @param return_meta_data
#' @param alread_unassigned
#' @param minimal_cell
#'
#' @return
#' @export
#'
#' @examples
unassign_cell<-function(sce, cluster, ignore_cell_prefix = NULL, min_cell_count = 30,
                        new_col = NULL, delete_unassigned = FALSE, return_meta_data = FALSE,
                        alread_unassigned = "unassigned", minimal_cell = 10){

  if(class(sce)!="Seurat") stop("sce must be a seurat object...")

  message(">>>--- Counts of celltypes before filtering: ")
  clu_com<-as.data.frame(table(as.character(sce@meta.data[,cluster])))
  print(table(as.character(sce@meta.data[,cluster])))
  message( paste0(">>>--- Counts of celltypes after filtering step-1:  Removing minimal cell counts: ", minimal_cell))
  clu_com<-clu_com[clu_com$Freq>=minimal_cell, ]
  # print(clu_com)
  clu_certein<-as.character(clu_com$Var1)
  if(is.null(new_col)){
    Idents(sce)<-cluster
    obj<-as.character(Idents(sce))
    sce@meta.data[,cluster] <- ifelse(!obj%in%clu_certein, alread_unassigned, obj)
  }else{

    sce@meta.data[, new_col]<-sce@meta.data[,cluster]
    Idents(sce)<- new_col
    obj<-as.character(Idents(sce))
    sce@meta.data[, new_col] <- ifelse(!obj%in%clu_certein, alread_unassigned, obj)
  }

  Idents(sce)<-cluster
  if(alread_unassigned%in%unique(as.character(Idents(sce)))){
    sce<-subset(sce, idents = alread_unassigned, invert = TRUE)
  }

  print(table(as.character(sce@meta.data[,cluster])))
  ###########################################################
  if(dim(clu_com)[1] >=2) {

    if(!is.null(ignore_cell_prefix)){

      for (i in 1:length(ignore_cell_prefix)) {

        cell_prefix<- ignore_cell_prefix[i]

        clu_com[grep(clu_com$Var1, pattern = cell_prefix), "Freq"] <- min_cell_count+100
      }

    }

    #delete-already-unassigned-cells
    if(alread_unassigned%in%clu_com$Var1) clu_com[clu_com$Var1==alread_unassigned, "Freq"] <- min_cell_count-5

    clu_certein<- as.character(clu_com[clu_com$Freq > min_cell_count,]$Var1)

    message(paste0(">>>--- Choose clusters with cell counts more than ", min_cell_count, " : "))
    message(" ")
    print(clu_certein)

    # obj<-as.character(sce@meta.data[,cluster])
    # print(table(obj))

    if(is.null(new_col)){
      Idents(sce)<-cluster
      obj<-as.character(Idents(sce))
      sce@meta.data[,cluster] <- ifelse(!obj%in%clu_certein, alread_unassigned, obj)
    }else{
      sce@meta.data[, new_col]<-sce@meta.data[,cluster]
      Idents(sce)<- new_col
      obj<-as.character(Idents(sce))
      sce@meta.data[, new_col] <- ifelse(!obj%in%clu_certein, alread_unassigned, obj)
    }
    # print(table(Idents(sce)))
    message(">>>--- Counts of celltypes after filtering step-2 : ")

  }
  ############################################################

  if(return_meta_data){
    input<-sce@meta.data

    if(delete_unassigned){
      if(is.null(new_col)) {
        index<-as.character(input[,cluster])%in%alread_unassigned
      }else{
        index<-as.character(input[,new_col])%in%alread_unassigned
      }
      input <- input[!index, ]
    }

    if(is.null(new_col)){
      print(table(as.character(input[, cluster])))
    }else{
      print(table(as.character(input[, new_col])))
    }

    return(input)

  }else{

    if(is.null(new_col)){
      Idents(sce)<-cluster
    }else{
      Idents(sce)<- new_col
    }
    #delete is true and unassign exit in clusters
    if(delete_unassigned & alread_unassigned%in%unique(Idents(sce))) {
      sce<-subset(sce, idents = alread_unassigned, invert = TRUE)
    }

    print(table(as.character(Idents(sce))))
    return(sce)
  }

}
