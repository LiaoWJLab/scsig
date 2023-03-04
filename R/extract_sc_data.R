







#' Title
#'
#' @param sce
#' @param vars
#' @param assay
#' @param combine_meta_data
#'
#' @return
#' @export
#'
#' @examples
extract_sc_data<-function(sce, vars = NULL, assay, slot = "scale.data", combine_meta_data = TRUE){


  exist<-Seurat::Assays(sce)
  message(paste0(">>>--- Assays of seurat object: "))
  message(paste(">>>---", exist, collapse = " "))
  assay<-assay[assay%in%exist]

  if(length(assay)==0) stop(">>>--- There no assay in object...")

  eset_cbind<-data.frame("ID" = rownames(sce@meta.data), "index" = 1:nrow(sce@meta.data))

  for(i in 1:length(assay)){

    method<-assay[i]
    DefaultAssay(sce)<-method

    eset<- SeuratObject:: GetAssayData(sce, assay = assay, slot = slot)

    feas_e<- rownames(eset)[rownames(eset)%in%unique(vars)]
    eset<- eset[feas_e, ]
    # print(head(eset))
    if(length(feas_e)==1){
     eset<- data.frame("ID" = as.character(names(eset)), vars = as.numeric(eset))
     colnames(eset)[2] <- feas_e
    }else{
      eset<- as.data.frame(t(eset))
      eset<- tibble:: rownames_to_column(eset, var = "ID")
    }
      # base::as.data.frame() %>%
      # tibble:: rownames_to_column(.,var = "id") %>%
      # dplyr:: filter(.$id%in%unique(vars)) %>%
      # column_to_rownames(.,var = "id") %>%
      # base:: t() %>%
      # base:: as.data.frame() %>%

    if(length(feas_e)==0) next

    # print(head(eset))

    if(length(assay)>1) colnames(eset)[2:ncol(eset)] <-paste0(colnames(eset)[2:ncol(eset)], "_", method)

    # colnames(eset)<-gsub(colnames(eset), pattern = "-", replacement =  "_")
    # colnames(eset)<-gsub(colnames(eset), pattern = "/", replacement = "_")
    # colnames(eset)<-gsub(colnames(eset), pattern = " ", replacement = "_")
    eset<-as.data.frame(eset)

    if(length(assay)==1){
      eset_cbind<-eset
    }else{
      eset_cbind<-inner_join(eset_cbind, eset, by ="ID")
    }

  }

  if(combine_meta_data){
    meta.data<-rownames_to_column(sce@meta.data, var = "ID")
    eset_cbind<- inner_join(meta.data, eset_cbind, by = "ID" )
  }

  # eset_cbind<-eset_cbind[,-which(colnames(eset_cbind)=="index")]
  # print(head(eset_cbind))
  return(eset_cbind)

}
