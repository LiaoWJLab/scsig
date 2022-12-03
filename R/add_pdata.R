





#' Title
#'
#' @param sce
#' @param pdata
#' @param id_sce
#' @param id_pdata
#'
#' @return
#' @export
#'
#' @examples
add_pdata<-function(sce, pdata, id_sce = "orig.ident", id_pdata = "ID"){

  pdata<-as.data.frame(pdata)
  if(id_pdata!="ID") colnames(pdata)[which(colnames(pdata)==id_pdata)]<-"ID"
  colnames(pdata)<-gsub(colnames(pdata), pattern = " ", replacement = "_")
  colnames(pdata)<-gsub(colnames(pdata), pattern = " ", replacement = "_")
  colnames(pdata)<-gsub(colnames(pdata), pattern = "__", replacement = "_")
  colnames(pdata)<-gsub(colnames(pdata), pattern = "-", replacement = "_")
  colnames(pdata)<-gsub(colnames(pdata), pattern = ",", replacement = "")
  colnames(pdata)

  #################################
  # help(AddMetaData)

  print(head(sce@meta.data))

  message("+++++++++++++++++++++++++++++++++++++++++++ ")

  message(">>>--- Propotion of pdata id in seurat metadata: ")
  print(summary(pdata$ID%in% as.character(sce@meta.data[,id_sce])))

  ###############################
  id<-data.frame( "rowname" = rownames(sce@meta.data), "id" = sce@meta.data[,id_sce] )

  pdata_meta<-merge(id, pdata, by.x = "id", by.y = "ID", all.x = T, all.y = F)
  # head(pdata_meta)
  pdata_meta<-pdata_meta[match(rownames(sce@meta.data), pdata_meta$rowname),]
  rownames(pdata_meta)<-NULL
  pdata_meta<-column_to_rownames(pdata_meta, var = "rowname")
  ################################
  sce<- AddMetaData(sce, pdata_meta)
  # head(sce@meta.data)
  return(sce)

}


