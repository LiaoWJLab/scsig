



#' enrich_cluster_degs
#'
#' @param sce.markers result of find_markers
#' @param cutoff_foldchange default is 0.5
#' @param cutoff_p_adj "enrichGO", "enrichKEGG"
#' @param methods
#' @param path
#' @param showCategory number of terms to show
#' @param cols
#' @param cluster
#' @param gene
#' @param p.adj
#' @param logfc
#' @param width
#' @param height
#' @param index default is NULL
#'
#' @return
#' @export
#'
#' @examples
enrich_cluster_degs<-function(sce.markers,
                              cutoff_foldchange = 0.5,
                              cutoff_p_adj      = 0.01,
                              methods           = c("enrichGO", "enrichKEGG"),
                              cluster           = "cluster",
                              gene              = "gene",
                              p.adj             = "p_val_adj",
                              logfc             = "avg_log2FC",
                              showCategory      = 12,
                              path              = NULL,
                              cols              = c('#b3eebe', "#46bac2", '#371ea3'),
                              index             = NULL,
                              width             = 9.36,
                              height            = 14.3){

  if(!is.null(path)){
    path<-creat_folder(path)
  }
  #########################################
  colnames(sce.markers)[which(colnames(sce.markers)== cluster )] <- "cluster"
  colnames(sce.markers)[which(colnames(sce.markers)== gene )] <- "gene"
  colnames(sce.markers)[which(colnames(sce.markers)== p.adj )] <- "p_val_adj"
  colnames(sce.markers)[which(colnames(sce.markers)== logfc )] <- "avg_log2FC"

  #########################################
  top_degs <- sce.markers %>% dplyr:: filter( .$p_val_adj <= cutoff_p_adj) %>%
    dplyr:: filter( .$avg_log2FC >= cutoff_foldchange) %>%
    dplyr:: group_by(cluster) %>%
    dplyr::select(cluster, gene)

  deg_list<-split(as.data.frame(top_degs), top_degs$cluster)
  deg_list<-lapply(deg_list, function(x) x[,"gene"])
  ####################################

  library(org.Hs.eg.db)
  for (y in 1:length(deg_list)) {
    #i=1
    ## 把SYMBOL改为ENTREZID
    deg_list[[y]]=as.character(na.omit(AnnotationDbi::select(org.Hs.eg.db,
                                                             keys = deg_list[[y]],
                                                             columns = 'ENTREZID',
                                                             keytype = 'SYMBOL')[,2]))}
  #####################################

  for (z in 1:length(methods)) {

    method<-methods[z]
    message(paste0(">>>>--- Processsing terms: ", method))

    if(method=="enrichGO"){
      xx <-clusterProfiler::compareCluster(deg_list,
                                           fun= method,
                                           OrgDb = org.Hs.eg.db,
                                           pvalueCutoff=0.99) # pvalueCutoff

      xx<- clusterProfiler::setReadable(xx, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
      # table(xx@compareClusterResult$Cluster)
      print(head(as.data.frame(xx)))

    }

    if(method=="enrichKEGG"){
      xx <-clusterProfiler::compareCluster(deg_list,
                                           fun= method,
                                           pvalueCutoff=0.99
                                           # OrgDb = org.Hs.eg.db,
                                           # count = 3,
      ) # pvalueCutoff
      xx<-clusterProfiler:: setReadable(xx, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
      # table(xx@compareClusterResult$Cluster)
      print(head(as.data.frame(xx)))

    }

    if(!is.null(index)) {
      writexl::write_xlsx(as.data.frame(xx), paste0(path$abspath, index,"-",z,"-",method, "-enrichment-result.xlsx"))
    }else{
      writexl::write_xlsx(as.data.frame(xx), paste0(path$abspath, z,"-",method, "-enrichment-result.xlsx"))
    }
   ##############################################################


    gg<-clusterProfiler:: dotplot(xx, showCategory = showCategory )+
      scale_color_gradientn(colours= cols,
                            guide=guide_colorbar(reverse=TRUE))+
      scale_y_discrete(name ="Description",labels=function(x) stringr:: str_wrap(x, width=40))+
      scale_size_continuous(range=c(2, 10)) +
      design_mytheme(legend.position = "right", legend.direction = "vertical", legend.size =0.25 )+
      labs(x="",
           y="",
           title = method)

    print(gg)
    if(!is.null(index)) {
      ggsave(gg, filename = paste0(index,"-",z,"-",method, "-enrichment-doplot.pdf"), width = width, height = height, path = path$folder_name )
    }else{
      ggsave(gg, filename = paste0(z,"-",method, "-enrichment-doplot.pdf"), width = width, height = height, path = path$folder_name )
    }


  }


}

