



#' Compare and show gene clusters functional profile for all identity classes
#'
#' This function combines `compareCluster()` and `dotplot()` both from package [clusterProfiler](http://127.0.0.1:60491/help/library/clusterProfiler/html/00Index.html)
<<<<<<< Updated upstream
#' to perform enrichment analysis and draw dotplot  for all identity classes.
#'
#' @param sce.markers A data frame of Gene expression markers for all identity classes(output from `FindAllMarkers`)
#' @param cutoff_foldchangeLimit Cutoff for filtering gene signatures which show on average, less than X-fold difference (log-scale) between the two groups of cells for enrichment analysis, default is 0.5
#' @param cutoff_p_adj Cutoff for filtering gene signatures which adjust P values are < n for enrichment  analysis, default is 0.01
#' @param methods A list of methods, choose from "groupGO", "enrichGO", "enrichKEGG", "enrichDO" or "enrichPathway"
#' @param path Path of the output saving directory
#' @param showCategory Number of terms to show, default is 5
#' @param cols Vector of colors to use for n-color gradient, default is c('#b3eebe', "#46bac2", '#371ea3'), users can define the colors manualy
#' @param cluster Name of the column where the clusters in the data frame, default is "cluster"
#' @param gene Name of the column where the genes in the data frame in the data framer, default is "gene"
#' @param p.adj Name of the column where P adjust value in the data framer, default is "p_val_adj"
#' @param logfc Name of the column where log-scaled fold change in the data frame, default is "avg_log2FC"
#' @param width Width of plot for saving
#' @param height Height of plot for saving
#' @param index Prefix of file name when saving
=======
#' to perform enrichment analysis and draw dot plot  for all identity classes.
>>>>>>> Stashed changes
#'
#' @param sce.markers A data frame of Gene expression markers for all identity classes(output from `FindAllMarkers`)
#' @param cutoff_foldchangeLimit Cutoff for filtering gene signatures which show on average, less than X-fold difference (log-scale) between the two groups of cells for enrichment analysis, default is 0.5
#' @param cutoff_p_adj Cutoff for filtering gene signatures which adjust P values are < n for enrichment  analysis, default is 0.01
#' @param methods A list of methods, choose from "groupGO", "enrichGO", "enrichKEGG", "enrichDO" or "enrichPathway"
#' @param path Path of the output saving directory
#' @param showCategory Number of terms to show, default is 5
#' @param cols Vector of colors to use for n-color gradient, default is c('#b3eebe', "#46bac2", '#371ea3'), users can define the colors manualy
#' @param cluster Name of the column where the clusters in the data frame, default is "cluster"
#' @param gene Name of the column where the genes in the data frame in the data framer, default is "gene"
#' @param p.adj Name of the column where P adjust value in the data framer, default is "p_val_adj"
#' @param logfc Name of the column where log-scaled fold change in the data frame, default is "avg_log2FC"
#' @param width Width of plot for saving
#' @param height Height of plot for saving
#' @param index Prefix of file name when saving
#'
#' @author Dongqiang Zeng
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
                              showCategory      = 5,
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
    # i=1
    # change SYMBOL ID into ENTREZID
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
                                           pvalueCutoff=0.99)
      # table(xx@compareClusterResult$Cluster)
      print(head(as_tibble(xx)))

    }

    if(method=="enrichKEGG"){
      xx <-clusterProfiler::compareCluster(deg_list,
                                           fun= method,
                                           pvalueCutoff=0.99
                                           # OrgDb = org.Hs.eg.db,
                                           # count = 3,
      )

      print(head(as_tibble(xx)))

    }

    if(!is.null(index)) {
      writexl::write_xlsx(as.data.frame(xx), paste0(path$abspath, index,"-",z,"-",method, "-enrichment-result.xlsx"))
    }else{
      writexl::write_xlsx(as.data.frame(xx), paste0(path$abspath, z,"-",method, "-enrichment-result.xlsx"))
    }
    ##############################################################
    gg<-clusterProfiler:: dotplot(xx, showCategory = showCategory,font.size =2.5)+
      scale_color_gradientn(colours= cols,
                            guide=guide_colorbar(reverse=TRUE))+
      scale_y_discrete(name ="Description",labels=function(x) stringr::str_wrap(x, width=40))+
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

