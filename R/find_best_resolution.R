




#' Find best resolution for clustering
#'
#' Find best resolution for clustering. First, for a given seurat object, identify clusters of cells by setting the resolution value from low to high.
#' Then create a plot of a clustering tree showing the relationship between clusters at different resolutions.
#' @param sce Suerat object
#' @param from Minimum resolution value for testing
#' @param to Maximum resolution value for testing
#' @param save_plot Whether to save plot
#' @param path Path of the output saving directory
#' @param palette If a string, will use that named palette. If a number, will index into the list of palettes of appropriate type. The list of available palettes can found in the Palettes section
#' @param index Prefix of file name for saving
#' @param alpha Either a numeric value giving the alpha of all nodes or the name of a metadata column to use for node transparency
#' @param verbose Whether to print output
#' @param return_sce Whether to return suerat object
#' @param graph.name Name of graph to use for the clustering algorithm
#' @param assay Assay to pull from, e.g. RNA, SCT, integrated
#' @param prefix String indicating columns containing clustering information
#'
#' @return Suerat object or a ggplot object
#' @export
#'
#' @examples
find_best_resolution<-function(sce, assay,  prefix, from = 0.2, to = 1.2, graph.name = NULL,  save_plot = T, path = NULL,
                               palette = "set2", index = NULL, alpha = 0.8, verbose = TRUE, return_sce = TRUE) {
  
  
  # creat direaction
  if(!is.null(path)){
    path<-creat_folder(path)
  }else{
    path<-"./"
  }
  ###########################################################
  DefaultAssay(sce)<- assay
  for (xx in seq(from, to, by = 0.2)) {
    message(paste0(">>>--- Resolution is ", xx))
    sce <- FindClusters(sce, resolution = xx, verbose = verbose, graph.name = graph.name)
    
    message(paste0(">>>-- Findclusters when resolution = "), xx, " :")
    print(summary(as.factor(Idents(sce))))
    message("------------------------------------------------------")
  }
  
  p<-clustree:: clustree(sce, alpha = alpha, prefix = prefix) +
    theme(legend.position = "bottom") +
    scale_color_brewer(palette = "Set2") +
    # scale_color_manual(values = palettes(category = "box", palette = palette)) +
    scale_edge_color_continuous(low = "grey80", high = "red")
  
  if(save_plot){
    if(is.null(index)) index<- 0
    ggsave(p, filename = paste0(index,"-finding-best-resolutions.pdf"),
           width = 10, height = 15, path = path$folder_name )
  }
  
  print(p)
  if(return_sce){
    return(sce)
  }else{
    return(p)
  }
  
}
