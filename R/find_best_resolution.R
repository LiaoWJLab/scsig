




#' Title
#'
#' @param sce
#' @param from
#' @param to
#' @param save_plot
#' @param path
#' @param palette
#' @param index
#' @param alpha
#' @param verbose
#' @param return_sce
#' @param graph.name
#' @param assay
#' @param prefix
#'
#' @return
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
