





#' Title
#'
#' @param sce 
#' @param group 
#' @param feas 
#' @param assay 
#' @param slot 
#' @param cols 
#' @param palette 
#' @param palette.heatmap 
#' @param path 
#' @param show_plot 
#' @param show_col 
#' @param fig.type 
#' @param prefix 
#'
#' @return
#' @export
#'
#' @examples
dong_heatmap<-function(sce, 
                       group, 
                       feas            = NULL,
                       assay           = NULL, 
                       slot            = "scale.data", 
                       cols            = "normal",
                       palette         = 1,
                       palette.heatmap = 1,
                       path            = NULL,
                       prefix          = 0,
                       show_plot       = TRUE,
                       show_col        = FALSE,
                       fig.type        = "pdf"){
  
  
  
  if(length(unique(sce@meta.data[, group]))==1) stop("Cell identity has only one level...")
  
  #define-path
  if(is.null(path)){
    file_name<-creat_folder(paste0("0-",group,"-DoHeatmap-plot"))
  }else{
    file_name<-creat_folder(paste0(path))
  }
  #############################################
  cols<- get_cols(cols = cols, palette = palette, show_col = show_col)
  
  show_vars<-feas
  Idents(sce)<- group
  ##################################
  wwidth<- 4 + length(unique(sce@meta.data[,group]))*0.4
  hheight<- 3 + length(unique(show_vars)) * 0.3
  ###################################
  
  # mapal <-rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256))
  mapal<- palettes(category = "heatmap", palette = palette.heatmap, show_col = show_col)
  ###################################
  pp<- DoHeatmap(object       =  sce,
                 features     = show_vars,
                 angle        = 60,
                 size         = 3.5,
                 group.colors = cols,
                 assay        = assay,
                 slot         = slot)+
    scale_fill_gradientn(colours = mapal)
  
  if(show_plot) print(pp)
  ###################################
  ggsave(plot =  pp, filename=paste0(prefix,"-DoHeatmap-",group,".",fig.type),
         path = file_name$folder_name,
         width = wwidth, height = hheight)
  return(pp)

}
