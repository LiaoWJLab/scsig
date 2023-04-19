





#' Modified feature expression heatmap
#'
#'The function is a modified version of `DoHeatmap{Seurat}`.
#'
#' @param sce Seurat object
#' @param group A vector of variables to group cells by; pass 'ident' to group by cell identity classes
#' @param feas A vector of features to plot, defaults to `VariableFeatures(object = object)`
#' @param assay Assay to pull from
#' @param slot Data slot to use, choose from 'raw.data', 'data', or 'scale.data'
#' @param disp.min Minimum display value (all values below are clipped)
#' @param disp.max Maximum display value (all values above are clipped); defaults to 2.5 if slot is 'scale.data', 6 otherwise
#' @param cols Vector of colors, each color corresponds to an identity class. This may also be a single character, such as normal and random,  to a palette as specified by `IOBR::palettes`. See `palettes` for details
#' @param palette Numeric value corresponding with colors to use for the color bar. Default is 1, other options: 2, 3, 4
#' @param palette.heatmap Numeric value corresponding with colors  to use for n-colour gradient, default is 1, other options: 2, 3, 4, 5, 6
#' @param path Path of the output saving directory
#' @param show_plot Whether to show the plot
#' @param show_col Whether to show color palettes
#' @param fig.type Format of plot saving, such as pdf and png
#' @param prefix Prefix of the file name
#'
#' @return a list of ggplot objects
#' @export
#'
#' @examples
#' data("pbmc_small")
#' scsig_heatmap(sce = pbmc_small,group= "RNA_snn_res.1")
scsig_heatmap<-function(sce,
                       group,
                       feas            = NULL,
                       assay           = NULL,
                       slot            = "scale.data",
                       disp.min        = -2.5,
                       disp.max        = 5,
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
  pp<- DoHeatmap(object       = sce,
                 features     = show_vars,
                 angle        = 60,
                 size         = 3.5,
                 group.colors = cols,
                 assay        = assay,
                 slot         = slot,
                 disp.min        = disp.min,
                 disp.max        = disp.max)+
    scale_fill_gradientn(colours = mapal)

  if(show_plot) print(pp)
  ###################################
  ggsave(plot =  pp, filename=paste0(prefix,"-DoHeatmap-",group,".",fig.type),
         path = file_name$folder_name,
         width = wwidth, height = hheight)
  return(pp)

}
