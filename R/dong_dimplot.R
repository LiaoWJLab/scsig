



#' Modified Dimplot
#'
#' @param sce Seurat object
#' @param reduction Which dimensionality reduction to use. If not specified, first searches for umap, then tsne, then pca
#' @param groups Name of one or more metadata columns to group (color) cells by (for example, orig.ident); pass 'ident' to group by identity class
#' @param split.by Name of a metadata column to split plot by
#' @param label logical.Whether to label the clusters
#' @param label.size size of label
#' @param pt.size size of point
#' @param cols selection of colors for plots, such as normal, random, default is normal. users can define the cols manualy
#' @param seed seed of the random number generator, defuat is 123
#' @param show_col logical. if TRUE, palettes will be displayed
#' @param width width of figure
#' @param height height of figure
#' @param w_index the parameter of width
#' @param w_add
#' @param max_category maximal number of groups
#' @param show_plot logical. if TRUE, plots will be displayed
#' @param path path of the output saving directory, default is null
#' @param index index number of folder name
#' @param legend.position position of legend
#' @param legend.direction direction of legend
#' @param legend.size font size of legend
#' @param plot_title_size title size
#' @param axis_title_size axis title size
#' @param axis_text_size axis text size
#' @param axis_angle axis angle
#' @param hjust horizontal justification
#' @param theme theme of figure, default is "classic"
#' @param save_plot logical. if TRUE, plots will be saved
#' @param palette color palette, default is 1, other options: 2,3,4
#' @param fig.type figure type, such as pdf and png
#'
#' @return
#' @export
#'
#' @examples
dong_dimplot<-function(sce,
                       reduction        = "umap",
                       groups           = "orig.ident",
                       split.by         = NULL,
                       label            = T,
                       label.size       = 5,
                       pt.size          = 0.5,
                       cols             = "normal",
                       seed             = 123,
                       palette          = 1,
                       show_col         = FALSE,
                       width            = 9.2,
                       height           = 8.2,
                       w_index          = 7,
                       w_add            = 2,
                       max_category     = 18,
                       show_plot        = T,
                       path             = NULL,
                       index            = 1,
                       theme            = "classic",
                       legend.position  = "right",
                       legend.direction = "horizontal",
                       legend.size      = 0.25,
                       plot_title_size  = 2,
                       axis_title_size  = 2,
                       axis_text_size   = 10,
                       axis_angle       = 0,
                       hjust            = 0.5,
                       save_plot        = TRUE,
                       fig.type         = "pdf"){


  if(!is.null(path)){
    file_store<-path
  }else{
    file_store<-paste0("Dimplot-dong")
  }

  if(!file.exists(file_store)) dir.create(file_store)
  abspath<-paste(getwd(),"/",file_store,"/",sep ="" )


  num_cols<-length(unique(as.character(sce@meta.data[,groups])))
  ##############################################
  if(length(cols)==1){
    if(cols=="random"){

      mycols<-palettes(category = "random", palette = palette, show_col = show_col)

      if(length(mycols)<num_cols){
        mycols<-c(mycols, palettes(category = "random", palette = 3, show_col = show_col))
      }

      mycols<-mycols[1:num_cols]
      message(">>>> Default seed is 123, you can change it by `seed`(parameter).")
      set.seed(seed)
      mycols<-mycols[sample(length(mycols), length(mycols))]
     ##########################################

      if(show_col) scales::show_col(mycols)

    }else if(cols == "normal"){

      mycols<-palettes(category = "random", palette = palette, show_col = show_col)
      if(length(mycols)<num_cols){
        mycols<-c(mycols, palettes(category = "random", palette = 3, show_col = show_col))
      }

      mycols<-mycols[1:num_cols]
    }
  }else{
    mycols<-cols
    if(show_col) scales::show_col(mycols)
  }
  ################################################


  for (i in 1:length(groups)) {

    group1<-groups[i]
    message(paste0(">>> Processing group:: ", group1))

    if(!group1%in%colnames(sce@meta.data)) stop(">>>>> groups not exist in meta.data of seurat object!")

    pp <- Seurat:: DimPlot(sce,
                           reduction  = reduction,
                           group.by   = group1,
                           label      = label,
                           split.by   = split.by,
                           label.size = label.size,
                           pt.size    = pt.size,
                           cols       = mycols)

    if(dim(sce)[2]>30000){
      limitsize<-TRUE
    }else{
      limitsize<-FALSE
    }

    if(length(unique(sce@meta.data[,group1]))> max_category & legend.position!="bottom"){

      w_add<-round(length(unique(sce@meta.data[,group1]))/max_category, 0) * w_add
      width<- width + w_add

    }

    if(!is.null(split.by)){
      width<-length(unique(sce@meta.data[,split.by]))*w_index
      fig.name<-paste0("-split-by-",split.by)

      if(length(unique(sce@meta.data[,split.by]))>3){
        limitsize<-FALSE
      }else{
        limitsize<-TRUE
      }
    }else{
      fig.name<-NULL
    }

    if(is.null(index)) index<-1
    #####################################

    mytheme<-design_mytheme(
      theme            = theme,
      plot_title_size  = plot_title_size,
      axis_title_size  = axis_title_size,
      axis_text_size   = axis_text_size,
      axis_angle       = axis_angle,
      hjust            = hjust,
      legend.position  = legend.position,
      legend.direction = legend.direction,
      legend.size      = legend.size
    )
    pp<-pp+mytheme #+guides(color = guide_legend(override.aes = list(size = legend.size)), size = "none")

    if(show_plot) print(pp)
    ####################################


    if(legend.position=="bottom") height<- height + 6

    if(save_plot){
      ggsave(pp,
             filename = paste0(index,"-",i,"-",reduction,"-",group1, fig.name,'.',fig.type),
             path = file_store,
             width = width,
             height = height,
             dpi = 300,
             limitsize = limitsize)

      message(paste0(">>>>>> Figure name is:: ", paste0(index,"-",i,"-",reduction,"-",group1, fig.name,'.',fig.type)))

    }

  }

  return(pp)
}



