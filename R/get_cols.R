



#' Title
#'
#' @param cols 
#' @param palette 
#' @param show_col 
#' @param seed 
#'
#' @return
#' @export
#'
#' @examples
get_cols<-function(cols = "normal", palette = 1, show_col = T, seed = 123){
  
  ##############################################
  if(length(cols)==1){
    if(cols=="random"){
      
      mycols<-palettes(category = "random", palette = palette, show_col = show_col)
      message(">>>> Default seed is 123, you can change it by `seed`(parameter).")
      set.seed(seed)
      mycols<-mycols[sample(length(mycols), length(mycols))]
      if(show_col) scales::show_col(mycols)
      
    }else if(cols == "normal"){
      
      mycols<-palettes(category = "random", palette = palette, show_col = show_col)
    }else{
      mycols<- palettes(category = "box", palette = cols, show_col = show_col)
    }
  }else{
    mycols<-cols
    if(show_col) scales::show_col(mycols)
  }
  ################################################
  return(mycols)
}
