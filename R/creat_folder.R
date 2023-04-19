



#'  Construct path to a file or directory
#'
#' @param f1 A character vector containing primary directory of path name
#' @param f2 A character vector containing second-level directory of path name. Default is NULL.
#' @param return Number of levels of the directory, choose from 1, 2 and 3. Default is NULL.
#' @param f3 A character vector containing third-level directory of path name. Default is NULL
#'
#' @return A list containing folder name and absolute path from the current working directory.
#' @export
#'
#' @examples
#' ##Construct a file path containing one level directory
#' path<-creat_folder("1")
#' ##Construct a file path containing three level directory
#' path<-creat_folder(f1="1",f2="2",f3="3")
#' ##Construct a file path containing two level directory
#' path<-creat_folder(f1="1",f2="2",f3="3", return=2)
creat_folder<-function(f1, f2 = NULL, f3 = NULL, return = NULL){


  if(!is.null(f3)){
    path<-file.path(getwd(), f1, f2, f3)
    if(!dir.exists(path)) dir.create(file.path(getwd(), f1, f2, f3), recursive = TRUE)
  }else if(!is.null(f2)){
    path<-file.path(getwd(), f1, f2)
    if(!dir.exists(path)) dir.create(file.path(getwd(), f1, f2), recursive = TRUE)
  }else{
    path<-file.path(getwd(), f1)
    if(!dir.exists(path)) dir.create(file.path(getwd(), f1), recursive = TRUE)
  }


  if(is.null(return)){

    if(!is.null(f3)){
      res<-list("folder_name" = paste0(f1, "/", f2, "/", f3),
                "abspath" = paste0(file.path(getwd(), f1, f2, f3), "/"))
    }else if(!is.null(f2)){
      res<-list("folder_name" = paste0(f1, "/", f2),
                "abspath" = paste0(file.path(getwd(), f1, f2), "/"))
    }else{
      res<-list("folder_name" = f1,
                "abspath" = paste0(file.path(getwd(), f1), "/"))
    }

  }else{

    if(return == 1){
      res<-list("folder_name" = f1,
                "abspath" = paste0(file.path(getwd(), f1), "/"))
    }else if(return == 2){
      res<-list("folder_name" = paste0(f1, "/", f2),
                "abspath" = paste0(file.path(getwd(), f1, f2), "/"))
    }else if(return == 3){
      res<-list("folder_name" = paste0(f1, "/", f2, "/", f3),
                "abspath" = paste0(file.path(getwd(), f1, f2, f3), "/"))
    }

  }

  return(res)
}
