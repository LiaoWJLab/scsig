





#' randomly strata cells
#'
#' @param input cell annotation with cell id
#' @param group annotation of cell type
#' @param propotion range between 0 to 1
#' @param minum_count minimum count of each cell type
#' @param minum_count_include
#' @param max_count
#' @param sub_cluster
#' @param cell_type
#'
#' @return
#' @export
#'
#' @examples
random_strata_cells<-function(input, group, propotion = 0.1, minum_count_include = 300, minum_count = 200, max_count = 1000, sub_cluster = NULL, cell_type = NULL){


  if(class(input)=="Seurat") input<- input@meta.data

  # group <- match.arg(group)
  input<-input %>%
    dplyr::filter(!is.na(group)) %>%
    dplyr::filter(group!="Undetermined") %>%
    dplyr::filter(group!="NA")%>%
    dplyr::filter(group!=" ")


  if(!is.null(sub_cluster)){
    input<-input[input[,sub_cluster]==cell_type, ]
  }

  #choose cell types with count larger than minum_count_include
  cell_freq<-as.data.frame(table(input[, group]))
  message(">>> Count of each cell type:")
  print(cell_freq)
  cell_id<-as.character(cell_freq[cell_freq$Freq>minum_count_include,"Var1"])
  message(">>> Cell types after sampling: ")

  if("Undetermined"%in%cell_id) cell_id<-cell_id[-c(which(cell_id=="Undetermined"))]

  print(cell_id)

  index<-input[,group]%in%cell_id
  input<-input[index,]
  index<-order(input[,group],decreasing = F)
  input<-input[index, ]

  input2<-input



  cell_freq<-as.data.frame(table(input2[, group]))
  if(minum_count > min(cell_freq$Freq)) stop(">>>-- `minum_count` must smaller than minimum count of data for sampling...")
  #####################################


  set.seed(1234)
  n<-round(table(input[,group])*propotion)
  print(n)
  size<-n
  names(size)<-NULL
  x <- sampling::strata(input, group, size = size, method = "srswor")

  input<-sampling:: getdata(input,x)
  input<-as.data.frame(input2[rownames(input2)%in%rownames(input), ])

  cell_freq<-as.data.frame(table(input[, group]))
  print(cell_freq)


  ###########################################

  #cell count
  cell_freq<-as.data.frame(table(input[, group]))
  cell_id_min<-as.character(cell_freq[cell_freq$Freq< minum_count, "Var1"])
  cell_id_max<-as.character(cell_freq[cell_freq$Freq> max_count, "Var1"])

  if(!is.null(cell_id_min)&length(cell_id_min)>0){

    message(paste0(">>>-- Adjustment for celltypes with counts smaller than ", minum_count, " : "))
    print(cell_id_min)

    for (dd in 1:length(cell_id_min)) {

      cell_loop<-cell_id_min[dd]
      input_loop<-input2[input2[,group]==cell_loop,]
      # print(dd)
      input_loop<-input_loop[sample(rownames(input_loop), minum_count), ]
      # print(dd+1)
      input<-input[!input[,group]==cell_loop,]
      input<-rbind(input, input_loop)
    }

  }



  if(!is.null(cell_id_max)&length(cell_id_max)>0){

    cell_id_max<-cell_id_max[!is.na(cell_id_max)]
    message(paste0(">>>-- Adjustment for celltypes with counts larger than ", max_count, " : "))
    print(cell_id_max)
    for (dd in 1:length(cell_id_max)) {

      cell_loop<-cell_id_max[dd]
      input_loop<-input2[input2[,group]==cell_loop,]
      # print(dd)
      input_loop<-input_loop[sample(rownames(input_loop), max_count), ]
      input<-input[!input[,group]==cell_loop,]
      input<-rbind(input, input_loop)
    }

  }

  message(">>> Cell types after adjustment: ")
  cell_freq<-as.data.frame(table(input[, group]))
  print(cell_freq)

  return(input)
}
