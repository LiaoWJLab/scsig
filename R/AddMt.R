



#' Title
#'
#' @param sce
#' @param species
#'
#' @return
#' @export
#'
#' @examples
AddMt<-function(sce, species='human'){

  message("species must be one of `human` or `mouse`.")

  if(species=='human'){

    sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
    sce[["percent.rp"]] <- PercentageFeatureSet(sce, pattern = "^RPL|^RPS")
  }
  if(species=='mouse'){
    sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^mt-")
    sce[["percent.rp"]] <- PercentageFeatureSet(sce, pattern = "^Rpl|^Rps")
  }
  return(sce)
}
