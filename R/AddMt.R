



#' Calculate the percentage of mitochondrial and ribosomal genes
#'
#' @param sce Seurat object
#' @param species Species name. Species must be one of human or mouse
#'
#' @return A Seurat object with the percentage of mitochondrial and ribosomal genes stored in metadata.
#' @export
#'
#' @examples
#' data(pbmc_small)
#' pbmc_small <- AddMt(pbmc_small , species= "human" )
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
