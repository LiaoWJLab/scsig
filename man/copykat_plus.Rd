% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/copykat_plus.R
\name{copykat_plus}
\alias{copykat_plus}
\title{modified_copykat}
\usage{
copykat_plus(
  sce,
  eset = NULL,
  file_type = "10X",
  data_path = NULL,
  project = NULL,
  index = 1,
  id_type = "symbol",
  normal_cell_id = "",
  check_data = TRUE,
  filter_data = TRUE,
  minFeature = 2000,
  minCount = 1000,
  percent.mt = 20
)
}
\arguments{
\item{sce}{seurat object}

\item{eset}{expression set, default is null}

\item{file_type}{if file is from 10X, data_path must be provide}

\item{data_path}{data_path of 10X data}

\item{project}{project name}

\item{index}{folder name of save_path}

\item{id_type}{symbol or ensembl}

\item{normal_cell_id}{}

\item{check_data}{if true}

\item{filter_data}{if TRUE, cells will be filtered by aforementioned criterion}

\item{minFeature}{minimum feature of cell}

\item{minCount}{minimum count of cell}

\item{percent.mt}{max of percentage of mitochondrial gens}
}
\description{
modified_copykat
}
\examples{
copykat_res<-copykat_plus(eset = sc_tnbc, file_type = "10X", data_path = NULL, index = 1, project = "TNBC", id_type = "symbol", check_data = TRUE, filter_data = FALSE, minFeature = 2000, minCount = 1000, percent.mt =20)

}
