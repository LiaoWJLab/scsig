% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/standard_sc.R
\name{standard_sc}
\alias{standard_sc}
\title{standard_Seurat_data_processing}
\usage{
standard_sc(
  eset = NULL,
  sce_object = NULL,
  assay = "RNA",
  file_type = "10X",
  data_path = NULL,
  project = "sc",
  qc_identity = NULL,
  nPCs = 30,
  res = 1,
  verbose = FALSE,
  index = 1,
  plot = FALSE,
  vars_dim = "Phase",
  minFeature = 2000,
  minCount = 1000,
  species = "human",
  cutoff_percent_mt = 25,
  filter_min_cells = FALSE,
  min_cell_count = 20,
  findmarkers = FALSE,
  cols = "normal",
  palette = 1,
  show_col = FALSE,
  seed = 123,
  already_normalized = FALSE,
  save_path = NULL,
  save_sce_object = TRUE
)
}
\arguments{
\item{eset}{Either a matrix-like object with unnormalized data with cells as columns and features as rows or an Assay-derived object}

\item{sce_object}{Seurat object}

\item{assay}{Assay to pull data from such as RNA,sct,integrated, default is RNA}

\item{file_type}{If file_type = "10X", data_path must be provided}

\item{data_path}{Absolute path of sparse data matrices provided by 10X genomics}

\item{project}{Project name for the Seurat object}

\item{qc_identity}{Default is Null}

\item{nPCs}{Maximum PC of input}

\item{res}{Value of the resolution parameter. The higher the value, the larger number of communities you will obtain.}

\item{verbose}{Whether to show progress bar for procedure}

\item{index}{Prefix of folder name}

\item{plot}{If TURE, relevant plot will be drawn}

\item{vars_dim}{A vector of variables to group cells by}

\item{minFeature}{Minimum features of cell}

\item{minCount}{Minimum feature counts of cell}

\item{species}{Species name, such as human or mouse}

\item{cutoff_percent_mt}{Cutoffs for filtering cells that have >n percent mitochondrial gene counts, default is 25}

\item{filter_min_cells}{Whether to filter out the clusters that contain < n cells, default is FALSE}

\item{min_cell_count}{Minimal cell counts of the clusters, default is 20}

\item{findmarkers}{If TRUE, top 6 marker genes of each cluster will be defined and drawn}

\item{cols}{Vector of colors, users can define the cols manually.  This may also be a single character, such as normal and random,  to a palette as specified by `palettes(){IOBR}`
See [palettes](http://127.0.0.1:60491/help/library/IOBR/html/palettes.html) for details}

\item{palette}{Numeric value corresponding with color palette. Default is 1, other options: 2, 3, 4}

\item{show_col}{Whether to show color palettes}

\item{seed}{Seed of the random number generator, default is 123. The parameter works when cols ="random"}

\item{already_normalized}{If TRUE, `sce` object will not be normalized}

\item{save_path}{Path of the output saving directory}

\item{save_sce_object}{If TRUE, Suerat object will be saved}
}
\value{
Returns the standard processed object
}
\description{
The standard data process includes quality control, dimension reduction, unsupervised clustering and differential expression
}
\examples{
## Not run:
data("pbmc_small")
pbmc_small<-standard_sc( sce_object=pbmc_small )
## End(Not run)
View the relationship between dimensionality reduction data and cell annotation
}
