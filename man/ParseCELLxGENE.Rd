% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cellxgene.R
\name{ParseCELLxGENE}
\alias{ParseCELLxGENE}
\title{Download CELLxGENE Datasets.}
\usage{
ParseCELLxGENE(
  meta,
  file.ext = c("rds", "h5ad"),
  out.folder = NULL,
  timeout = 3600,
  quiet = FALSE,
  parallel = TRUE
)
}
\arguments{
\item{meta}{Metadata used to download, can be from \code{ExtractCELLxGENEMeta},
should contain dataset_id, rds_id/h5ad_id (depend on \code{file.ext}) and name columns.}

\item{file.ext}{The valid file extension for download. When NULL, use "rds" and "h5ad". Default: c("rds", "h5ad").}

\item{out.folder}{The output folder. Default: NULL (current working directory).}

\item{timeout}{Maximum request time. Default: 3600.}

\item{quiet}{Logical value, whether to show downloading progress. Default: FALSE (show).}

\item{parallel}{Logical value, whether to download parallelly. Default: TRUE. When "libcurl" is available for \code{download.file},
the parallel is done by default (\code{parallel} can be FALSE).}
}
\value{
Dataframe contains failed datasets or NULL.
}
\description{
Download CELLxGENE Datasets.
}
\examples{
\dontrun{
# all available datasets
all.cellxgene.datasets <- ShowCELLxGENEDatasets()
# human 10x v2 and v3 datasets
human.10x.cellxgene.meta <- ExtractCELLxGENEMeta(
  all.samples.df = all.cellxgene.datasets,
  assay = c("10x 3' v2", "10x 3' v3"),
  organism = "Homo sapiens"
)
# download, need to provide the output folder
ParseCELLxGENE(meta = human.10x.cellxgene.meta, out.folder = "/path/to/output")
}
}
\references{
https://gist.github.com/ivirshup/f1a1603db69de3888eacb4bdb6a9317a
}
