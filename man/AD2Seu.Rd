% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IObenchmark.R
\name{AD2Seu}
\alias{AD2Seu}
\title{Convert AnnData to SeuratObject using SeuratDisk/sceasy/scDIOR/schard/SeuratDisk+scDIOR.}
\usage{
AD2Seu(
  anndata.file,
  method = c("SeuratDisk", "sceasy", "scDIOR", "schard", "SeuratDisk+scDIOR"),
  assay = "RNA",
  load.assays = "RNA",
  slot = "counts",
  use.raw = TRUE
)
}
\arguments{
\item{anndata.file}{The file contains AnnData.}

\item{method}{Method used to perform conversion, choose from "SeuratDisk", "sceasy",    "scDIOR", "schard",
"SeuratDisk+scDIOR". Default: "SeuratDisk".}

\item{assay}{Name to store assay data as. Default: RNA.}

\item{load.assays}{Which assays to load, used when \code{method} is "SeuratDisk". Default: RNA.}

\item{slot}{Slot to store adata.X, used when \code{method} is "sceasy". Default: counts.}

\item{use.raw}{Logical value, whether to use adata.raw, used when \code{method} is "schard". Default: TRUE.}
}
\value{
A SeuratObject.
}
\description{
Convert AnnData to SeuratObject using SeuratDisk/sceasy/scDIOR/schard/SeuratDisk+scDIOR.
}
\examples{
# SeuratDisk
ann.seu <- AD2Seu(anndata.file = "pbmc3k.h5ad", method = "SeuratDisk", assay = "RNA", load.assays = c("RNA"))
# sceasy
ann.sceasy <- AD2Seu(anndata.file = "pbmc3k.h5ad", method = "sceasy", assay = "RNA", slot = "scale.data")
# scDIOR
ann.scdior <- AD2Seu(anndata.file = "pbmc3k.h5ad", method = "scDIOR", assay = "RNA")
# schard
ann.schard <- AD2Seu(anndata.file = "pbmc3k.h5ad", method = "schard", assay = "RNA", use.raw = T)
# SeuratDisk+scDIOR
ann.seuscdior <- AD2Seu(anndata.file = "pbmc3k.h5ad", method = "SeuratDisk+scDIOR", assay = "RNA", load.assays = c("RNA"))
}
