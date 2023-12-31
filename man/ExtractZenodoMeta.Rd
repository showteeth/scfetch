% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Zenodo.R
\name{ExtractZenodoMeta}
\alias{ExtractZenodoMeta}
\title{Prepare Dataframe with Zenodo DOIs.}
\usage{
ExtractZenodoMeta(doi, file.ext = c("rdata", "h5ad"))
}
\arguments{
\item{doi}{A vector of Zenodo DOIs, should start with "10.5281/zenodo.".}

\item{file.ext}{The valid file extension for download. When NULL, use all files. Default: c("rdata", "h5ad").}
}
\value{
Dataframe contains files with valid extension in given Zenodo DOI.
}
\description{
Prepare Dataframe with Zenodo DOIs.
}
\examples{
\donttest{
zebrafish.df <- ExtractZenodoMeta(doi = "10.5281/zenodo.7243603")
ExtractZenodoMeta(doi = "10.5281/zenodo.48065") # Restricted Access
# vector of dois
multi.dois <- ExtractZenodoMeta(doi = c(
  "1111", "10.5281/zenodo.7243603",
  "10.5281/zenodo.7244441"
))
}
}
