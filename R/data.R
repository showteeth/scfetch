#' All Sample Metadata of PanglaoDB Datasets
#'
#' A dataset contains all sample metadata and cell types of PanglaoDB datasets.
#'
#' @format A data frame with 1,368 rows and 7 variables:
#' \describe{
#'   \item{SRA}{The SRA identifier}
#'   \item{SRS}{The SRS identifier}
#'   \item{Tissue}{The tissue of the sample}
#'   \item{Protocol}{The protocol used to generate this sample}
#'   \item{Species}{The species of this sample}
#'   \item{Cells}{Total cell number of this sample}
#'   \item{CellType}{Predicted cell types, separated by comma}
#' }
"PanglaoDBMeta"


#' All Sample Composition of PanglaoDB Datasets
#'
#' A dataset contains all sample composition of PanglaoDB datasets.
#'
#' @format A data frame with 19,449 rows and 8 variables:
#' \describe{
#'   \item{SRA}{The SRA identifier}
#'   \item{SRS}{The SRS identifier}
#'   \item{Tissue}{The tissue of the sample}
#'   \item{Protocol}{The protocol used to generate this sample}
#'   \item{Species}{The species of this sample}
#'   \item{Cluster}{Seurat cluster}
#'   \item{Cells}{Cluster cell number}
#'   \item{Cell Type}{Predicted cluster cell types}
#' }
"PanglaoDBComposition"
