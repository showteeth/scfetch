#' Stat Database Attributes.
#'
#' @param df All metadata, can be \code{PanglaoDBMeta} and obtained with \code{ShowCBDatasets}, \code{ShowCELLxGENEDatasets},
#' and \code{ShowHCAProjects}.
#' @param filter Vector of attributes.
#' @param database Database name, choose from "PanglaoDB", "UCSC", "CELLxGENE", "HCA". Default: "PanglaoDB".
#'
#' @return List of attributes information, including attribute, value and number.
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' \dontrun{
#' # PanglaoDB
#' StatDBAttribute(df = PanglaoDBMeta, filter = c("species", "protocol"), database = "PanglaoDB")
#' # UCSC Cell Browser, need users to provide the json folder
#' ucsc.cb.samples <- ShowCBDatasets(lazy = TRUE, json.folder = NULL, update = FALSE)
#' StatDBAttribute(df = ucsc.cb.samples, filter = c("organism", "organ"), database = "UCSC")
#' # CELLxGENE
#' all.cellxgene.datasets <- ShowCELLxGENEDatasets()
#' StatDBAttribute(
#'   df = all.cellxgene.datasets, filter = c("organism", "sex"),
#'   database = "CELLxGENE"
#' )
#' # HCA
#' all.hca.projects <- ShowHCAProjects()
#' StatDBAttribute(df = all.hca.projects, filter = c("organism", "sex"), database = "HCA")
#' }
StatDBAttribute <- function(df, filter, database = c("PanglaoDB", "UCSC", "CELLxGENE", "HCA")) {
  # check parameters
  database <- match.arg(arg = database)
  # prepare filter vector
  if (database == "PanglaoDB") {
    filter.vec <- c("Species", "Protocol", "Tissue", "SRA", "SRS", "CellType")
    names(filter.vec) <- c("species", "protocol", "tissue", "sra", "srs", "cell.type")
  } else if (database == "UCSC") {
    filter.vec <- c("shortLabel", "subLabel", "body_parts", "diseases", "organisms", "projects")
    names(filter.vec) <- c("collection", "sub.collection", "organ", "disease", "organism", "project")
  } else if (database == "CELLxGENE") {
    filter.vec <- c("organism", "self_reported_ethnicity", "sex", "tissue", "disease", "assay", "suspension_type", "cell_type")
    names(filter.vec) <- c("organism", "ethnicity", "sex", "tissue", "disease", "assay", "suspension.type", "cell.type")
  } else if (database == "HCA") {
    filter.vec <- c(
      "genusSpecies", "biologicalSex", "organ", "organPart", "disease", "sampleEntityType",
      "preservationMethod", "libraryConstructionApproach", "nucleicAcidSource", "selectedCellType",
      "instrumentManufacturerModel"
    )
    names(filter.vec) <- c(
      "organism", "sex", "organ", "organ.part", "disease", "sample.type", "preservation.method",
      "protocol", "suspension.type", "cell.type", "sequencing.type"
    )
  }
  # get valid filter
  valid.filter <- intersect(filter, names(filter.vec))
  # check valid filter
  if (length(valid.filter) == 0) {
    stop("The filter you provided is not valid, choose from: ", paste0(names(filter.vec), collapse = ", "))
  }
  # get filter values
  valid.filter.list <- CheckFilter(df = df, filter = valid.filter, all.filter = filter.vec, database = database)
  return(valid.filter.list)
}
