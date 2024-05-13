#' Stat Database Attributes.
#'
#' @param df All metadata, can be \code{PanglaoDBMeta} and obtained with \code{ShowCBDatasets}, \code{ShowCELLxGENEDatasets},
#' and \code{ShowHCAProjects}. Skip when \code{use.census} is TRUE. Default: NULL.
#' @param filter Vector of attributes.
#' @param database Database name, choose from "PanglaoDB", "UCSC", "CELLxGENE", "HCA". Default: "PanglaoDB".
#' @param combine Logical value, whether to combine all attributes in \code{filter} for summarising.
#' Default: FALSE.
#' @param use.census Logical value, whether to use CZ CELLxGENE Census to summarise metadata. Default: FALSE.
#' @param census.version The version of the Census, e.g., "2024-05-13", or "latest" or "stable". Default: stable.
#' @param organism Organism, should be in lower case and replace space with '_'. Default: FALSE (human).
#'
#' @return List of attributes information (attribute, value and number) when \code{combine} is FALSE,
#' dataframe when \code{combine} is TRUE.
#' @importFrom magrittr %>%
#' @importFrom dplyr filter if_all everything group_by_all summarise n arrange desc
#' @importFrom tidyr separate_rows all_of
#' @import cellxgene.census
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
StatDBAttribute <- function(df = NULL, filter, database = c("PanglaoDB", "UCSC", "CELLxGENE", "HCA"), combine = FALSE,
                            use.census = FALSE, census.version = "stable", organism = NULL) {
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
  # check cell_type
  if (use.census && database == "CELLxGENE") {
    message(
      "The use.census is true, ",
      "we will use R package cellxgene.census to access CELLxGENE. Please note the Census release date, ",
      "there may be a delay (control by census.version)!"
    )
    # require cellxgene.census
    if (!require("cellxgene.census", quietly = TRUE, character.only = TRUE)) {
      stop(
        "Please install cellxgene.census package! You can try install.packages('cellxgene.census', ",
        "repos=c('https://chanzuckerberg.r-universe.dev', 'https://cloud.r-project.org'))"
      )
    }
    suppressWarnings(suppressMessages(library("cellxgene.census", character.only = TRUE)))
    if (is.null(organism)) {
      message("'cellxgene.census' requires organism, the default is human")
      organism <- "homo_sapiens"
    } else {
      if (grepl(pattern = "^[A-Z]+", x = organism)) {
        message("Detect upper case in ", organism, ". Convert to lower case!")
        organism <- tolower(organism)
      }
      if (grepl(pattern = "[[:space:]]+", x = organism)) {
        message("Detect space in ", organism, ". Replace with '_'!")
        organism <- gsub(pattern = "[[:space:]]+", replacement = "_", x = organism)
      }
    }
    census <- cellxgene.census::open_soma(census_version = census.version)
    sp <- census$get("census_data")$get(organism)
    all.cols <- sp$obs$colnames()
    valid.census.filter <- intersect(filter, all.cols)
    invalid.census.filter <- setdiff(filter, all.cols)
    if (length(invalid.census.filter) > 0) {
      message(
        "Detect invalid filter: ", paste(invalid.census.filter, collapse = ", "),
        ", please choose from: ", paste(all.cols, collapse = ", ")
      )
    }
    census.obs.df <- sp$obs$read(column_names = valid.census.filter)
    census.obs.df <- as.data.frame(census.obs.df$concat())
    # close census to release memory and other resources
    census$close()
    census.obs.df.stat <- census.obs.df %>%
      dplyr::group_by_all() %>%
      dplyr::summarise(Num = dplyr::n()) %>%
      dplyr::arrange(dplyr::desc(Num))
    return(census.obs.df.stat)
  } else {
    if (is.null(df)) {
      stop("Please provide metadata!")
    }
    # get valid filter
    valid.filter <- intersect(filter, names(filter.vec))
    # check valid filter
    if (length(valid.filter) == 0) {
      stop("The filter you provided is not valid, choose from: ", paste0(names(filter.vec), collapse = ", "))
    }
    # get filter values
    valid.filter.res <- CheckFilter(df = df, filter = valid.filter, all.filter = filter.vec, database = database, combine = combine)
    return(valid.filter.res)
  }
}
