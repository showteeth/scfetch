# extract all projects
ExtractHCAProjects <- function(catalog = NULL) {
  # base urls
  hca.base.url <- "https://service.azul.data.humancellatlas.org"
  # get all catalogs
  catalog.url <- paste0(hca.base.url, "/index/catalogs")
  catalog.content <- URLRetrieval(catalog.url)
  catalog.vec <- sapply(names(catalog.content$catalogs), function(x) {
    if (!catalog.content$catalogs[[x]]$internal) x
  }) %>% unlist()
  if (!is.null(catalog)) {
    catalog.vec <- intersect(catalog, catalog.vec)
    if (length(catalog.vec) == 0) {
      stop("Please check catalog you proided, choose from: dcp29, dcp30, dcp1, lm2, lm3.")
    }
  }
  # get catalog project list
  hca.projects.list <- lapply(catalog.vec, function(x) {
    cat.prj.url <- paste0(hca.base.url, "/index/projects?catalog=", x, "&size=100")
    cat.prj <- RecurURLRetrieval(cat.prj.url) %>% as.data.frame()
    cat.prj$catalog <- x
    return(cat.prj)
  })
  # get catalog project df
  hca.projects.df <- data.table::rbindlist(hca.projects.list, fill = TRUE) %>% as.data.frame()
  # remove duplicated projects in different catalogs
  hca.projects.df <- hca.projects.df %>% dplyr::distinct(.data[["entryId"]], .keep_all = TRUE)
  return(hca.projects.df)
}

#' Show All Available Projects in Human Cell Atlas.
#'
#' @param catalog The catalog of the projects. Different catalogs may share some projects. Choose from "dcp29",
#' "dcp30", "dcp1", "lm2", "lm3", one or multiple values. Default: NULL (all catalogs, remove duplicated projects).
#'
#' @return Dataframe contains all available projects.
#' @importFrom magrittr %>%
#' @importFrom curl curl_fetch_memory
#' @importFrom jsonlite fromJSON
#' @importFrom data.table rbindlist
#' @importFrom dplyr distinct
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' # # all available projects
#' # all.hca.projects = ShowHCAProjects()
ShowHCAProjects <- function(catalog = NULL) {
  # get all projects information
  hca.projects.df <- ExtractHCAProjects(catalog = catalog)
  # get project detail information
  hca.projects.detail.list <- lapply(1:nrow(hca.projects.df), function(x) {
    x.df <- hca.projects.df[x, ]

    # entryid and catalog
    entryId <- x.df$entryId
    catalog <- x.df$catalog

    # procotol information
    x.df.protocol <- x.df$protocols[[1]]
    workflow <- HCAPasteCol(x.df.protocol, col = "workflow")
    libraryConstructionApproach <- HCAPasteCol(x.df.protocol, col = "libraryConstructionApproach")
    nucleicAcidSource <- HCAPasteCol(x.df.protocol, col = "nucleicAcidSource")
    instrumentManufacturerModel <- HCAPasteCol(x.df.protocol, col = "instrumentManufacturerModel")
    pairedEnd <- HCAPasteCol(x.df.protocol, col = "pairedEnd")

    # source
    x.df.source <- x.df$sources[[1]]
    sourceId <- HCAPasteColdf(x.df.source, col = "sourceId")
    sourceSpec <- HCAPasteColdf(x.df.source, col = "sourceSpec")

    # project
    x.df.projects <- x.df$projects[[1]]
    projectId <- HCAPasteColdf(x.df.projects, col = "projectId")
    projectTitle <- HCAPasteColdf(x.df.projects, col = "projectTitle")
    projectShortname <- HCAPasteColdf(x.df.projects, col = "projectShortname")
    laboratory <- HCAPasteCol(x.df.projects, col = "laboratory")
    estimatedCellCount <- HCAPasteColdf(x.df.projects, col = "estimatedCellCount")
    projectDescription <- HCAPasteColdf(x.df.projects, col = "projectDescription")
    publications <- HCAPasteColdf(x.df.projects$publications[[1]], col = "publicationTitle")
    accessions <- HCAPasteColdf(x.df.projects$accessions[[1]], col = "accession")
    accessible <- HCAPasteColdf(x.df.projects, col = "accessible")

    # sample
    x.df.samples <- x.df$samples[[1]]
    sampleEntityType <- HCAPasteCol(df = x.df.samples, col = "sampleEntityType")
    organ <- HCAPasteCol(df = x.df.samples, col = "effectiveOrgan")
    sampleID <- HCAPasteCol(df = x.df.samples, col = "id")
    organPart <- HCAPasteCol(df = x.df.samples, col = "organPart")
    disease <- HCAPasteCol(df = x.df.samples, col = "disease")
    preservationMethod <- HCAPasteCol(df = x.df.samples, col = "preservationMethod")

    # cell line
    x.df.cellLines <- x.df$cellLines[[1]]
    cellLineID <- HCAPasteCol(df = x.df.cellLines, col = "id")
    cellLineType <- HCAPasteCol(df = x.df.cellLines, col = "cellLineType")
    cellLinemodelOrgan <- HCAPasteCol(df = x.df.cellLines, col = "modelOrgan")

    # Organism
    x.df.organisms <- x.df$donorOrganisms[[1]]
    donorCount <- ifelse(is.null(x.df.organisms$donorCount), "", x.df.organisms$donorCount)
    developmentStage <- HCAPasteCol(df = x.df.organisms, col = "developmentStage")
    genusSpecies <- HCAPasteCol(df = x.df.organisms, col = "genusSpecies")
    biologicalSex <- HCAPasteCol(df = x.df.organisms, col = "biologicalSex")

    # organoids
    x.df.organoids <- x.df$organoids[[1]]
    organoidsID <- HCAPasteCol(df = x.df.organoids, col = "id")
    organoidsmodelOrgan <- HCAPasteCol(df = x.df.organoids, col = "modelOrgan")
    organoidsmodelOrganPart <- HCAPasteCol(df = x.df.organoids, col = "modelOrganPart")

    # cellSuspensions
    x.df.cellSuspensions <- x.df$cellSuspensions[[1]]
    selectedCellType <- HCAPasteCol(df = x.df.cellSuspensions, col = "selectedCellType")

    # date
    x.df.date <- x.df$dates[[1]]
    lastModifiedDate <- HCAPasteColdf(df = x.df.date, col = "lastModifiedDate")

    # return final dataframe
    data.frame(
      projectTitle = projectTitle, projectId = projectId, projectShortname = projectShortname,
      projectDescription = projectDescription, publications = publications, laboratory = laboratory,
      accessions = accessions, accessible = accessible, estimatedCellCount = estimatedCellCount,
      sampleEntityType = sampleEntityType, organ = organ, organPart = organPart, sampleID = sampleID,
      disease = disease, preservationMethod = preservationMethod, donorCount = donorCount,
      developmentStage = developmentStage, genusSpecies = genusSpecies, biologicalSex = biologicalSex,
      selectedCellType = selectedCellType, catalog = catalog, entryId = entryId, sourceId = sourceId, sourceSpec = sourceSpec,
      workflow = workflow, libraryConstructionApproach = libraryConstructionApproach, nucleicAcidSource = nucleicAcidSource,
      instrumentManufacturerModel = instrumentManufacturerModel, pairedEnd = pairedEnd, cellLineID = cellLineID,
      cellLineType = cellLineType, cellLinemodelOrgan = cellLinemodelOrgan, organoidsID = organoidsID,
      organoidsmodelOrgan = organoidsmodelOrgan, organoidsmodelOrganPart = organoidsmodelOrganPart,
      lastModifiedDate = lastModifiedDate
    )
  })
  hca.projects.detail.df <- data.table::rbindlist(hca.projects.detail.list, fill = TRUE) %>% as.data.frame()
  return(hca.projects.detail.df)
}

#' Extract Metadata of Human Cell Atlas Projects with Attributes.
#'
#' @param all.projects.df All detail information of HCA projects, obtained with \code{ShowHCAProjects}.
#' @param organism The organism of the projects, choose from "Homo sapiens", "Mus musculus",
#' "Macaca mulatta", "canis lupus familiaris", one or multiple values. Default: NULL (All).
#' @param sex The sex of the projects, choose from "female", "male", "mixed", "unknown",
#' one or multiple values. Default: NULL (All).
#' @param organ The organ of the projects (e.g. brain), obtain available values with \code{StatDBAttribute},
#' one or multiple values. Default: NULL (All).
#' @param organ.part The organ part of the projects (e.g. cortex), obtain available values with \code{StatDBAttribute},
#' one or multiple values. Default: NULL (All).
#' @param disease The disease of the projects (e.g. normal), obtain available values with \code{StatDBAttribute},
#' one or multiple values. Default: NULL (All).
#' @param sample.type The sex of the projects, choose from "specimens", "organoids", "cellLines",
#' one or multiple values. Default: NULL (All).
#' @param preservation.method The preservation method of the projects (e.g. fresh), obtain available values with \code{StatDBAttribute},
#' one or multiple values. Default: NULL (All).
#' @param protocol The protocol of the projects (e.g. 10x 3' v2), obtain available values with \code{StatDBAttribute},
#' one or multiple values. Default: NULL (All).
#' @param suspension.type The suspension type of the projects, choose from "single cell", "single nucleus", "bulk cell", "bulk nuclei",
#' one or multiple values. Default: NULL (All).
#' @param cell.type The cell type of the projects (e.g. neuron), obtain available values with \code{StatDBAttribute},
#' one or multiple values. Default: NULL (All).
#' @param cell.num Cell number filter. If NULL, no filter; if one value, lower filter; if two values, low and high filter.
#' Deault: NULL(without filtering).
#' @param sequencing.type The sequencing instrument type of the projects (e.g. illumina hiseq 2500),
#' obtain available values with \code{StatDBAttribute}, one or multiple values. Default: NULL (All).
#'
#' @return Dataframe contains filtered projects.
#' @importFrom magrittr %>%
#' @importFrom curl curl_fetch_memory
#' @importFrom jsonlite fromJSON
#' @importFrom data.table rbindlist
#' @export
#' @references https://bioconductor.org/packages/release/bioc/html/hca.html
#'
#' @examples
#' # # all available projects
#' # all.hca.projects = ShowHCAProjects()
#' # # all human projects
#' # all.human.projects = ExtractHCAMeta(all.projects.df = all.hca.projects, organism = "Homo sapiens")
#' # # all human and 10x 3' v2
#' # all.human.10x.projects = ExtractHCAMeta(all.projects.df = all.hca.projects,  organism = "Homo sapiens",
#' #                                         protocol = c("10x 3' v2", "10x 3' v3"))
ExtractHCAMeta <- function(all.projects.df, organism = NULL, sex = NULL, organ = NULL, organ.part = NULL, disease = NULL,
                           sample.type = NULL, preservation.method = NULL, protocol = NULL,
                           suspension.type = NULL, cell.type = NULL, cell.num = NULL, sequencing.type = NULL) {
  # all projects detail dataframe
  hca.projects.detail.df <- all.projects.df
  # extract row index under different filter
  organism.idx <- HCAAttrFilter(df = hca.projects.detail.df, attr = "genusSpecies", attr.value = organism)
  sex.idx <- HCAAttrFilter(df = hca.projects.detail.df, attr = "biologicalSex", attr.value = sex)
  organ.idx <- HCAAttrFilter(df = hca.projects.detail.df, attr = "organ", attr.value = organ)
  organ.part.idx <- HCAAttrFilter(df = hca.projects.detail.df, attr = "organPart", attr.value = organ.part)
  disease.idx <- HCAAttrFilter(df = hca.projects.detail.df, attr = "disease", attr.value = disease)
  sample.type.idx <- HCAAttrFilter(df = hca.projects.detail.df, attr = "sampleEntityType", attr.value = sample.type)
  preservation.method.idx <- HCAAttrFilter(df = hca.projects.detail.df, attr = "preservationMethod", attr.value = preservation.method)
  protocol.idx <- HCAAttrFilter(df = hca.projects.detail.df, attr = "libraryConstructionApproach", attr.value = protocol)
  suspension.type.idx <- HCAAttrFilter(df = hca.projects.detail.df, attr = "nucleicAcidSource", attr.value = suspension.type)
  cell.type.idx <- HCAAttrFilter(df = hca.projects.detail.df, attr = "selectedCellType", attr.value = cell.type)
  sequencing.type.idx <- HCAAttrFilter(df = hca.projects.detail.df, attr = "instrumentManufacturerModel", attr.value = sequencing.type)
  if (is.null(cell.num)) {
    cnum.idx <- 1:nrow(hca.projects.detail.df)
  } else if (length(cell.num) == 1) {
    cnum.idx <- which(hca.projects.detail.df$estimatedCellCount > as.numeric(cell.num))
  } else {
    cnum.idx <- which(hca.projects.detail.df$estimatedCellCount > as.numeric(cell.num[1]) &
      hca.projects.detail.df$estimatedCellCount < as.numeric(cell.num[2]))
  }
  # filter on the whole dataset
  valid.idx <- Reduce(intersect, list(
    organism.idx, sex.idx, organ.idx, organ.part.idx, disease.idx, sample.type.idx,
    preservation.method.idx, protocol.idx, suspension.type.idx, cell.type.idx, sequencing.type.idx, cnum.idx
  ))
  used.sample.df <- hca.projects.detail.df[valid.idx, ]
  rownames(used.sample.df) <- NULL
  return(used.sample.df)
}
