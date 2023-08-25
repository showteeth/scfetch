#' Extract Metadata of Human Cell Atlas Projects.
#'
#' @return Dataframe contains all available projects.
#' @importFrom magrittr %>%
#' @importFrom curl curl_fetch_memory
#' @importFrom jsonlite fromJSON
#' @importFrom data.table rbindlist
#' @export
#' @references https://bioconductor.org/packages/release/bioc/html/hca.html
#'
#' @examples
#' # # all available projects
#' # all.hca.projects = ExtractHCAMeta()
ExtractHCAMeta <- function() {
  # base urls
  hca.base.url <- "https://service.azul.data.humancellatlas.org"
  # get all catalogs
  catalog.url <- paste0(hca.base.url, "/index/catalogs")
  catalog.content <- URLRetrieval(catalog.url)
  catalog.vec <- sapply(names(catalog.content$catalogs), function(x) {
    if (!catalog.content$catalogs[[x]]$internal) x
  }) %>% unlist()
  # get catalog project list
  hca.projects.list <- lapply(catalog.vec, function(x) {
    cat.prj.url <- paste0(hca.base.url, "/index/projects?catalog=", x, "&size=100")
    RecurURLRetrieval(cat.prj.url)
  })
  # get catalog project df
  hca.projects.df <- data.table::rbindlist(hca.projects.list, fill = TRUE) %>% as.data.frame()
  # get project detail information
  hca.projects.detail.list <- lapply(1:nrow(hca.projects.df), function(x) {
    message(x)
    x.df <- hca.projects.df[x, ]
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
    # 下载链接可以直接从 x.df.projects$matrices 以及 x.df.projects$contributedAnalyses 中提取

    # samples和specimens中的信息重复
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
    x.df.cellSuspensions <- x.df$cellSuspensions
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
      selectedCellType = selectedCellType, sourceId = sourceId, sourceSpec = sourceSpec, workflow = workflow,
      libraryConstructionApproach = libraryConstructionApproach, nucleicAcidSource = nucleicAcidSource,
      instrumentManufacturerModel = instrumentManufacturerModel, pairedEnd = pairedEnd, cellLineID = cellLineID,
      cellLineType = cellLineType, cellLinemodelOrgan = cellLinemodelOrgan, organoidsID = organoidsID,
      organoidsmodelOrgan = organoidsmodelOrgan, organoidsmodelOrganPart = organoidsmodelOrganPart,
      lastModifiedDate = lastModifiedDate
    )
  })
  hca.projects.detail.df <- data.table::rbindlist(hca.projects.detail.list, fill = TRUE) %>% as.data.frame()
  return(hca.projects.detail.df)
}
