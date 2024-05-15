#' Show All Available Datasets in CELLxGENE.
#'
#' @return Dataframe contains all available datasets.
#' @importFrom magrittr %>%
#' @importFrom curl curl_fetch_memory
#' @importFrom jsonlite fromJSON flatten
#' @importFrom data.table rbindlist
#' @export
#' @references https://gist.github.com/ivirshup/f1a1603db69de3888eacb4bdb6a9317a
#'
#' @examples
#' \donttest{
#' # all available datasets
#' all.cellxgene.datasets <- ShowCELLxGENEDatasets()
#' }
ShowCELLxGENEDatasets <- function() {
  # urls
  cellxgene.base.url <- "https://api.cellxgene.cziscience.com/curation/v1/"
  cellxgene.collections.url <- paste0(cellxgene.base.url, "collections/")
  # extract all collections
  cellxgene.collections.content <- URLRetrieval(cellxgene.collections.url)

  # extract all datasets
  cellxgene.collections.list <- lapply(1:nrow(cellxgene.collections.content), function(x) {
    cellxgene.collection.content <- cellxgene.collections.content[x, ]
    rownames(cellxgene.collection.content) <- NULL
    cellxgene.sc.url <- paste0(cellxgene.collections.url, cellxgene.collections.content[x, "collection_id"])
    cellxgene.sc.content <- URLRetrieval(cellxgene.sc.url)
    cellxgene.sc.datasets <- jsonlite::flatten(cellxgene.sc.content$datasets)
    colnames(cellxgene.sc.datasets) <- gsub(pattern = "^title$", replacement = "dataset_description", x = colnames(cellxgene.sc.datasets))
    # add metadata
    cellxgene.collection.content$title <- cellxgene.sc.content$name
    cellxgene.collection.content$description <- cellxgene.sc.content$description
    cellxgene.collection.content$contact <- cellxgene.sc.content$contact_name
    cellxgene.collection.content$contact_email <- cellxgene.sc.content$contact_email
    cellxgene.collection.content <- cellxgene.collection.content[c(
      "title", "description", "doi", "contact", "contact_email",
      "collection_id", "collection_url", "consortia",
      "curator_name", "visibility"
    )]
    cellxgene.collection.final <- cbind(cellxgene.collection.content, cellxgene.sc.datasets) %>% as.data.frame()
    return(cellxgene.collection.final)
  })
  # create all datasets dataframe
  cellxgene.collections.datasets.df <- data.table::rbindlist(cellxgene.collections.list, fill = TRUE) %>%
    as.data.frame()
  # modify columns
  cellxgene.collections.datasets.df <-
    PasteAttrCXG(
      df = cellxgene.collections.datasets.df,
      attr = c(
        "assay", "cell_type", "organism", "self_reported_ethnicity", "sex", "tissue",
        "disease", "development_stage"
      ), col = "label"
    )
  cellxgene.collections.datasets.df <-
    PasteAttr(df = cellxgene.collections.datasets.df, attr = c("batch_condition", "suspension_type", "donor_id"))
  # add h5ad and rds information
  cellxgene.collections.datasets.list <- lapply(1:nrow(cellxgene.collections.datasets.df), function(x) {
    x.df <- cellxgene.collections.datasets.df[x, ]
    x.df.dataset <- x.df$assets[[1]]
    # x.df$dataset_id <- unique(x.df.dataset$dataset_id)
    if ("RDS" %in% unique(x.df.dataset$filetype)) {
      x.rds.idx <- which(x.df.dataset$filetype == "RDS")
      x.df$rds_id <- x.df.dataset$url[x.rds.idx]
    } else {
      x.df$rds_id <- NA
    }
    if ("H5AD" %in% unique(x.df.dataset$filetype)) {
      x.h5ad.idx <- which(x.df.dataset$filetype == "H5AD")
      x.df$h5ad_id <- x.df.dataset$url[x.h5ad.idx]
    } else {
      x.df$h5ad_id <- NA
    }
    x.df
  })
  # final dataframe
  cellxgene.collections.datasets.final <- data.table::rbindlist(cellxgene.collections.datasets.list, fill = TRUE) %>%
    as.data.frame()
  return(cellxgene.collections.datasets.final)
}

#' Extract Metadata of CELLxGENE Datasets with Attributes.
#'
#' @param all.samples.df All detail information of CELLxGENE datasets, obtained with \code{ShowCELLxGENEDatasets}.
#' @param organism The organism of the datasets, choose from "Homo sapiens", "Mus musculus", "Callithrix jacchus",
#' "Macaca mulatta", "Sus scrofa domesticus", one or multiple values. Default: NULL (All).
#' @param ethnicity The ethnicity of the datasets, choose from "Asian", "European", "unknown", "na", "African", "Bangladeshi",
#' "British", "Irish", "East Asian", "African American", "Hispanic or Latin American", "African American or Afro-Caribbean",
#' "European American", "Finnish", "Indian", "Japanese", "Korean", "Malaysian", "Singaporean Chinese", "American", "Pacific Islander",
#' "admixed ancestry", "Eskimo", "Han Chinese", "Greater Middle Eastern  (Middle Eastern, North African or Persian)", "multiethnic",
#' "Jewish Israeli", "South Asian", "Oceanian", "Chinese", one or multiple values. Default: NULL (All).
#' @param sex The sex of the datasets, choose from "female", "male", "unknown", one or multiple values. Default: NULL (All).
#' @param tissue The tissue of the datasets, obtain available values with \code{StatDBAttribute}. One or multiple values. Default: NULL (All).
#' @param disease The disease of the datasets, obtain available values with \code{StatDBAttribute}. One or multiple values. Default: NULL (All).
#' @param assay The assay of the datasets, choose from "10x 3' v1", "10x 3' v2", "10x 3' v3", "10x 3' transcription profiling",
#' "10x 5' v1", "10x 5' v2", "10x 5' transcription profiling", "10x multiome", "10x scATAC-seq", "sci-RNA-seq", "Drop-seq",
#' "Smart-seq", "Smart-seq2", "Smart-seq v4", "snmC-Seq2", "Visium Spatial Gene Expression", "Seq-Well", "Seq-Well S3", "Patch-seq",
#' "sci-Plex", "BD Rhapsody Targeted mRNA", "BD Rhapsody Whole Transcriptome Analysis", "Slide-seqV2", "GEXSCOPE technology", "inDrop",
#' "microwell-seq", "CEL-seq2", "STRT-seq", "DroNc-seq", "MERFISH", "scATAC-seq", "MARS-seq", "TruDrop", one or multiple values. Default: NULL (All).
#' @param suspension.type The suspension type of the datasets, choose from "nucleus", "cell", "na", one or multiple values. Default: NULL (All).
#' @param cell.type The cell type of the datasets, obtain available values with \code{StatDBAttribute}. One or multiple values. Default: NULL (All).
#' @param cell.num Cell number filter. If NULL, no filter; if one value, lower filter; if two values, low and high filter.
#' Deault: NULL(without filtering).
#'
#' @return Dataframe contains filtered datasets.
#' @importFrom magrittr %>%
#' @importFrom curl curl_fetch_memory
#' @importFrom jsonlite fromJSON flatten
#' @importFrom data.table rbindlist
#' @export
#' @references https://gist.github.com/ivirshup/f1a1603db69de3888eacb4bdb6a9317a
#'
#' @examples
#' \donttest{
#' # all available datasets
#' all.cellxgene.datasets <- ShowCELLxGENEDatasets()
#' # human 10x v2 and v3 datasets
#' human.10x.cellxgene.meta <- ExtractCELLxGENEMeta(
#'   all.samples.df = all.cellxgene.datasets,
#'   assay = c("10x 3' v2", "10x 3' v3"),
#'   organism = "Homo sapiens"
#' )
#' }
ExtractCELLxGENEMeta <- function(all.samples.df, organism = NULL, ethnicity = NULL, sex = NULL, tissue = NULL, disease = NULL,
                                 assay = NULL, suspension.type = NULL, cell.type = NULL, cell.num = NULL) {
  # all datasets information
  cellxgene.collections.datasets.final <- all.samples.df
  # extract row index under different filter
  organism.idx <- cellxgeneAttrFilter(df = cellxgene.collections.datasets.final, attr = "organism", attr.value = organism)
  ethnicity.idx <- cellxgeneAttrFilter(df = cellxgene.collections.datasets.final, attr = "self_reported_ethnicity", attr.value = ethnicity)
  sex.idx <- cellxgeneAttrFilter(df = cellxgene.collections.datasets.final, attr = "sex", attr.value = sex)
  tissue.idx <- cellxgeneAttrFilter(df = cellxgene.collections.datasets.final, attr = "tissue", attr.value = tissue)
  disease.idx <- cellxgeneAttrFilter(df = cellxgene.collections.datasets.final, attr = "disease", attr.value = disease)
  assay.idx <- cellxgeneAttrFilter(df = cellxgene.collections.datasets.final, attr = "assay", attr.value = assay)
  suspension.type.idx <- cellxgeneAttrFilter(df = cellxgene.collections.datasets.final, attr = "suspension_type", attr.value = suspension.type)
  cell.type.idx <- cellxgeneAttrFilter(df = cellxgene.collections.datasets.final, attr = "cell_type", attr.value = cell.type)
  if (is.null(cell.num)) {
    cnum.idx <- 1:nrow(cellxgene.collections.datasets.final)
  } else if (length(cell.num) == 1) {
    cnum.idx <- which(cellxgene.collections.datasets.final$cell_count > as.numeric(cell.num))
  } else {
    cnum.idx <- which(cellxgene.collections.datasets.final$cell_count > as.numeric(cell.num[1]) &
      cellxgene.collections.datasets.final$cell_count < as.numeric(cell.num[2]))
  }
  # filter on the whole dataset
  valid.idx <- Reduce(intersect, list(
    organism.idx, ethnicity.idx, sex.idx, tissue.idx, disease.idx, assay.idx,
    suspension.type.idx, cell.type.idx, cnum.idx
  ))
  used.sample.df <- cellxgene.collections.datasets.final[valid.idx, ]
  rownames(used.sample.df) <- NULL
  return(used.sample.df)
}

#' Download CELLxGENE Datasets and Return SeuratObject.
#'
#' @param meta Metadata used to download, can be from \code{ExtractCELLxGENEMeta},
#' should contain dataset_id, rds_id/h5ad_id (depend on \code{file.ext}) and name columns.
#' Skip when \code{use.census} is TRUE. Default: NULL.
#' @param file.ext The valid file extension for download. When NULL, use "rds" and "h5ad". Default: c("rds", "h5ad").
#' @param out.folder The output folder. Default: NULL (current working directory).
#' @param timeout Maximum request time. Default: 3600.
#' @param quiet Logical value, whether to show downloading progress. Default: FALSE (show).
#' @param parallel Logical value, whether to download parallelly. Default: TRUE. When "libcurl" is available for \code{download.file},
#' the parallel is done by default (\code{parallel} can be FALSE).
#' @param return.seu Logical value, whether to load downloaded datasets to Seurat. Valid when rds in \code{file.ext} and all
#' datasets download successfully. Default: FALSE.
#' @param merge Logical value, whether to merge Seurat list when there are multiple rds files,
#' used when \code{return.seu} is TRUE. Default: FALSE.
#' @param use.census Logical value, whether to use CZ CELLxGENE Census to download and subset datasets. Default: FALSE.
#' @param census.version The version of the Census, e.g., "2024-05-13", or "latest" or "stable". Default: stable.
#' @param organism Organism, should be in lower case and replace space with '_'. Default: FALSE (human).
#' @param obs.value.filter Filter expression for cell's metadata,
#' e.g., cell_type == 'B cell' & tissue_general == 'lung' & disease == 'COVID-19'. Default: NULL.
#' @param obs.keys Columns to fetch for the cell's metadata. e.g., c("cell_type", "tissue_general", "disease", "sex").
#' @param include.genes Genes to include, e.g, \code{include.genes} = c('ENSG00000161798', 'ENSG00000188229')
#' same as \code{var_value_filter} = "feature_id %in% c('ENSG00000161798', 'ENSG00000188229')" in \code{\link{get_seurat}}.
#' Default: NULL.
#' @param ... Parameters for \code{\link{get_seurat}}, used when \code{use.census} is TRUE.
#'
#' @return Dataframe contains failed datasets, SeuratObject (\code{return.seu} is TRUE, rds in \code{file.ext}) and
#' NULL (\code{return.seu} is FALSE or rds not in \code{file.ext}).
#' @importFrom httr POST stop_for_status content
#' @importFrom jsonlite fromJSON
#' @importFrom parallel detectCores mclapply
#' @importFrom utils download.file
#' @importFrom rlang parse_expr
#' @importFrom dplyr filter
#' @import cellxgene.census
#' @export
#' @references https://gist.github.com/ivirshup/f1a1603db69de3888eacb4bdb6a9317a
#'
#' @examples
#' \dontrun{
#' # all available datasets
#' all.cellxgene.datasets <- ShowCELLxGENEDatasets()
#' # human 10x v2 and v3 datasets
#' human.10x.cellxgene.meta <- ExtractCELLxGENEMeta(
#'   all.samples.df = all.cellxgene.datasets,
#'   assay = c("10x 3' v2", "10x 3' v3"),
#'   organism = "Homo sapiens"
#' )
#' # download, need to provide the output folder
#' ParseCELLxGENE(meta = human.10x.cellxgene.meta, out.folder = "/path/to/output")
#' }
ParseCELLxGENE <- function(meta = NULL, file.ext = c("rds", "h5ad"), out.folder = NULL, timeout = 3600, quiet = FALSE,
                           parallel = TRUE, return.seu = FALSE, merge = TRUE, use.census = FALSE, census.version = "stable",
                           organism = NULL, obs.value.filter = NULL, obs.keys = NULL, include.genes = NULL, ...) {
  if (use.census) {
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
      organism <- "Homo sapiens"
    } else {
      if (grepl(pattern = "^[a-z]+", x = organism)) {
        message("Detect lower case in ", organism, ". Convert to title case!")
        organism <- paste(toupper(substring(organism, 1, 1)),
          tolower(substring(organism, 2, nchar(organism))),
          sep = ""
        )
      }
      if (grepl(pattern = "_", x = organism)) {
        message("Detect '_' in ", organism, ". Replace with spacce!")
        organism <- gsub(pattern = "_", replacement = " ", x = organism)
      }
    }
    census <- cellxgene.census::open_soma(census_version = census.version)
    message("Access CELLxGENE data and create SeuratObject!")
    # create gene filter
    if (!is.null(include.genes)) {
      include.genes.filter <- paste0("feature_id %in% c( '", paste(include.genes, collapse = "', '"), "' )")
      seu.obj <- cellxgene.census::get_seurat(
        census = census, organism = organism, obs_column_names = obs.keys,
        var_value_filter = include.genes.filter, obs_value_filter = obs.value.filter, ...
      )
    } else {
      seu.obj <- cellxgene.census::get_seurat(
        census = census, organism = organism, obs_column_names = obs.keys,
        obs_value_filter = obs.value.filter, ...
      )
    }
    # close census to release memory and other resources
    census$close()
    return(seu.obj)
  } else {
    # check file extension
    if (is.null(file.ext)) {
      warning("There is no file extension provided, use rds and h5ad.")
      file.ext <- c("rds", "h5ad")
    }
    file.ext <- intersect(file.ext, c("rds", "h5ad"))
    if (length(file.ext) == 0) {
      stop("Please provide valid file extension: rds, h5ad.")
    }
    # prepare download urls
    download.urls <- c()
    download.df.list <- list()
    # prepare rds
    if ("rds" %in% file.ext) {
      if (!"rds_id" %in% colnames(meta)) {
        stop("The meta dataframe you provided doesn't contain rds!")
      }
      rds.urls.list <- PrepareCELLxGENEUrls(df = meta, fe = "rds")
      rds.urls <- rds.urls.list$urls
      download.df.list <- c(download.df.list, list(rds.urls.list$df))
      download.urls <- c(download.urls, rds.urls)
    }
    # prepare h5ad
    if ("h5ad" %in% file.ext) {
      if (!"h5ad_id" %in% colnames(meta)) {
        stop("The meta dataframe you provided doesn't contain h5ad!")
      }
      h5ad.urls.list <- PrepareCELLxGENEUrls(df = meta, fe = "h5ad")
      h5ad.urls <- h5ad.urls.list$urls
      download.df.list <- c(download.df.list, list(h5ad.urls.list$df))
      download.urls <- c(download.urls, h5ad.urls)
    }
    download.df <- data.table::rbindlist(download.df.list, fill = TRUE) %>% as.data.frame()
    # prepare output folder
    if (is.null(out.folder)) {
      out.folder <- getwd()
    }
    if (!dir.exists(out.folder)) {
      message(out.folder, " does not exist, create automatically!")
      dir.create(out.folder, recursive = TRUE)
    }
    names(download.urls) <- file.path(out.folder, names(download.urls))
    # download urls
    # set timeout
    env.timeout <- getOption("timeout")
    on.exit(options(timeout = env.timeout)) # restore timeout
    options(timeout = timeout)
    message("Start downloading!")
    if (isTRUE(parallel)) {
      # prepare cores
      cores.used <- min(parallel::detectCores(), length(download.urls))
      down.status <- parallel::mclapply(X = 1:length(download.urls), FUN = function(x) {
        utils::download.file(url = download.urls, destfile = names(download.urls), quiet = quiet, mode = "wb")
      }, mc.cores = cores.used)
    } else {
      down.status <- utils::download.file(url = download.urls, destfile = names(download.urls), quiet = quiet, mode = "wb")
    }
    # process failed datasets
    down.status <- unlist(down.status)
    fail.status <- which(down.status != 0)
    if (length(fail.status) == 0) {
      message("All datasets downloaded successfully!")
      if (isTRUE(return.seu)) {
        seu.obj <- LoadRDS2Seurat(
          out.folder = out.folder, merge = merge, obs.value.filter = obs.value.filter,
          obs.keys = obs.keys, include.genes = include.genes
        )
        return(seu.obj)
      } else {
        return(NULL)
      }
    } else {
      fail.datasets.id <- download.df[fail.status, "dataset_id"] %>% unique()
      fail.meta <- meta[meta$dataset_id %in% fail.datasets.id, ]
      return(fail.meta)
    }
  }
}
