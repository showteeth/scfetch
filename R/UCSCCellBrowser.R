#' Show All Available Datasets in UCSC Cell Browser.
#'
#' @return Dataframe contains all available datasets.
#' @importFrom magrittr %>%
#' @importFrom jsonlite fromJSON flatten
#' @importFrom data.table rbindlist
#' @importFrom dplyr select filter mutate starts_with
#' @export
#'
#' @examples
#' # ucsc.cb.samples = ShowCBDatasets()
ShowCBDatasets <- function() {
  # parse all datasets json
  all.datasets.json <- jsonlite::fromJSON(txt = "https://cells.ucsc.edu/dataset.json")
  all.datasets.df <- jsonlite::flatten(all.datasets.json$datasets)
  colnames(all.datasets.df) <- gsub(pattern = ".*\\.", replacement = "", x = colnames(all.datasets.df))
  # modify elements
  all.datasets.df <- PasteAttr(df = all.datasets.df, attr = c(
    "tags", "diseases", "organisms", "body_parts",
    "projects", "life_stages", "domains", "sources", "assays"
  ))
  # extract all samples
  all.samples.df <- ExtractSample(all.datasets.df) %>% as.data.frame()
  # modify columns
  all.samples.df <- PasteAttr(df = all.samples.df, attr = c(
    "tags", "diseases", "organisms", "body_parts",
    "projects", "life_stages", "domains", "sources", "assays"
  ))
  # remove unused columns
  all.samples.df <- all.samples.df %>% dplyr::select(-c("md5", "hasFiles", "isCollection", "datasetCount", "collectionCount"))
  # add label
  all.samples.df.single <- all.samples.df %>% dplyr::filter(!grepl(x = name, pattern = "/"))
  all.samples.df.multi <- all.samples.df %>%
    dplyr::filter(grepl(x = name, pattern = "/")) %>%
    dplyr::mutate(parent = gsub(pattern = "([^/]*).*", replacement = "\\1", x = name))
  all.datasets.df.multi <- all.datasets.df[!is.na(all.datasets.df$isCollection), c("shortLabel", "name", "tags", "organisms", "body_parts", "diseases", "projects")]
  colnames(all.datasets.df.multi) <- gsub(pattern = "^", replacement = "parent_", x = colnames(all.datasets.df.multi))
  all.samples.df.multi <- merge(all.samples.df.multi, all.datasets.df.multi, by.y = "parent_name", by.x = "parent")
  # modify elements
  all.samples.df.multi <- InheritParient(df = all.samples.df.multi, attr = c("tags", "organisms", "body_parts", "diseases", "projects"))
  # modify label
  colnames(all.samples.df.multi) <- gsub(pattern = "^shortLabel$", replacement = "subLabel", x = colnames(all.samples.df.multi))
  all.samples.df.multi$subLabel <- as.character(all.samples.df.multi$subLabel)
  colnames(all.samples.df.multi) <- gsub(pattern = "^parent_shortLabel$", replacement = "shortLabel", x = colnames(all.samples.df.multi))
  all.samples.df.multi$shortLabel <- as.character(all.samples.df.multi$shortLabel)
  # remove parent columns
  all.samples.df.multi <- all.samples.df.multi %>% dplyr::select(!dplyr::starts_with("parent"))
  # final dataframe
  all.samples.final <- data.table::rbindlist(list(all.samples.df.single, all.samples.df.multi), fill = TRUE) %>% as.data.frame()
  all.samples.final <- all.samples.final[c(
    "shortLabel", "subLabel", "name", "tags", "body_parts", "diseases", "organisms",
    "projects", "life_stages", "domains", "sources", "sampleCount", "assays"
  )]
  # add sample information
  base.url <- "https://cells.ucsc.edu/"
  desc.files <- paste0(base.url, all.samples.final$name, "/desc.json")
  names(desc.files) <- all.samples.final$name
  desc.list <- lapply(names(desc.files), function(x) {
    sd.json <- jsonlite::fromJSON(txt = desc.files[x])
    used.attr <- c(
      "title" = "title", "paper" = "paper_url", "abstract" = "abstract", "matrix" = "matrixFile", "unit" = "unitDesc",
      "coords" = "coordFiles", "methods" = "methods", "geo" = "geo_series"
    )
    sd.df <- ExtractDesc(lst = sd.json, attr = used.attr)
    sd.df$name <- x
  })
  desc.df <- do.call(rbind, desc.list)
  all.samples.final <- merge(all.samples.final, desc.df, by = "name")
  all.samples.final$matrix <- basename(all.samples.final$matrix)
  # return value
  return(all.samples.final)
}


#' Extract Datasets with Attributes.
#'
#' @param collection The collection of the datasets, corresponds to \code{shortLabel} column
#' of \code{ShowCBDatasets}. Default: NULL (without filtering).
#' @param sub.collection The sub-collection of the datasets, corresponds to \code{subLabel} column
#' of \code{ShowCBDatasets}. Default: NULL (without filtering).
#' @param organ The organ of the datasets, corresponds to \code{body_parts} column
#' of \code{ShowCBDatasets}. Default: NULL (without filtering).
#' @param disease The disease of the datasets, corresponds to \code{diseases} column
#' of \code{ShowCBDatasets}. Default: NULL (without filtering).
#' @param organism The specie of the datasets, corresponds to \code{organisms} column
#' of \code{ShowCBDatasets}. Default: NULL (without filtering).
#' @param project The project of the datasets, corresponds to \code{projects} column
#' of \code{ShowCBDatasets}. Default: NULL (without filtering).
#' @param fuzzy.match Logical value, whether to perform fuzzy match with provided attribute values. Default: TRUE.
#'
#' @return Dataframe contains filtered datasets.
#' @importFrom magrittr %>%
#' @importFrom jsonlite fromJSON flatten
#' @importFrom data.table rbindlist
#' @importFrom dplyr select filter mutate starts_with
#' @export
#'
#' @examples
#' #
ExtractCBDatasets <- function(collection = NULL, sub.collection = NULL, organ = NULL, disease = NULL, organism = NULL,
                              project = NULL, fuzzy.match = TRUE) {
  # all sample dataframe
  all.samples.df <- ShowCBDatasets()
  # extract row index under different filter
  collection.idx <- CheckParas(df = all.samples.df, column = "shortLabel", para.value = collection, fuzzy.match = fuzzy.match)
  sub.collection.idx <- CheckParas(df = all.samples.df, column = "subLabel", para.value = sub.collection, fuzzy.match = fuzzy.match)
  organ.idx <- CheckParas(df = all.samples.df, column = "body_parts", para.value = organ, fuzzy.match = fuzzy.match)
  disease.idx <- CheckParas(df = all.samples.df, column = "diseases", para.value = disease, fuzzy.match = fuzzy.match)
  organism.idx <- CheckParas(df = all.samples.df, column = "organisms", para.value = organism, fuzzy.match = fuzzy.match)
  project.idx <- CheckParas(df = all.samples.df, column = "projects", para.value = project, fuzzy.match = fuzzy.match)
  # filter on the whole dataset
  valid.idx <- Reduce(intersect, list(collection.idx, sub.collection.idx, organ.idx, disease.idx, organism.idx, project.idx))
  used.sample.df <- all.samples.df[valid.idx, ]
  return(used.sample.df)
}

#' Download UCSC Cell Browser Datasets.
#'
#' @param sample.df Dataframe contains used datasets. Default: NULL.
#' @param collection The collection of the datasets, corresponds to \code{shortLabel} column
#' of \code{ShowCBDatasets}. Default: NULL (without filtering).
#' @param sub.collection The sub-collection of the datasets, corresponds to \code{subLabel} column
#' of \code{ShowCBDatasets}. Default: NULL (without filtering).
#' @param organ The organ of the datasets, corresponds to \code{body_parts} column
#' of \code{ShowCBDatasets}. Default: NULL (without filtering).
#' @param disease The disease of the datasets, corresponds to \code{diseases} column
#' of \code{ShowCBDatasets}. Default: NULL (without filtering).
#' @param organism The specie of the datasets, corresponds to \code{organisms} column
#' of \code{ShowCBDatasets}. Default: NULL (without filtering).
#' @param project The project of the datasets, corresponds to \code{projects} column
#' of \code{ShowCBDatasets}. Default: NULL (without filtering).
#' @param fuzzy.match Logical value, whether to perform fuzzy match with provided attribute values. Default: TRUE.
#' @param merge Logical value, whether to merge Seurat list. Default: FALSE.
#'
#' @return Seurat object (if \code{merge} is TRUE) or list of Seurat objects (if \code{merge} is FALSE).
#' @importFrom magrittr %>%
#' @importFrom jsonlite fromJSON flatten
#' @importFrom data.table rbindlist fread
#' @importFrom dplyr select filter mutate starts_with full_join
#' @importFrom Seurat CreateSeuratObject
#' @importFrom purrr reduce
#' @importFrom tibble column_to_rownames
#' @export
#'
#' @examples
#' #
DownloadCBDatasets <- function(sample.df = NULL, collection = NULL, sub.collection = NULL, organ = NULL, disease = NULL, organism = NULL,
                               project = NULL, fuzzy.match = TRUE, merge = TRUE) {
  # prepare samples for download
  if (!is.null(sample.df)) {
    # use provided dataframe to download data
    used.sample.df <- sample.df
  } else {
    # all sample dataframe
    all.samples.df <- ShowCBDatasets()
    # extract row index under different filter
    collection.idx <- CheckParas(df = all.samples.df, column = "shortLabel", para.value = collection, fuzzy.match = fuzzy.match)
    sub.collection.idx <- CheckParas(df = all.samples.df, column = "subLabel", para.value = sub.collection, fuzzy.match = fuzzy.match)
    organ.idx <- CheckParas(df = all.samples.df, column = "body_parts", para.value = organ, fuzzy.match = fuzzy.match)
    disease.idx <- CheckParas(df = all.samples.df, column = "diseases", para.value = disease, fuzzy.match = fuzzy.match)
    organism.idx <- CheckParas(df = all.samples.df, column = "organisms", para.value = organism, fuzzy.match = fuzzy.match)
    project.idx <- CheckParas(df = all.samples.df, column = "projects", para.value = project, fuzzy.match = fuzzy.match)
    # filter on the whole dataset
    valid.idx <- Reduce(intersect, list(collection.idx, sub.collection.idx, organ.idx, disease.idx, organism.idx, project.idx))
    used.sample.df <- all.samples.df[valid.idx, ]
  }
  # parepare exp
  base.url <- "https://cells.ucsc.edu/"
  sample.names <- gsub(pattern = "/", replacement = "_", x = used.sample.df$name)
  # prepare exp matrix
  exp.urls <- paste0(base.url, used.sample.df$name, "/exprMatrix.tsv.gz")
  names(exp.urls) <- sample.names
  # out.exp.files = file.path(out.folder, paste(sample.names, "exprMatrix.tsv.gz", sep = "_"))
  # names(out.exp.files) = sample.names
  # prepare metadata
  meta.urls <- paste0(base.url, used.sample.df$name, "/meta.tsv")
  names(meta.urls) <- sample.names
  # out.meta.files = file.path(out.folder, paste(sample.names, "meta.tsv", sep = "_"))
  # names(out.meta.files) = sample.names
  # prepare coord
  coord.files <- sapply(1:nrow(used.sample.df), function(x) {
    paste0(base.url, used.sample.df$name[x], "/", strsplit(x = used.sample.df$coords[x], split = ", ")[[1]])
  })
  names(coord.files) <- sample.names
  # create seurat object
  seu.obj.list <- lapply(sample.names, function(x) {
    message(x)
    Load2Seurat(exp.file = exp.urls[x], meta.file = meta.urls[x], coord.file = coord.files[[x]], name = x)
  })
  # merge or not
  if (isTRUE(merge)) {
    seu.obj <- mergeExperiments(seu.obj.list)
    return(seu.obj)
  } else {
    return(seu.obj.list)
  }
}
