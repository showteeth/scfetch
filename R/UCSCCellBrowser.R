#' Show All Available Datasets in UCSC Cell Browser.
#'
#' @param lazy Logical value, whether to always load datasets online or locally. Default: TRUE (load locally).
#' @param json.folder Folder used to save or load datasets json files. Default: NULL (current working directory).
#' @param update Logical value, whther to update local datasets json. Default: FALSE. For the first time, you should set \code{lazy} TRUE and
#' \code{update} TRUE to save json files.
#' @param quiet Logical value, whether to show downloading progress. Default: FALSE (show).
#'
#' @return Dataframe contains all available datasets.
#' @importFrom magrittr %>%
#' @importFrom jsonlite fromJSON flatten
#' @importFrom data.table rbindlist
#' @importFrom dplyr select filter mutate starts_with
#' @importFrom utils download.file
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' \dontrun{
#' # first time run (lazy mode), need users to provide json folder
#' ucsc.cb.samples <- ShowCBDatasets(lazy = TRUE, update = TRUE)
#' # second time run (lazy mode), need users to provide json folder
#' ucsc.cb.samples <- ShowCBDatasets(lazy = TRUE, update = FALSE)
#' }
ShowCBDatasets <- function(lazy = TRUE, json.folder = NULL, update = FALSE, quiet = FALSE) {
  # parse all datasets json
  if (lazy) {
    if (is.null(json.folder)) {
      json.folder <- getwd()
    }
    base.url <- json.folder
    if (update) {
      message("Lazy mode is on (save time for next time), save all json file to ", json.folder)
      # download top-level json
      utils::download.file(url = "https://cells.ucsc.edu/dataset.json", destfile = file.path(base.url, "dataset.json"), quiet = quiet, mode = "wb")
      # top-level info
      all.datasets.json <- jsonlite::fromJSON(txt = file.path(base.url, "dataset.json"))
      base.url <- "https://cells.ucsc.edu"
    } else {
      message("Lazy mode is on, load downloaded json from ", json.folder)
      json.files <- list.files(path = json.folder, pattern = "json$")
      if (length(json.files) == 0) {
        stop("There is no json file in ", json.folder, " please re-run ShowCBDatasets with update = TRUE to download the json files!")
      }
      all.datasets.json <- jsonlite::fromJSON(txt = file.path(json.folder, "dataset.json"))
      base.url <- json.folder
    }
    all.datasets.df <- jsonlite::flatten(all.datasets.json$datasets)
    colnames(all.datasets.df) <- gsub(pattern = ".*\\.", replacement = "", x = colnames(all.datasets.df))
    # modify elements
    all.datasets.df <- PasteAttr(df = all.datasets.df, attr = c(
      "tags", "diseases", "organisms", "body_parts",
      "projects", "life_stages", "domains", "sources", "assays"
    ))
    # extract all samples
    all.samples.df <- ExtractSample(df = all.datasets.df, base.url = base.url, json.folder = json.folder, quiet = quiet) %>%
      as.data.frame()
    # desc folder
    desc.folder <- json.folder
    datasets.folder <- json.folder
  } else {
    message("Lazy mode is off, always load json online (always up-to-date), this is also always time-consuming.")
    all.datasets.json <- jsonlite::fromJSON(txt = "https://cells.ucsc.edu/dataset.json")
    all.datasets.df <- jsonlite::flatten(all.datasets.json$datasets)
    colnames(all.datasets.df) <- gsub(pattern = ".*\\.", replacement = "", x = colnames(all.datasets.df))
    # modify elements
    all.datasets.df <- PasteAttr(df = all.datasets.df, attr = c(
      "tags", "diseases", "organisms", "body_parts",
      "projects", "life_stages", "domains", "sources", "assays"
    ))
    # extract all samples
    all.samples.df <- ExtractSampleOnline(all.datasets.df) %>% as.data.frame()
    # desc folder
    desc.folder <- "https://cells.ucsc.edu"
    datasets.folder <- "https://cells.ucsc.edu"
  }
  # split columns
  all.sample.dup.index <- duplicated(colnames(all.samples.df))
  # all.sample.dup.cols = colnames(all.samples.df)[all.sample.dup.index]
  # all.samples.df.unique = all.samples.df[!all.sample.dup.index]
  # all.samples.df.dup = all.samples.df[all.sample.dup.index]
  # # modify columns
  # attr.cols = c(
  #   "tags", "diseases", "organisms", "body_parts",
  #   "projects", "life_stages", "domains", "sources", "assays"
  # )
  # unique.valid = intersect(attr.cols, colnames(all.samples.df.unique))
  # all.samples.df.unique <- PasteAttr(df = all.samples.df.unique, attr = unique.valid)
  # dup.valid = intersect(attr.cols, colnames(all.samples.df.dup))
  # all.samples.df.dup <- PasteAttr(df = all.samples.df.dup, attr = dup.valid)
  # # deal with duplicated columns
  # for (col in all.sample.dup.cols){
  #   col.df = data.frame(col1=all.samples.df.unique[, col], col2 = all.samples.df.dup[, col])
  #   col.value = apply(col.df, 1, function(x){
  #     if(x["col1"] == ""){
  #       if(x["col2"] == ""){
  #         x["col1"]
  #       }else{
  #         x["col2"]
  #       }
  #     }else{
  #       x["col1"]
  #     }
  #   })
  #   all.samples.df.unique[, col] = col.value
  # }
  all.samples.df <- all.samples.df[!all.sample.dup.index]
  # modify columns
  all.samples.df <- PasteAttr(df = all.samples.df, attr = c(
    "tags", "diseases", "organisms", "body_parts",
    "projects", "life_stages", "domains", "sources", "assays"
  ))
  # remove unused columns
  all.samples.df <- all.samples.df %>% dplyr::select(-c("md5", "hasFiles", "isCollection", "datasetCount", "collectionCount"))
  # add label
  all.samples.df.single <- all.samples.df %>% dplyr::filter(!grepl(x = .data[["name"]], pattern = "/"))
  all.samples.df.multi <- all.samples.df %>%
    dplyr::filter(grepl(x = .data[["name"]], pattern = "/")) %>%
    dplyr::mutate(parent = gsub(pattern = "([^/]*).*", replacement = "\\1", x = .data[["name"]]))
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
  # add matrix information
  datasets.files <- file.path(datasets.folder, all.samples.final$name, "dataset.json")
  names(datasets.files) <- all.samples.final$name
  datasets.list <- lapply(names(datasets.files), function(x) {
    x.json <- jsonlite::fromJSON(txt = datasets.files[x])
    x.file.df <- jsonlite::flatten(as.data.frame(x.json$fileVersions))
    mat.name <- ifelse("outMatrix.fname" %in% colnames(x.file.df), basename(x.file.df[, "outMatrix.fname"]), "")
    barcode.name <- ifelse("barcodes.fname" %in% colnames(x.file.df), basename(x.file.df[, "barcodes.fname"]), "")
    feature.name <- ifelse("features.fname" %in% colnames(x.file.df), basename(x.file.df[, "features.fname"]), "")
    x.file.df <- data.frame(matrix = mat.name, barcode = barcode.name, feature = feature.name)
    x.file.df$name <- x
    x.file.df
  })
  dts.df <- do.call(rbind, datasets.list)
  all.samples.final <- merge(all.samples.final, dts.df, by = "name") %>% as.data.frame()
  all.samples.final$matrixType <- ifelse(all.samples.final$barcode == "", "matrix", "10x")
  # add sample information
  desc.files <- file.path(desc.folder, all.samples.final$name, "desc.json")
  names(desc.files) <- all.samples.final$name
  desc.list <- lapply(names(desc.files), function(x) {
    sd.json <- jsonlite::fromJSON(txt = desc.files[x])
    used.attr <- c(
      "title" = "title", "paper" = "paper_url", "abstract" = "abstract", "unit" = "unitDesc",
      "coords" = "coordFiles", "methods" = "methods", "geo" = "geo_series"
    )
    sd.df <- ExtractDesc(lst = sd.json, attr = used.attr)
    sd.df$name <- x
    sd.df
  })
  desc.df <- do.call(rbind, desc.list)
  all.samples.final <- merge(all.samples.final, desc.df, by = "name") %>% as.data.frame()
  # return value
  return(all.samples.final)
}


#' Extract UCSC Cell Browser Datasets with Attributes.
#'
#' @param all.samples.df Dataframe contains all samples metadata, obtained with \code{ShowCBDatasets}.
#' @param collection The collection of the datasets, corresponds to \code{shortLabel} column
#' of \code{ShowCBDatasets}, obtain available values with \code{StatDBAttribute}. Default: NULL (without filtering).
#' @param sub.collection The sub-collection of the datasets, corresponds to \code{subLabel} column
#' of \code{ShowCBDatasets}, obtain available values with \code{StatDBAttribute}. Default: NULL (without filtering).
#' @param organ The organ of the datasets, corresponds to \code{body_parts} column
#' of \code{ShowCBDatasets}, obtain available values with \code{StatDBAttribute}. Default: NULL (without filtering).
#' @param disease The disease of the datasets, corresponds to \code{diseases} column
#' of \code{ShowCBDatasets}, obtain available values with \code{StatDBAttribute}. Default: NULL (without filtering).
#' @param organism The specie of the datasets, corresponds to \code{organisms} column
#' of \code{ShowCBDatasets}, obtain available values with \code{StatDBAttribute}. Default: NULL (without filtering).
#' @param project The project of the datasets, corresponds to \code{projects} column
#' of \code{ShowCBDatasets}, obtain available values with \code{StatDBAttribute}. Default: NULL (without filtering).
#' @param fuzzy.match Logical value, whether to perform fuzzy match with provided attribute values. Default: TRUE.
#' @param cell.num Cell number filter. If NULL, no filter; if one value, lower filter; if two values, low and high filter.
#' Deault: NULL(without filtering).
#'
#' @return Dataframe contains filtered datasets.
#' @export
#'
#' @examples
#' \dontrun{
#' # lazy mode, load datasets json files locally, need users to provide json folder
#' ucsc.cb.samples <- ShowCBDatasets(lazy = TRUE, json.folder = NULL, update = FALSE)
#' # cell number is between 1000 and 2000
#' hbb.sample.df <- ExtractCBDatasets(
#'   all.samples.df = ucsc.cb.samples, organ = c("brain", "blood"),
#'   organism = "Human (H. sapiens)", cell.num = c(1000, 2000)
#' )
#' }
ExtractCBDatasets <- function(all.samples.df, collection = NULL, sub.collection = NULL, organ = NULL, disease = NULL, organism = NULL,
                              project = NULL, fuzzy.match = TRUE, cell.num = NULL) {
  # extract row index under different filter
  collection.idx <- CheckParas(df = all.samples.df, column = "shortLabel", para.value = collection, fuzzy.match = fuzzy.match)
  sub.collection.idx <- CheckParas(df = all.samples.df, column = "subLabel", para.value = sub.collection, fuzzy.match = fuzzy.match)
  organ.idx <- CheckParas(df = all.samples.df, column = "body_parts", para.value = organ, fuzzy.match = fuzzy.match)
  disease.idx <- CheckParas(df = all.samples.df, column = "diseases", para.value = disease, fuzzy.match = fuzzy.match)
  organism.idx <- CheckParas(df = all.samples.df, column = "organisms", para.value = organism, fuzzy.match = fuzzy.match)
  project.idx <- CheckParas(df = all.samples.df, column = "projects", para.value = project, fuzzy.match = fuzzy.match)
  if (is.null(cell.num)) {
    cnum.idx <- 1:nrow(all.samples.df)
  } else if (length(cell.num) == 1) {
    cnum.idx <- which(all.samples.df$sampleCount > as.numeric(cell.num))
  } else {
    cnum.idx <- which(all.samples.df$sampleCount > as.numeric(cell.num[1]) & all.samples.df$sampleCount < as.numeric(cell.num[2]))
  }
  # filter on the whole dataset
  valid.idx <- Reduce(intersect, list(collection.idx, sub.collection.idx, organ.idx, disease.idx, organism.idx, project.idx, cnum.idx))
  used.sample.df <- all.samples.df[valid.idx, ]
  rownames(used.sample.df) <- NULL
  return(used.sample.df)
}

#' Extract Cell Type Composition of UCSC Cell Browser Datasets.
#'
#' @param json.folder Folder contains datasets json files, same as \code{json.folder} of \code{ShowCBDatasets}.
#' Default: NULL (current working directory).
#' @param sample.df Dataframe contains used datasets. Default: NULL.
#' @param all.samples.df Dataframe contains all samples metadata, obtained with \code{ShowCBDatasets}. Default: NULL.
#' \code{sample.df} and \code{all.samples.df} cannot be both NULL.
#' @param collection The collection of the datasets, corresponds to \code{shortLabel} column
#' of \code{all.samples.df}, obtain available values with \code{StatDBAttribute}. Default: NULL (without filtering).
#' @param sub.collection The sub-collection of the datasets, corresponds to \code{subLabel} column
#' of \code{all.samples.df}, obtain available values with \code{StatDBAttribute}. Default: NULL (without filtering).
#' @param organ The organ of the datasets, corresponds to \code{body_parts} column
#' of \code{all.samples.df}, obtain available values with \code{StatDBAttribute}. Default: NULL (without filtering).
#' @param disease The disease of the datasets, corresponds to \code{diseases} column
#' of \code{all.samples.df}, obtain available values with \code{StatDBAttribute}. Default: NULL (without filtering).
#' @param organism The specie of the datasets, corresponds to \code{organisms} column
#' of \code{all.samples.df}, obtain available values with \code{StatDBAttribute}. Default: NULL (without filtering).
#' @param project The project of the datasets, corresponds to \code{projects} column
#' of \code{all.samples.df}, obtain available values with \code{StatDBAttribute}. Default: NULL (without filtering).
#' @param fuzzy.match Logical value, whether to perform fuzzy match with provided attribute values. Default: TRUE.
#' @param cell.num Cell number filter. If NULL, no filter; if one value, lower filter; if two values, low and high filter. Deault: NULL.
#'
#' @return Dataframe contains sample information and cell type composition.
#' @importFrom jsonlite fromJSON
#' @importFrom data.table rbindlist
#' @importFrom dplyr select
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' \dontrun{
#' # lazy mode, load datasets json files locally, need users to provide json folder
#' ucsc.cb.samples <- ShowCBDatasets(lazy = TRUE, json.folder = NULL, update = FALSE)
#' # cell number is between 1000 and 2000
#' hbb.sample.df <- ExtractCBDatasets(
#'   all.samples.df = ucsc.cb.samples, organ = c("brain", "blood"),
#'   organism = "Human (H. sapiens)", cell.num = c(1000, 2000)
#' )
#' hbb.sample.ct <- ExtractCBComposition(json.folder = NULL, sample.df = hbb.sample.df)
#' }
ExtractCBComposition <- function(json.folder = NULL, sample.df = NULL, all.samples.df = NULL, collection = NULL, sub.collection = NULL, organ = NULL,
                                 disease = NULL, organism = NULL, project = NULL, fuzzy.match = TRUE, cell.num = NULL) {
  # prepare samples for download
  if (!is.null(sample.df)) {
    # use provided dataframe to download data
    used.sample.df <- sample.df
  } else {
    if (is.null(all.samples.df)) {
      stop("Please provide all.samples.df, obtained with ShowCBDatasets.")
    }
    used.sample.df <- ExtractCBDatasets(
      all.samples.df = all.samples.df, collection = collection, sub.collection = sub.collection, organ = organ,
      disease = disease, organism = organism, project = project, fuzzy.match = fuzzy.match, cell.num = cell.num
    )
  }
  # prepare json folder
  if (is.null(json.folder)) {
    json.folder <- getwd()
  }
  # check json files
  expect.json.files <- file.path(json.folder, used.sample.df$name, "dataset.json")
  available.json.files <- list.files(path = json.folder, pattern = "dataset.json", recursive = TRUE, full.names = TRUE)
  diff.json.files <- setdiff(expect.json.files, available.json.files)
  valid.json.files <- intersect(expect.json.files, available.json.files)
  if (length(diff.json.files) > 0) {
    warning(paste0(diff.json.files, collapse = ","), " json files are not available, please check!")
  }
  if (length(valid.json.files) > 0) {
    cell.type.list <- lapply(valid.json.files, function(x) {
      ct.name <- gsub(pattern = paste0(json.folder, "/", "(.*?)", "/dataset.json"), replacement = "\\1", x = x)
      ct.json <- jsonlite::fromJSON(txt = x)
      valid.label.id <- intersect(c("labelField", "clusterField"), names(ct.json))
      cid <- ct.json[[valid.label.id[1]]]
      meta.fields <- ct.json$metaFields
      cid.idx <- which(meta.fields$label == cid)
      ct.list <- ct.json$metaFields$valCounts[cid.idx][[1]]
      if (is.null(ct.list)) {
        ct.df <- data.frame(CellType = NA, Num = NA)
      } else {
        ct.df <- ct.list %>% as.data.frame()
        colnames(ct.df) <- c("CellType", "Num")
        ct.df$Num <- as.numeric(ct.df$Num)
      }
      ct.df$name <- ct.name
      ct.df
    })
    cell.type.df <- data.table::rbindlist(cell.type.list, fill = TRUE) %>% as.data.frame()
    cell.type.df <- merge(used.sample.df, cell.type.df, by = "name", all.y = TRUE) %>% as.data.frame()
    cell.type.df <- cell.type.df %>%
      dplyr::select(-c("matrix", "barcode", "feature", "matrixType", "unit", "coords", "name")) %>%
      dplyr::relocate(c("CellType", "Num"), .after = "subLabel")
  } else {
    cell.type.df <- data.frame()
  }
  return(cell.type.df)
}

#' Download UCSC Cell Browser Datasets.
#'
#' @param sample.df Dataframe contains used datasets. Default: NULL.
#' @param all.samples.df Dataframe contains all samples metadata, obtained with \code{ShowCBDatasets}. Default: NULL.
#' \code{sample.df} and \code{all.samples.df} cannot be both NULL.
#' @param collection The collection of the datasets, corresponds to \code{shortLabel} column
#' of \code{all.samples.df}, obtain available values with \code{StatDBAttribute}. Default: NULL (without filtering).
#' @param sub.collection The sub-collection of the datasets, corresponds to \code{subLabel} column
#' of \code{all.samples.df}, obtain available values with \code{StatDBAttribute}. Default: NULL (without filtering).
#' @param organ The organ of the datasets, corresponds to \code{body_parts} column
#' of \code{all.samples.df}, obtain available values with \code{StatDBAttribute}. Default: NULL (without filtering).
#' @param disease The disease of the datasets, corresponds to \code{diseases} column
#' of \code{all.samples.df}, obtain available values with \code{StatDBAttribute}. Default: NULL (without filtering).
#' @param organism The specie of the datasets, corresponds to \code{organisms} column
#' of \code{all.samples.df}, obtain available values with \code{StatDBAttribute}. Default: NULL (without filtering).
#' @param project The project of the datasets, corresponds to \code{projects} column
#' of \code{all.samples.df}, obtain available values with \code{StatDBAttribute}. Default: NULL (without filtering).
#' @param fuzzy.match Logical value, whether to perform fuzzy match with provided attribute values. Default: TRUE.
#' @param cell.num Cell number filter. If NULL, no filter; if one value, lower filter; if two values, low and high filter. Deault: NULL.
#' @param timeout Maximum request time when loading data online. Default: 1000.
#' @param merge Logical value, whether to merge Seurat list. Default: FALSE.
#'
#' @return Seurat object (if \code{merge} is TRUE) or list of Seurat objects (if \code{merge} is FALSE).
#' @importFrom magrittr %>%
#' @importFrom data.table fread
#' @importFrom dplyr full_join
#' @importFrom Seurat CreateSeuratObject as.sparse
#' @importFrom purrr reduce
#' @importFrom tibble column_to_rownames
#' @importFrom Matrix readMM
#' @importFrom methods new
#' @export
#'
#' @examples
#' \dontrun{
#' # lazy mode, load datasets json files locally, need users to provide json folder
#' ucsc.cb.samples <- ShowCBDatasets(lazy = TRUE, json.folder = NULL, update = FALSE)
#' # cell number is between 1000 and 2000
#' hbb.sample.df <- ExtractCBDatasets(
#'   all.samples.df = ucsc.cb.samples, organ = c("brain", "blood"),
#'   organism = "Human (H. sapiens)", cell.num = c(1000, 2000)
#' )
#' hbb.sample.seu <- ParseCBDatasets(sample.df = hbb.sample.df)
#' # test 10x and matrix load
#' complex.df <- ucsc.cb.samples[c(1, 927, 379), ] # two 10x and one matrix
#' complex.seu.list <- ParseCBDatasets(sample.df = test.df, merge = F)
#' }
ParseCBDatasets <- function(sample.df = NULL, all.samples.df = NULL, collection = NULL, sub.collection = NULL, organ = NULL,
                            disease = NULL, organism = NULL, project = NULL, fuzzy.match = TRUE, cell.num = NULL,
                            timeout = 1000, merge = TRUE) {
  # prepare samples for download
  if (!is.null(sample.df)) {
    # use provided dataframe to download data
    used.sample.df <- sample.df
  } else {
    if (is.null(all.samples.df)) {
      stop("Please provide all.samples.df, obtained with ShowCBDatasets.")
    }
    used.sample.df <- ExtractCBDatasets(
      all.samples.df = all.samples.df, collection = collection, sub.collection = sub.collection, organ = organ,
      disease = disease, organism = organism, project = project, fuzzy.match = fuzzy.match, cell.num = cell.num
    )
  }
  # set timeout
  env.timeout <- getOption("timeout")
  on.exit(options(timeout = env.timeout)) # restore timeout
  options(timeout = timeout)
  # base url
  base.url <- "https://cells.ucsc.edu/"
  seu.obj.list <- list()
  for (rn in 1:nrow(used.sample.df)) {
    sample.name <- gsub(pattern = "/", replacement = "_", x = used.sample.df[rn, "name"])
    matrix.type <- used.sample.df[rn, "matrixType"]
    # prepare exp matrix
    exp.url <- paste0(base.url, used.sample.df[rn, "name"], "/", used.sample.df[rn, "matrix"])
    # prepare barcode and feature
    if (matrix.type == "matrix") {
      barcode.url <- NULL
      feature.url <- NULL
    } else if (matrix.type == "10x") {
      barcode.url <- paste0(base.url, used.sample.df[rn, "name"], "/", used.sample.df[rn, "barcode"])
      feature.url <- paste0(base.url, used.sample.df[rn, "name"], "/", used.sample.df[rn, "feature"])
    }
    # prepare metadata
    meta.url <- paste0(base.url, used.sample.df[rn, "name"], "/meta.tsv")
    # prepare coord
    coord.file <- paste0(base.url, used.sample.df[rn, "name"], "/", strsplit(x = used.sample.df[rn, "coords"], split = ", ")[[1]])
    seu.obj.list[[rn]] <- Load2Seurat(
      exp.file = exp.url, barcode.url = barcode.url, feature.url = feature.url,
      meta.file = meta.url, coord.file = coord.file, name = sample.name
    )
  }
  # merge or not
  if (isTRUE(merge)) {
    seu.obj <- mergeExperiments(seu.obj.list)
    return(seu.obj)
  } else {
    return(seu.obj.list)
  }
}
