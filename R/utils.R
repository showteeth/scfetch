# merge Seurat object: https://github.com/dosorio/rPanglaoDB/blob/master/R/mergeExperiments.R
mergeExperiments <- function(experimentList) {
  for (i in seq_along(experimentList)[-1]) {
    experimentList[[1]] <- suppressWarnings(merge(experimentList[[1]], experimentList[[i]]))
    experimentList[[i]] <- new("Seurat")
  }
  experimentList <- experimentList[[1]]
  return(experimentList)
}

# check columns existence
CheckColumns <- function(df, columns) {
  if (!all(columns %in% colnames(df))) {
    miss.cols <- setdiff(columns, colnames(df))
    stop(paste0(paste(miss.cols, collapse = ", "), " does not exist, Please Check!"))
  }
}

# used in UCSCCellBrowser, merge multiple attributes
PasteAttr <- function(df, attr) {
  for (at in attr) {
    df[[at]] <- sapply(df[[at]], function(x) {
      paste0(x, collapse = ", ")
    })
  }
  return(df)
}

# used in UCSCCellBrowser, recursively extract samples
ExtractSample <- function(df, base.url, json.folder, quiet) {
  # prepare json
  if (base.url != json.folder) {
    df.json <- file.path(base.url, df$name, "dataset.json")
    names(df.json) <- df$name
    df.json.folder <- file.path(json.folder, df$name)
    names(df.json.folder) <- df$name
    df.desc <- file.path(base.url, df$name, "desc.json")
    names(df.desc) <- df$name
    dird <- sapply(df.json.folder, function(x) {
      dir.create(x, showWarnings = FALSE, recursive = TRUE)
    })
    down.status <- lapply(df$name, function(x) {
      utils::download.file(url = df.json[x], destfile = file.path(df.json.folder[x], "dataset.json"), quiet = quiet, mode = "wb", method = "wget")
      utils::download.file(url = df.desc[x], destfile = file.path(df.json.folder[x], "desc.json"), quiet = quiet, mode = "wb", method = "wget")
    })
  }
  # process
  if (!"isCollection" %in% colnames(df)) {
    return(df)
  } else {
    cf <- df[!is.na(df$isCollection), ]
    sf <- df[is.na(df$isCollection), ]
    cu.json.folder <- file.path(json.folder, cf$name)
    cul <- lapply(file.path(cu.json.folder, "dataset.json"), function(x) {
      x.json <- jsonlite::fromJSON(txt = x)
      x.df <- jsonlite::flatten(x.json$datasets)
      colnames(x.df) <- gsub(pattern = ".*\\.", replacement = "", x = colnames(x.df))
      x.df
    })
    cu.df <- data.table::rbindlist(cul, fill = TRUE)
    # df = data.table::rbindlist(list(sf, cu.df), fill = TRUE)
    # return(list(sf, ExtractSample(cu.df)))
    return(data.table::rbindlist(list(sf, ExtractSample(df = cu.df, base.url = base.url, json.folder = json.folder, quiet = quiet)), fill = TRUE))
  }
}

# used in UCSCCellBrowser, recursively extract samples online
ExtractSampleOnline <- function(df) {
  base.url <- "https://cells.ucsc.edu/"
  if (!"isCollection" %in% colnames(df)) {
    return(df)
  } else {
    cf <- df[!is.na(df$isCollection), ]
    sf <- df[is.na(df$isCollection), ]
    cu <- file.path(base.url, cf$name, "dataset.json")
    cul <- lapply(cu, function(x) {
      x.json <- jsonlite::fromJSON(txt = x)
      x.df <- jsonlite::flatten(x.json$datasets)
      colnames(x.df) <- gsub(pattern = ".*\\.", replacement = "", x = colnames(x.df))
      x.df
    })
    cu.df <- data.table::rbindlist(cul, fill = TRUE)
    # df = data.table::rbindlist(list(sf, cu.df), fill = TRUE)
    # return(list(sf, ExtractSample(cu.df)))
    return(data.table::rbindlist(list(sf, ExtractSampleOnline(cu.df)), fill = TRUE))
  }
}

# used in UCSCCellBrowser, inherit attributes from parents
InheritParient <- function(df, attr) {
  for (at in attr) {
    df[[at]] <- ifelse(df[[at]] == "", ifelse(df[[paste0("parent_", at)]] == "", "", paste0(df[[paste0("parent_", at)]], "|parent")), df[[at]])
  }
  return(df)
}

# used in UCSCCellBrowser, extract sample attribute
ExtractDesc <- function(lst, attr) {
  at.list <- list()
  for (atn in names(attr)) {
    at <- attr[atn]
    if (at %in% names(lst)) {
      at.value <- paste0(lst[[at]], collapse = ", ")
    } else {
      at.value <- ""
    }
    at.list[[atn]] <- at.value
  }
  return(as.data.frame(at.list))
}

# used in UCSCCellBrowser, filter attributes and return index
CheckParas <- function(df, column, para.value, fuzzy.match = TRUE) {
  # convert to lower case to avoid case-sensitive
  all.values <- gsub("\\|parent$", replacement = "", x = tolower(df[[column]]))
  if (is.null(para.value)) {
    message("Use all ", column, " as input!")
    # return row index
    value <- 1:nrow(df)
  } else {
    para.value <- tolower(para.value)
    # deal with fuzzy match
    if (fuzzy.match) {
      value.list <- sapply(para.value, function(x) grep(pattern = x, x = all.values, fixed = TRUE))
      value.list.len <- sapply(value.list, function(x) length(x))
      invalid.value <- names(value.list.len[value.list.len == 0])
      value <- unique(unlist(value.list))
    } else {
      if (column %in% c("body_parts", "diseases", "organisms", "projects")) {
        # value contains dot
        value.list <- list()
        for (pv in para.value) {
          pv.vec <- c()
          for (avi in 1:length(all.values)) {
            av <- all.values[avi]
            av.vec <- strsplit(x = av, split = ", ")[[1]]
            if (pv %in% av.vec) {
              pv.vec <- c(pv.vec, avi)
            }
          }
          value.list[[pv]] <- pv.vec
        }
        # get invalid value
        invalid.value <- setdiff(para.value, names(value.list))
        value <- unique(unlist(value.list))
      } else {
        # value doesn't contain dot
        value.list <- lapply(para.value, function(x) {
          which(all.values %in% x)
        })
        names(value.list) <- para.value
        value.list.len <- sapply(value.list, function(x) length(x))
        invalid.value <- names(value.list.len[value.list.len == 0])
        value <- unique(unlist(value.list))
      }
    }
    # print invalid value
    if (length(invalid.value) > 0) {
      message(paste0(invalid.value, collapse = ", "), " are not valid for ", column)
    }
    # deal with empty value
    if (length(value) == 0) {
      warning("There is no valid value under ", paste(para.value, collapse = ", "), " in ", column, ". This filter is invalid!")
      value <- 1:nrow(df)
    }
  }
  return(value)
}

# used in UCSCCellBrowser, create seurat object (add coord to metadata)
Load2Seurat <- function(exp.file, barcode.url = NULL, feature.url = NULL,
                        meta.file, coord.file = NULL, name = NULL) {
  # source: https://cellbrowser.readthedocs.io/en/master/load.html
  # read matrix
  if (is.null(barcode.url)) {
    mat <- data.table::fread(exp.file, check.names = FALSE)
    # get genes
    genes <- gsub(".+[|]", "", mat[, 1][[1]])
    # modify mat genes
    mat <- data.frame(mat[, -1], row.names = genes, check.names = FALSE)
  } else {
    # with default parameters
    mat <- Read10XOnline(matrix.url = exp.file, barcode.url = barcode.url, feature.url = feature.url)
  }
  # read metadata
  meta <- data.frame(data.table::fread(meta.file, check.names = FALSE), row.names = 1)
  if (is.null(coord.file)) {
    seu.obj <- Seurat::CreateSeuratObject(counts = mat, project = name, meta.data = meta)
  } else {
    # prepare coord file
    coord.list <- lapply(1:length(coord.file), function(x) {
      coord.name <- gsub(pattern = ".coords.tsv.gz", replacement = "", x = basename(coord.file[x]))
      coord.df <- data.frame(data.table::fread(coord.file[x], check.names = FALSE), row.names = 1)
      colnames(coord.df) <- paste(coord.name, 1:ncol(coord.df), sep = "_")
      coord.df$Barcode <- rownames(coord.df)
      return(coord.df)
    })
    if (length(coord.file) == 1) {
      all.coord.df <- coord.list[[1]]
    } else {
      all.coord.df <- coord.list %>% purrr::reduce(dplyr::full_join, by = "Barcode")
    }
    # merge metadata
    meta <- merge(meta, all.coord.df, by.x = 0, by.y = "Barcode", all.x = TRUE) %>%
      tibble::column_to_rownames(var = "Row.names")
    seu.obj <- Seurat::CreateSeuratObject(counts = mat, project = name, meta.data = meta)
  }
  return(seu.obj)
}

# used in UCSCCellBrowser
# source: https://github.com/satijalab/seurat/blob/master/R/utilities.R#L1949
ExtractField <- function(string, field = 1, delim = "_") {
  fields <- as.numeric(x = unlist(x = strsplit(x = as.character(x = field), split = ",")))
  if (length(x = fields) == 1) {
    return(strsplit(x = string, split = delim)[[1]][field])
  }
  return(paste(strsplit(x = string, split = delim)[[1]][fields], collapse = delim))
}

# used in UCSCCellBrowser, load cellranger output to matrix
Read10XOnline <- function(matrix.url, barcode.url, feature.url, gene.column = 2,
                          cell.column = 1, unique.features = TRUE, strip.suffix = FALSE) {
  # load matrix
  data <- Matrix::readMM(file = gzcon(url(matrix.url)))
  # load barcode
  cell.barcodes <- as.data.frame(data.table::fread(barcode.url, header = FALSE))
  cn <- ifelse(ncol(x = cell.barcodes) > 1, cell.column, 1)
  cell.names <- cell.barcodes[, cn]
  if (all(grepl(pattern = "\\-1$", x = cell.names)) & strip.suffix) {
    cell.names <- as.vector(x = as.character(x = sapply(
      X = cell.names,
      FUN = ExtractField,
      field = 1,
      delim = "-"
    )))
  }
  # matrix colnames
  colnames(x = data) <- cell.names
  # load feature
  feature.names <- as.data.frame(data.table::fread(feature.url, header = FALSE))
  # modify gene column
  gene.column <- min(ncol(feature.names), gene.column)
  if (any(is.na(x = feature.names[, gene.column]))) {
    warning(
      "Some features names are NA. Replacing NA names with ID from the opposite column requested",
      call. = FALSE,
      immediate. = TRUE
    )
    na.features <- which(x = is.na(x = feature.names[, gene.column]))
    replacement.column <- ifelse(test = gene.column == 2, yes = 1, no = 2)
    feature.names[na.features, gene.column] <- feature.names[na.features, replacement.column]
  }

  # modify matrix rownames
  if (unique.features) {
    rownames(x = data) <- make.unique(names = feature.names[, gene.column])
  } else {
    rownames(x = data) <- feature.names[, gene.column]
  }
  # In cell ranger 3.0, a third column specifying the type of data was added
  # and we will return each type of data as a separate matrix
  if (ncol(x = feature.names) > 2) {
    data_types <- factor(x = feature.names$V3)
    lvls <- levels(x = data_types)
    if (length(x = lvls) > 1 && length(x = full.data) == 0) {
      message("10X data contains more than one type and is being returned as a list containing matrices of each type.")
    }
    expr_name <- "Gene Expression"
    if (expr_name %in% lvls) { # Return Gene Expression first
      lvls <- c(expr_name, lvls[-which(x = lvls == expr_name)])
    }
    data <- lapply(
      X = lvls,
      FUN = function(l) {
        return(data[data_types == l, , drop = FALSE])
      }
    )
    names(x = data) <- lvls
  } else {
    data <- list(data)
  }
  # convert to dgCMatrix
  final.data <- Seurat::as.sparse(data[[1]])
  return(final.data)
}
