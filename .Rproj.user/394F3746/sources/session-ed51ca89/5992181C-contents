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
ExtractSample <- function(df) {
  base.url <- "https://cells.ucsc.edu/"
  if (!"isCollection" %in% colnames(df)) {
    return(df)
  } else {
    cf <- df[!is.na(df$isCollection), ]
    sf <- df[is.na(df$isCollection), ]
    cu <- paste0(base.url, cf$name, "/dataset.json")
    cul <- lapply(cu, function(x) {
      x.json <- jsonlite::fromJSON(txt = x)
      x.df <- jsonlite::flatten(x.json$datasets)
      colnames(x.df) <- gsub(pattern = ".*\\.", replacement = "", x = colnames(x.df))
      x.df
    })
    cu.df <- data.table::rbindlist(cul, fill = TRUE)
    # df = data.table::rbindlist(list(sf, cu.df), fill = TRUE)
    # return(list(sf, ExtractSample(cu.df)))
    return(data.table::rbindlist(list(sf, ExtractSample(cu.df)), fill = TRUE))
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
Load2Seurat <- function(exp.file, meta.file, coord.file = NULL, name = NULL) {
  # source: https://cellbrowser.readthedocs.io/en/master/load.html
  # read matrix
  mat <- data.table::fread(exp.file, check.names = FALSE)
  # read metadata
  meta <- data.frame(data.table::fread(meta.file, check.names = FALSE), row.names = 1)
  # get genes
  genes <- gsub(".+[|]", "", mat[, 1][[1]])
  # modify mat genes
  mat <- data.frame(mat[, -1], row.names = genes, check.names = FALSE)
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
