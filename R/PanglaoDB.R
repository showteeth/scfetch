#' Show Metadata of scRNA-seq Datasets in PanglaoDB.
#'
#' @param specie The specie of the datasets, choose from "Homo sapiens", "Mus musculus", one or multiple value. Default: NULL (All).
#' @param protocol Protocol used to generate the datasets, choose from "10x chromium", "drop-seq", "microwell-seq",
#' "C1 Fluidigm", "inDrops", "Smart-seq2", "CEL-seq", one or multiple value. Default: NULL (All).
#' @param tissue The tissue of the datasets. Default: NULL (All).
#' @param show.cell.type Logical value, whether to show inferred cell type. Default: TRUE.
#'
#' @return Dataframe contains SRA, SRS, Tissue, Protocol, Specie, Cells, CellType (inferred).
#' @importFrom magrittr %>%
#' @importFrom rPanglaoDB getSampleList getSampleComposition
#' @importFrom dplyr filter
#' @export
#'
#' @examples
#' # human.meta = ShowPanglaoDBMeta(specie = "Homo sapiens", protocol = c("Smart-seq2", "10x chromium"))
ShowPanglaoDBMeta <- function(specie = NULL, protocol = NULL, tissue = NULL, show.cell.type = TRUE) {
  # get all sample metadata
  all.meta <- rPanglaoDB::getSampleList()
  # modify SMART-seq2 to Smart-seq2
  all.meta$Protocol <- gsub(pattern = "SMART-seq2", replacement = "Smart-seq2", x = all.meta$Protocol)

  # no attribute filter
  if (is.null(specie) && is.null(protocol) && is.null(tissue)) {
    used.meta <- all.meta
  }
  # test specie
  if (!is.null(specie)) {
    valid.specie <- intersect(specie, unique(all.meta$Species))
    if (length(valid.specie) == 0) {
      warning("Please provide valid specie! The returned dataframe contains metadata without specie filtering, please check!")
      used.meta <- all.meta
    } else {
      used.meta <- all.meta %>% dplyr::filter(Species %in% valid.specie)
    }
  } else {
    used.meta <- all.meta
  }
  # test protocol
  if (!is.null(protocol)) {
    valid.protocol <- intersect(protocol, unique(used.meta$Protocol))
    if (length(valid.protocol) == 0) {
      warning("Please provide valid protocol! The returned dataframe contains metadata without protocol filtering, please check!")
    } else {
      used.meta <- used.meta %>% dplyr::filter(Protocol %in% valid.protocol)
    }
  }
  # test tissue
  if (!is.null(tissue)) {
    valid.tissue <- intersect(tissue, unique(used.meta$Tissue))
    if (length(valid.tissue) == 0) {
      warning("Please provide valid tissue! The returned dataframe contains metadata without tissue filtering, please check!")
    } else {
      used.meta <- used.meta %>% dplyr::filter(Tissue %in% valid.tissue)
    }
  }
  # get sample cell type
  if (show.cell.type) {
    if (nrow(used.meta) > 0) {
      used.meta.ct <- sapply(used.meta$SRS, function(x) {
        tryCatch(
          {
            x.com <- rPanglaoDB::getSampleComposition(srs = x, verbose = FALSE)
            paste(unique(x.com[["Cell Type"]]), collapse = ", ")
          },
          error = function(e) {
            message("This is an error when obtaining inferred cell types of ", x)
            "None"
          }
        )
      })
      used.meta.ct.df <- as.data.frame(used.meta.ct)
      colnames(used.meta.ct.df) <- "CellType"
      used.meta <- merge(used.meta, used.meta.ct.df, by.x = "SRS", by.y = 0, all.x = TRUE)
      used.meta <- used.meta[c("SRA", "SRS", "Tissue", "Protocol", "Species", "Cells", "CellType")]
    }
  }
  # replace srs with '=' with notused
  used.meta$SRS <- gsub(pattern = "nSRS=[0-9]*", replacement = "notused", x = used.meta$SRS)
  return(used.meta)
}

#' Parse PanglaoDB Data.
#'
#' @param meta Metadata contains "SRA", "SRS", "Tissue", "Protocol", "Species", can be obtained with \code{ShowPanglaoDBMeta}.
#' @param cell.type Extract samples with specified cell types. For samples without SRS (notused), this value can only be "All" or "None", or
#' these samples will be filtered. Default: "All".
#' @param include.gene Include cells expressing the genes. Default: NA.
#' @param exclude.gene Exclude cells expressing the genes. Default: NA.
#' @param merge Logical value, whether to merge Seurat list. Default: FALSE.
#'
#' @return Seurat object (if \code{merge} is TRUE) or list of Seurat objects (if \code{merge} is FALSE).
#' @importFrom magrittr %>%
#' @importFrom rPanglaoDB getSampleList getSampleComposition getSamples
#' @importFrom dplyr filter
#' @importFrom pbapply pbapply
#' @importFrom Matrix rowSums colMeans
#' @importFrom Seurat CreateSeuratObject
#' @importFrom utils read.table
#' @export
#'
#' @examples
#' # hsa.meta = ShowPanglaoDBMeta(specie = "Homo sapiens", protocol = c("Smart-seq2", "10x chromium"), show.cell.type = TRUE)
#' # hsa.seu = ParsePanglaoDB(hsa.meta, merge = TRUE)
ParsePanglaoDB <- function(meta, cell.type = "All", include.gene = NA, exclude.gene = NA, merge = FALSE) {
  # check columns
  CheckColumns(df = meta, columns = c("SRA", "SRS", "Tissue", "Protocol", "Species"))

  # split meta to save time
  meta.normal <- meta %>% dplyr::filter(SRS != "notused")
  meta.abnormal <- meta %>% dplyr::filter(SRS == "notused")

  # get Seurat object
  ## normal
  if (nrow(meta.normal) > 0) {
    message("Processing ", nrow(meta.normal), " samples with SRS!")
    normal.seu <- rPanglaoDB::getSamples(
      sra = meta.normal$SRA, srs = meta.normal$SRS, tissue = meta.normal$Tissue,
      protocol = meta.normal$Protocol, specie = meta.normal$Species, celltype = cell.type,
      include = include.gene, exclude = exclude.gene, merge = FALSE
    )
  } else {
    normal.seu <- list()
  }
  ## abnormal
  if (nrow(meta.abnormal) > 0) {
    if (cell.type != "All" && cell.type != "None") {
      message("There is no ", cell.type, " in samples without SRS (notused)!")
      abnormal.seu <- list()
    } else {
      message("Processing ", nrow(meta.abnormal), " samples without SRS (notused)!")
      abnormal.seu <- getSamples_internal(
        sra = meta.abnormal$SRA, srs = meta.abnormal$SRS, tissue = meta.abnormal$Tissue,
        protocol = meta.abnormal$Protocol, specie = meta.abnormal$Species, celltype = cell.type,
        include = include.gene, exclude = exclude.gene, merge = FALSE
      )
    }
  } else {
    abnormal.seu <- list()
  }

  # merge list
  meta.seu <- c(normal.seu, abnormal.seu)
  # merge or not
  if (merge) {
    meta.seu <- mergeExperiments(meta.seu)
  }
  return(meta.seu)
}

# get samples without SRS
getSamples_internal <- function(sra = "All", srs = "All", tissue = "All", protocol = "All", specie = "All", celltype = "All", include = NA, exclude = NA, merge = TRUE) {
  # SampleList
  sampleList <- rPanglaoDB::getSampleList()
  sampleList$SRS <- gsub(pattern = "nSRS=[0-9]*", replacement = "notused", x = sampleList$SRS)

  # Filters
  SRA <- match.arg(arg = sra, choices = unique(c("All", sampleList$SRA)), several.ok = TRUE)
  if (isTRUE("All" %in% SRA)) {
    SRA <- unique(sampleList$SRA)
  }
  SRS <- match.arg(arg = srs, choices = unique(c("All", "notused", sampleList$SRS)), several.ok = TRUE)
  if (isTRUE("All" %in% SRS)) {
    SRS <- unique(sampleList$SRS)
  }

  Tissue <- match.arg(arg = tissue, choices = unique(c("All", sampleList$Tissue)), several.ok = TRUE)
  if (isTRUE("All" %in% Tissue)) {
    Tissue <- unique(sampleList$Tissue)
  }
  Protocol <- match.arg(arg = protocol, choices = unique(c("All", sampleList$Protocol)), several.ok = TRUE)
  if (isTRUE("All" %in% Protocol)) {
    Protocol <- unique(sampleList$Protocol)
  }
  Specie <- match.arg(arg = specie, choices = unique(c("All", sampleList$Species)), several.ok = TRUE)
  if (isTRUE("All" %in% Specie)) {
    Specie <- unique(sampleList$Species)
  }

  # Applying filter
  F1 <- sampleList$SRA %in% SRA
  F2 <- sampleList$SRS %in% SRS
  F3 <- sampleList$Tissue %in% Tissue
  F4 <- sampleList$Protocol %in% Protocol
  F5 <- sampleList$Species %in% Specie
  sampleList <- sampleList[F1 & F2 & F3 & F4 & F5, ]

  # Error
  if (nrow(sampleList) == 0) {
    message("0 Samples Found")
    return()
  }

  # Filtering by cell-type
  ctList.raw <- lapply(1:nrow(sampleList), function(x) {
    if (sampleList[x, "SRS"] == "notused") {
      sc <- sampleList[x, ]
      sc[["Cluster"]] <- "None"
      sc[["Cell Type"]] <- "None"
      sc[c("SRA", "SRS", "Tissue", "Protocol", "Species", "Cluster", "Cells", "Cell Type")]
    } else {
      rPanglaoDB::getSampleComposition(srs = sampleList[x, "SRS"], verbose = FALSE)
    }
  })
  ctList <- do.call(rbind, ctList.raw)
  rownames(ctList) <- NULL

  CellType <- match.arg(arg = celltype, choices = unique(c("All", ctList$`Cell Type`)), several.ok = TRUE)
  if (isTRUE("All" %in% CellType)) {
    CellType <- unique(ctList$`Cell Type`)
  }
  ctList <- ctList[ctList$`Cell Type` %in% CellType, ]
  sampleList <- sampleList[sampleList$SRS %in% ctList$SRS, ]

  # Error
  if (nrow(sampleList) == 0) {
    message("0 Samples Found")
    return()
  }

  dataSets <- pbapply::pbapply(sampleList, 1, function(X) {
    oConnection <- paste0("https://panglaodb.se/data_dl.php?sra=", X[1], "&srs=", X[2], "&filetype=R&datatype=readcounts")
    oConnection <- url(oConnection, headers = list(
      `Connection` = "keep-alive",
      `User-Agent` = "Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:93.0) Gecko/20100101 Firefox/93.0"
    ))
    try(load(oConnection), silent = TRUE)
    if (exists("sm")) {
      rownames(sm) <- unlist(lapply(strsplit(rownames(sm), "-ENS|_ENS"), function(X) {
        X[1]
      }))
      sm <- sm[Matrix::rowSums(sm) > 0, ]
      if (X[2] == "notused") {
        sm <- suppressWarnings(Seurat::CreateSeuratObject(sm, project = paste(as.character(X[1]), as.character(X[2]), sep = "_")))
        cellTypes <- "None"
        cClusters <- NA
      } else {
        cNames <- rPanglaoDB::getSampleComposition(srs = as.character(X[2]), verbose = FALSE)
        rownames(cNames) <- cNames$Cluster
        tempFile <- tempfile()
        download.file(paste0("https://panglaodb.se/data_dl.php?sra=", X[1], "&srs=", X[2], "&datatype=clusters&filetype=txt"),
          destfile = tempFile, quiet = TRUE,
          headers = list(
            `Connection` = "keep-alive",
            `User-Agent` = "Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:93.0) Gecko/20100101 Firefox/93.0"
          )
        )
        cClusters <- utils::read.table(tempFile, row.names = 1)
        sm <- sm[, colnames(sm) %in% rownames(cClusters)]
        cClusters <- cClusters[colnames(sm), ]
        names(cClusters) <- colnames(sm)
        # Capital gene names to allow integration across Human and Mice
        rownames(sm) <- toupper(rownames(sm))
        sm <- suppressWarnings(Seurat::CreateSeuratObject(sm, project = as.character(X[2])))
        cellTypes <- cNames[as.character(cClusters), ]$`Cell Type`
        names(cellTypes) <- colnames(sm)
      }
      sm$Sample <- sm$orig.ident
      sm$CellTypes <- cellTypes
      sm$panglaoCluster <- as.character(cClusters)
      sm$Tissue <- X[["Tissue"]]
      sm <- subset(sm, cells = colnames(sm)[sm$CellTypes %in% CellType])
      sm$CellTypes[sm$CellTypes %in% "Unknown"] <- NA
      sm$CellTypes[sm$CellTypes %in% "None"] <- NA
      sm$Specie <- X[["Species"]]

      # Filtering by genes
      include <- include[include %in% rownames(sm@assays$RNA@counts)]
      exclude <- exclude[exclude %in% rownames(sm@assays$RNA@counts)]
      filterCells <- FALSE
      if (length(include) > 0) {
        include <- Matrix::colMeans(sm@assays$RNA@counts[include, , drop = FALSE] != 0) == 1
        filterCells <- TRUE
      } else {
        include <- rep(TRUE, ncol(sm))
      }
      if (length(exclude) > 0) {
        exclude <- Matrix::colMeans(sm@assays$RNA@counts[exclude, , drop = FALSE] != 0) != 0
        filterCells <- TRUE
      } else {
        exclude <- rep(FALSE, ncol(sm))
      }
      if (isTRUE(filterCells)) {
        filterCells <- (include & !exclude)
        if (any(filterCells)) {
          sm <- subset(sm, cells = colnames(sm)[filterCells])
        } else {
          sm <- list()
        }
      }
      close.connection(oConnection)
    } else {
      sm <- list()
    }
    return(sm)
  })
  names(dataSets) <- ifelse(sampleList$SRS == "notused", paste(sampleList$SRA, sampleList$SRS, sep = "_"), sampleList$SRS)

  dataSets <- dataSets[unlist(lapply(dataSets, class)) %in% "Seurat"]

  if (isTRUE(merge)) {
    dataSets <- mergeExperiments(dataSets)
  }
  return(dataSets)
}
