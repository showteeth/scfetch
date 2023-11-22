#' Download Matrix from GEO and Load to Seurat.
#'
#' @param acce GEO accession number.
#' @param platform Platform information/field. Disable when \code{down.supp} is TRUE. Default: NULL (disable).
#' @param down.supp Logical value, whether to download supplementary files to create count matrix. If TRUE, always
#' download supplementary files. If FALSE, use \code{ExpressionSet} (If contains non-integer or empty,
#' download supplementary files automatically). Default: FALSE.
#' @param supp.idx The index of supplementary files to download. This should be consistent with \code{platform}. Default: 1.
#' @param timeout Timeout for \code{\link{download.file}}. Default: 3600.
#' @param data.type The data type of the dataset, choose from "sc" (single-cell) and "bulk" (bulk). Default: "sc".
#' @param supp.type The type of downloaded supplementary files, choose from count (count matrix file or single count matrix file)
#' and 10x (cellranger output files, contains barcodes, genes/features and matrix). Default: count.
#' @param out.folder Output folder to save 10x files. Default: NULL (current working directory).
#' @param gene2feature Logical value, whether to rename \code{genes.tsv.gz} to \code{features.tsv.gz}. Default: TRUE.
#' @param merge Logical value, whether to merge Seurat list when there are multiple 10x files (\code{supp.type} is 10x). Default: FALSE.
#' @param ... Parameters for \code{\link{getGEO}}.
#'
#' @return If \code{data.type} is "sc", return Seurat object (if \code{merge} is TRUE) or Seurat object list (if \code{merge} is FALSE).
#' If \code{data.type} is "bulk", return count matrix.
#' @importFrom magrittr %>%
#' @importFrom GEOquery getGEO getGEOSuppFiles gunzip
#' @importFrom Biobase annotation experimentData pData phenoData notes sampleNames exprs
#' @importFrom tools file_ext
#' @importFrom utils untar
#' @importFrom data.table fread
#' @importFrom openxlsx read.xlsx
#' @importFrom Seurat Read10X CreateSeuratObject
#' @importFrom methods new
#' @export
#'
#' @examples
#' \dontrun{
#' # the supp files are count matrix
#' GSE94820.seu <- ParseGEO(acce = "GSE94820", down.supp = TRUE, supp.idx = 1, supp.type = "count")
#' # the supp files are cellranger output files: barcodes, genes/features and matrix
#' # need users to provide the output folder
#' GSE200257.seu <- ParseGEO(
#'   acce = "GSE200257", down.supp = TRUE, supp.idx = 1, supp.type = "10x",
#'   out.folder = "/path/to/output/folder"
#' )
#' }
ParseGEO <- function(acce, platform = NULL, down.supp = FALSE, supp.idx = 1, timeout = 3600, data.type = c("sc", "bulk"),
                     supp.type = c("count", "10x"), out.folder = NULL, gene2feature = TRUE, merge = TRUE, ...) {
  # check parameters
  data.type <- match.arg(arg = data.type)
  supp.type <- match.arg(arg = supp.type)
  # check platform
  if (down.supp) {
    message("Download supplementary files to generate matrix!")
    pf.obj <- NULL
  } else {
    message("Extract expression data from eSets!")
    if (is.null(platform)) {
      stop("Platform is required to extract expression data!")
    }
    # get GEO object
    pf.obj <- GEOobj(acce = acce, platform = platform, ...)
  }
  # change supp type to count when bulk
  if (data.type == "bulk") {
    supp.type <- "count"
  }
  # extract counts matrix
  pf.count <- ExtractGEOExp(
    pf.obj = pf.obj, acce = acce, supp.idx = supp.idx, down.supp = down.supp,
    timeout = timeout, supp.type = supp.type, out.folder = out.folder, gene2feature = gene2feature
  )
  if (data.type == "bulk") {
    return(pf.count)
  } else if (data.type == "sc") {
    # load seurat
    if (is.null(pf.count) && supp.type == "10x") {
      message("Loading data to Seurat!")
      all.samples.folder <- dir(out.folder, full.names = TRUE)
      # check file
      valid.samples.folder <- Check10XFiles(folders = all.samples.folder, gene2feature = gene2feature)
      if (length(valid.samples.folder) == 0) {
        stop("No valid sample folder detected under ", out.folder, ". Please check!")
      }
      # load to seurat
      seu.list <- sapply(valid.samples.folder, function(x) {
        x.mat <- Seurat::Read10X(data.dir = x)
        seu.obj <- Seurat::CreateSeuratObject(counts = x.mat, project = basename(x))
        seu.obj
      })
      if (isTRUE(merge)) {
        seu.obj <- mergeExperiments(seu.list)
      } else {
        seu.obj <- seu.list
      }
    } else if (!is.null(pf.count) && supp.type == "count") {
      seu.obj <- Seurat::CreateSeuratObject(counts = pf.count, project = acce)
    }
    return(seu.obj)
  }
}

#' Extract Sample Metadata from GEO.
#'
#' @param acce GEO accession number.
#' @param platform Platform information/field. Default: NULL (all platforms).
#' @param ... Parameters for \code{\link{getGEO}}.
#'
#' @return Dataframe contains all metadata of provided GEO accession number.
#' @importFrom magrittr %>%
#' @importFrom GEOquery getGEO
#' @importFrom Biobase annotation experimentData pData phenoData notes sampleNames exprs
#' @export
#'
#' @examples
#' \donttest{
#' # users may need to set the size of the connection buffer
#' # Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 60)
#' # extract metadata of specified platform
#' GSE200257.meta <- ExtractGEOMeta(acce = "GSE200257", platform = "GPL24676")
#' }
ExtractGEOMeta <- function(acce, platform = NULL, ...) {
  # get GEO object
  if (is.null(platform)) {
    geo.obj <- GEOquery::getGEO(GEO = acce, ...)
    pfs <- sapply(geo.obj, function(x) {
      Biobase::annotation(x)
    })
    names(geo.obj) <- pfs
  } else {
    geo.obj <- list()
    pf.obj <- GEOobj(acce = acce, platform = platform, ...)
    geo.obj[[platform]] <- pf.obj
  }
  # extract metadata
  geo.meta.list <- lapply(names(geo.obj), function(x) {
    # extract general information
    pf.info <- ExtractGEOInfo(pf.obj = geo.obj[[x]], sample.wise = FALSE)
    pf.info$Platform <- x
    # select meta data
    pf.meta <- ExtractGEOSubMeta(pf.obj = geo.obj[[x]])
    pf.meta$Platform <- x
    # merge all dataframe
    pf.all <- merge(pf.meta, pf.info, by = "Platform", all.x = TRUE) %>% as.data.frame()
    pf.all
  })
  # all metadta
  if (length(geo.meta.list) > 1) {
    geo.meta.df <- data.table::rbindlist(geo.meta.list, fill = TRUE) %>% as.data.frame()
  } else {
    geo.meta.df <- geo.meta.list[[1]]
  }
  return(geo.meta.df)
}

# connect to GEO, extract GEO object, extract platform object
GEOobj <- function(acce, platform, ...) {
  # obtain GEO object
  geo.obj <- GEOquery::getGEO(GEO = acce, ...)

  # extract platform
  pfs <- sapply(geo.obj, function(x) {
    Biobase::annotation(x)
  })
  if (!platform %in% pfs) {
    stop(paste("The platform you provides is not valid!", paste(pfs, collapse = ", ")))
  }
  # extract platform data
  pf.idx <- which(pfs == platform[1])
  pf.obj <- geo.obj[[pf.idx]]
  return(pf.obj)
}

# merge cols into one
SimpleCol <- function(df, col) {
  cols <- grep(pattern = col, x = colnames(df), value = T)
  df.cols <- df[cols]
  value <- paste(c(t(df.cols)), collapse = ". ")
  return(value)
}

#' Extract GEO Study Information.
#'
#' @param pf.obj GEO object of platform.
#' @param sample.wise Logical value, whether to extract sample-wise information. Default: FALSE.
#'
#' @return A dataframe.
#'
ExtractGEOInfo <- function(pf.obj, sample.wise = FALSE) {
  # platform information
  pf.info <- Biobase::experimentData(pf.obj)
  # additional information
  pf.info.add.df <- as.data.frame(Biobase::pData(Biobase::phenoData(pf.obj)))
  if (sample.wise) {
    pf.info.final <- cbind(data.frame(
      Title = pf.info@title,
      Type = Biobase::notes(pf.info)$type,
      Abstract = pf.info@abstract,
      Design = pf.info@other$overall_design,
      SampleCount = length(Biobase::sampleNames(pf.obj)),
      SupplementaryFile = gsub(pattern = "\n", replacement = ", ", pf.info@other$supplementary_file),
      PMID = gsub(pattern = "\n", replacement = ", ", pf.info@pubMedIds), pf.info.add.df
    )) %>%
      as.data.frame()
  } else {
    used.cols <- c("organism", "molecule", "strategy", "extract_protocol", "data_processing")
    pf.info.add.used <- pf.info.add.df[, grep(pattern = paste(used.cols, collapse = "|"), colnames(pf.info.add.df), value = T)]
    pf.info.add.sim <- apply(pf.info.add.used, 2, function(x) {
      paste(unique(x), collapse = ", ")
    }) %>%
      t() %>%
      as.data.frame()
    # final information
    pf.info.final <- data.frame(
      Title = pf.info@title,
      Type = Biobase::notes(pf.info)$type,
      Organism = SimpleCol(df = pf.info.add.sim, col = "organism"),
      Abstract = pf.info@abstract,
      Design = pf.info@other$overall_design,
      SampleCount = length(Biobase::sampleNames(pf.obj)),
      Molecule = SimpleCol(df = pf.info.add.sim, col = "molecule"),
      ExtractProtocol = SimpleCol(df = pf.info.add.sim, col = "extract_protocol"),
      LibraryStrategy = SimpleCol(df = pf.info.add.sim, col = "strategy"),
      DataProcessing = SimpleCol(df = pf.info.add.sim, col = "data_processing"),
      SupplementaryFile = gsub(pattern = "\n", replacement = ", ", pf.info@other$supplementary_file),
      Contact = paste(pf.info@name, pf.info@contact, sep = "; "),
      PMID = gsub(pattern = "\n", replacement = ", ", pf.info@pubMedIds)
    )
  }
  return(pf.info.final)
}

#' Extract Sample Metadata.
#'
#' @param pf.obj GEO object of platform.
#'
#' @return A dataframe.
#'
ExtractGEOSubMeta <- function(pf.obj) {
  # extract sample detail information
  pf.info <- as.data.frame(Biobase::pData(Biobase::phenoData(pf.obj)))
  # select used basic cols
  valid.cols <- intersect(colnames(pf.info), c(c("title", "geo_accession", "source_name_ch1", "description")))
  pf.info.used <- pf.info[valid.cols]
  # process characteristics
  pf.info.charac <- pf.info[grep(pattern = "^characteristics", x = colnames(pf.info))]
  ## modify colnames
  pf.info.charac.colnames <-
    unique(apply(pf.info.charac, 2, function(x) {
      gsub(pattern = "(.*?): (.*)", replacement = "\\1", x = x)
    }))
  colnames(pf.info.charac) <- pf.info.charac.colnames
  ## modify values
  pf.info.charac <- apply(pf.info.charac, 2, function(x) {
    gsub(pattern = "(.*?): (.*)", replacement = "\\2", x = x)
  })
  ## final meta
  pf.meta <- cbind(pf.info.used, pf.info.charac) %>% as.data.frame()

  return(pf.meta)
}


#' Extract Raw Count Matrix from Supplementary Files.
#'
#' @param acce GEO accession number.
#' @param timeout Timeout for \code{\link{download.file}}. Default: 3600.
#' @param supp.idx The index of supplementary files to download. Default: 1.
#'
#' @return A dataframe.
#'
ExtractGEOExpSupp <- function(acce, timeout = 3600, supp.idx = 1) {
  # create tmp folder
  tmp.folder <- tempdir()
  # get current timeout
  if (!is.null(timeout)) {
    message("Change Timeout to: ", timeout)
    env.timeout <- getOption("timeout")
    on.exit(options(timeout = env.timeout)) # restore timeout
    options(timeout = timeout)
  }
  # download supplementary file
  # supp.down.log <- GEOquery::getGEOSuppFiles(GEO = acce, baseDir = tmp.folder)
  supp.down.log <- tryCatch(
    expr = {
      GEOquery::getGEOSuppFiles(GEO = acce, baseDir = tmp.folder)
    },
    error = function(e) {
      print(e)
      stop("You can change the timeout with a larger value.")
    }
  )
  if (supp.idx > nrow(supp.down.log)) {
    stop("Please provide valid supplementary file index.")
  }
  supp.file.path <- rownames(supp.down.log)[supp.idx]
  # remove unused supplementary file
  unused.supp <- setdiff(rownames(supp.down.log), supp.file.path)
  unused.supp.remove <- file.remove(unused.supp)
  # file unzip
  file.ext <- tools::file_ext(supp.file.path)
  if (file.ext == "gz") {
    # gunzip file
    GEOquery::gunzip(supp.file.path, overwrite = TRUE)
    supp.file.path <- gsub(pattern = "\\.gz", replacement = "", x = supp.file.path)
    # update file extension
    file.ext <- tools::file_ext(supp.file.path)
  }
  if (file.ext == "tar") {
    # untar
    utils::untar(supp.file.path, exdir = file.path(tmp.folder, acce, "sample"))
    # unzip
    unzip.log <- sapply(
      list.files(file.path(tmp.folder, acce, "sample"), full.names = TRUE, pattern = "gz$"),
      function(x) {
        GEOquery::gunzip(x, overwrite = TRUE)
      }
    )
    # read files
    count.list <- lapply(
      list.files(file.path(tmp.folder, acce, "sample"), full.names = TRUE),
      function(x) {
        sample.count <- data.table::fread(file = x) %>% as.data.frame()
        colnames(sample.count) <- c("GeneName", gsub(pattern = "(GSM[0-9]*).*", replacement = "\\1", x = basename(x)))
        sample.count
      }
    )
    # create count matrix
    count.mat <- Reduce(f = function(x, y) {
      merge.mat <- merge(x, y, by = "GeneName", all = T)
    }, x = count.list)
    rownames(count.mat) <- count.mat$GeneName
    count.mat$GeneName <- NULL
  } else {
    if (file.ext %in% c("xlsx", "xls")) {
      # read excel file
      count.mat <- openxlsx::read.xlsx(xlsxFile = supp.file.path, rowNames = TRUE)
    } else if (file.ext %in% c("csv", "tsv", "txt")) {
      # read text file
      count.mat <- data.table::fread(file = supp.file.path) %>% as.data.frame()
      # the first column must be gene
      rownames(count.mat) <- count.mat[, 1]
      count.mat[, 1] <- NULL
    }
  }
  return(count.mat)
}

#' Fortmat Supplementary Files to 10x.
#'
#' @param acce GEO accession number.
#' @param supp.idx The index of supplementary files to download. Default: 1.
#' @param timeout Timeout for \code{\link{download.file}}. Default: 3600.
#' @param out.folder Output folder to save 10x files. Default: NULL (current working directory).
#' @param gene2feature Logical value, whether to rename \code{genes.tsv.gz} to \code{features.tsv.gz}.
#' Default: TURE.
#'
#' @return NULL
#'
ExtractGEOExpSupp10x <- function(acce, supp.idx = 1, timeout = 3600,
                                 out.folder = NULL, gene2feature = TRUE) {
  # create tmp folder
  tmp.folder <- tempdir()
  # get current timeout
  if (!is.null(timeout)) {
    message("Change Timeout to: ", timeout)
    env.timeout <- getOption("timeout")
    on.exit(options(timeout = env.timeout)) # restore timeout
    options(timeout = timeout)
  }
  # download supp file
  supp.down.log <- tryCatch(
    expr = {
      GEOquery::getGEOSuppFiles(GEO = acce, baseDir = tmp.folder)
    },
    error = function(e) {
      print(e)
      stop("You can change the timeout with a larger value.")
    }
  )
  # check supp.idx
  if (supp.idx > nrow(supp.down.log)) {
    stop("Please provide valid supplementary file index.")
  }
  # get used supplementary file
  supp.file.path <- rownames(supp.down.log)[supp.idx]
  # remove unused supplementary file
  unused.supp <- setdiff(rownames(supp.down.log), supp.file.path)
  unused.supp.remove <- file.remove(unused.supp)
  # file unzip
  file.ext <- tools::file_ext(supp.file.path)
  if (file.ext == "gz") {
    # gunzip file
    GEOquery::gunzip(supp.file.path, overwrite = TRUE)
    supp.file.path <- gsub(pattern = "\\.gz", replacement = "", x = supp.file.path)
    # update file extension
    file.ext <- tools::file_ext(supp.file.path)
  }
  if (file.ext == "tar") {
    # untar
    utils::untar(supp.file.path, exdir = file.path(tmp.folder, acce, "sample"))
    # recognize valid files: barcodes.tsv.gz, genes.tsv.gz, matrix.mtx.gz and features.tsv.gz
    valid.pat <- "barcodes.tsv.gz$|genes.tsv.gz$|matrix.mtx.gz$|features.tsv.gz$"
    all.files <- list.files(file.path(tmp.folder, acce, "sample"), full.names = TRUE, pattern = valid.pat)
    # change file name
    if (gene2feature) {
      change.name.log <- sapply(all.files, function(x) {
        if (grepl(pattern = "genes.tsv.gz$", x = x)) {
          new.name <- gsub(pattern = "genes.tsv.gz$", replacement = "features.tsv.gz", x = x)
          file.rename(from = x, to = new.name)
        }
      })
      all.files <- list.files(file.path(tmp.folder, acce, "sample"), full.names = TRUE, pattern = valid.pat)
    }
    # prepare out folder
    if (is.null(out.folder)) {
      out.folder <- getwd()
    }
    # get folder
    all.sample.folder <- sapply(all.files, function(x) {
      # get basename and dirname
      file.name <- basename(x)
      dir.name <- dirname(x)
      # remove file type tag
      file.name <- gsub(pattern = valid.pat, replacement = "", x = file.name)
      # remove possible _ and .
      file.name <- gsub(pattern = "[_.]$", replacement = "", x = file.name)
      file.folder <- file.path(out.folder, file.name)
    })
    # create folder and move file
    move.file.log <- sapply(all.files, function(x) {
      # get folder name
      folder.name <- all.sample.folder[x]
      # create folder
      if (!dir.exists(folder.name)) {
        dir.create(path = folder.name, recursive = TRUE)
      }
      new.file.name <- gsub(pattern = paste0(".*(barcodes.tsv.gz$|genes.tsv.gz$|matrix.mtx.gz$|features.tsv.gz$)"), replacement = "\\1", x = x)
      # move file
      copy.tag <- file.copy(from = x, to = file.path(folder.name, new.file.name))
      # remove the original file
      remove.tag <- file.remove(x)
      copy.tag
    })
    message("Process 10x fiels done! All files are in ", out.folder)
  } else {
    stop("Does not support non-tar file for 10x mode!")
  }
}


#' Extract Raw Count Matrix from Supplementary Files or Fortmat Supplementary Files to 10x.
#'
#' @param acce GEO accession number.
#' @param supp.idx The index of supplementary files to download. Default: 1.
#' @param timeout Timeout for \code{\link{download.file}}. Default: 3600.
#' @param supp.type The type of downloaded supplementary files, choose from count (count matrix file or single count matrix file)
#' and 10x (cellranger output files, contains barcodes, genes/features and matrix). Default: count.
#' @param out.folder Output folder to save 10x files. Default: NULL (current working directory).
#' @param gene2feature Logical value, whether to rename \code{genes.tsv.gz} to \code{features.tsv.gz}. Default: TRUE.
#' Default: TURE.
#'
#' @return Count matrix (\code{supp.type} is count) or NULL (\code{supp.type} is 10x).
#'
ExtractGEOExpSuppAll <- function(acce, supp.idx = 1, timeout = 3600,
                                 supp.type = c("count", "10x"), out.folder = NULL, gene2feature = TRUE) {
  if (supp.type == "count") {
    count.mat <- ExtractGEOExpSupp(acce = acce, supp.idx = supp.idx, timeout = timeout)
    return(count.mat)
  } else if (supp.type == "10x") {
    ExtractGEOExpSupp10x(acce = acce, supp.idx = supp.idx, timeout = timeout, out.folder = out.folder, gene2feature = gene2feature)
    return(NULL)
  }
}

#' Extract Raw Count Matrix or Fortmat Supplementary Files to 10x.
#'
#' @param pf.obj GEO object of platform.
#' @param acce GEO accession number.
#' @param supp.idx The index of supplementary files to download. Default: 1.
#' @param down.supp Logical value, whether to download supplementary files to create count matrix. If TRUE, always
#' download supplementary files. If FALSE, use \code{ExpressionSet} (If contains non-integer or emoty,
#' download supplementary files automatically). Default: FALSE.
#' @param timeout Timeout for \code{\link{download.file}}. Default: 3600.
#' @param supp.type The type of downloaded supplementary files, choose from count (count matrix file or single count matrix file)
#' and 10x (cellranger output files, contains barcodes, genes/features and matrix). Default: count.
#' @param out.folder Output folder to save 10x files. Default: NULL (current working directory).
#' @param gene2feature Logical value, whether to rename \code{genes.tsv.gz} to \code{features.tsv.gz}. Default: TRUE.
#'
#' @return Count matrix (\code{supp.type} is count) or NULL (\code{supp.type} is 10x).
#'
ExtractGEOExp <- function(pf.obj, acce, supp.idx = 1, down.supp = FALSE, timeout = 3600,
                          supp.type = c("count", "10x"), out.folder = NULL, gene2feature = TRUE) {
  # check parameters
  supp.type <- match.arg(arg = supp.type)
  # download supplementary files
  if (down.supp) {
    exp.data <- ExtractGEOExpSuppAll(
      acce = acce, supp.idx = supp.idx, timeout = timeout,
      supp.type = supp.type, out.folder = out.folder, gene2feature = gene2feature
    )
  } else {
    expr.mat <- Biobase::exprs(pf.obj)
    if (nrow(expr.mat) == 0) {
      message("Matrix not available! Downloading supplementary files.")
      exp.data <- ExtractGEOExpSuppAll(
        acce = acce, supp.idx = supp.idx, timeout = timeout,
        supp.type = supp.type, out.folder = out.folder, gene2feature = gene2feature
      )
    } else {
      if (all(expr.mat %% 1 == 0)) {
        exp.data <- expr.mat
      } else {
        message("Matrix contains non-integer values! Downloading supplementary files.")
        exp.data <- ExtractGEOExpSuppAll(
          acce = acce, supp.idx = supp.idx, timeout = timeout,
          supp.type = supp.type, out.folder = out.folder, gene2feature = gene2feature
        )
      }
    }
  }
  return(exp.data)
}
