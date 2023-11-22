#' Export SeuratObject to Other Formats.
#'
#' @param seu.obj A seurat object.
#' @param assay Which assay to use. Default: NULL (get with \code{\link{DefaultAssay}}).
#' @param reduction Name of DimReduc to set to main reducedDim in cds.
#' @param to The target format, chosen from "SCE" (SingleCellExperiment), "AnnData", "CellDataSet", "cell_data_set", "loom".
#' Default: "SCE".
#' @param anndata.file File used to save AnnData results. Default: NULL.
#' @param loom.file File used to save loom results. Default: NULL.
#' @param conda.path Conda environment path, used when \code{to} is "AnnData". Default: NULL.
#' @param ... Parameter for \code{\link{as.SingleCellExperiment}}, \code{sceasy::convertFormat}, \code{\link{as.CellDataSet}},
#' \code{as.cell_data_set}, \code{SaveLoom}, corresponding to \code{to}.
#'
#' @return Object corresponding to \code{to}.
#' @export
#' @importFrom Seurat DefaultAssay as.SingleCellExperiment as.CellDataSet
#' @importFrom reticulate use_condaenv
#' @importFrom methods is
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' # export to SingleCellExperiment
#' sce.obj <- ExportSeurat(seu.obj = pbmc_small, assay = "RNA", to = "SCE")
#' # export to CellDataSet
#' cds.obj <- ExportSeurat(seu.obj = pbmc_small, assay = "RNA", reduction = "tsne", to = "CellDataSet")
#' # export to cell_data_set
#' cds3.obj <- ExportSeurat(seu.obj = pbmc_small, assay = "RNA", to = "cell_data_set")
#' # export to AnnData, need users to provide the conda path and the output file
#' ExportSeurat(
#'   seu.obj = pbmc_small, assay = "RNA", to = "AnnData", conda.path = "/path/to/anaconda3",
#'   anndata.file = "/path/to/pbmc_small.h5ad"
#' )
#' # export to loom, need users to provide the output file
#' ExportSeurat(
#'   seu.obj = pbmc_small, assay = "RNA", to = "loom",
#'   loom.file = "/path/to/pbmc_small.loom"
#' )
#' }
ExportSeurat <- function(seu.obj, assay = NULL, reduction = NULL,
                         to = c("SCE", "AnnData", "CellDataSet", "cell_data_set", "loom"),
                         anndata.file = NULL, loom.file = NULL, conda.path = NULL, ...) {
  # check parameters
  to <- match.arg(arg = to)
  # check object
  if (!methods::is(seu.obj, "Seurat")) {
    stop("Please provide valid Seurat object!")
  }

  # get default assay
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(object = seu.obj)
  }

  # convert SeuratObject to other formats
  if (to == "SCE") {
    message("Convert SeuratObject to SingleCellExperiment (suitable for scater)!")
    result.obj <- Seurat::as.SingleCellExperiment(seu.obj, assay = assay, ...)
    return(result.obj)
  } else if (to == "AnnData") {
    message("Convert SeuratObject to AnnData (suitable for scanpy)!")
    # prepare anndata file
    if (is.null(anndata.file)) {
      warning("You do not provide anndata file, generate automatically in current working directory!")
      seu.name <- deparse(substitute(seu.obj))
      anndata.file <- file.path(getwd(), paste0(seu.name, ".h5ad"))
    }
    # h5Seurat.file <- gsub(pattern = ".h5ad$", replacement = ".h5seurat", x = anndata.file)
    # SeuratDisk::SaveH5Seurat(seu.obj, overwrite = overwrite, filename = h5Seurat.file)
    # SeuratDisk::Convert(source = h5Seurat.file, assay = assay, dest = "h5ad", overwrite = overwrite, ...)
    # # remove h5Seurat file
    # file.remove(h5Seurat.file)
    # set conda env
    # python environment
    if (!is.null(conda.path)) {
      reticulate::use_condaenv(conda.path, required = TRUE)
    }
    if (!requireNamespace("sceasy", quietly = TRUE)) {
      stop("Can not find sceasy package, install with devtools::install_github('cellgeni/sceasy')!")
    }
    sceasy::convertFormat(seu.obj,
      from = "seurat", to = "anndata", drop_single_values = F,
      assay = assay, outFile = anndata.file, ...
    )
  } else if (to == "CellDataSet") {
    message("Convert SeuratObject to CellDataSet (suitable for Monocle)!")
    result.obj <- Seurat::as.CellDataSet(x = seu.obj, assay = assay, reduction = reduction, ...)
    if (tolower(reduction) != "tsne") {
      # tsne stored in reducedDimA, others stored in reducedDimS and transposed
      # bug of as.CellDataSet: the dimension is transposed
      # https://github.com/satijalab/seurat/blob/57564e282a1ba54c655192d907d83df21164dc78/R/objects.R#L851
      result.obj@reducedDimS <- t(result.obj@reducedDimS)
    }
    return(result.obj)
  } else if (to == "cell_data_set") {
    message("Convert SeuratObject to cell_data_set (suitable for Monocle3)!")
    if (!requireNamespace("SeuratWrappers", quietly = TRUE)) {
      stop("Can not find SeuratWrappers package, install with devtools::install_github('satijalab/seurat-wrappers')!")
    }
    if (is.null(reduction)) {
      result.obj <- SeuratWrappers::as.cell_data_set(x = seu.obj, assay = assay, ...)
    } else {
      result.obj <- SeuratWrappers::as.cell_data_set(x = seu.obj, assay = assay, reduction = reduction, ...)
    }
    return(result.obj)
  } else if (to == "loom") {
    message("Convert SeuratObject to loom!")
    # prepare loom file
    if (is.null(loom.file)) {
      warning("You do not provide loom file, generate automatically in current working directory!")
      seu.name <- deparse(substitute(seu.obj))
      loom.file <- file.path(getwd(), paste0(seu.name, ".loom"))
    }
    # Convert SingleCellExperiment to loom
    if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
      stop("Can not find SeuratDisk package, install with devtools::install_github('mojaveazure/seurat-disk')!")
    }
    SeuratDisk::SaveLoom(object = seu.obj, filename = loom.file, overwrite = TRUE, ...)
  }
}

#' Convert Other Formats to SeuratObject.
#'
#' @param obj Other formats object (eg: SingleCellExperiment, CellDataSet).
#' Default: NULL (used when \code{from} is "AnnData").
#' @param assay Assay name to store expression matrices in SeuratObject. Default: RNA.
#' @param from The source formats, chosen from "SCE" (SingleCellExperiment), "AnnData", "CellDataSet", "cell_data_set". Default: "SCE".
#' @param count.assay The assay of source formats to save raw counts,
#' used when \code{from} is "SCE" or cell_data_set. Default: counts.
#' @param data.assay The assay of source formats to save log transformed counts,
#' used when \code{from} is "SCE" or cell_data_set. Default: logcounts.
#' @param slot Slot to store expression data as, used when \code{from} is "CellDataSet". Default: counts.
#' @param anndata.file The file contains AnnData. Default: NULL.
#' @param loom.file The file contains loom. Default: NULL.
#' @param conda.path Conda environment path, used when \code{from} is "AnnData". Default: NULL.
#' @param ... Parameter for \code{\link{as.Seurat}}, \code{sceasy::convertFormat}, \code{\link{as.Seurat}}, \code{\link{as.Seurat}},
#' \code{\link{as.Seurat}}, corresponding to \code{from}.
#'
#' @return A Seurat object.
#' @importFrom Seurat as.Seurat
#' @importFrom reticulate use_condaenv
#' @importFrom SummarizedExperiment assayNames
#' @importFrom scater logNormCounts
#' @importFrom methods is
#' @export
#'
#' @examples
#' \dontrun{
#' # import data from SingleCellExperiment
#' seu.obj <- ImportSeurat(
#'   obj = sce.obj, from = "SCE", count.assay = "counts",
#'   data.assay = "logcounts", assay = "RNA"
#' )
#' # import data from CellDataSet
#' seu.obj <- ImportSeurat(obj = cds.obj, from = "CellDataSet", count.assay = "counts", assay = "RNA")
#' # import data from cell_data_set
#' seu.obj <- ImportSeurat(
#'   obj = sce.obj, from = "cell_data_set", count.assay = "counts",
#'   data.assay = "logcounts", assay = "RNA"
#' )
#' # import data from AnnData, need users to provide the file for conversion
#' seu.obj <- ImportSeurat(anndata.file = "path/to/h5ad", from = "AnnData", assay = "RNA")
#' # import data from loom, need users to provide the file for conversion
#' seu.obj <- ImportSeurat(loom.file = "path/to/loom", from = "loom")
#' }
ImportSeurat <- function(obj = NULL, assay = "RNA", from = c("SCE", "AnnData", "CellDataSet", "cell_data_set", "loom"),
                         count.assay = "counts", data.assay = "logcounts", slot = "counts",
                         anndata.file = NULL, loom.file = NULL, conda.path = NULL, ...) {
  # check parameters
  from <- match.arg(arg = from)

  # convert other formats to SeuratObject
  if (from == "SCE") {
    message("Convert SingleCellExperiment to SeuratObject!")
    # check object
    if (is.null(obj)) {
      stop("Please provide SingleCellExperiment with obj!")
    } else if (!methods::is(obj, "SingleCellExperiment")) {
      stop("Please provide valid SingleCellExperiment object!")
    }
    # check assays
    ava.assays <- SummarizedExperiment::assayNames(obj)
    if (!count.assay %in% ava.assays) {
      stop("There is no count assay, check object or count.assay!")
    }
    if (!data.assay %in% ava.assays) {
      message("There is no data assay (data.assay), create with logNormCounts!")
      obj <- scater::logNormCounts(obj)
      data.assay <- "logcounts"
    }
    # convert
    seu.obj <- Seurat::as.Seurat(
      x = obj, counts = count.assay, data = data.assay, assay = assay, ...
    )
    # change reductions names
    names(seu.obj@reductions) <- tolower(names(seu.obj@reductions))
  } else if (from == "AnnData") {
    message("Convert AnnData to SeuratObject!")
    if (is.null(anndata.file)) {
      stop("Please provide a file a to AnnData results.")
    } else {
      # SeuratDisk::Convert(source = anndata.file, assay = assay, dest = "h5seurat", overwrite = overwrite, ...)
      # h5Seurat.file <- gsub(pattern = ".h5ad$", replacement = ".h5seurat", x = anndata.file)
      # seu.obj <- SeuratDisk::LoadH5Seurat(h5Seurat.file, assay = assay)
      # python environment
      if (!is.null(conda.path)) {
        reticulate::use_condaenv(conda.path, required = TRUE)
      }
      if (!requireNamespace("sceasy", quietly = TRUE)) {
        stop("Can not find sceasy package, install with devtools::install_github('cellgeni/sceasy')!")
      }
      seu.obj <- sceasy::convertFormat(anndata.file,
        from = "anndata", to = "seurat",
        assay = assay, outFile = NULL, ...
      )
    }
  } else if (from == "CellDataSet") {
    message("Convert CellDataSet (Monocle) to SeuratObject!")
    # check object
    if (is.null(obj)) {
      stop("Please provide CellDataSet with obj!")
    } else if (!methods::is(obj, "CellDataSet")) {
      stop("Please provide valid CellDataSet object!")
    }
    # convert
    seu.obj <- Seurat::as.Seurat(x = obj, slot = slot, assay = assay, ...)
    # change reductions names
    names(seu.obj@reductions) <- tolower(names(seu.obj@reductions))
  } else if (from == "cell_data_set") {
    message("Convert cell_data_set (Monocle3) to SeuratObject!")
    # check object
    if (is.null(obj)) {
      stop("Please provide cell_data_set with obj!")
    } else if (!methods::is(obj, "cell_data_set")) {
      stop("Please provide valid cell_data_set object!")
    }
    # convert
    seu.obj <- Seurat::as.Seurat(
      x = obj, counts = count.assay, data = data.assay,
      assay = assay, ...
    )
    # change reductions names
    names(seu.obj@reductions) <- tolower(names(seu.obj@reductions))
  } else if (from == "loom") {
    message("Convert loom to SeuratObject!")
    if (is.null(loom.file)) {
      stop("Please provide a file a to loom results.")
    } else {
      if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
        stop("Can not find SeuratDisk package, install with devtools::install_github('mojaveazure/seurat-disk')!")
      }
      loom.info <- SeuratDisk::Connect(filename = loom.file, mode = "r")
      seu.obj <- Seurat::as.Seurat(loom.info, ...)
    }
  }
  return(seu.obj)
}

#' Data Format Conversion between SingleCellExperiment and AnnData.
#'
#' @param from The source data format to convert, chosen from SingleCellExperiment and AnnData.
#' Default: SingleCellExperiment.
#' @param to The target data format to convert, chosen from AnnData and SingleCellExperiment.
#' Default: AnnData.
#' @param sce The SingleCellExperiment object to convert. Default: NULL.
#' @param anndata.file File used to save or contains AnnData results. Default: NULL.
#' @param slot.name Slot name used to save count matrix, used when converting from AnnData to SingleCellExperiment.
#' Default: counts.
#' @param ... Parameters for \code{writeH5AD} and \code{readH5AD}.
#'
#' @return NULL or SingleCellExperiment.
#' @importFrom SingleCellExperiment reducedDimNames reducedDimNames<-
#' @importFrom methods is
#' @export
#'
#' @examples
#' \dontrun{
#' library(scRNAseq)
#' seger <- SegerstolpePancreasData()
#' SCEAnnData(from = "SingleCellExperiment", to = "AnnData", sce = seger, X_name = "counts")
#' # need users to provide the output file
#' sce <- SCEAnnData(
#'   from = "AnnData", to = "SingleCellExperiment",
#'   anndata.file = "path/to/seger.h5ad"
#' )
#' }
SCEAnnData <- function(from = c("SingleCellExperiment", "AnnData"),
                       to = c("AnnData", "SingleCellExperiment"),
                       sce = NULL, anndata.file = NULL, slot.name = "counts", ...) {
  # check parameters
  from <- match.arg(arg = from)
  to <- match.arg(arg = to)
  if (!requireNamespace("zellkonverter", quietly = TRUE)) {
    stop("Can not find zellkonverter package, install with BiocManager::install('zellkonverter')!")
  }
  # conversion
  if (from == "SingleCellExperiment" & to == "AnnData") {
    message("Convert SingleCellExperiment to AnnData.")
    # check SingleCellExperiment
    if (is.null(sce)) {
      stop("Please provide SingleCellExperiment with sce!")
    } else if (!methods::is(sce, "SingleCellExperiment")) {
      stop("Please provide valid SingleCellExperiment object!")
    }
    # check h5ad file
    if (is.null(anndata.file)) {
      warning("You do not provide h5ad file, generate automatically in current working directory!")
      sce.name <- deparse(substitute(sce))
      anndata.file <- file.path(getwd(), paste0(sce.name, ".h5ad"))
    }
    # Convert SingleCellExperiment to AnnData
    zellkonverter::writeH5AD(sce = sce, file = anndata.file, ...)
    return(NULL)
  } else if (from == "AnnData" & to == "SingleCellExperiment") {
    message("Convert AnnData to SingleCellExperiment.")
    if (is.null(anndata.file)) {
      stop("Please provide AnnData with anndata.file!")
    }
    # Convert AnnData to SingleCellExperiment
    sce <- zellkonverter::readH5AD(file = anndata.file, ...)
    # change slot name
    names(sce@assays) <- slot.name
    SingleCellExperiment::reducedDimNames(sce) <- toupper(gsub(pattern = "X_", replacement = "", x = SingleCellExperiment::reducedDimNames(sce)))
    return(sce)
  } else {
    stop(paste0("Invalid conversion from ", from, " to ", to, "."))
  }
}

#' Data Format Conversion between SingleCellExperiment and loom.
#'
#' @param from The source data format to convert, chosen from SingleCellExperiment and AnnData.
#' Default: SingleCellExperiment.
#' @param to The target data format to convert, chosen from AnnData and SingleCellExperiment.
#' Default: loom.
#' @param sce The SingleCellExperiment object to convert. Default: NULL.
#' @param loom.file File used to save or contains loom results. Default: NULL.
#' @param ... Parameters for \code{sceasy::convertFormat} and \code{sceasy::convertFormat}.
#'
#' @return NULL or SingleCellExperiment.
#' @importFrom LoomExperiment SingleCellLoomExperiment export import
#' @importFrom SummarizedExperiment assayNames
#' @importFrom methods as is
#' @export
#'
#' @examples
#' \dontrun{
#' # convert from loom to SingleCellExperiment, need users to provide the loom file
#' sce.obj <- SCELoom(
#'   from = "loom", to = "SingleCellExperiment",
#'   loom.file = "path/to/loom"
#' )
#' # convert from SingleCellExperiment to loom, need users to provide the loom file
#' SCELoom(
#'   from = "SingleCellExperiment", to = "loom", sce = sce.obj,
#'   loom.file = "path/to/loom"
#' )
#' }
SCELoom <- function(from = c("SingleCellExperiment", "loom"),
                    to = c("loom", "SingleCellExperiment"),
                    sce = NULL, loom.file = NULL, ...) {
  # check parameters
  from <- match.arg(arg = from)
  to <- match.arg(arg = to)
  # convension
  if (from == "SingleCellExperiment" & to == "loom") {
    message("Convert SingleCellExperiment to loom.")
    # check SingleCellExperiment
    if (is.null(sce)) {
      stop("Please provide SingleCellExperiment with sce!")
    } else if (!methods::is(sce, "SingleCellExperiment")) {
      stop("Please provide valid SingleCellExperiment object!")
    }
    # check loom file
    if (is.null(loom.file)) {
      warning("You do not provide loom file, generate automatically in current working directory!")
      sce.name <- deparse(substitute(sce))
      loom.file <- file.path(getwd(), paste0(sce.name, ".loom"))
    }
    if (file.exists(loom.file)) {
      message(loom.file, " exists!")
    } else {
      # Convert SingleCellExperiment to loom
      # sceasy::convertFormat(sce, from = "sce", to = "loom", outFile = loom.file, ...)
      suppressMessages(suppressWarnings(sce2loom_internal(obj = sce, outFile = loom.file, ...)))
    }
  } else if (from == "loom" & to == "SingleCellExperiment") {
    message("Convert loom to SingleCellExperiment.")
    if (is.null(loom.file)) {
      stop("Please provide loom with loom.file!")
    }
    # Convert loom to SingleCellExperiment
    # sce <- sceasy::convertFormat(loom.file, from = "loom", to = "sce", ...)
    sce <- loom2sce_internal(inFile = loom.file, ...)
    return(sce)
  } else {
    stop(paste0("Invalid conversion from ", from, " to ", to, "."))
  }
}


# modified from https://github.com/cellgeni/sceasy/blob/master/R/functions.R to save rownames and colnames
sce2loom_internal <- function(obj, outFile, ...) {
  if (!requireNamespace("LoomExperiment")) {
    stop("This function requires the 'LoomExperiment' package.")
  }
  scle <- LoomExperiment::SingleCellLoomExperiment(obj)
  if (!is.null(outFile)) {
    LoomExperiment::export(
      scle, outFile,
      matrix = SummarizedExperiment::assayNames(scle)[1],
      colnames_attr = "CellID", rownames_attr = "Gene", ...
    )
  }
}

# modified from https://github.com/cellgeni/sceasy/blob/master/R/functions.R to identify rownames and colnames
loom2sce_internal <- function(inFile, ...) {
  if (!requireNamespace("LoomExperiment")) {
    stop("This function requires the 'LoomExperiment' package.")
  }
  if (!requireNamespace("SingleCellExperiment")) {
    stop("This function requires the 'SingleCellExperiment' package.")
  }
  scle <- LoomExperiment::import(inFile, rownames_attr = "Gene", colnames_attr = "CellID", ...)
  sce <- methods::as(scle, "SingleCellExperiment")
  return(sce)
}
