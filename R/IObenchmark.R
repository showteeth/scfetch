#' Convert SeuratObject to AnnData using SeuratDisk/sceasy/scDIOR.
#'
#' @param seu.obj A seurat object.
#' @param method Method used to perform conversion, choose from "SeuratDisk", "sceasy",	"scDIOR". Default: "SeuratDisk".
#' @param out.folder Output folder. Default: NULL (current working directory).
#' @param out.filename Output file name, e.g., test.h5ad. Default: NULL (\code{seu.obj} name, method, h5ad.).
#' @param assay Assay name to store expression matrices in SeuratObject. Default: RNA.
#' @param slot Slot for adata.X, used when \code{method} is "sceasy". Default: counts.
#' @param save.scale Logical value, whether to save scaled count matrix, used when \code{method} is "SeuratDisk", "scDIOR".
#' Default: FALSE. When \code{method} is "SeuratDisk" and \code{save.scale} is TRUE (FALSE), scaled data is stored in X，
#' log-normalized data is in raw.X (log-normalized data is stored in X，raw count matrix is in raw.X); When \code{method}
#' is "scDIOR" and \code{save.scale} is TRUE (FALSE), scaled data is stored in X，log-normalized data is in raw.X,
#' raw count matrix is in layers (log-normalized data is stored in X，raw count matrix is in layers).
#' @param conda.path Conda environment path, used when \code{method} is "sceasy". Default: NULL.
#'
#' @return Run log.
#' @export
#' @importFrom Seurat GetAssayData
#' @importFrom SeuratDisk SaveH5Seurat Convert
#' @importFrom sceasy convertFormat
#' @importFrom dior write_h5
#'
#' @examples
#' # SeuratDisk
#' Seu2AD(seu.obj = pbmc3k.final, method = "SeuratDisk", out.folder = "benchmark", assay = "RNA", save.scale = TRUE)
#' # sceasy
#' Seu2AD(seu.obj = pbmc3k.final, method = "sceasy", out.folder = "benchmark", assay = "RNA", slot = "counts", conda.path = "/path/to/conda")
#' # scDIOR
#' Seu2AD(seu.obj = pbmc3k.final, method = "scDIOR", out.folder = "benchmark", assay = "RNA", save.scale = TRUE)
Seu2AD <- function(seu.obj, method = c("SeuratDisk", "sceasy", "scDIOR"), out.folder = NULL,
                   out.filename = NULL, assay = "RNA", slot = "counts", save.scale = FALSE, conda.path = NULL) {
  # check parameters
  method <- match.arg(arg = method)
  # check conda env
  if (!is.null(conda.path)) {
    reticulate::use_condaenv(conda.path, required = TRUE)
  }
  # check folder
  if (is.null(out.folder)) {
    out.folder <- getwd()
  }
  if (!dir.exists(out.folder)) {
    message(out.folder, " does not exist, create automatically!")
    dir.create(path = out.folder, showWarnings = FALSE)
  }
  # out name
  out.name <- deparse(substitute(seu.obj))
  # conversion
  if (method == "SeuratDisk") {
    if (is.null(out.filename)) {
      seu.out.name <- file.path(out.folder, paste0(out.name, "_SeuratDisk.h5Seurat"))
    } else {
      seu.out.name <- file.path(out.folder, out.filename)
    }
    if (save.scale) {
      seu.scale <- Seurat::GetAssayData(object = seu.obj, slot = "scale.data", assay = assay)
      if (nrow(seu.scale) == 0) {
        message("There is no scale.data in seu.obj!")
      }
    } else {
      seu.obj[["RNA"]]@scale.data <- matrix(numeric(0), 0, 0)
    }
    seu.log <- tryCatch(
      {
        SeuratDisk::SaveH5Seurat(seu.obj, filename = seu.out.name, overwrite = TRUE)
        SeuratDisk::Convert(seu.out.name, dest = "h5ad", assay = assay, overwrite = TRUE)
      },
      error = function(cond) {
        message("There is an error when using SeuratDisk: ", cond)
      }
    )
    return(seu.log)
  } else if (method == "sceasy") {
    if (is.null(out.filename)) {
      sceasy.out.name <- file.path(out.folder, paste0(out.name, "_sceasy.h5ad"))
    } else {
      sceasy.out.name <- file.path(out.folder, out.filename)
    }
    sceasy.log <- tryCatch(
      {
        # reticulate::use_condaenv("/Applications/anaconda3", required = TRUE)
        # or set RETICULATE_PYTHON = "/Applications/anaconda3/bin/python" in Renvion
        sceasy::convertFormat(seu.obj,
          from = "seurat", to = "anndata", drop_single_values = FALSE,
          outFile = sceasy.out.name, main_layer = slot, assay = assay
        )
      },
      error = function(cond) {
        message("There is an error when using sceasy: ", cond)
      }
    )
    return(sceasy.log)
  } else if (method == "scDIOR") {
    if (is.null(out.filename)) {
      scdior.out.name <- file.path(out.folder, paste0(out.name, "_scDIOR.h5"))
    } else {
      scdior.out.name <- file.path(out.folder, out.filename)
    }
    scdior.log <- tryCatch(
      {
        dior::write_h5(
          data = seu.obj, object.type = "seurat", file = scdior.out.name,
          assay.name = assay, save.scale = save.scale
        )
        # adata = diopy.input.read_h5(file = 'pbmc3k.h5') # require diopy to load h5 to AnnData
      },
      error = function(cond) {
        message("There is an error when using scDIOR: ", cond)
      }
    )
    return(scdior.log)
  }
}

#' Convert SingleCellExperiemnt to AnnData using sceasy/scDIOR/zellkonverter.
#'
#' @param sce.obj A SingleCellExperiment object.
#' @param method Method used to perform conversion, choose from "sceasy",	"scDIOR", "zellkonverter". Default: "sceasy".
#' @param out.folder Output folder. Default: NULL (current working directory).
#' @param out.filename Output file name, e.g., test.h5ad. Default: NULL (\code{seu.obj} name, method, h5ad.).
#' @param slot Slot for adata.X, used when \code{method} is "sceasy", "zellkonverter". Default: counts.
#' @param conda.path Conda environment path, used when \code{method} is "zellkonverter" or "sceasy". Default: NULL.
#'
#' @return Run log.
#' @export
#' @importFrom sceasy convertFormat
#' @importFrom dior write_h5
#' @importFrom reticulate import
#' @importFrom zellkonverter SCE2AnnData
#'
#' @examples
#' # sceasy
#' SCE2AD(sce.obj = pbmc3k.sce, method = "sceasy", out.folder = "benchmark", slot = "rawcounts", conda.path = "/path/to/conda")
#' # scDIOR
#' pbmc3k.sce.scdior <- pbmc3k.sce
#' library(SingleCellExperiment)
#' # scDIOR does not support varm in rowData
#' rowData(pbmc3k.sce.scdior)$varm <- NULL
#' SCE2AD(sce.obj = pbmc3k.sce.scdior, method = "scDIOR", out.folder = "benchmark")
#' # zellkonverter
#' SCE2AD(sce.obj = pbmc3k.sce, method = "zellkonverter", out.folder = "benchmark", slot = "rawcounts", conda.path = "/path/to/conda")
SCE2AD <- function(sce.obj, method = c("sceasy", "scDIOR", "zellkonverter"), out.folder = NULL,
                   out.filename = NULL, slot = "counts", conda.path = NULL) {
  # check parameters
  method <- match.arg(arg = method)
  # check conda env
  if (!is.null(conda.path)) {
    reticulate::use_condaenv(conda.path, required = TRUE)
  }
  # check folder
  if (is.null(out.folder)) {
    out.folder <- getwd()
  }
  if (!dir.exists(out.folder)) {
    message(out.folder, " does not exist, create automatically!")
    dir.create(path = out.folder, showWarnings = FALSE)
  }
  # out name
  out.name <- deparse(substitute(sce.obj))
  # conversion
  if (method == "sceasy") {
    if (is.null(out.filename)) {
      sceasy.out.name <- file.path(out.folder, paste0(out.name, "_sceasy.h5ad"))
    } else {
      sceasy.out.name <- file.path(out.folder, out.filename)
    }
    sceasy.log <- tryCatch(
      {
        # reticulate::use_condaenv("/Applications/anaconda3", required = TRUE)
        # or set RETICULATE_PYTHON = "/Applications/anaconda3/bin/python" in Renvion
        sceasy::convertFormat(sce.obj,
          from = "sce", to = "anndata", drop_single_values = FALSE,
          outFile = sceasy.out.name, main_layer = slot
        )
      },
      error = function(cond) {
        message("There is an error when using sceasy: ", cond)
      }
    )
    return(sceasy.log)
  } else if (method == "scDIOR") {
    if (is.null(out.filename)) {
      scdior.out.name <- file.path(out.folder, paste0(out.name, "_scDIOR.h5"))
    } else {
      scdior.out.name <- file.path(out.folder, out.filename)
    }
    scdior.log <- tryCatch(
      {
        dior::write_h5(
          data = sce.obj, object.type = "singlecellexperiment",
          file = scdior.out.name
        )
        # adata = diopy.input.read_h5(file = 'pbmc3k.sce.h5') # require diopy to load h5 to AnnData
      },
      error = function(cond) {
        message("There is an error when using scDIOR: ", cond)
      }
    )
    return(scdior.log)
  } else if (method == "zellkonverter") {
    if (is.null(out.filename)) {
      zell.out.name <- file.path(out.folder, paste0(out.name, "_zellkonverter.h5ad"))
    } else {
      zell.out.name <- file.path(out.folder, out.filename)
    }
    zell.log <- tryCatch(
      {
        anndata <- reticulate::import("anndata")
        adata <- zellkonverter::SCE2AnnData(sce.obj, X_name = slot)
        adata$write_h5ad(zell.out.name)
      },
      error = function(cond) {
        message("There is an error when using zellkonverter: ", cond)
      }
    )
    return(zell.log)
  }
}

#' Convert AnnData to SeuratObject using SeuratDisk/sceasy/scDIOR/schard/SeuratDisk+scDIOR.
#'
#' @param anndata.file The file contains AnnData.
#' @param method Method used to perform conversion, choose from "SeuratDisk", "sceasy",	"scDIOR", "schard",
#' "SeuratDisk+scDIOR". Default: "SeuratDisk".
#' @param assay Name to store assay data as. Default: RNA.
#' @param load.assays Which assays to load, used when \code{method} is "SeuratDisk". Default: RNA.
#' @param slot Slot to store adata.X, used when \code{method} is "sceasy". Default: counts.
#' @param use.raw Logical value, whether to use adata.raw, used when \code{method} is "schard". Default: TRUE.
#'
#' @return A SeuratObject.
#' @export
#' @importFrom SeuratDisk Convert LoadH5Seurat
#' @importFrom hdf5r H5File
#' @importFrom dior read_h5ad
#' @importFrom Seurat Assays GetAssayData CreateAssayObject
#' @importFrom sceasy convertFormat
#' @importFrom schard h5ad2seurat
#'
#' @examples
#' # SeuratDisk
#' ann.seu <- AD2Seu(anndata.file = "pbmc3k.h5ad", method = "SeuratDisk", assay = "RNA", load.assays = c("RNA"))
#' # sceasy
#' ann.sceasy <- AD2Seu(anndata.file = "pbmc3k.h5ad", method = "sceasy", assay = "RNA", slot = "scale.data")
#' # scDIOR
#' ann.scdior <- AD2Seu(anndata.file = "pbmc3k.h5ad", method = "scDIOR", assay = "RNA")
#' # schard
#' ann.schard <- AD2Seu(anndata.file = "pbmc3k.h5ad", method = "schard", assay = "RNA", use.raw = T)
#' # SeuratDisk+scDIOR
#' ann.seuscdior <- AD2Seu(anndata.file = "pbmc3k.h5ad", method = "SeuratDisk+scDIOR", assay = "RNA", load.assays = c("RNA"))
AD2Seu <- function(anndata.file, method = c("SeuratDisk", "sceasy", "scDIOR", "schard", "SeuratDisk+scDIOR"), assay = "RNA",
                   load.assays = "RNA", slot = "counts", use.raw = TRUE) {
  # check parameters
  method <- match.arg(arg = method)

  # check file
  if (!file.exists(anndata.file)) {
    stop(anndata.file, " does not exist, please check!")
  }
  # conversion
  if (grepl(pattern = "SeuratDisk", x = method)) {
    # SeuratDisk
    seu <- tryCatch(
      {
        SeuratDisk::Convert(anndata.file, dest = "h5seurat", overwrite = TRUE, assay = assay)
        h5seurat.file <- gsub(pattern = "h5ad$", replacement = "h5seurat", x = anndata.file)
        # https://github.com/mojaveazure/seurat-disk/issues/109
        f <- hdf5r::H5File$new(h5seurat.file, "r+")
        groups <- f$ls(recursive = TRUE)
        for (name in groups$name[grepl("categories", groups$name)]) {
          names <- strsplit(name, "/")[[1]]
          names <- c(names[1:length(names) - 1], "levels")
          new_name <- paste(names, collapse = "/")
          f[[new_name]] <- f[[name]]
        }
        for (name in groups$name[grepl("codes", groups$name)]) {
          names <- strsplit(name, "/")[[1]]
          names <- c(names[1:length(names) - 1], "values")
          new_name <- paste(names, collapse = "/")
          f[[new_name]] <- f[[name]]
          grp <- f[[new_name]]
          grp$write(args = list(1:grp$dims), value = grp$read() + 1)
        }
        f$close_all()
        SeuratDisk::LoadH5Seurat(h5seurat.file, assays = load.assays)
      },
      error = function(cond) {
        message("There is an error when using SeuratDisk: ", cond)
      }
    )
    if (grepl(pattern = "scDIOR", x = method)) {
      # scDIOR
      seu.scdior <- tryCatch(
        {
          dior::read_h5ad(file = anndata.file, assay_name = assay, target.object = "seurat")
        },
        error = function(cond) {
          message("There is an error when using scDIOR: ", cond)
        }
      )
      # add additional assays
      all.assays <- Seurat::Assays(seu.scdior)
      unused.assays <- setdiff(all.assays, assay)
      if (length(unused.assays) > 0) {
        for (ay in unused.assays) {
          # https://github.com/JiekaiLab/dior/blob/2b1ea47b6661c8a10d9455f3baeeccb8f12be2f0/R/seuratIO.R#L65
          # https://github.com/satijalab/seurat-object/blob/58bf437fe058dd78913d9ef7b48008a3e24a306a/R/assay.R#L157
          assay.data <- Seurat::GetAssayData(object = seu.scdior[[ay]], slot = "counts")
          seu[[ay]] <- Seurat::CreateAssayObject(counts = assay.data)
        }
      }
      # add graphs
      seu@graphs <- seu.scdior@graphs
    }
  } else if (method == "sceasy") {
    seu <- tryCatch(
      {
        sceasy::convertFormat(anndata.file,
          from = "anndata", to = "seurat",
          main_layer = slot, assay = assay
        )
      },
      error = function(cond) {
        message("There is an error when using sceasy: ", cond)
      }
    )
  } else if (method == "scDIOR") {
    seu <- tryCatch(
      {
        dior::read_h5ad(file = anndata.file, assay_name = assay, target.object = "seurat")
      },
      error = function(cond) {
        message("There is an error when using scDIOR: ", cond)
      }
    )
  } else if (method == "schard") {
    seu <- tryCatch(
      {
        schard::h5ad2seurat(file = anndata.file, use.raw = use.raw, assay = assay)
      },
      error = function(cond) {
        message("There is an error when using schard: ", cond)
      }
    )
  }
  return(seu)
}

#' Convert AnnData to SingleCellExperiemnt using scDIOR/zellkonverter/schard.
#'
#' @param anndata.file The file contains AnnData.
#' @param method Method used to perform conversion, choose from "scDIOR", "zellkonverter", "schard". Default: "scDIOR".
#' @param assay The type of data, used when \code{method} is "scDIOR". Default: "RNA" (scRNA-seq data).
#' @param slot Name used when saving adata.X as an assay, used when \code{method} is "zellkonverter". Default: "counts".
#' @param use.raw Logical value, whether to use adata.raw. Default: TRUE. When \code{method} is "scDIOR" and \code{use.raw}
#' is TRUE (FALSE), raw.X -> assays (X/layers -> assays); When \code{method} is "zellkonverter" and \code{use.raw} is TRUE (FALSE),
#' raw.X -> altExp, X and layers -> assays (X and layers -> assays); When \code{method} is
#' "schard" and \code{use.raw} is TRUE (FALSE), raw.X -> assays (X -> assays).
#' @param conda.path Conda environment path, used when \code{method} is "scDIOR" or "zellkonverter". Default: NULL.
#'
#' @return A SingleCellExperiment object.
#' @export
#' @importFrom reticulate import
#' @importFrom dior read_h5 read_h5ad
#' @importFrom zellkonverter AnnData2SCE
#' @importFrom schard h5ad2sce
#'
#' @examples
#' # scDIOR
#' sce.scdior <- AD2SCE(anndata.file = "pbmc3k.h5ad", method = "scDIOR", assay = "RNA", use.raw = TRUE, conda.path = "/path/to/conda")
#' # zellkonverter
#' sce.zell <- AD2SCE(anndata.file = "pbmc3k.h5ad", method = "zellkonverter", slot = "scale.data", use.raw = TRUE, conda.path = "/path/to/conda")
#' # schard
#' sce.schard <- AD2SCE(anndata.file = "pbmc3k.h5ad", method = "schard", use.raw = TRUE)
AD2SCE <- function(anndata.file, method = c("scDIOR", "zellkonverter", "schard"), assay = "RNA",
                   slot = "counts", use.raw = TRUE, conda.path = NULL) {
  # check parameters
  method <- match.arg(arg = method)
  # check conda env
  if (!is.null(conda.path)) {
    reticulate::use_condaenv(conda.path, required = TRUE)
  }
  # check file
  if (!file.exists(anndata.file)) {
    stop(anndata.file, " does not exist, please check!")
  }
  # conversion
  if (method == "scDIOR") {
    sce <- tryCatch(
      {
        if (use.raw) {
          anndata <- reticulate::import("anndata")
          adata <- anndata$read_h5ad(anndata.file)
          diopy <- reticulate::import("diopy")
          h5.file <- gsub(pattern = ".h5ad$", replacement = "_tmp.h5", x = anndata.file)
          diopy$output$write_h5(adata = adata$raw$to_adata(), file = h5.file, assay_name = assay, save_X = TRUE)
          dior::read_h5(file = h5.file, target.object = "singlecellexperiment")
        } else {
          dior::read_h5ad(file = anndata.file, assay_name = assay, target.object = "singlecellexperiment")
        }
      },
      error = function(cond) {
        message("There is an error when using scDIOR: ", cond)
      }
    )
  } else if (method == "zellkonverter") {
    sce <- tryCatch(
      {
        anndata <- reticulate::import("anndata")
        adata <- anndata$read_h5ad(anndata.file)
        zellkonverter::AnnData2SCE(adata, X_name = slot, raw = use.raw)
      },
      error = function(cond) {
        message("There is an error when using zellkonverter: ", cond)
      }
    )
  } else if (method == "schard") {
    sce <- tryCatch(
      {
        schard::h5ad2sce(anndata.file, use.raw = use.raw)
      },
      error = function(cond) {
        message("There is an error when using schard: ", cond)
      }
    )
  }
}
