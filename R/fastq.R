#' Extract Runs with GEO Accession Number or GSM Number.
#'
#' @param gsm GSM number. Default: NULL (use \code{acce}).
#' @param acce GEO accession number. Default: NULL (use \code{gsm}).
#' \code{acce} and \code{gsm} cannot be both NULL.
#' @param platform Platform information/field. Default: NULL (all platforms).
#' @param parallel Logical value, whether to process GSM parallelly. Default: TRUE.
#' @param ... Parameters for \code{\link{getGEO}}. Used when \code{acce} is not NULL.
#'
#' @return Dataframe contains GSM and Runs.
#' @importFrom magrittr %>%
#' @importFrom GEOquery getGEO
#' @importFrom Biobase annotation experimentData pData phenoData notes sampleNames exprs
#' @importFrom parallel detectCores
#' @export
#'
#' @examples
#' \dontrun{
#' GSE186003.runs <- ExtractRun(acce = "GSE186003", platform = "GPL24247", parallel = FALSE)
#' }
ExtractRun <- function(gsm = NULL, acce = NULL, platform = NULL, parallel = TRUE, ...) {
  # get GSM
  if (is.null(gsm)) {
    message("Extract all GSM with acce: ", acce, " and platform: ", platform)
    gsm.meta <- ExtractGEOMeta(acce = acce, platform = platform, ...)
    gsm <- gsm.meta$geo_accession
  }
  # prepare core
  if (parallel) {
    cores.used <- min(parallel::detectCores(), length(gsm))
  } else {
    cores.used <- 1
  }
  # extract srr
  rnsgeofastq <- requireNamespace("GEOfastq", quietly = TRUE)
  if (!rnsgeofastq) {
    stop("Can not find GEOfastq package, install with devtools::install_github('alexvpickering/GEOfastq')!")
  }
  gsm.run.df <- GEOfastq::crawl_gsms(gsm, max.workers = cores.used)
  if (is.null(gsm.run.df)) {
    stop("There is no valid srr numbers available, please check the raw data available under ", paste0(gsm, collapse = ", "))
  } else {
    gsm.run.df <- gsm.run.df[c("title", "gsm_name", "experiment", "run")]
  }
  return(gsm.run.df)
}

#' Download SRA.
#'
#' @param gsm.df Dataframe contains GSM and Run numbers, obtained from \code{ExtractRun}.
#' @param prefetch.path Path to prefetch. Default: NULL (conduct automatic detection).
#' @param out.folder Output folder. Default: NULL (current working directory).
#' @param prefetch.paras Parameters for \code{prefetch}. Default: "-X 100G".
#'
#' @return Dataframe contains failed runs or NULL.
#' @export
#'
#' @examples
#' \dontrun{
#' # need users to provide the prefetch.path and out.folder
#' GSE186003.runs <- ExtractRun(acce = "GSE186003", platform = "GPL24247")
#' GSE186003.down <- DownloadSRA(
#'   gsm.df = GSE186003.runs, prefetch.path = "/path/to/prefetch",
#'   out.folder = "/path/to/output"
#' )
#' }
DownloadSRA <- function(gsm.df, prefetch.path = NULL, out.folder = NULL, prefetch.paras = "-X 100G") {
  # check dataframe
  if (nrow(gsm.df) == 0) {
    stop("Please provide valid gsm.df!")
  }
  CheckColumns(df = gsm.df, columns = c("gsm_name", "run"))
  # all samples
  all.samples <- unique(as.character(gsm.df$gsm_name))
  # prepare output folder
  if (is.null(out.folder)) {
    out.folder <- getwd()
  }
  # get prefetch path
  if (is.null(prefetch.path)) {
    # specify prefetch path
    prefetch.path <- Sys.which("prefetch")
    if (prefetch.path == "") {
      stop("Can not find prefetch automatically, please specify the path!")
    }
  } else {
    prefetch.path <- prefetch.path
  }
  # prepare sample output folder
  samples.folder <- file.path(out.folder, gsm.df$gsm_name)
  names(samples.folder) <- gsm.df$run
  all.runs.down <- sapply(names(samples.folder), function(x) {
    sf <- samples.folder[x]
    RunPrefetch(sra = x, prefetch.path = prefetch.path, out.folder = sf, prefetch.paras = prefetch.paras)
  })
  # select fail samples
  fail.flag <- sapply(names(all.runs.down), function(x) {
    !is.null(all.runs.down[[x]])
  })
  fail.run <- names(all.runs.down)[fail.flag]
  if (length(fail.run) > 0) {
    fail.df <- gsm.df[gsm.df$run %in% fail.run, ]
    message(
      length(fail.run), " samples failed to download: ",
      paste0(fail.run, collapse = ","), ". Return dataframe contains failed samples, you can re-run!"
    )
    return(fail.df)
  } else {
    message("All samples downloaded successfully!")
    return(NULL)
  }
}

# run prefetch per run
RunPrefetch <- function(sra, prefetch.path, out.folder, prefetch.paras) {
  # prepare output folder
  if (!dir.exists(out.folder)) {
    dir.create(path = out.folder, recursive = TRUE)
  }
  cwd <- getwd()
  on.exit(setwd(cwd))
  # change directory to avoid bam download bug
  setwd(out.folder)
  # prefetch command
  prefetch.cmd <- paste(prefetch.path, prefetch.paras, "-O", out.folder, sra)
  # run command
  message(paste("Calling Prefetch: ", prefetch.cmd))
  prefetch.status <- system(prefetch.cmd, intern = TRUE)
  prefetch.status.code <- attr(prefetch.status, "status")
  if (!is.null(prefetch.status.code)) {
    warning("Run prefetch error on: ", sra, ", remove partial files!")
    do.call(file.remove, list(list.files(out.folder, full.names = TRUE)))
    return(sra)
  } else {
    message("Download successful: ", sra)
    return(NULL)
  }
}

#' Split SRA to fastq Files and Format to 10x Standard Style.
#'
#' @param sra.folder Folder contains all sras, obtained from \code{DownloadSRA}. Default: NULL.
#' @param sra.path Paths of sras. \code{sra.folder} and \code{sra.path} cannot be both NULL. Default: NULL.
#' @param fastq.type The source of fastq files, choose from 10x (use \code{--split-files} to split sra) or
#' other (use \code{--split-3} to split sra). Default: 10x.
#' @param split.cmd.path The full command path used to split, can be path to parallel-fastq-dump,
#' fasterq-dump and fastq-dump. Default: NULL (conduct automatic detection).
#' @param sratools.path Path to sratoolkit bin. When \code{split.cmd.path} is path to parallel-fastq-dump,
#' it requires sratoolkit. Default: NULL (conduct automatic detection).
#' @param split.cmd.paras Parameters for \code{split.cmd.path}. Default: NULL.
#' @param split.cmd.threads Threads, used when \code{split.cmd.path} is path to parallel-fastq-dump or fasterq-dump.
#' Default: NULL (1).
#' @param format.10x Logical value, whether to format split fastqs to 10x standard format. Default: TRUE.
#' @param remove.raw Logical value, whether to remove old split fastqs (unformatted), used when \code{format.10x} is TRUE. Default: FALSE.
#'
#' @return NULL or paths of failed sras.
#' @importFrom utils read.table
#' @export
#'
#' @examples
#' \dontrun{
#' # need users to provide the prefetch.path, sra.folder, split.cmd.path, sratools.path and out.folder
#' GSE186003.runs <- ExtractRun(acce = "GSE186003", platform = "GPL24247")
#' GSE186003.down <- DownloadSRA(
#'   gsm.df = GSE186003.runs, prefetch.path = "/path/to/prefetch",
#'   out.folder = "/path/to/output"
#' )
#' GSE186003.split <- SplitSRA(
#'   sra.folder = "/path/to/output",
#'   split.cmd.path = "/path/to/parallel-fastq-dump",
#'   sratools.path = "/path/to/sra/bin", fastq.type = "10x",
#'   split.cmd.threads = 4
#' )
#' }
SplitSRA <- function(sra.folder = NULL, sra.path = NULL, fastq.type = c("10x", "other"), split.cmd.path = NULL, sratools.path = NULL,
                     split.cmd.paras = NULL, split.cmd.threads = NULL, format.10x = TRUE, remove.raw = FALSE) {
  # check parameters
  fastq.type <- match.arg(arg = fastq.type)

  # check split.cmd.paras
  if (!is.null(split.cmd.paras)) {
    split.cmd.paras.sep <- unlist(strsplit(split.cmd.paras, split = "[[:space:]]+"))
    # invalid scp
    invalid.scp <- c(
      "--gzip", "--bzip2", "--split-files", "--split-3", "--split-spot",
      "-s", "-S", "-3", "--outdir", "-O", "-Z", "--stdout", "-e", "--threads",
      "--sra-id"
    )
    # removed split.cmd.paras
    split.cmd.paras.remove <- intersect(split.cmd.paras.sep, invalid.scp)
    if (length(split.cmd.paras.remove) > 0) {
      message("The split.cmd.paras you proided should not contain: ", paste(invalid.scp, collapse = ", "), ". Please check!")
    }
    # valid split.cmd.paras
    split.cmd.paras.sep <- setdiff(split.cmd.paras.sep, invalid.scp)
    split.cmd.paras <- paste(split.cmd.paras.sep, collapse = " ")
  }

  # get fastq-dump/fasterq-dump/parallel-fastq-dump path
  if (is.null(split.cmd.path)) {
    # specify prefetch path
    all.split.cmd.path <- Sys.which(c("parallel-fastq-dump", "fasterq-dump", "fastq-dump"))
    valid.split.cmd.path <- all.split.cmd.path[which(all.split.cmd.path != "")]
    if (length(valid.split.cmd.path) == 0) {
      stop("Can not find parallel-fastq-dump, fasterq-dump or fastq-dump automatically, please specify the path!")
    }
    split.cmd.path <- valid.split.cmd.path[1]
  } else {
    split.cmd.path <- split.cmd.path
  }
  # get cmd used
  split.name <- basename(split.cmd.path)
  if (split.name == "parallel-fastq-dump") {
    if (is.null(sratools.path)) {
      sratools.path <- Sys.which("fastq-dump")
      if (sratools.path == "") {
        stop("parallel-fastq-dump requires sratools, please specify the sratools.path!")
      }
      sratools.path <- basename(sratools.path)
    } else {
      sratools.path <- sratools.path
    }
  }
  # prepare thread
  if (split.name %in% c("parallel-fastq-dump", "fasterq-dump")) {
    if (is.null(split.cmd.threads)) {
      split.cmd.threads <- 1
    }
  }
  # get all sra
  if (!is.null(sra.folder)) {
    all.sras <- list.files(path = sra.folder, pattern = "sra$", full.names = TRUE, recursive = TRUE)
  } else if (!is.null(sra.path)) {
    all.sras <- sra.path
  } else {
    stop("Please provide either valid sra.folder or sra.path!")
  }
  # check sras
  if (length(all.sras) == 0) {
    stop("There is no sra files. Please check!")
  }
  # run split
  all.sras.split <- sapply(all.sras, function(x) {
    RunSplit(
      sra.path = x, fastq.type = fastq.type, split.cmd.path = split.cmd.path, sratools.path = sratools.path,
      split.cmd.paras = split.cmd.paras, split.cmd.threads = split.cmd.threads
    )
  })
  # select fail sras
  fail.flag <- sapply(names(all.sras.split), function(x) {
    !is.null(all.sras.split[[x]])
  })
  fail.sra <- names(all.sras.split)[fail.flag]
  # format fastqs to 10x standard names
  if (format.10x) {
    success.sra <- setdiff(names(all.sras.split), fail.sra)
    if (length(success.sra) > 0) {
      success.sra.folders <- dirname(success.sra)
      format.flag <- sapply(success.sra.folders, function(x) {
        IdentifyReads(fastq.folder = x, remove = remove.raw)
      })
    }
  }
  if (length(fail.sra) > 0) {
    message(
      length(fail.sra), " runs failed to split: ",
      paste0(fail.sra, collapse = ","), ". Return paths, you can re-run with sra.path!"
    )
    return(fail.sra)
  } else {
    message("All runs split successfully!")
    return(NULL)
  }
}

# split single sra file
RunSplit <- function(sra.path, fastq.type, split.cmd.path, sratools.path, split.cmd.paras, split.cmd.threads) {
  # prepare output folder
  out.folder <- dirname(sra.path)
  # check the fastq type
  if (fastq.type == "10x") {
    sp <- "--split-files"
  } else if (fastq.type == "other") {
    sp <- "--split-3"
  }
  # get cmd used
  split.name <- basename(split.cmd.path)
  # split command
  if (split.name == "parallel-fastq-dump") {
    split.cmd <- paste(
      split.cmd.path, split.cmd.paras, sp, "--gzip", "-O", out.folder,
      "--sra-id", sra.path, "--threads", split.cmd.threads
    )
    split.cmd <- paste0("export PATH=", sratools.path, ":$PATH;", split.cmd)
  } else if (split.name == "fasterq-dump") {
    split.cmd <- paste(
      split.cmd.path, split.cmd.paras, sp, "-O", out.folder,
      "--threads", split.cmd.threads, sra.path
    )
  } else {
    split.cmd <- paste(split.cmd.path, split.cmd.paras, sp, "--gzip", "-O", out.folder, sra.path)
  }
  # run command
  message(paste("Calling", split.name, ": ", split.cmd))
  split.status <- system(split.cmd, intern = TRUE)
  split.status.code <- attr(split.status, "status")
  if (!is.null(split.status.code)) {
    warning("Run", split.name, " error on: ", sra.path, ", remove partial files!")
    do.call(file.remove, list(list.files(out.folder, full.names = TRUE, pattern = "fastq.gz$|fastq$")))
    return(sra.path)
  } else {
    message("Split successfully: ", sra.path)
    return(NULL)
  }
}

# identify read1, read2, index and rename to 10x style
IdentifyReads <- function(fastq.folder, remove = FALSE) {
  all.fastqs <- list.files(path = fastq.folder, pattern = "fastq$|fastq.gz$|fq$|fq.gz$", full.names = TRUE)
  if (length(all.fastqs) == 1) {
    warning("There is only one fastq detected under: ", fastq.folder, ". Escape!")
  } else if (length(all.fastqs) >= 2) {
    all.fastqs.len <- sapply(all.fastqs, GetFastqLen)
    possible.index <- names(all.fastqs.len)[all.fastqs.len %in% c(8, 10)]
    valid.fastqs <- setdiff(all.fastqs, possible.index)
    if (length(valid.fastqs) >= 2) {
      # get possible read1
      possible.r1 <- names(all.fastqs.len)[all.fastqs.len %in% c(26, 28)]
      if (length(possible.r1) == 0) {
        possible.r1 <- valid.fastqs[1]
        message("There is no fastq with read length 26 or 28, choose the first (non-index): ", possible.r1, " as read1! Please check manually.")
        remove <- FALSE
      } else if (length(possible.r1) >= 2) {
        possible.r1 <- valid.fastqs[1]
        message("There are multiple fastqs with read length 26 or 28, choose the first (non-index): ", possible.r1, " as read1! Please check manually.")
        remove <- FALSE
      }
      # get possible read2
      possible.r2 <- setdiff(valid.fastqs, possible.r1)
      if (length(possible.r2) >= 2) {
        possible.r2 <- possible.r2[1]
        message("There are multiple possible read2 fastqs, choose the first (non-index): ", possible.r2, " as read2! Please check manually.")
        remove <- FALSE
      }
      # get folder and new name
      fastq.dir <- dirname(possible.r1)
      fastq.name <- gsub(pattern = "(.*)_[1-9]*.*", replacement = "\\1", x = basename(possible.r1))
      fq1.new.name <- file.path(fastq.dir, paste0(fastq.name, "_S1_L001_R1_001.fastq.gz"))
      fq2.new.name <- file.path(fastq.dir, paste0(fastq.name, "_S1_L001_R2_001.fastq.gz"))
      # rename
      fq1.copy.tag <- file.copy(from = possible.r1, to = fq1.new.name)
      fq2.copy.tag <- file.copy(from = possible.r2, to = fq2.new.name)
      # remove old files
      if (remove) {
        fq1.remove.tag <- file.remove(possible.r1)
        fq2.remove.tag <- file.remove(possible.r2)
      }
    } else {
      warning(
        "Valid fastq file is less than 2 after filtering possible index fastqs: ",
        paste0(possible.index, collapse = ", ")
      )
    }
  }
}

# get fastq max length, used in IdentifyReads
GetFastqLen <- function(fq.file) {
  fq.head <- utils::read.table(
    file = fq.file, nrows = 10000, sep = "\t",
    comment.char = "", stringsAsFactors = FALSE
  )
  fq.valid <- fq.head[1:nrow(fq.head) %% 4 == 2, ]
  fq.len <- max(sapply(fq.valid, nchar))
  return(fq.len)
}
