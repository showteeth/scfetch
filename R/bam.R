#' Download bam.
#'
#' @param gsm.df Dataframe contains GSM and Run numbers, obtained from \code{ExtractRun}.
#' @param bam.type The source of bam files to download, choose from 10x (e.g. CellRanger) or other. Default: 10x.
#' @param prefetch.path Path to prefetch. Default: NULL (conduct automatic detection).
#' @param samdump.path Path to sam-dump, used when \code{bam.type} is other. Default: NULL (conduct automatic detection).
#' @param out.folder Output folder. Default: NULL (current working directory).
#' @param prefetch.paras Parameters for \code{prefetch}. This should not contain --type or -T values. Default: "-X 100G".
#' @param samdump.paras Parameters for \code{sam-dump}. Default: "".
#'
#' @return Dataframe contains failed runs or NULL.
#' @export
#'
#' @examples
#' \dontrun{
#' # need users to provide prefetch.path
#' GSE138266.runs <- ExtractRun(acce = "GSE138266", platform = "GPL18573")
#' GSE138266.down <- DownloadBam(
#'   gsm.df = GSE138266.runs, bam.type = "10x",
#'   prefetch.path = "/path/to/prefetch",
#'   out.folder = "/path/to/output"
#' )
#' }
DownloadBam <- function(gsm.df, bam.type = c("10x", "other"), prefetch.path = NULL, samdump.path = NULL,
                        out.folder = NULL, prefetch.paras = "-X 100G", samdump.paras = "") {
  # check parameters
  bam.type <- match.arg(arg = bam.type)
  # check prefetch paras
  if (grepl(pattern = "--type|-T", x = prefetch.paras)) {
    stop("Please remomve --type|-T parameter in prefetch.paras!")
  }
  if (bam.type == "10x") {
    message("Download 10x original bam files (keep custom tags generated by 10x softwares)!")
    # add paras to download 10x original bam
    prefetch.paras <- paste0(prefetch.paras, " --type TenX ")
    bam.down <- DownloadSRA(gsm.df = gsm.df, prefetch.path = prefetch.path, out.folder = out.folder, prefetch.paras = prefetch.paras)
    return(bam.down)
  } else if (bam.type == "other") {
    message("Download normal bam files (non-10x), download sra first and convert to bam!")
    sra.down <- DownloadSRA(gsm.df = gsm.df, prefetch.path = prefetch.path, out.folder = out.folder, prefetch.paras = prefetch.paras)
    # convert sra to bam
    all.sras <- list.files(path = out.folder, pattern = "sra$", full.names = TRUE, recursive = TRUE)
    # get sam-dump path
    if (is.null(samdump.path)) {
      # specify sam-dump path
      samdump.path <- Sys.which("sam-dump")
      if (samdump.path == "") {
        stop("Can not find sam-dump automatically, please specify the path!")
      }
    } else {
      samdump.path <- samdump.path
    }
    # run sam-dump
    samdump.run <- sapply(all.sras, function(x) {
      RunSamdump(sra = x, samdump.path = samdump.path, samdump.paras = samdump.paras)
    })
    # select fail sras
    fail.flag <- sapply(names(samdump.run), function(x) {
      !is.null(samdump.run[[x]])
    })
    fail.sra <- names(samdump.run)[fail.flag]
    if (length(fail.sra) > 0) {
      fail.run <- gsub(pattern = ".sra$", replacement = "", x = basename(fail.sra))
      fail.df <- gsm.df[gsm.df$run %in% fail.run, ]
      if (!is.null(sra.down)) {
        fail.df <- as.data.frame(rbind(fail.df, sra.down))
        message(
          length(sra.down$run), " runs failed download and ", length(fail.sra), " sras failed to convert to bam.",
          "Return gsm.df, you can re-run with gsm.df!"
        )
      } else {
        message(
          length(fail.sra), " sras failed to convert to bam: ",
          paste0(fail.sra, collapse = ","), ". Return gsm.df, you can re-run with gsm.df!"
        )
      }
      return(fail.df)
    } else {
      if (!is.null(sra.down)) {
        fail.df <- sra.down
        message(
          length(sra.down$run), " runs failed download: ", paste0(sra.down$run, collapse = ","),
          ". Return gsm.df, you can re-run with gsm.df!"
        )
        return(fail.df)
      } else {
        message("All sra download and convert to bam successfully!")
        return(NULL)
      }
    }
  }
}

# convert sra to bam files
RunSamdump <- function(sra, samdump.path, samdump.paras) {
  out.bam <- file.path(dirname(sra), gsub(pattern = "sra$", replacement = "bam", x = basename(sra)))
  # sam-dump command
  samdump.cmd <- paste(samdump.path, samdump.paras, sra, "| samtools view -bS - -o ", out.bam)
  # run command
  message(paste("Calling sam-dump:", samdump.cmd))
  samdump.status <- system(samdump.cmd, intern = TRUE)
  samdump.status.code <- attr(samdump.status, "status")
  if (!is.null(samdump.status.code)) {
    warning("Run sam-dump error on: ", sra, ", remove partial files!")
    file.remove(out.bam)
    return(sra)
  } else {
    message("Run sam-dump successful on: ", sra)
    return(NULL)
  }
}

#' Convert bam files to fastq files.
#'
#' @param bam.folder Folder contains bam files, obtained from \code{DownloadSRA}. Default: NULL.
#' @param bam.path Paths of bams. \code{bam.folder} and \code{bam.path} cannot be both NULL. Default: NULL.
#' @param bam.type The source of bam files, choose from 10x (e.g. CellRanger) or other. Default: 10x.
#' @param pair.end The bam files are pair-end or single-end, used when \code{bam.type} is other. Default: NULL (identify with flag).
#' @param bamtofastq.path Path to 10x bamtofastq (\code{bam.type} is 10x) or samtools (\code{bam.type} is other).
#' Default: NULL (conduct automatic detection with bamtofastq_linux/samtools).
#' @param bamtofastq.paras Parameters for \code{bamtofastq.path}. Default: "--nthreads 4".
#' @param sort.name Logical value, whether the bam files are sorted by name, required when \code{bam.type} is other. Default: FALSE.
#' @param sort.thread The number of threads for bam sorting, used when \code{bam.type} is other. Default: 4.
#'
#' @return NULL or paths of failed bams.
#' @export
#'
#' @examples
#' \dontrun{
#' # need users to provide prefetch.path and bamtofastq.path
#' GSE138266.runs <- ExtractRun(acce = "GSE138266", platform = "GPL18573")
#' GSE138266.down <- DownloadBam(
#'   gsm.df = GSE138266.runs, prefetch.path = "/path/to/prefetch",
#'   out.folder = "/path/to/output"
#' )
#' GSE138266.convert <- Bam2Fastq(
#'   bam.folder = "/path/to/output",
#'   bamtofastq.path = "/path/to/bamtofastq_linux or samtools",
#'   bamtofastq.paras = "--nthreads 4"
#' )
#' }
Bam2Fastq <- function(bam.folder = NULL, bam.path = NULL, bam.type = c("10x", "other"), pair.end = NULL,
                      bamtofastq.path = NULL, bamtofastq.paras = "--nthreads 4", sort.name = FALSE, sort.thread = 4) {
  # check parameters
  bam.type <- match.arg(arg = bam.type)

  # get bamtofastq/samtools path
  if (is.null(bamtofastq.path)) {
    if (bam.type == "10x") {
      # specify bamtofastq path
      bamtofastq.path <- Sys.which("bamtofastq_linux")
    } else {
      # specify samtools path
      bamtofastq.path <- Sys.which("samtools")
    }
    if (bamtofastq.path == "") {
      stop("Can not find bamtofastq_linux/samtools automatically, please specify the path!")
    }
  } else {
    bamtofastq.path <- bamtofastq.path
  }
  # prepare bam files
  if (!is.null(bam.folder)) {
    all.bams <- list.files(path = bam.folder, pattern = "bam$", full.names = TRUE, recursive = TRUE)
  } else if (!is.null(bam.path)) {
    all.bams <- bam.path
  } else {
    stop("Please provide either valid bam.folder or bam.path!")
  }
  # check bams
  if (length(all.bams) == 0) {
    stop("There is no bam files detected. Please check!")
  }
  # run conversion
  all.bams.convert <- sapply(all.bams, function(x) {
    Runbamtofastq(
      bam.path = x, bam.type = bam.type, pair.end = pair.end, bamtofastq.path = bamtofastq.path,
      bamtofastq.paras = bamtofastq.paras, sort.name = sort.name, sort.thread = sort.thread
    )
  })
  # select fail bams
  fail.flag <- sapply(names(all.bams.convert), function(x) {
    !is.null(all.bams.convert[[x]])
  })
  fail.bam <- names(all.bams.convert)[fail.flag]
  # return
  if (length(fail.bam) > 0) {
    message(
      length(fail.bam), " runs failed to convert to fastqs: ",
      paste0(fail.bam, collapse = ","), ". Return paths, you can re-run with bam.path!"
    )
    return(fail.bam)
  } else {
    message("All bams convert successfully!")
    return(NULL)
  }
}

# convert bam to fastq: for 10x bam or general bam files
Runbamtofastq <- function(bam.path, bam.type, pair.end, bamtofastq.path, bamtofastq.paras, sort.name, sort.thread) {
  out.folder <- file.path(dirname(bam.path), "bam2fastq")
  # bamtofastq command
  if (bam.type == "10x") {
    bamtofastq.cmd <- paste(bamtofastq.path, bamtofastq.paras, bam.path, out.folder)
  } else if (bam.type == "other") {
    if (isFALSE(sort.name)) {
      sortn.bam <- gsub(pattern = ".bam$", replacement = ".sortname.bam", x = bam.path)
      bam.prefix <- gsub(pattern = ".bam$", replacement = "", x = basename(bam.path))
      sort.cmd <- paste(bamtofastq.path, "sort -@", sort.thread, "-T", bam.prefix, "-o", sortn.bam, bam.path)
      sort.status <- system(sort.cmd, intern = TRUE)
      bam.path <- sortn.bam
    }
    # check pair or single
    if (is.null(pair.end)) {
      bam.end <- CheckBam(bam = bam.path, samtools.path = bamtofastq.path)
    } else {
      bam.end <- pair.end
    }
    # create output folder
    if (!dir.exists(out.folder)) {
      dir.create(path = out.folder, recursive = TRUE)
    }
    if (bam.end) {
      fq1.name <- file.path(out.folder, gsub(pattern = ".bam$", replacement = "_1.fastq.gz", x = basename(bam.path)))
      fq2.name <- file.path(out.folder, gsub(pattern = ".bam$", replacement = "_2.fastq.gz", x = basename(bam.path)))
      bamtofastq.cmd <- paste(bamtofastq.path, "fastq", bamtofastq.paras, "-1", fq1.name, "-2", fq2.name, "-0 /dev/null -s /dev/null", bam.path)
    } else {
      fq.name <- file.path(out.folder, gsub(pattern = ".bam$", replacement = ".fastq.gz", x = basename(bam.path)))
      bamtofastq.cmd <- paste(bamtofastq.path, "fastq", bamtofastq.paras, bam.path, "|gzip - >", fq.name)
    }
  }
  # run command
  message(paste("Calling bamtofastq:", bamtofastq.cmd))
  bamtofastq.status <- system(bamtofastq.cmd, intern = TRUE)
  bamtofastq.status.code <- attr(bamtofastq.status, "status")
  if (!is.null(bamtofastq.status.code)) {
    warning("Run bamtofastq error on: ", bam.path, ", please remove ", out.folder, " and re-run!")
    return(bam.path)
  } else {
    message("Conversion successful: ", bam.path)
    return(NULL)
  }
}

# # convert bam to fastq: only for 10x
# Runbamtofastq <- function(bam.path, bamtofastq.path, bamtofastq.paras) {
#   out.folder <- file.path(dirname(bam.path), "bam2fastq")
#   # bamtofastq command
#   bamtofastq.cmd <- paste(bamtofastq.path, bamtofastq.paras, bam.path, out.folder)
#   # run command
#   message(paste("Calling bamtofastq:", bamtofastq.cmd))
#   bamtofastq.status <- system(bamtofastq.cmd, intern = TRUE)
#   bamtofastq.status.code <- attr(bamtofastq.status, "status")
#   if (!is.null(bamtofastq.status.code)) {
#     warning("Run bamtofastq error on: ", bam.path, ", please remove ", out.folder, " and re-run!")
#     return(bam.path)
#   } else {
#     message("Conversion successful: ", bam.path)
#     return(NULL)
#   }
# }
