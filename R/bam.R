#' Download bam.
#'
#' @param gsm.df Dataframe contains GSM and Run numbers, obtained from \code{ExtractRun}.
#' @param prefetch.path Path to prefetch. Default: NULL (conduct automatic detection).
#' @param out.folder Output folder. Default: NULL (current working directory).
#' @param prefetch.paras Parameters for \code{prefetch}. This should not contain --type or -T values. Default: "-X 100G".
#'
#' @return Dataframe contains failed runs or NULL.
#' @export
#'
#' @examples
#' # GSE138266.runs = ExtractRun(acce = "GSE138266", platform = "GPL18573")
#' # GSE138266.down = DownloadBam(gsm.df = GSE138266.runs, prefetch.path = "/path/to/prefetch",
#' #                              out.folder = "/path/to/output")
DownloadBam <- function(gsm.df, prefetch.path = NULL, out.folder = NULL, prefetch.paras = "-X 100G") {
  # check prefetch paras
  if (grepl(pattern = "--type|-T", x = prefetch.paras)) {
    stop("Please remomve --type|-T parameter in prefetch.paras!")
  }
  # add paras to download 10x original bam
  prefetch.paras <- paste0(prefetch.paras, " --type TenX ")
  bam.down <- DownloadSRA(gsm.df = gsm.df, prefetch.path = prefetch.path, out.folder = out.folder, prefetch.paras = prefetch.paras)
  return(bam.down)
}

#' Convert 10x bam to fastqs.
#'
#' @param bam.folder Folder contains bam files, obtained from \code{DownloadSRA}. Default: NULL.
#' @param bam.path Paths of bams. \code{bam.folder} and \code{bam.path} cannot be both NULL. Default: NULL.
#' @param bamtofastq.path Path to 10x bamtofastq. Default: NULL (conduct automatic detection with bamtofastq_linux).
#' @param bamtofastq.paras Parameters for \code{bamtofastq.path}. Default: "--nthreads 4".
#'
#' @return NULL or paths of failed bams.
#' @export
#'
#' @examples
#' # GSE138266.runs = ExtractRun(acce = "GSE138266", platform = "GPL18573")
#' # GSE138266.down = DownloadBam(gsm.df = GSE138266.runs, prefetch.path = "/path/to/prefetch",
#' #                              out.folder = "/path/to/output")
#' # GSE138266.convert = Bam2Fastq(bam.folder = "/path/to/output", bamtofastq.path = "/path/to/bamtofastq_linux",
#' #                               bamtofastq.paras = "--nthreads 4")
Bam2Fastq <- function(bam.folder = NULL, bam.path = NULL, bamtofastq.path = NULL, bamtofastq.paras = "--nthreads 4") {
  # get bamtofastq path
  if (is.null(bamtofastq.path)) {
    # specify bamtofastq path
    bamtofastq.path <- Sys.which("bamtofastq_linux")
    if (bamtofastq.path == "") {
      stop("Can not find bamtofastq_linux automatically, please specify the path!")
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
    Runbamtofastq(bam.path = x, bamtofastq.path = bamtofastq.path, bamtofastq.paras = bamtofastq.paras)
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

# convert bam to fastq
Runbamtofastq <- function(bam.path, bamtofastq.path, bamtofastq.paras) {
  out.folder <- file.path(dirname(bam.path), "bam2fastq")
  # bamtofastq command
  bamtofastq.cmd <- paste(bamtofastq.path, bamtofastq.paras, bam.path, out.folder)
  # run command
  message(paste("Calling bamtofastq: ", bamtofastq.cmd))
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
