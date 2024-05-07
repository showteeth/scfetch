#' Run CellRanger on Downloaded FASTQ Files.
#'
#' @param sample.dir Directory contains all samples.
#' @param ref Path of folder containing 10x-compatible
#' transcriptome reference (\code{--transcriptome}).
#' @param localcores Set max cores the pipeline may request at one time (\code{--localcores}).
#' Only applies to local jobs (\code{--jobmode=local}). Default: 4.
#' @param localmem Set max GB the pipeline may request at one time (\code{--localmem}).
#' Only applies to local jobs (\code{--jobmode=local}). Default: 16.
#' @param out.folder Output folder. Default: NULL (current working directory).
#' @param cr.path Path to cellranger. Default: NULL (conduct automatic detection).
#' @param cr.paras Parameters for \code{cellranger}.
#' Default: "--chemistry=auto --jobmode=local".
#'
#' @return Vector contains failed samples or NULL.
#' @export
#'
#' @examples
#' \dontrun{
#' RunCellRanger(
#'   sample.dir = "/path/to/fastq",
#'   ref = "/path/to/cellranger_tiny_ref/3.0.0",
#'   out.folder = "/path/to/results",
#'   cr.path = "/path/to/cellranger-x.x.x/cellranger"
#' )
#' }
RunCellRanger <- function(sample.dir, ref, localcores = 4, localmem = 16, out.folder = NULL,
                          cr.path = NULL, cr.paras = "--chemistry=auto --jobmode=local") {
  # check the path
  if (!dir.exists(sample.dir)) {
    stop(sample.dir, " doesn't exist, please check and re-run!")
  }
  sample.fq.dir <- dir(path = sample.dir, full.names = TRUE)
  all.sample.cr <- sapply(sample.fq.dir, function(x) {
    RunCellRangerSingle(
      fq.dir = x, transcriptome = ref, localcores = localcores,
      localmem = localmem, out.folder = out.folder,
      cr.path = cr.path, cr.paras = cr.paras
    )
  })
  # select fail samples
  fail.flag <- sapply(names(all.sample.cr), function(x) {
    !is.null(all.sample.cr[[x]])
  })
  fail.sample <- all.sample.cr[fail.flag]
  if (length(fail.sample) > 0) {
    message(
      "CellRanger failed on: ",
      paste0(fail.sample, collapse = ","), ". Please check and re-run!"
    )
    return(fail.sample)
  } else {
    message("CellRanger run successfully on all samples!")
    return(NULL)
  }
}

# run CellRanger on single sample
RunCellRangerSingle <- function(fq.dir, transcriptome, localcores = 4, localmem = 16, out.folder = NULL,
                                cr.path = NULL, cr.paras = "--chemistry=auto --jobmode=local") {
  # get cellranger path
  if (is.null(cr.path)) {
    # specify cellranger path
    cr.path <- Sys.which("cellranger")
    if (cr.path == "") {
      stop("Can not find cellranger automatically, please specify the path!")
    }
  } else {
    cr.path <- cr.path
  }
  # prepare output folder
  if (is.null(out.folder)) {
    out.folder <- getwd()
  }
  if (!dir.exists(out.folder)) {
    message(out.folder, " does not exist, create automatically!")
    dir.create(out.folder, recursive = TRUE)
  }
  # check the path
  if (!dir.exists(fq.dir)) {
    message(fq.dir, " doesn't exist, please check and re-run!")
    return(basename(fq.dir))
  } else {
    fq.files <- list.files(path = fq.dir, pattern = "fastq.gz$")
    if (length(fq.files) == 0) {
      message("There is no fastq.gz files under ", fq.dir, " , please check and re-run!")
      return(basename(fq.dir))
    } else {
      fq1.file <- grep(pattern = ".*_S[0-9]_L[0-9]{3}_R1_[0-9]{3}.fastq.gz", x = fq.files, value = T)
      fq2.file <- grep(pattern = ".*_S[0-9]_L[0-9]{3}_R2_[0-9]{3}.fastq.gz", x = fq.files, value = T)
      if (length(fq1.file) > 0 && length(fq2.file) > 0) {
        sample.id <- basename(fq.dir)
        sample.name <- basename(fq.dir)
        # check additional paras
        if (grepl(pattern = "--id|--transcriptome|--fastqs|--sample|--localcores|--localmem", cr.paras)) {
          message(
            cr.paras, " overlaps with built-in paras (--id, --transcriptome, --fastqs, --sample,",
            " --localcores, --localmem), please check and re-run!"
          )
        }
        # cellranger command
        cr.cmd <- paste0(
          "cd ", out.folder, " && ", cr.path, " count --id=", sample.id, " --transcriptome=", transcriptome,
          " --fastqs=", fq.dir, " --sample=", sample.name, " --localcores=", localcores,
          " --localmem=", localmem, " ", cr.paras
        )
        # run command
        message(paste("Calling CellRanger:", cr.cmd))
        cr.status <- system(cr.cmd, intern = TRUE)
        cr.status.code <- attr(cr.status, "status")
        if (!is.null(cr.status.code)) {
          cr.status.msg <- paste(cr.status, collapse = " ")
          warning(
            "Run CellRanger error on: ", sample.name, ". Error message :",
            cr.status.msg, " .Please check and re-run!"
          )
          return(sample.name)
        } else {
          message("Finish CellRanger: ", sample.name)
          return(NULL)
        }
      } else {
        message("There is no R1 or R2 fastq.gz files under ", fq.dir, " , please check and re-run!")
        return(basename(fq.dir))
      }
    }
  }
}

#' Run STAR on Downloaded FASTQ Files.
#'
#' @param sample.dir Directory contains all samples.
#' @param ref Path of folder containing STAR version-compatible reference.
#' @param out.folder Output folder. Default: NULL (current working directory).
#' @param thread The number of threads to use. Default: 4.
#' @param star.path Path to \code{STAR}. Default: NULL (conduct automatic detection).
#' @param star.paras Parameters for \code{STAR}.
#' Default: "--outBAMsortingThreadN 4 --twopassMode None".
#'
#' @return Vector contains failed samples or NULL.
#' @export
#'
#' @examples
#' \dontrun{
#' RunSTAR(sample.dir = "/path/to/fastq",
#'         ref = "/path/to/star/index",
#'         out.folder = "/path/to/star/mapping",
#'         star.path = "/path/to/bin/STAR")
#' }
RunSTAR <- function(sample.dir, ref, out.folder = NULL, thread = 4, star.path = NULL,
                    star.paras = "--outBAMsortingThreadN 4 --twopassMode None") {
  # check the path
  if (!dir.exists(sample.dir)) {
    stop(sample.dir, " doesn't exist, please check and re-run!")
  }
  sample.fq.dir <- dir(path = sample.dir, full.names = TRUE)
  all.sample.star <- sapply(sample.fq.dir, function(x) {
    RunSTARSingle(
      fq.dir = x, ref = ref, out.folder = out.folder, thread = thread,
      star.path = star.path, star.paras = star.paras
    )
  })
  # select fail samples
  fail.flag <- sapply(names(all.sample.star), function(x) {
    !is.null(all.sample.star[[x]])
  })
  fail.sample <- all.sample.star[fail.flag]
  if (length(fail.sample) > 0) {
    message(
      "STAR failed on: ",
      paste0(fail.sample, collapse = ","), ". Please check and re-run!"
    )
    return(fail.sample)
  } else {
    message("STAR run successfully on all samples!")
    return(NULL)
  }
}

# run STAR on single sample
RunSTARSingle <- function(fq.dir, ref, out.folder = NULL, thread = 4, star.path = NULL,
                          star.paras = "--outBAMsortingThreadN 4 --twopassMode None") {
  # get STAR path
  if (is.null(star.path)) {
    # specify STAR path
    star.path <- Sys.which("STAR")
    if (star.path == "") {
      stop("Can not find STAR automatically, please specify the path!")
    }
  } else {
    star.path <- star.path
  }
  # prepare output folder
  if (is.null(out.folder)) {
    out.folder <- getwd()
  }
  if (!dir.exists(out.folder)) {
    message(out.folder, " does not exist, create automatically!")
    dir.create(out.folder, recursive = TRUE)
  }
  # check the path
  if (!dir.exists(fq.dir)) {
    message(fq.dir, " doesn't exist, please check and re-run!")
    return(basename(fq.dir))
  } else {
    fq.files <- list.files(path = fq.dir, pattern = "fastq.gz$")
    if (length(fq.files) == 0) {
      message("There is no fastq.gz files under ", fq.dir, " , please check and re-run!")
      return(basename(fq.dir))
    } else {
      sample.name <- basename(fq.dir)
      out.folder <- file.path(out.folder, sample.name, "")
      # prepare reads
      pair.r1 <- grep(pattern = paste0(sample.name, "_1.fastq.gz"), x = fq.files, value = T)
      pair.r2 <- grep(pattern = paste0(sample.name, "_2.fastq.gz"), x = fq.files, value = T)
      single.read <- grep(pattern = paste0(sample.name, ".fastq.gz"), x = fq.files, value = T)
      # check additional paras
      if (grepl(pattern = "--runThreadN|--genomeDir|--readFilesIn|--outSAMtype|--readFilesCommand|--quantMode|--outFileNamePrefix", star.paras)) {
        message(
          star.paras, " overlaps with built-in paras (--runThreadN, --genomeDir, --readFilesIn, --outSAMtype,",
          " --readFilesCommand, --quantMode, --outFileNamePrefix), please check and re-run!"
        )
      }
      # prepare command
      if (length(pair.r1) == 1 && length(pair.r2) == 1) {
        message("Detected pair-end fastqs!")
        # STAR command
        star.cmd <- paste(
          star.path, "--runThreadN", thread, "--genomeDir", ref,
          "--readFilesIn", file.path(fq.dir, pair.r1), file.path(fq.dir, pair.r2),
          "--outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --quantMode GeneCounts",
          "--outFileNamePrefix", out.folder, star.paras, "&& mv",
          paste0(out.folder, "ReadsPerGene.out.tab"), paste0(out.folder, sample.name, ".txt")
        )
      } else if (length(single.read) == 1) {
        message("Detected single-end fastq!")
        # STAR command
        star.cmd <- paste(
          star.path, "--runThreadN", thread, "--genomeDir", ref,
          "--readFilesIn", file.path(fq.dir, single.read),
          "--outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --quantMode GeneCounts",
          "--outFileNamePrefix", out.folder, star.paras, "&& mv",
          paste0(out.folder, "ReadsPerGene.out.tab"), paste0(out.folder, sample.name, ".txt")
        )
      } else {
        message(
          "There is no valid fastq.gz files under ", fq.dir,
          ". For pair-end,the fastq files should be: ", paste0(sample.name, "_[12].fastq.gz"),
          ". For single-end, the fastq file should be: ", paste0(sample.name, ".fastq.gz"),
          ". Please check and re-run!"
        )
        return(sample.name)
      }
      # run command
      message(paste("Calling STAR:", star.cmd))
      star.status <- system(star.cmd, intern = TRUE)
      star.status.code <- attr(star.status, "status")
      if (!is.null(star.status.code)) {
        star.status.msg <- paste(star.status, collapse = " ")
        warning(
          "Run STAR error on: ", sample.name, ". Error message :",
          star.status.msg, " .Please check and re-run!"
        )
        return(sample.name)
      } else {
        message("Finish STAR: ", sample.name)
        return(NULL)
      }
    }
  }
}
