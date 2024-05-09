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
#' RunSTAR(
#'   sample.dir = "/path/to/fastq",
#'   ref = "/path/to/star/index",
#'   out.folder = "/path/to/star/mapping",
#'   star.path = "/path/to/bin/STAR"
#' )
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

#' Pipe FASTQ files to SeuratObject and DESeqDataSet.
#'
#' @param sample.dir Directory contains all samples.
#' @param ref Path of folder containing 10x-compatible transcriptome \code{\link{RunCellRanger}}
#' STAR \code{\link{RunSTAR}} reference.
#' @param method Mapping methods, choose from CellRanger (10x Genomics) and STAR (Smart-seq2 or bulk RNA-seq).
#' Default: CellRanger.
#' @param localcores The max cores used \code{\link{RunCellRanger}}. Default: 4.
#' @param localmem The max memory (GB) used \code{\link{RunCellRanger}}. Default: 16.
#' @param thread The number of threads to use \code{\link{RunSTAR}}. Default: 4.
#' @param out.folder Output folder. Default: NULL (current working directory).
#' @param st.path Path to \code{STAR} or \code{cellranger}. Default: NULL (conduct automatic detection).
#' @param st.paras Parameters for \code{STAR} or \code{cellranger}.
#' Default: "--chemistry=auto --jobmode=local".
#' @param merge Logical, whether to merge the SeuratObjects, use when \code{method} is CellRanger. Default: TRUE.
#' @param count.col Column contains used count data (2: unstranded; 3: \code{stranded=yes}; 4: \code{stranded=reverse}),
#' use when \code{method} is STAR. Default: 2.
#' @param meta.data Dataframe contains sample information for DESeqDataSet, use when \code{method} is STAR. Default: NULL.
#' @param fmu Column of \code{meta.data} contains group information. Default: NULL.
#'
#' @return SeuratObject, DESeqDataSet or NULL.
#' @importFrom magrittr %>%
#' @importFrom Seurat Read10X CreateSeuratObject
#' @importFrom methods new
#' @importFrom data.table fread
#' @importFrom magrittr %>%
#' @importFrom stats formula
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @export
#'
#' @examples
#' \dontrun{
#' # run CellRanger (10x Genomics)
#' seu <- Fastq2R(
#'   sample.dir = "/path/to/fastq",
#'   ref = "/path/to/10x/ref",
#'   method = "CellRanger",
#'   out.folder = "/path/to/results",
#'   st.path = "/path/to/cellranger"
#' )
#' # run STAR (Smart-seq2 or bulk RNA-seq)
#' deobj <- Fastq2R(
#'   sample.dir = "/path/to/fastq",
#'   ref = "/path/to/star/ref",
#'   method = "STAR",
#'   out.folder = "/path/to/results",
#'   st.path = "/path/to/STAR",
#'   st.paras = "--outBAMsortingThreadN 4 --twopassMode None"
#' )
#' }
Fastq2R <- function(sample.dir, ref, method = c("CellRanger", "STAR"), localcores = 4, localmem = 16, thread = 4,
                    out.folder = NULL, st.path = NULL, st.paras = "--chemistry=auto --jobmode=local",
                    merge = TRUE, count.col = 2, meta.data = NULL, fmu = NULL) {
  # check parameters
  method <- match.arg(arg = method)
  if (method == "CellRanger") {
    res <- RunCellRanger(
      sample.dir = sample.dir, ref = ref, localcores = localcores,
      localmem = localmem, out.folder = out.folder,
      cr.path = st.path, cr.paras = st.paras
    )
    if (is.null(res)) {
      all.cr.folder <- dir(out.folder, full.names = TRUE)
      all.samples.folder <- file.path(all.cr.folder, "outs", "filtered_feature_bc_matrix")
      # check file
      valid.samples.folder <- Check10XFiles(folders = all.samples.folder, gene2feature = TRUE)
      if (length(valid.samples.folder) == 0) {
        stop("No valid sample folder detected under ", out.folder, ". Please check!")
      }
      # load to seurat
      seu.list <- sapply(valid.samples.folder, function(x) {
        x.mat <- Seurat::Read10X(data.dir = x)
        seu.obj <- Seurat::CreateSeuratObject(counts = x.mat, project = basename(x))
        seu.obj
      })
      # merge SeuratObject
      if (isTRUE(merge)) {
        out.obj <- mergeExperiments(seu.list)
      } else {
        out.obj <- seu.list
      }
      return(out.obj)
    } else {
      message("Some samples failed to run, skipping loading into Seurat!")
      return(NULL)
    }
  } else if (method == "STAR") {
    res <- RunSTAR(
      sample.dir = sample.dir, ref = ref, out.folder = out.folder, thread = thread,
      star.path = st.path, star.paras = st.paras
    )
    if (is.null(res)) {
      all.txt <- list.files(path = out.folder, pattern = ".txt$", recursive = TRUE, full.names = TRUE)
      # read files
      count.list <- lapply(
        all.txt,
        function(x) {
          sample.count <- data.table::fread(file = x, select = c(1, count.col)) %>% as.data.frame()
          sample.count <- sample.count[!sample.count[[1]] %in% c("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous"), ]
          colnames(sample.count) <- c("GeneName", gsub(pattern = "(.*).txt", replacement = "\\1", x = basename(x)))
          sample.count
        }
      )
      # create count matrix
      count.mat <- Reduce(f = function(x, y) {
        merge.mat <- merge(x, y, by = "GeneName", all = T)
      }, x = count.list)
      rownames(count.mat) <- count.mat$GeneName
      count.mat$GeneName <- NULL
      # loadding into DESeq2
      de.obj <- Loading2DESeq2(mat = count.mat, meta = meta.data, fmu = fmu)
      return(de.obj)
    } else {
      message("Some samples failed to run, skipping loading into DESeq2!")
      return(NULL)
    }
  }
}
