#' Prepare Dataframe with Zenodo DOI.
#'
#' @param doi Zenodo DOI, should start with "10.5281/zenodo.".
#' @param file.ext The valid file extension for download. When NULL, use all files. Default: c("rdata", "h5ad").
#'
#' @return Dataframe contains files with valid extension in given Zenodo DOI.
#' @importFrom magrittr %>%
#' @importFrom curl curl_fetch_memory
#' @importFrom jsonlite fromJSON
#' @importFrom dplyr filter
#' @export
#'
#' @examples
#' # zebrafish.df = ExtractZenodoMeta(doi = "10.5281/zenodo.7243603")
#' # ExtractZenodoMeta(doi = "10.5281/zenodo.48065") # Restricted Access
ExtractZenodoMeta <- function(doi, file.ext = c("rdata", "h5ad")) {
  # prepare link
  record.id <- gsub(pattern = "10.5281/zenodo.", replacement = "", x = doi, fixed = TRUE)
  record.api <- paste0("https://zenodo.org/api/records/", record.id)
  # download page info
  record.page <- curl::curl_fetch_memory(record.api)
  # convert to string
  record.content <- jsonlite::fromJSON(rawToChar(record.page$content))
  # access status
  record.access <- record.content$metadata$access_right
  if (record.access == "open") {
    # extract files
    record.files <- record.content$files
    # filter files
    if (is.null(file.ext)) {
      record.files.used <- record.files
    } else {
      record.files.used <- record.files %>% dplyr::filter(type %in% file.ext)
    }
    # check the data
    if (nrow(record.files.used) == 0) {
      return(NULL)
    }
    # prepare md5sum
    record.files.used$checksum <- gsub(pattern = "md5:", replacement = "", record.files.used$checksum)
    record.files.used.final <- data.frame(
      title = record.content$metadata$title, description = record.content$metadata$description,
      url = record.files.used$links$self, filename = basename(record.files.used$links$self),
      md5 = record.files.used$checksum, license = record.content$metadata$license$id
    )
  } else {
    message(doi, " is not open access!")
    record.files.used.final <- data.frame()
  }
  return(record.files.used.final)
}

#' Download Data with Zenodo DOI.
#'
#' @param doi A vector of Zenodo DOIs to download. Default: NULL.
#' @param file.ext The valid file extension for download. When NULL, use all files. Default: c("rdata", "rds", "h5ad").
#' @param doi.df DOI dataframe for download. This is useful when something wrong happens in downloading
#' (e.g. MD5 verification failure, \code{DownloadZenodo} will return a dataframe contains failure terms.). Default: NULL.
#' It is required to provide either \code{doi} or \code{doi.df}.
#' @param out.folder The output folder. Default: NULL (current working directory).
#' @param timeout Maximum request time. Default: 1000.
#' @param quiet Logical value, whether to show downloading progress. Default: FALSE (show).
#' @param parallel Logical value, whether to download parallelly. Default: TRUE. When "libcurl" is available for \code{download.file},
#' the parallel is done by default (\code{parallel} can be FALSE).
#'
#' @return When successful, NULL. When MD5 verification failure, a dataframe contains failure terms.
#' @importFrom magrittr %>%
#' @importFrom curl curl_fetch_memory
#' @importFrom jsonlite fromJSON
#' @importFrom dplyr filter
#' @importFrom parallel detectCores mclapply
#' @importFrom utils download.file
#' @importFrom tools md5sum
#' @export
#'
#' @examples
#' # zebrafish.dois = ParseZenodo(doi = c("1111", "10.5281/zenodo.7243603", "10.5281/zenodo.7244441"), file.ext = c("rdata", "rds"),
#' #                              out.folder = "/path/to/outfoder")
ParseZenodo <- function(doi = NULL, file.ext = c("rdata", "rds", "h5ad"), doi.df = NULL, out.folder = NULL, timeout = 1000,
                        quiet = FALSE, parallel = TRUE) {
  if (!is.null(doi.df)) {
    doi.df <- doi.df
  } else if (!is.null(doi)) {
    # check doi
    doi.status <- startsWith(x = doi, prefix = "10.5281/zenodo.")
    if (!all(doi.status)) {
      wrong.doi <- doi[!doi.status]
      message(paste0(wrong.doi, collapse = ", "), " are not valid dois, please check!")
      doi <- doi[doi.status]
    }
    # file extension to lower
    file.ext <- tolower(file.ext)
    # prepare data frame
    doi.list <- lapply(doi, function(x) {
      ExtractZenodoMeta(doi = x, file.ext = file.ext)
    })
    doi.df <- do.call(rbind, doi.list)
  } else {
    stop("Please provide either doi or doi.df!")
  }
  # prepare output folder
  if (is.null(out.folder)) {
    out.folder <- getwd()
  }
  doi.df$filename <- file.path(out.folder, doi.df$filename)
  # set timeout
  env.timeout <- getOption("timeout")
  options(timeout = timeout)
  message("Start downloading!")
  if (isTRUE(parallel)) {
    # prepare cores
    cores.used <- min(parallel::detectCores(), nrow(doi.df))
    down.status <- parallel::mclapply(X = 1:nrow(doi.df), FUN = function(x) {
      utils::download.file(url = doi.df[x, "url"], destfile = doi.df[x, "filename"], quiet = quiet, mode = "wb")
    }, mc.cores = cores.used)
  } else {
    down.status <- utils::download.file(url = doi.df$url, destfile = doi.df$filename, quiet = quiet, mode = "wb")
  }
  message("Finish downloading!")
  # restore timeout
  options(timeout = env.timeout)
  # check the md5sum
  down.md5 <- tools::md5sum(doi.df$filename)
  raw.md5 <- doi.df$md5
  names(raw.md5) <- doi.df$filename
  md5.check <- down.md5 == raw.md5
  if (all(md5.check)) {
    message("Download and MD5 verification successful!")
    return(NULL)
  } else {
    wrong.md5.df <- doi.df[!md5.check, ]
    # restore filename
    wrong.md5.df$filename <- basename(wrong.md5.df$filename)
    message("MD5 verification failure for: ", paste0(wrong.md5.df$filename, collapse = ", "), " . You can re-run ParseZenodo with doi.df!")
    return(wrong.md5.df)
  }
}
