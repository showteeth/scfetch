#' Extract Metadata of CELLxGENE Datasets.
#'
#' @return Dataframe contains all available datasets.
#' @importFrom magrittr %>%
#' @importFrom curl curl_fetch_memory
#' @importFrom jsonlite fromJSON flatten
#' @importFrom data.table rbindlist
#' @export
#'
#' @examples
#' # # all available datasets
#' # cellxgene.meta = ExtractCELLxGENEMeta()
ExtractCELLxGENEMeta <- function() {
  # urls
  cellxgene.base.url <- "https://api.cellxgene.cziscience.com/dp/v1/"
  cellxgene.collections.url <- paste0(cellxgene.base.url, "collections/")
  # extract all collections
  cellxgene.collections.content <- URLRetrieval(cellxgene.collections.url)
  cellxgene.collections.df <- cellxgene.collections.content$collections
  colnames(cellxgene.collections.df) <- c(
    "collection_created_at", "collection_id",
    "collection_owner", "collection_visibility"
  )
  # extract all datasets
  cellxgene.collections.list <- lapply(1:nrow(cellxgene.collections.df), function(x) {
    cellxgene.collection.df <- cellxgene.collections.df[x, ]
    rownames(cellxgene.collection.df) <- NULL
    cellxgene.sc.url <- paste0(cellxgene.collections.url, cellxgene.collections.df[x, "collection_id"])
    cellxgene.sc.content <- URLRetrieval(cellxgene.sc.url)
    cellxgene.sc.datasets <- jsonlite::flatten(cellxgene.sc.content$datasets)
    # add metadata
    cellxgene.collection.df$title <- cellxgene.sc.content$name
    cellxgene.collection.df$description <- cellxgene.sc.content$description
    cellxgene.collection.df$contact <- cellxgene.sc.content$contact_name
    cellxgene.collection.df$contact_email <- cellxgene.sc.content$contact_email
    cellxgene.collection.df <- cellxgene.collection.df[c(
      "title", "description", "contact", "contact_email",
      "collection_created_at", "collection_id",
      "collection_owner", "collection_visibility"
    )]
    cellxgene.collection.final <- cbind(cellxgene.collection.df, cellxgene.sc.datasets) %>% as.data.frame()
    return(cellxgene.collection.final)
  })
  # create all datasets dataframe
  cellxgene.collections.datasets.df <- data.table::rbindlist(cellxgene.collections.list, fill = TRUE) %>%
    as.data.frame()
  # modify columns
  cellxgene.collections.datasets.df <-
    PasteAttrCXG(
      df = cellxgene.collections.datasets.df,
      attr = c(
        "assay", "cell_type", "organism", "self_reported_ethnicity", "sex", "tissue",
        "disease", "development_stage"
      ), col = "label"
    )
  cellxgene.collections.datasets.df <-
    PasteAttrCXG(
      df = cellxgene.collections.datasets.df,
      attr = c("dataset_deployments"), col = "url"
    )
  cellxgene.collections.datasets.df <-
    PasteAttr(df = cellxgene.collections.datasets.df, attr = c("batch_condition", "suspension_type", "donor_id"))
  # add h5ad and rds information
  cellxgene.collections.datasets.list <- lapply(1:nrow(cellxgene.collections.datasets.df), function(x) {
    x.df <- cellxgene.collections.datasets.df[x, ]
    x.df$dataset_id <- unique(x.df$dataset_assets[[1]]$dataset_id)
    if ("RDS" %in% unique(x.df$dataset_assets[[1]]$filetype)) {
      x.rds.idx <- which(x.df$dataset_assets[[1]]$filetype == "RDS")
      x.df$rds_id <- x.df$dataset_assets[[1]]$id[x.rds.idx]
      x.df$rds_s3_uri <- x.df$dataset_assets[[1]]$s3_uri[x.rds.idx]
      x.df$rds_user_submitted <- x.df$dataset_assets[[1]]$user_submitted[x.rds.idx]
    } else {
      x.df$rds_id <- NA
      x.df$rds_s3_uri <- NA
      x.df$rds_user_submitted <- NA
    }
    if ("H5AD" %in% unique(x.df$dataset_assets[[1]]$filetype)) {
      x.h5ad.idx <- which(x.df$dataset_assets[[1]]$filetype == "H5AD")
      x.df$h5ad_id <- x.df$dataset_assets[[1]]$id[x.h5ad.idx]
      x.df$h5ad_s3_uri <- x.df$dataset_assets[[1]]$s3_uri[x.h5ad.idx]
      x.df$h5ad_user_submitted <- x.df$dataset_assets[[1]]$user_submitted[x.h5ad.idx]
    } else {
      x.df$h5ad_id <- NA
      x.df$h5ad_s3_uri <- NA
      x.df$h5ad_user_submitted <- NA
    }
    x.df
  })
  # final dataframe
  cellxgene.collections.datasets.final <- data.table::rbindlist(cellxgene.collections.datasets.list, fill = TRUE) %>%
    as.data.frame()
  return(cellxgene.collections.datasets.final)
}
