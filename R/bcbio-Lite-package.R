#' bcbioLite
#'
#' @import SummarizedExperiment
#' @import tximport
#' @import janitor
#' @import tibble
#' @import dplyr
#' @importFrom dplyr "%>%"
#' @import rio
#' @import DESeq2
#' @importFrom S4Vectors DataFrame
#'
"_PACKAGE"

# Alias required for help links in downstream packages
#' @importFrom SummarizedExperiment colData
#' @export
SummarizedExperiment::colData

#' @importFrom SummarizedExperiment rowData
#' @export
SummarizedExperiment::rowData