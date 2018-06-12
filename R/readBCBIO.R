#' readBCBIO
#'
#' This function is .
#' @param projectDir
#' @keywords dataLoading
#' @param projectDir path to final bcbio project folder. Normally
#'   named with date format.
#' @export
#' @examples
#' readBCBIO()

readBCBIO <- function(projectDir){
    library(dplyr)
    projectDir <- normalizePath(projectDir)
    sampleMetadata <- read.csv(file.path(projectDir,"metadata.csv"),row.names = 1)
    metrics <- import(file.path(projectDir,"multiqc","multiqc_data","multiqc_bcbio_metrics.txt"))
    metrics <- janitor::clean_names(metrics,case = "small_camel")
    sampleDirs = file.path(projectDir,"..", rownames(sampleMetadata))
    tx2genes_file <- file.path(projectDir,"tx2gene.csv")
    salmon_files = file.path(sampleDirs, "salmon", "quant.sf")
    sf_files = salmon_files
    names(sf_files) = rownames(sampleMetadata)
    # FIX check files exists and match sampleMetadata
    tx2gene = read.table(tx2genes_file, sep=",", row.names=NULL, header=FALSE)
    txi.salmon = tximport(sf_files,
                          type="salmon",
                          tx2gene=tx2gene,
                          countsFromAbundance="lengthScaledTPM")
    rawCounts = round(data.frame(txi.salmon$counts, check.names=FALSE),0) %>%
        as.matrix()
    colData <- sampleMetadata %>% as.data.frame() %>%
        .[colnames(rawCounts),, drop = FALSE] %>%
        rownames_to_column() %>%
        mutate_all(as.factor) %>%
        mutate_all(droplevels) %>%
        column_to_rownames() %>%
        as("data.frame")

    metadata <- list(metrics = metrics,
                     countsFromAbundance = txi.salmon$countsFromAbundance)
    se <- SummarizedExperiment(assays = list(raw = rawCounts,
                                             tpm = txi.salmon$abundance,
                                             length = txi.salmon$length),
                               colData = colData,
                               metadata = metadata)
    invisible(se)
}

