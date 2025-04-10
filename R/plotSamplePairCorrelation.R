#' Plot sample pairwise correlation in a Heatmap
#'
#' `plotSamplePairCorrelation()` visualizes the pairwise correlation of
#'  CPM values across samples in a `barbieQ` object,
#'  using a heatmap with a checkerboard-like pattern.
#'
#' @param barbieQ A `barbieQ` object created by the [createBarbieQ] function.
#' @param sampleOrder A character vector of names of the factors in the
#'  sample conditions provided by `barbieQ` or `sampleMetadata`,
#'  specifying the order in which samples should be arranged.
#'  Defaults to the original order of factors in the sample conditions.
#' @param sampleMetadata A `matrix`, `data.frame` or `DataFrame` of sample conditions,
#'  where each factor is represented in a separate column. Defaults to NULL,
#'  in which case sample conditions are inherited from `colData(barbieQ)$sampleMetadata`.
#' @param sampleGroup A string representing the name of a factor from the
#'  sample conditions passed by `barbieQ` or `sampleMetadata`, or a vector of
#'  sample conditions, indicating the primary factor to split sample slices.
#' @param method A string specifying the correlation method to use.
#'  Defaults to 'pearson'. Options include: 'pearson', 'spearman'.
#'
#' @return A 'Heatmap' S4 object displaying the pairwise correlation between
#'  samples in a Heatmap.
#'
#' @export
#'
#' @import ComplexHeatmap
#' @importFrom circlize colorRamp2
#' @importFrom grid gpar
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr pull
#' @importFrom dplyr arrange
#' @importFrom dplyr across
#' @importFrom dplyr row_number
#' @importFrom dplyr all_of
#' @importFrom stats cor
#' @import data.table
#' @importClassesFrom S4Vectors DataFrame
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment assays
#' @importFrom S4Vectors metadata
#'
#' @examples
#' ## sample conditions and color palettes
#' sampleConditions <- data.frame(
#'   Treat = factor(rep(c('ctrl', 'drug'), each = 6)),
#'   Time = rep(rep(seq_len(2), each = 3), 2)
#' )
#' conditionColor <- list(
#'   Treat = c(ctrl = '#999999', drug = '#112233'),
#'   Time = c('1' = '#778899', '2' = '#998877')
#' )
#' ## Barcode count data
#' nbarcodes <- 50
#' nsamples <- 12
#' barcodeCount <- abs(matrix(10, nbarcodes, nsamples))
#' barcodeCount[seq(21, 50), ] <- 0.0001
#' rownames(barcodeCount) <- paste0('Barcode', seq_len(nbarcodes))
#' ## create a `barbieQ` object
#' myBarbieQ <- createBarbieQ(barcodeCount, sampleConditions, conditionColor)
#' plotSamplePairCorrelation(myBarbieQ)
plotSamplePairCorrelation <- function(barbieQ, sampleOrder = NULL, sampleMetadata = NULL,
    sampleGroup = NULL, method = "pearson") {
    ## extract sampleMetadata and primary effector based on arguments
    sampleMetadata <- extractSampleMetadataAndPrimaryFactor(barbieQ = barbieQ, sampleMetadata = sampleMetadata,
        sampleGroup = sampleGroup)
    primaryFactor <- S4Vectors::metadata(sampleMetadata)$primaryFactor

    ## set the primary effector as sample splitter displaying at bottom
    topsampleMetadata <- sampleMetadata[, primaryFactor, drop = FALSE]
    ## the rest of effectors displayed at top
    bottomsampleMetadata <- sampleMetadata[, dplyr::setdiff(colnames(sampleMetadata), primaryFactor),
        drop = FALSE]
    ## check sampleOrder
    if (is.null(sampleOrder)) {
        sampleOrder <- c(colnames(topsampleMetadata), colnames(bottomsampleMetadata))
    }
    ## check method
    method <- match.arg(method, c("pearson", "spearman"))

    ## reorder the columns of the data.frame based on the specified order
    sampleMetadata <- sampleMetadata[, sampleOrder, drop = FALSE]
    ## extract sample order
    rowOrder <- as.data.frame(sampleMetadata) %>%
        dplyr::mutate(rowNumber = dplyr::row_number()) %>%
        dplyr::arrange(dplyr::across(dplyr::all_of(sampleOrder))) %>%
        dplyr::pull(rowNumber)

    ## compute annotation obejct
    sampleAnnotationColumn <- HeatmapAnnotation(df = sampleMetadata, annotation_name_side = "right",
        annotation_name_gp = grid::gpar(fontsize = 10), col = S4Vectors::metadata(barbieQ)$factorColors)
    sampleAnnotationRow <- rowAnnotation(df = sampleMetadata, annotation_name_side = "bottom",
        annotation_name_gp = grid::gpar(fontsize = 10), col = S4Vectors::metadata(barbieQ)$factorColors,
        show_legend = FALSE)
    ## calculate correlation
    corMat <- stats::cor(log2(SummarizedExperiment::assays(barbieQ)$CPM + 1), method = method)
    message("displaying ", method, " correlation coefficient between samples on Barcode log2 CPM+1.")

    hp <- Heatmap(corMat, name = "corCoef", width = unit(6, "cm"), height = unit(6, "cm"),
        col = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red")), heatmap_legend_param = list(at = c(-1,
            0, 1)), cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE,
        show_column_names = FALSE, top_annotation = sampleAnnotationColumn, left_annotation = sampleAnnotationRow,
        column_order = rowOrder, row_order = rowOrder)

    return(hp)
}
