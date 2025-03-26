#' Plot Barcode CPM (or occurrence) across samples in a Heatmap
#'
#' `plotBarcodeHeatmap()` visualizes Barcode output across samples with
#'  the option to include sample annotations. The Heatmap can display either:
#'   * CPM: log2(CPM+1)
#'   * occurrence: 0 or 1
#'
#' @param barbieQ A `barbieQ` object created by the [createBarbieQ] function.
#' @param barcodeMetric A string indicating what to visualize.
#'  Defaults to 'CPM'. Options include: 'CPM' and 'occurrence'.
#' @param splitSamples A logical barcodeMetric deciding whether to split samples
#'  into slices. Defaults to FALSE.
#' @param sampleMetadata A `matrix`, `data.frame` or `DataFrame` of sample conditions,
#'  where each factor is represented in a separate column. Defaults to NULL,
#'  in which case sample conditions are inherited from `colData(barbieQ)$sampleMetadata`.
#' @param sampleGroup A string representing the name of a factor from the
#'  sample conditions passed by `barbieQ` or `sampleMetadata`, or a vector of
#'  sample conditions, indicating the primary factor to split sample slices.
#' @param barcodeAnnotation A row Annotation object created by the
#'  [ComplexHeatmap::rowAnnotation] function. Defaults to NULL, which means
#'  no Barcode annotation will be displayed in the Heatmap.
#' @param sampleAnnotation A column Annotation object created by the
#'  [ComplexHeatmap::HeatmapAnnotation] function. Defaults to NULL, which means
#'  the sample annotations are generated from the sample conditions provided by
#'  `barbieQ` and `sampleMetadata`.
#'
#' @return A `Heatmap` S4 object displaying the heatmap of Barcode
#'  data across samples, optionally annotated with sample and Barcode
#'  information.
#'
#' @export
#'
#' @import ComplexHeatmap
#' @importFrom circlize colorRamp2
#' @importFrom grid gpar
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
#' plotBarcodeHeatmap(myBarbieQ)
plotBarcodeHeatmap <- function(barbieQ, barcodeMetric = "CPM", splitSamples = FALSE, sampleMetadata = NULL,
    sampleGroup = NULL, barcodeAnnotation = NULL, sampleAnnotation = NULL) {

    ## check which barcodeMetric to visualize
    barcodeMetric <- match.arg(barcodeMetric, c("CPM", "occurrence"))

    ## extract sampleMetadata and primary effector based on arguments
    sampleMetadata <- extractSampleMetadataAndPrimaryFactor(barbieQ = barbieQ, sampleMetadata = sampleMetadata,
        sampleGroup = sampleGroup)
    primaryFactor <- S4Vectors::metadata(sampleMetadata)$primaryFactor

    ## set the primary effector as sample splitter displaying at bottom
    topFactor <- sampleMetadata[, primaryFactor, drop = FALSE]
    ## the rest of effectors displayed at top
    bottomFactor <- sampleMetadata[, setdiff(colnames(sampleMetadata), primaryFactor), drop = FALSE]

    groupAnnotation <- NULL

    if (is.null(sampleAnnotation)) {
        sampleAnnotation <- HeatmapAnnotation(df = bottomFactor, annotation_name_side = "left",
            annotation_name_gp = grid::gpar(fontsize = 10), col = S4Vectors::metadata(barbieQ)$factorColors)
        groupAnnotation <- HeatmapAnnotation(df = topFactor, annotation_name_side = "left",
            annotation_name_gp = grid::gpar(fontsize = 10), col = S4Vectors::metadata(barbieQ)$factorColors)
    }


    if (splitSamples) {
        splitBy <- topFactor
    } else {
        splitBy <- NULL
    }

    ## choose barcodeMetric to be visualised
    if (barcodeMetric == "CPM") {
        mat <- log2(SummarizedExperiment::assays(barbieQ)$CPM + 1) %>%
            as.matrix()
        matTitle <- "log2 CPM+1"
        colorFun <- circlize::colorRamp2(c(min(mat), mean(mat), max(mat)), c("blue", "white",
            "red"))
    } else {
        mat <- (SummarizedExperiment::assays(barbieQ)$occurrence + 1 - 1) %>%
            as.matrix()
        matTitle <- "occurrence"
        colorFun <- structure(c(2, 4), names = c("1", "0"))
    }

    hp <- Heatmap(mat, name = matTitle, width = unit(6, "cm"), height = unit(6, "cm"), cluster_rows = TRUE,
        cluster_columns = TRUE, show_row_names = FALSE, show_column_names = FALSE, column_title = paste0(ncol(mat),
            " Samples"), row_title = paste0(nrow(mat), " Barcodes"), col = colorFun, right_annotation = barcodeAnnotation,
        top_annotation = groupAnnotation, bottom_annotation = sampleAnnotation, column_split = splitBy,
        cluster_column_slices = FALSE)

    return(hp)
}
