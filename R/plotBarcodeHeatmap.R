#' Plot Barcode CPM (or occurrence) across samples in a Heatmap
#'
#' `plotBarcodeHeatmap()` visualizes Barcode output across samples with
#'  the option to include sample annotations. The Heatmap can display either:
#'   * CPM: log2(CPM+1)
#'   * occurrence: 0 or 1
#'
#' @param barbieQ A `barbieQ` object created by the [createBarbieQ] function.
#' @param colorMapTo A string indicating mapping matrix color to barcode
#'  proportion or occurrence, with or without transformation.
#'  Defaults to 'asin-sqrt proportion'. Options include: "asin-sqrt proportion", 
#'  "logit proportion", "proportion", and "occurrence".
#' @param showRawProportion A logical value, deciding whether to label matrix color
#'  legend by raw proportion when color is mapped to transformed data specified by
#'  `colorMapTo`. Default to TRUE.
#' @param splitSamples A logical value, deciding whether to split samples
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
#' @importFrom stats setNames
#' @importClassesFrom S4Vectors DataFrame
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment assays
#' @importFrom SummarizedExperiment assay
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
plotBarcodeHeatmap <- function(barbieQ, colorMapTo = "asin-sqrt proportion", showRawProportion = TRUE,
    splitSamples = FALSE, sampleMetadata = NULL, sampleGroup = NULL, barcodeAnnotation = NULL, sampleAnnotation = NULL) {

    ## check which colorMapTo to visualize
    colorMapTo <- match.arg(colorMapTo, c("asin-sqrt proportion", "logit proportion", "proportion", "occurrence"))
    
    ## ---- Part1: Annotate columns (samples) by conditions ----
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
    
    ## ---- Part2: map the main matrix color to the input data ----
    ## get colorMapTo alias
    aliasMetric <- stats::setNames(c("asin(sqrt(prop.))", "logit(prop.)", "proportion", "occurrence"),
      c("asin-sqrt proportion", "logit proportion", "proportion", "occurrence"))
    ## set mat name
    matName <- aliasMetric[colorMapTo]

    ## visualize binary values when colorMapTo == `occurrence`
    if (colorMapTo == "occurrence") {
      mat <- (SummarizedExperiment::assays(barbieQ)$occurrence + 1 - 1) %>% as.matrix()
      colorFun <- structure(c(2, 4), names = c("1", "0"))
    } else if(colorMapTo == "proportion"){
      ## alternatively visualize raw `proportion`, ...
      mat <- SummarizedExperiment::assays(barbieQ)$proportion %>% as.matrix()
    } else if(colorMapTo == "asin-sqrt proportion") {
      mat <- SummarizedExperiment::assays(barbieQ)$proportion %>% as.matrix()
      mat <- asin(sqrt(mat))
    } else if(colorMapTo == "logit proportion") {
      countPlus <- SummarizedExperiment::assay(barbieQ) + 0.5
      proportionPlus <- countPlus / colSums(countPlus)
      mat <- log((proportionPlus)/(1 - proportionPlus)) %>% as.matrix()
    }
    
    ## choose color mapping function for the continuous proportion, ...
    if (colorMapTo != "occurrence") {
      colorFun <- circlize::colorRamp2(c(min(mat), mean(mat), max(mat)), c("blue", "white","red"))
    } 
    
    ## draw the heatmap to get the legend breaks
    ht <- ComplexHeatmap::Heatmap(mat, name = matName, width = unit(6, "cm"), height = unit(6, "cm"), cluster_rows = TRUE,
        cluster_columns = TRUE, show_row_names = FALSE, show_column_names = FALSE, column_title = paste0(ncol(mat),
            " Samples"), row_title = paste0(nrow(mat), " Barcodes"), col = colorFun, right_annotation = barcodeAnnotation,
        top_annotation = groupAnnotation, bottom_annotation = sampleAnnotation, column_split = splitBy,
        cluster_column_slices = FALSE)
    
    ## ---- Part 3: Reverse the legend breaks back to the original proportion ----
    ## ---- when colorMapto == the transformed, and showRawProp is true ----
    if((colorMapTo == "asin-sqrt proportion" || colorMapTo == "logit proportion") && showRawProportion) {
      # Build the heatmap object to access breaks
      ht_built <- draw(ht)
      # Get the auto-generated breaks from the legend
      lgd <- ht_built@ht_list[[matName]]@matrix_color_mapping
      breaks_trans <- lgd@levels  # transformed break positions
      
      # Convert breaks back to original scale using inverse: sin^2(x)
      if(colorMapTo == "asin-sqrt proportion") {
        breaks_orig <- sin(breaks_trans)^2
      } else {
        ## inverse: e^x/(1+e^x) ---- but this is reversed to the proportion based on countPlus
        breaks_orig <- exp(breaks_trans)/(1+ exp(breaks_trans))
      }
      
      # re-draw heatmap with custom legend labels but same break positions
      ht <- ComplexHeatmap::Heatmap(
        mat, name = "proportion", width = unit(6, "cm"), height = unit(6, "cm"), col = colorFun, 
        cluster_rows = TRUE, cluster_columns = TRUE, show_row_names = FALSE, show_column_names = FALSE, 
        column_title = paste0(ncol(mat), " Samples"), row_title = paste0(nrow(mat), " Barcodes"), 
        right_annotation = barcodeAnnotation, top_annotation = groupAnnotation, 
        bottom_annotation = sampleAnnotation, column_split = splitBy, cluster_column_slices = FALSE,
        ##same transformed positions, but with original scale labels
        heatmap_legend_param = list(at = breaks_trans, labels = format(round(breaks_orig, 4), scientific = FALSE))
        )
      message(paste0("matrix color is mapped to `", colorMapTo, "` but labeled by raw proportion."))
    }

    return(ht)
}
