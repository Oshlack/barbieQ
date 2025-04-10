#' Plot Barcodes categorized by significant change using a dot plot
#'
#' `plotBarcodePValue()` use a dot plot to visualize the significance
#'  level of each Barcode in differential proportion or occurrence,
#'  as determined by the [testBarcodeSignif] function.
#'  P.values are plotted against other properties of Barcodes
#'  specified by `xAxis`.
#'
#' @param barbieQ A `barbieQ` object created by the [createBarbieQ] function,
#'  updated with Barcode test results by calling the [testBarcodeSignif] function.
#' @param xAxis A string indicating what to visualise on the x scale of the
#'  dot plot. Options include: 'avgRank', 'totalOcc', 'avgLogCPM',
#'  and 'avgProportion'. Defaults to 'avgRank', representing the average
#'  rank of Barcodes across samples.
#'
#' @return A `ggplot` S3 class object displaying the significance level
#'  against other properties of Barcodes in a dot plot.
#'
#' @export
#'
#' @import ggplot2
#' @importFrom stats setNames
#' @import data.table
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment rowData
#' @importFrom SummarizedExperiment assays
#' @importFrom S4Vectors metadata
#'
#'
#' @examples
#' Block <- c(1, 1, 2, 3, 3, 4, 1, 1, 2, 3, 3, 4)
#' Treat <- factor(rep(c('ctrl', 'drug'), each = 6))
#' Time <- rep(rep(seq_len(2), each = 3), 2)
#' nbarcodes <- 50
#' nsamples <- 12
#' count <- abs(matrix(rnorm(nbarcodes * nsamples), nbarcodes, nsamples))
#' rownames(count) <- paste0('Barcode', seq_len(nbarcodes))
#' barbieQ <- createBarbieQ(count, data.frame(Treat = Treat, Time = Time))
#' testBB <- testBarcodeSignif(barbieQ, sampleGroup = 'Treat')
#' plotBarcodePValue(barbieQ = testBB)
plotBarcodePValue <- function(barbieQ, xAxis = "avgRank") {
    ## extract testing results and information
    statsDf <- SummarizedExperiment::rowData(barbieQ)$testingBarcode
    design <- S4Vectors::metadata(statsDf)$design
    contrastGroups <- S4Vectors::metadata(statsDf)$contrastGroups

    ## extract color code
    colorCode <- S4Vectors::metadata(barbieQ)$factorColors$testingBarcode

    ## define a custom color/shape palette for test results
    customShape <- stats::setNames(c(21, 24, 23), c(contrastGroups, "n.s."))
    ## y axis will be p.value x axis will be optional: total occurrence, average rank,

    ## check xAxis
    xOptions <- c("avgRank", "totalOcc", "avgLogCPM", "avgProp")
    xAxis <- match.arg(xAxis, xOptions)
    xTitle <- stats::setNames(c("Average rank of Barcode across samples", "N. samples in which Barcode occurs",
        "Average Barcode Log2 CPM+1 across samples", "Average Barcode proportion across samples"),
        xOptions)

    ## data.frame for ggplot
    mydata <- data.frame(statsDf, totalOcc = rowSums(SummarizedExperiment::assays(barbieQ)$occurrence),
        avgRank = rowMeans(SummarizedExperiment::assays(barbieQ)$rank), avgLogCPM = log2(rowMeans((SummarizedExperiment::assays(barbieQ)$CPM +
            1))), avgProp = rowMeans(SummarizedExperiment::assays(barbieQ)$proportion))
    ## visualize by ggplot
    p <- ggplot(mydata, aes(x = mydata[, xAxis], y = -log10(adj.P.Val))) + geom_point(aes(color = tendencyTo,
        shape = tendencyTo, fill = tendencyTo), size = 2, stroke = 1) + theme_classic() +
        theme(aspect.ratio = 1) + labs(title = paste0(S4Vectors::metadata(statsDf)$method),
        y = "-log10(adj.P.Value)", x = xTitle[xAxis]) + geom_hline(yintercept = -log10(0.05),
        linetype = "dashed", color = "grey") + scale_color_manual(values = colorCode) +
        scale_shape_manual(values = customShape) + scale_fill_manual(values = alpha(colorCode,
        0.2))
    ## reverse x scale if displaying Barcode rank
    if (xAxis == "avgRank") {
        p <- p + annotate("text", x = min(mydata[, xAxis]) * 1.1, y = -log10(0.05), label = "P.Value = 0.05",
            vjust = 1.5, hjust = 1) + scale_x_reverse()
    } else {
        p <- p + annotate("text", x = max(mydata[, xAxis]) * 0.9, y = -log10(0.05), label = "P.Value = 0.05",
            vjust = 1.5, hjust = 1)
    }

    return(p)
}

#' Plot Barcodes catogorized by significant change using a dot plot
#'
#' `plotBarcodeMA()` use a dot plot to visualize the variation and mean of
#'  each Barcode in differential proportion or occurrence,
#'  as determined by the [testBarcodeSignif] function.
#'
#' @param barbieQ A `barbieQ` object created by the [createBarbieQ] function,
#'  updated with Barcode test results by calling the [testBarcodeSignif] function.
#'
#' @return A `ggplot` S3 class object displaying the significance level
#'  against other properties of Barcodes in a dot plot.
#'
#' @export
#'
#' @import ggplot2
#' @importFrom stats setNames
#' @import data.table
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment rowData
#' @importFrom SummarizedExperiment assays
#' @importFrom S4Vectors metadata
#'
#'
#' @examples
#' Block <- c(1, 1, 2, 3, 3, 4, 1, 1, 2, 3, 3, 4)
#' Treat <- factor(rep(c('ctrl', 'drug'), each = 6))
#' Time <- rep(rep(seq_len(2), each = 3), 2)
#' nbarcodes <- 50
#' nsamples <- 12
#' count <- abs(matrix(rnorm(nbarcodes * nsamples), nbarcodes, nsamples))
#' rownames(count) <- paste0('Barcode', seq_len(nbarcodes))
#' barbieQ <- createBarbieQ(count, data.frame(Treat = Treat, Time = Time))
#' testBB <- testBarcodeSignif(barbieQ, sampleGroup = 'Treat')
#' plotBarcodeMA(barbieQ = testBB)
plotBarcodeMA <- function(barbieQ) {
    ## extract testing results and information
    statsDf <- SummarizedExperiment::rowData(barbieQ)$testingBarcode
    design <- S4Vectors::metadata(statsDf)$design
    contrastGroups <- S4Vectors::metadata(statsDf)$contrastGroups
    method <- S4Vectors::metadata(statsDf)$method
    if(method == "diffProp") {
      transformation <- S4Vectors::metadata(statsDf)$transformation
    }

    ## extract color code
    colorCode <- S4Vectors::metadata(barbieQ)$factorColors$testingBarcode

    ## define a custom color/shape palette for test results
    customShape <- stats::setNames(c(21, 24, 23), c(contrastGroups, "n.s."))

    ## data.frame for ggplot
    mydata <- data.frame(statsDf, totalOcc = rowSums(SummarizedExperiment::assays(barbieQ)$occurrence))

    ## visualize by ggplot
    if (method == "diffProp") {
        p <- ggplot(mydata, aes(x = Amean, y = meanDiff)) + geom_point(aes(color = tendencyTo,
            shape = tendencyTo, fill = tendencyTo), size = 2, stroke = 1) + theme_classic() +
            theme(aspect.ratio = 1) + labs(title = paste0(S4Vectors::metadata(statsDf)$method,
            " MA plot"), y = paste0("Prop. Diff. (", transformation, " trans.)"), 
            x = paste0("Mean Prop. (", transformation, " trans.)")) +
            scale_color_manual(values = colorCode) + scale_shape_manual(values = customShape) +
            scale_fill_manual(values = alpha(colorCode, 0.2))
    } else if (method == "diffOcc") {
        p <- ggplot(mydata, aes(x = totalOcc, y = logOR)) + geom_point(aes(color = tendencyTo,
            shape = tendencyTo, fill = tendencyTo), size = 4, stroke = 1) + theme_classic() +
            theme(aspect.ratio = 1) + labs(title = paste0(S4Vectors::metadata(statsDf)$method,
            " 'MA' plot"), y = "LogOR", x = "Total Occurrence Freq.") + scale_color_manual(values = colorCode) +
            scale_shape_manual(values = customShape) + scale_fill_manual(values = alpha(colorCode,
            0.2))
    } 

    return(p)
}


#' Plot Barcodes categorized by significant change using a Heatmap
#'
#' `plotSignifBarcodeHeatmap()` uses the Heatmap annotations to visualize the
#'  significance level of each Barcode in differential proportion or occurrence,
#'  as determined by the [testBarcodeSignif] function.
#'
#' @param barbieQ A `barbieQ` object created by the [createBarbieQ] function,
#'  updated with Barcode test results by calling the [testBarcodeSignif] function.
#' @param barcodeMetric A string indicating what to visualize.
#'  Defaults to `CPM`. Options include: 'CPM' and 'occurrence'.
#' @param sampleAnnotation A column Annotation object created by the
#'  [ComplexHeatmap::HeatmapAnnotation] function. Defaults to samples annotated
#'  by the groups to be compared.
#'
#' @return A `Heatmap` S4 class object displaying the significance level
#'  of Barcodes in Heatmap annotations.
#'
#' @export
#'
#' @import ComplexHeatmap
#' @importFrom grid gpar
#' @importFrom dplyr setdiff
#' @importFrom stats setNames
#' @import data.table
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment rowData
#' @importFrom SummarizedExperiment assays
#' @importFrom S4Vectors metadata
#'
#' @examples
#' Block <- c(1, 1, 2, 3, 3, 4, 1, 1, 2, 3, 3, 4)
#' Treat <- factor(rep(c('ctrl', 'drug'), each = 6))
#' Time <- rep(rep(seq_len(2), each = 3), 2)
#' nbarcodes <- 50
#' nsamples <- 12
#' count <- abs(matrix(rnorm(nbarcodes * nsamples), nbarcodes, nsamples))
#' rownames(count) <- paste0('Barcode', seq_len(nbarcodes))
#' barbieQ <- createBarbieQ(count, data.frame(Treat = Treat, Time = Time))
#' testBB <- testBarcodeSignif(barbieQ, sampleGroup = 'Treat')
#' plotSignifBarcodeHeatmap(barbieQ = testBB)
plotSignifBarcodeHeatmap <- function(barbieQ, barcodeMetric = "CPM", sampleAnnotation = NULL) {

    ## extract testing results and information
    statsDf <- SummarizedExperiment::rowData(barbieQ)$testingBarcode
    contrastGroups <- S4Vectors::metadata(statsDf)$contrastGroups
    method <- S4Vectors::metadata(statsDf)$method
    ## extract design based on tests
    if(method == "diffProp") {
      design <- S4Vectors::metadata(statsDf)$design
    } else if (method == "diffOcc") {
      design <- S4Vectors::metadata(statsDf)$pseudo.design
    }

    ## extract color code
    colorCode <- S4Vectors::metadata(barbieQ)$factorColors$testingBarcode

    ## customize row annotation
    barcodeAnnotation <- rowAnnotation(TendencyTo = statsDf$tendencyTo, annotation_name_side = "top",
        annotation_name_gp = grid::gpar(fontsize = 10), col = list(TendencyTo = colorCode),
        show_legend = TRUE, show_annotation_name = TRUE)
    
    ## extract relavent conditions
    if(all(names(contrastGroups) %in% c("levelLow", "levelHigh"))) {
      splitGroupHigh <- strsplit(contrastGroups["levelHigh"], " \\+ ")[[1]]
      splitGroupLow <- strsplit(contrastGroups["levelLow"], " \\+ ")[[1]]
    }
    ## extract sample groups
    if(all(splitGroupHigh %in% colnames(design)) ||
       all(splitGroupLow %in% colnames(design))) {
      ## when contrast is factor: two levels
      GroupHigh <- design[, splitGroupHigh, drop = FALSE] |> rowSums()
      GroupLow <- design[, splitGroupLow, drop = FALSE] |> rowSums()
      GroupVec <- rep("others", ncol(barbieQ))
      names(GroupVec) <- rownames(design)
      GroupVec[GroupHigh == 1] <- contrastGroups["levelHigh"]
      GroupVec[GroupLow == 1] <- contrastGroups["levelLow"]
      splitSamples <- TRUE
      sampleAnnotation <- HeatmapAnnotation(testingGroups = GroupVec, annotation_name_side = "left",
        annotation_name_gp = grid::gpar(fontsize = 10), col = list(testingGroups = c(colorCode, "others" = "grey")))
    } else {
      ## when contrast is numeric
      GroupVec <- NULL
      splitSamples <- FALSE
      sampleAnnotation <- NULL
    }

    hp <- plotBarcodeHeatmap(barbieQ = barbieQ, barcodeMetric = barcodeMetric, splitSamples = splitSamples,
        sampleMetadata = data.frame(testingGroups = GroupVec), sampleGroup = "testingGroups",
        barcodeAnnotation = barcodeAnnotation, sampleAnnotation = sampleAnnotation)

    return(hp)
}
