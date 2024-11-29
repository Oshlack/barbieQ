#' Plot Barcodes catogorized by significant change using a dot plot
#'
#' `plotBarcodeBiasScatterPlot()` use a dot plot to visualize the significance
#'  level of each Barcode in differential proportion or occurrence,
#'  as determined by the [testBarcodeBias] function.
#'  P.values are plotted against other properties of Barcodes
#'  specified by `xAxis`.
#'
#' @param barbieQ A `barbieQ` object created by the [createBarbieQ] function,
#'  updated with Barcode test results by calling the [testBarcodeBias] function.
#' @param elementName A string indicating name of the test
#'  conducted and stored in `barbieQ$testBarcode`. Default to the last
#'  test conducted and stored.
#' @param reorderRank A logical value deciding whether to reorder
#'  Barcode ranks within each sample. Defaults to FALSE.
#' @param pValuesAdjusted A logical value indicating which p.value to present.
#'  TRUE to present adjusted p.values. FALSE to present original p.values.
#'  Defaults to TRUE.
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
#'
#' @examples
#' Block <- c(1, 1, 2, 3, 3, 4, 1, 1, 2, 3, 3, 4)
#' Treat <- factor(rep(c("ctrl", "drug"), each = 6))
#' Time <- rep(rep(seq_len(2), each = 3), 2)
#' nbarcodes <- 50
#' nsamples <- 12
#' count <- abs(matrix(rnorm(nbarcodes * nsamples), nbarcodes, nsamples))
#' rownames(count) <- paste0("Barcode", seq_len(nbarcodes))
#' barbieQ <- createBarbieQ(count, data.frame(Treat = Treat, Time = Time))
#' testBB <- testBarcodeBias(barbieQ, sampleGroups = "Treat")
#' plotBarcodeBiasScatterPlot(barbieQ = testBB, elementName = "diffProp_Treat")
plotBarcodeBiasScatterPlot <- function(
    barbieQ, elementName = NULL, reorderRank = FALSE,
    pValuesAdjusted = TRUE, xAxis = "avgRank") {
  ## extract test resilts and information
  if (is.null(elementName)) {
    elementName <- names(barbieQ$testBarcodes)[length(names(barbieQ$testBarcodes))]
  }
  if (is.null(barbieQ$testBarcodes[[elementName]])) {
    stop("test results not specified or not found")
  }
  testInfo <- barbieQ$testBarcodes[[elementName]]
  statMat <- testInfo$results
  methodLs <- testInfo$methods
  ## define a custom color/shape palette for test results
  customShape <- stats::setNames(c(21, 24, 23), c(methodLs$contrastLevels, "n.s."))
  customColor <- barbieQ$factorColors[[elementName]]
  ## y axis will be p.value x axis will be optional: total occurrence, average
  ## rank, average log CPM reorder Barcode rank within samples.
  if (reorderRank) {
    rank <- apply(barbieQ$rank, 2, rank)
  } else {
    rank <- barbieQ$rank
  }
  ## choose p.values
  if (pValuesAdjusted) {
    p.value <- statMat$adj.p.value
  } else {
    p.value <- statMat$p.value
  }
  ## check xAxis
  xOptions <- c("avgRank", "totalOcc", "avgLogCPM", "avgProportion")
  xAxis <- match.arg(xAxis, xOptions)
  xTitle <- stats::setNames(
    c(
      "Average rank of Barcode across samples", "Number of samples in which Barcode occurs",
      "Average Barcode Log2 CPM+1 across samples", "Average Barcode proportion across samples"
    ),
    xOptions
  )
  ## data.frame for ggplot
  mydata <- data.frame(
    direction = statMat$direction, minusLogP = -log10(p.value), totalOcc = rowSums(barbieQ$occurrence),
    avgRank = rowMeans(barbieQ$rank), BarcodeID = rownames(barbieQ$assay), avgLogCPM = log2(rowMeans((barbieQ$CPM +
      1))), avgProportion = rowMeans(barbieQ$proportion)
  )
  ## visualize by ggplot
  p <- ggplot(mydata, aes(x = mydata[, xAxis], y = minusLogP, text = BarcodeID)) +
    geom_point(aes(
      color = direction,
      shape = direction, fill = direction
    ), size = 4, stroke = 1) +
    theme_classic() +
    theme(aspect.ratio = 1) +
    labs(
      title = paste0(methodLs$aim, " : ", methodLs$contrastVector),
      y = "-log10(p.value)", x = xTitle[xAxis]
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed", color = "grey"
    ) +
    scale_color_manual(values = customColor) +
    scale_shape_manual(values = customShape) +
    scale_fill_manual(values = alpha(
      customColor,
      0.2
    ))
  ## reverse x scale if displaying Barcode rank
  if (xAxis == "avgRank") {
    p <- p + annotate("text",
      x = min(mydata[, xAxis]) * 1.1, y = -log10(0.05), label = "p.value = 0.05",
      vjust = 1.5, hjust = 1
    ) + scale_x_reverse()
  } else {
    p <- p + annotate("text",
      x = max(mydata[, xAxis]) * 0.9, y = -log10(0.05), label = "p.value = 0.05",
      vjust = 1.5, hjust = 1
    )
  }

  return(p)
}


#' Plot Barcodes categorized by significant change using a Heatmap
#'
#' `plotBarcodeBiasHeatmap()` uses the Heatmap annotations to visualize the
#'  significance level of each Barcode in differential proportion or occurrence,
#'  as determined by the [testBarcodeBias] function.
#'
#' @param barbieQ A `barbieQ` object created by the [createBarbieQ] function,
#'  updated with Barcode test results by calling the [testBarcodeBias] function.
#' @param value A string indicating what to visualize.
#'  Defaults to 'CPM'. Options include: 'CPM' and 'occurrence'.
#' @param elementName A string indicating name of the test
#'  conducted and stored in `barbieQ$testBarcode`. Default to the last
#'  test conducted and stored.
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
#'
#' @examples
#' Block <- c(1, 1, 2, 3, 3, 4, 1, 1, 2, 3, 3, 4)
#' Treat <- factor(rep(c("ctrl", "drug"), each = 6))
#' Time <- rep(rep(seq_len(2), each = 3), 2)
#' nbarcodes <- 50
#' nsamples <- 12
#' count <- abs(matrix(rnorm(nbarcodes * nsamples), nbarcodes, nsamples))
#' rownames(count) <- paste0("Barcode", seq_len(nbarcodes))
#' barbieQ <- createBarbieQ(count, data.frame(Treat = Treat, Time = Time))
#' testBB <- testBarcodeBias(barbieQ, sampleGroups = "Treat")
#' plotBarcodeBiasHeatmap(barbieQ = testBB, elementName = "diffProp_Treat")
plotBarcodeBiasHeatmap <- function(barbieQ, value = "CPM", elementName = NULL, sampleAnnotation = NULL) {
  ## extract test resilts and information
  if (is.null(elementName)) {
    elementName <- names(barbieQ$testBarcodes)[length(names(barbieQ$testBarcodes))]
  }
  if (is.null(barbieQ$testBarcodes[[elementName]])) {
    stop("test results not specified or not found")
  }
  testInfo <- barbieQ$testBarcodes[[elementName]]
  statMat <- testInfo$results
  methodLs <- testInfo$methods
  modelTargets <- testInfo$targets
  ## define a custom color/shape palette for test results
  customShape <- stats::setNames(c(21, 24, 23), c(methodLs$contrastLevels, "n.s."))
  customColor <- barbieQ$factorColors[[elementName]]

  ## customize row annotation
  barcodeAnnotation <- rowAnnotation(
    Direction = statMat$direction, annotation_name_side = "top",
    annotation_name_gp = grid::gpar(fontsize = 10), col = list(Direction = customColor),
    show_legend = TRUE, show_annotation_name = TRUE
  )

  ## adjust the order of slices based on contrast levels in the test
  restLevels <- dplyr::setdiff(levels(modelTargets[, methodLs$contrastVector]), methodLs$contrastLevels)
  levels(modelTargets[, methodLs$contrastVector]) <- c(methodLs$contrastLevels, restLevels)

  hp <- plotBarcodeHeatmap(
    barbieQ = barbieQ, value = value, splitSamples = TRUE, targets = modelTargets,
    sampleGroups = methodLs$contrastVector, barcodeAnnotation = barcodeAnnotation,
    sampleAnnotation = sampleAnnotation
  )

  return(hp)
}
