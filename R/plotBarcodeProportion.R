#' Plot Barcode contributions as mean Barcode proportion across samples
#'
#' `plotBarcodeProportion()` visualizes the mean proportion of each Barcode
#'  across samples using a bar plot, allowing for an easy comparison of
#'  Barcode contributions.
#'
#' @param barbieQ A `SummarizedExperiment` object created by the [createBarbieQ] function.
#' @param absoluteProportion A logical value indicating whether to present
#'  absolute Barcode mean proportion (across samples) or relative values across Barcodes, 
#'  Defaults to FALSE,
#'  which means it will present the percentage of (Barcode mean proportion) across Barcodes.
#' @param coordFixed A logical value indicating whether to coordinate the x
#'  and y scales. Defaults to FALSE, meaning the scales are not coordinated.
#' @param colorGradient A logical value indicating whether to apply a gradient
#'  to the bar colors. Defaults to FALSE, which means colors will be
#'  consistent across all bars.
#'
#' @return A `ggplot` S3 class object displaying Barcode contributions in a
#'  bar plot.
#'
#' @export
#'
#' @import ggplot2
#' @importFrom circlize colorRamp2
#' @import data.table
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment assays
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
#' plotBarcodeProportion(myBarbieQ)
plotBarcodeProportion <- function(barbieQ, absoluteProportion = FALSE, coordFixed = FALSE,
    colorGradient = FALSE) {
    ## compute mean Barcode proportion across samples as contribution
    proportion <- SummarizedExperiment::assays(barbieQ)$proportion
    contribution <- rowMeans(proportion)
    relativeContribution <- contribution/sum(contribution) * 100
    ## rank the relative contribution of each Barcode
    rank <- rank(-relativeContribution, ties.method = "first")

    ## choose y axis to plot relative or absolute contribution
    if (absoluteProportion) {
        contributionToPlot <- contribution
        yTitle <- "Barcode Mean Proportion"
    } else {
        contributionToPlot <- relativeContribution
        yTitle <- "Relative Barcode Mean Proportion (%)"
    }

    ## extract plotting data
    data <- data.frame(individual = rownames(barbieQ), contributionToPlot = contributionToPlot,
        rank = as.numeric(rank))
    ## define color function and assign bar colors
    colorFun <- circlize::colorRamp2(c(min(relativeContribution), max(relativeContribution)),
        c("#FFC1E0", "#FF3399"))
    if (colorGradient) {
        barColor <- colorFun(relativeContribution)
    } else {
        barColor <- "#FF3399"
    }

    p <- ggplot(data, aes(x = rank, y = contributionToPlot)) + geom_bar(stat = "identity",
        alpha = 1, color = alpha(barColor, 0.2), fill = barColor, linewidth = 1) + labs(title = "Barcode Contribution",
        x = paste0(nrow(data), " Barcodes"), y = yTitle) + theme_minimal() + theme(axis.ticks = element_blank(),
        panel.background = element_blank(), panel.grid = element_blank(), axis.text.x = element_blank())

    if (coordFixed) {
        p <- p + coord_fixed()
    } else {
        p <- p + theme(aspect.ratio = 1)
    }

    return(p)
}
