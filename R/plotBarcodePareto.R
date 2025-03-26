#' Plot contributions of Barcodes split by \emph{top} v.s. \emph{bottom}
#'  in a circular bar plot.
#'
#' After the [tagTopBarcodes] function tags Barcodes as either
#'  \emph{top} or \emph{bottom}, `plotBarcodePareto()` visualizes the
#'  proportion of each Barcode and separates them by these tags
#'  in a circular bar plot, also known as a Pareto plot.
#'
#' @param barbieQ A `SummarizedExperiment` object created by the [createBarbieQ] function.
#' @param absoluteProportion A logical value indicating whether to present
#'  absolute Barcode mean proportion (across samples) or relative values across Barcodes, 
#'  Defaults to FALSE,
#'  which means it will present the percentage of (Barcode mean proportion) across Barcodes.
#'
#' @return A `ggplot` S3 class object displaying a circular bar plot,
#'  highlighting the relative total proportion of each Barcode across samples.
#'
#' @note To save the plot with its original aspect ratio, use the
#'  `ggplot2::ggsave` function and set `width = 8` and `height = 6`,
#'  or an equivalent size.
#'
#' @export
#'
#' @importFrom dplyr arrange
#' @importFrom dplyr group_by
#' @importFrom dplyr summarize
#' @importFrom dplyr rowwise
#' @importFrom dplyr mutate
#' @importFrom magrittr %>%
#' @importFrom graphics title
#' @importFrom stats end
#' @importFrom stats start
#' @import ggplot2
#' @import data.table
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment assays
#' @importFrom SummarizedExperiment rowData
#'
#' @examples
#' ## sample conditions and color palettes
#' sampleConditions <- data.frame(
#'   Treat = factor(rep(c("ctrl", "drug"), each = 6)),
#'   Time = rep(rep(seq_len(2), each = 3), 2)
#' )
#' conditionColor <- list(
#'   Treat = c(ctrl = "#999999", drug = "#112233"),
#'   Time = c("1" = "#778899", "2" = "#998877")
#' )
#' ## Barcode count data
#' nbarcodes <- 50
#' nsamples <- 12
#' barcodeCount <- abs(matrix(10, nbarcodes, nsamples))
#' barcodeCount[seq(21, 50), ] <- 0.0001
#' rownames(barcodeCount) <- paste0("Barcode", seq_len(nbarcodes))
#' ## create a `barbieQ` object
#' myBarbieQ <- createBarbieQ(barcodeCount, sampleConditions, conditionColor)
#' myBarbieQ <- tagTopBarcodes(myBarbieQ)
#' plotBarcodePareto(myBarbieQ)
plotBarcodePareto <- function(barbieQ, absoluteProportion = FALSE) {
  ## extract the top tag for each barcode.
  flag <- SummarizedExperiment::rowData(barbieQ)$isTopBarcode$isTop
  proportion <- SummarizedExperiment::assays(barbieQ)$proportion
  contribution <- rowMeans(proportion)
  relativeContribution <- contribution/sum(contribution) * 100
  ## create dataset
  ntop <- sum(flag)
  nbottom <- sum((!flag))
  topTag <- factor(flag, levels = c(TRUE, FALSE))
  
  ## choose y axis to plot relative or absolute contribution
  if(absoluteProportion) {
    contributionToPlot <- contribution
    yTitle <- "Barcode Mean Proportion"
  } else {
    contributionToPlot <- relativeContribution
    yTitle <- "Relative Barcode Mean Proportion (%)"
  }
  
  data <- data.frame(
    individual = rownames(barbieQ),
    group = topTag,
    value = contributionToPlot
  )
  ## set a number of 'empty bar' to add at the end of each group
  emptyBar <- 5
  toAdd <- data.frame(matrix(NA, emptyBar * nlevels(data$group), ncol(data)))
  colnames(toAdd) <- colnames(data)
  toAdd$group <- rep(levels(data$group), each = emptyBar)
  data <- rbind(data, toAdd)
  ## order data
  data <- data %>% dplyr::arrange(group, -value)
  data$id <- nrow(data) |> seq()
  ## add label
  labelData <- data
  numberBar <- nrow(labelData)
  angle <- 90 - 360 * (labelData$id - 0.5) / numberBar
  labelData$hjust <- ifelse(angle < -90, 1, 0)
  labelData$angle <- ifelse(angle < -90, angle + 180, angle)

  ## prepare a data frame for base lines
  baseData <- data %>%
    dplyr::group_by(group) %>%
    dplyr::summarize(start = min(id), end = max(id) - emptyBar) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(title = mean(c(start, end)))
  ## calculate the number of clones for each group
  baseData$labs <- ifelse(baseData$group == "TRUE", ntop, nbottom)
  ## prepare a data frame for grid (scales)
  gridData <- baseData
  gridData$end <- gridData$end[c(nrow(gridData), seq_len(nrow(gridData)) - 1)] + 1
  gridData$start <- gridData$start - 1
  gridData <- gridData[-1, ]
  ## total length of all bars
  pi <- max(data$id)

    p <- ggplot(data, aes(x = as.factor(id), y = value)) +
      ## add the bars with a blue color
      geom_bar(stat = "identity", alpha = 1, color = "#999966", linewidth = 0.5) +
      ## add a val=100/75/50/25 lines - compute it first to make sure barplots are OVER it.
      geom_segment(
        data = gridData,
        aes(x = -pi * 0.005, y = 20, xend = -pi * 0.001, yend = 20),
        colour = "grey", alpha = 1, linewidth = 0.3, inherit.aes = FALSE
      ) +
      geom_segment(
        data = gridData,
        aes(x = -pi * 0.005, y = 15, xend = -pi * 0.001, yend = 15),
        colour = "grey", alpha = 1, linewidth = 0.3, inherit.aes = FALSE
      ) +
      geom_segment(
        data = gridData,
        aes(x = -pi * 0.005, y = 10, xend = -pi * 0.001, yend = 10),
        colour = "grey", alpha = 1, linewidth = 0.3, inherit.aes = FALSE
      ) +
      geom_segment(
        data = gridData,
        aes(x = -pi * 0.005, y = 5, xend = -pi * 0.001, yend = 5),
        colour = "grey", alpha = 1, linewidth = 0.3, inherit.aes = FALSE
      ) +
      ## add text showing the value of each lines
      annotate("text",
        x = c(rep(pi * 0.995, 4), pi * 0.05),
        y = c(5, 10, 15, 20, 30),
        label = c("5%", "10%", "15%", "20%", yTitle),
        color = "#999966", size = 3, angle = 0, fontface = "bold", hjust = 1
      ) +
      ## the negative value controls the size of the central circle, the positive to add size over each bar
      ylim(-8, 30) +
      ## custom the theme: no axis title and no cartesian grid
      theme_minimal() +
      theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        ## remove unnecessary margin around plot
        plot.margin = unit(rep(0, 4), "cm"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.5)
      ) +
      ## make the coordinate polar instead of cartesian.
      coord_polar(start = 0) +
      ## add base line information
      geom_segment(
        data = baseData,
        aes(x = start, y = -1, xend = end, yend = -1, color = group),
        alpha = 0.8, linewidth = 2, inherit.aes = FALSE
      ) +
      geom_text(
        data = baseData,
        aes(x = title * 1.05, y = c(-4, -6), label = labs, color = group),
        hjust = c(1, 1), alpha = 0.8, size = 3, fontface = "bold",
        inherit.aes = FALSE, show.legend = FALSE
      ) +
      labs(color = "Num of Barcodes") +
      scale_color_manual(
        values = c("TRUE" = "#FF3399", "FALSE" = "#0066FF"),
        labels = c("Top Barcodes", "Others")
      )

  return(p)
}
