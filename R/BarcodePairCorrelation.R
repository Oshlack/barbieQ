#' Plot Barcode pairwise correlation
#'
#' `plotBarcodePairCorrelation()` visualizes the correlation of each pair of
#'  Barcodes in the `barbieQ` object using a dot plot.
#'  "Co-existing" Barcodes are identified as showing high correlation in
#'  Barcode proportions across samples using the [clusterCorrelatingBarcodes]
#'  function.
#'  Visualizing the pair wise correlation assists in determining the threshold
#'  for tagging the highly correlated co-exisitng Barcode clusters.
#'
#' @param barbieQ A `barbieQ` object created by the [createBarbieQ] function.
#' @param method A string specifying the correlation method to use.
#'  Defaults to "pearson". Options include: "pearson", "kendall", "spearman".
#' @param dataVisual A string indicating what to present against
#'  correlation in the dot plot. Defaults to `mean`, representing the mean CPM
#'  for each pair. The alternative option is "max".
#' @param corCutoff A numeric value that sets the threshold for high correlation
#'  Defaults to 0.95
#' @param dataCutoff A numeric value that sets the minimum level of
#'  Barcode pair's log2(mean CPM) for a Barcode pair to be considered as
#'  highly correlated co-existing Barcodes. Defaults to 0
#' @param BarcodeClusters A `list` of known groups containing different Barcodes
#'  or a `vector`/`array` indicating Barcode groups. Defaults to NULL.
#'
#' @return A `ggplot` S3 object displaying a dot plot of the correlation in CPM
#'  between each pair of Barcodes, plotted against the mean or max of their CPM.
#'
#' @export
#'
#' @import ggplot2
#' @importFrom magrittr  %>%
#' @import data.table
#'
#' @examples
#' nbarcodes <- 50
#' nsamples <- 12
#' count <- abs(matrix(rnorm(nbarcodes * nsamples), nbarcodes, nsamples))
#' rownames(count) <- paste0("Barcode", seq_len(nbarcodes))
#' barbieQ <- createBarbieQ(count)
#' plotBarcodePairCorrelation(barbieQ, BarcodeClusters = c(rep(seq_len(10), 5)))
plotBarcodePairCorrelation <- function(
    barbieQ, method = "pearson", dataVisual = "mean",
    corCutoff = 0.95, dataCutoff = 0, BarcodeClusters = NULL) {
  ## check dataVisual
  dataVisual <- match.arg(dataVisual, c("mean", "max"))

  ## call preprocessing function to extract pair wise information
  processedInfo <- extractBarcodePairs(
    barbieQ,
    method = method, BarcodeClusters = BarcodeClusters
  )
  corTestResults <- processedInfo$corTestResults
  knownPairDf <- processedInfo$knownPairDf

  ## choose what to present on y axis
  yAxis <- log2(corTestResults[, dataVisual] + 1)
  yTitle <- paste0("log2 (", dataVisual, " CPM +1)")

  ## ggplot2 has a bug that gives wrong warning of 'Removed rows containing missing values'
  suppressWarnings({
    ## plotting correlations
    p <- ggplot(corTestResults, aes(x = cor)) +
      geom_histogram(
        aes(y = (after_stat(count)) / max(after_stat(count)) * max(yAxis)),
        binwidth = 0.05, alpha = 0.3, fill = "grey"
      ) +
      geom_point(aes(y = yAxis, color = signif, shape = knownCorrelating)) +
      stat_ecdf(geom = "step", aes(y = ..y.. * max(yAxis), color = signif), alpha = 0.5) +
      scale_color_manual(values = c("n.s." = "#00BFC4", "*" = "#F8766D")) +
      scale_shape_manual(values = c("TRUE" = 2, "FALSE" = 1)) +
      theme_classic() +
      theme(aspect.ratio = 1) +
      labs(
        x = paste0(method, " correlation of pairwise Barcode CPM"),
        y = yTitle
      ) +
      scale_x_continuous(expand = c(0.05, 0), limits = c(-1, 1)) +
      scale_y_continuous(
        sec.axis = sec_axis(~ . / max(yAxis), name = "cummulative freq.")
      ) +
      geom_hline(yintercept = dataCutoff, linetype = "dashed", color = "#7CAE00", alpha = 0.8) +
      geom_vline(xintercept = corCutoff, linetype = "dashed", color = "#C77CFF", alpha = 0.8) +
      annotate("text",
        x = corCutoff - 0.1, y = max(yAxis) * 1.03,
        color = "#C77CFF", alpha = 1, size = 3,
        label = paste0("Cor=", corCutoff)
      ) +
      annotate("text",
        x = corCutoff - 0.4, y = dataCutoff - max(yAxis) * 0.02,
        color = "#7CAE00", alpha = 1, size = 3,
        label = paste0("log2(", dataVisual, "+1)=", dataCutoff)
      )
  })

  return(p)
}


#' Cluster "co-existing" Barcodes based on pairwise correlation
#'
#' `clusterCorrelatingBarcodes()` groups "co-existing" Barcode into clusters by
#'  identifying "co-existing" Barcode pairs based on two criteria:
#'   * The proportions of Barcodes across samples exhibit high correlation,
#'    exceeding the threshold specified by `corCutoff`.
#'   * The log2(mean CPM) of the Barcode pair exceeds the threshold specified
#'    by `dataCutoff`.
#'  These two parameters can be manually optimized using the visualization
#'  function [plotBarcodePairCorrelation].
#'
#' @param barbieQ A `barbieQ` object created by the [createBarbieQ] function.
#' @param method A string specifying the correlation method to use.
#'  Defaults to "pearson". Options include: "pearson", "kendall", "spearman".
#' @param corCutoff A numeric value that sets the threshold for high correlation
#'  Defaults to 0.95
#' @param dataCutoff A numeric value that sets the minimum level of
#'  Barcode pair's log2(mean CPM) for a Barcode pair to be considered as
#'  highly correlated co-existing Barcodes. Defaults to 0
#' @param BarcodeClusters A `list` of known groups containing different Barcodes
#'  or a `vector`/`array` indicating Barcode groups. Defaults to NULL.
#' @param plotClusters A logical value. TRUE returns a plot showing
#'  the predicted Barcode clusters. FALSE returns an updated `barbieQ` object.
#'  Defaults to FALSE.
#'
#' @return A `barbieQ` object updated with adding a `data.frame` component of
#'  Barcode cluster, called `BarcodeCluster`.
#'
#' @export
#'
#' @importFrom dplyr filter
#' @importFrom igraph graph_from_edgelist
#' @importFrom igraph clusters
#' @importFrom magrittr %>%
#'
#' @examples
#' nbarcodes <- 50
#' nsamples <- 12
#' count <- abs(matrix(rnorm(nbarcodes * nsamples), nbarcodes, nsamples))
#' rownames(count) <- paste0("Barcode", seq_len(nbarcodes))
#' barbieQ <- createBarbieQ(count)
#' clusterCorrelatingBarcodes(barbieQ, BarcodeClusters = c(rep(seq_len(10), 5)))
clusterCorrelatingBarcodes <- function(
    barbieQ, method = "pearson", corCutoff = 0.95, dataCutoff = 0,
    BarcodeClusters = NULL, plotClusters = FALSE) {
  ## call preprocessing function to extract pair wise information
  processedInfo <- extractBarcodePairs(
    barbieQ,
    method = method, BarcodeClusters = BarcodeClusters
  )
  corTestResults <- processedInfo$corTestResults
  knownPairDf <- processedInfo$knownPairDf

  ## determine cluster based on high pair wise correlation on CPM
  corTestResultsHigh <- corTestResults %>%
    dplyr::filter(signif == "*") %>%
    dplyr::filter(cor >= corCutoff) %>%
    dplyr::filter(mean >= dataCutoff)

  ## combine known correlating pairs and predicted pairs
  totalPairs <- rbind(
    knownPairDf,
    corTestResultsHigh[, c("BarcodeX1", "BarcodeX2")]
  )
  ## plotting highly correlating Barcode pairs
  g <- igraph::graph_from_edgelist(
    totalPairs %>% as.matrix(),
    directed = FALSE
  )
  p <- plot(g, vertex.label.cex = 0.4)
  ## cluster Barcode groups based on correlating pairs
  groups <- igraph::components(g)$membership
  predictList <- list()
  if (length(groups) > 0L) predictList <- split(names(groups), groups)
  ## create an array indicating which group Barcodes belong to
  BarcodeGroupArray <- setNames(
    numeric(nrow(barbieQ$assay)), rownames(barbieQ$assay)
  )
  BarcodeGroupArray[names(groups)] <- groups
  ## assign single Barcodes by unique group names
  tagNoGroup <- which(BarcodeGroupArray == 0)
  BarcodeGroupArray[tagNoGroup] <- -seq_along(tagNoGroup)
  ## save Barcode groups in barbieQ object
  barbieQ$BarcodeCluster <- data.frame(corCluster = BarcodeGroupArray)
  ## message discovered clusters
  message(
    "predicting ", length(predictList), " clusters, including ", length(groups), " Barcodes."
  )

  if (plotClusters) {
    returnWhat <- p
  } else {
    returnWhat <- barbieQ
  }

  return(returnWhat)
}

#' Generate pairwise Barcode correlation and standardize
#'  known Barcode cluster information
#'
#' @param barbieQ A `barbieQ` object created by the [createBarbieQ] function.
#' @param method A string specifying the correlation method to use.
#'  Defaults to "pearson". Options include: "pearson", "kendall", "spearman".
#' @param BarcodeClusters A `list` of known groups containing different Barcodes
#'  or a `vector`/`array` indicating Barcode groups. Defaults to NULL.
#'
#' @return A `list` containing two `data.frame`s:
#'  * The correlation of all pairs of Barcodes
#'  * The names of known correlated Barcodes
#'
#' @importFrom utils combn
#' @importFrom dplyr setequal
#' @importFrom magrittr %>%
#' @importFrom stats cor.test
#' @importFrom stats p.adjust
#' @import data.table
#'
#' @noRd
#'
#' @examples \donttest{
#' nbarcodes <- 50
#' nsamples <- 12
#' count <- abs(matrix(rnorm(nbarcodes * nsamples), nbarcodes, nsamples))
#' rownames(count) <- paste0("Barcode", seq_len(nbarcodes))
#' barbieQ <- createBarbieQ(count)
#' barbieQ:::extractBarcodePairs(
#'   barbieQ,
#'   BarcodeClusters = c(rep(seq_len(10), 5))
#' )
#' }
extractBarcodePairs <- function(
    barbieQ, method = "pearson", BarcodeClusters = NULL) {
  ## check barbieQ
  checkBarbieQDimensions(barbieQ)
  ## check method
  method <- match.arg(method, c("pearson", "kendall", "spearman"))
  ## extract data
  mat <- barbieQ$CPM
  ## confirm Barcodde IDs are provided
  if (is.null(rownames(mat))) rownames(mat) <- rownames(barbieQ$assay)
  if (is.null(rownames(mat))) {
    rownames(mat) <- paste0(
      "Barcode", seq_len(nrow(mat))
    )
  }
  rownames(mat) <- make.unique(rownames(mat))
  ## check BarcodeClusters
  if (is.null(BarcodeClusters)) {
    BarcodeClusters <- seq_along(rownames(barbieQ$assay))
  }
  ## if BarcodeClusters is a vector or array
  if ((is.vector(BarcodeClusters) || is.factor(BarcodeClusters)) &&
    !is.list(BarcodeClusters)) {
    if (length(BarcodeClusters) != nrow(barbieQ$assay)) {
      stop("'BarcodeCluster' must be of the same dimension as Barcodes.")
    }
    if (is.null(names(BarcodeClusters))) {
      names(BarcodeClusters) <- rownames(barbieQ$assay)
    }
    ## convert the BarcodeCluster array into a list of groups
    BarcodeClustersList <- split(names(BarcodeClusters), BarcodeClusters)
  } else if (is.list(BarcodeClusters)) {
    BarcodeClustersList <- BarcodeClusters
  } else {
    stop("'BarcodeCluster' must be a vector/facter or a list.")
  }
  ## convert cluster list into data.frame of each pair
  knownPairList <- lapply(BarcodeClustersList, function(xList) {
    if (length(xList) >= 2) utils::combn(xList, 2) %>% t()
  })
  knownPairDf <- do.call(rbind, knownPairList)
  # if knownPairDf is null, make it an empty dataframe of two columns
  if (is.null(knownPairDf)) knownPairDf <- data.frame()

  ## transpose barcode into columns for the convenience of calculating correlation matrix
  mat <- t(mat)
  ## compute the correlation test on each pair of Barcodes
  ## omitting NA values by setting use="complete.obs"
  corTestList <- utils::combn(ncol(mat), 2, function(idx) {
    corTest <- stats::cor.test(
      mat[, idx[1]], mat[, idx[2]],
      alternative = "greater", use = "complete.obs", method = method
    )
    data.frame(
      pair = paste0(colnames(mat)[idx], collapse = "."),
      BarcodeX1 = colnames(mat)[idx[1]],
      BarcodeX2 = colnames(mat)[idx[2]],
      cor = corTest$estimate,
      p.value = corTest$p.value,
      mean = mean(mat[, idx], na.rm = TRUE),
      max = max(mat[, idx], na.rm = TRUE),
      ## assess if the current pair exist in the known pair list
      knownCorrelating = apply(
        knownPairDf, 1, function(row) {
          dplyr::setequal(colnames(mat)[idx], row)
        }
      ) %>% any()
    )
  }, simplify = FALSE)
  ## binding the result of each pair into a long data.frame
  corTestResults <- do.call(rbind, corTestList)
  ## adding an adjusted pvalue column
  corTestResults <- corTestResults %>%
    mutate(adj.p.value = stats::p.adjust(p.value, method = "BH")) %>%
    mutate(signif = ifelse(adj.p.value < 0.05, "*", "n.s."))

  return(list(
    corTestResults = corTestResults,
    knownPairDf = knownPairDf
  ))
}
