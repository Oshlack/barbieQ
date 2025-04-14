#' Plot Barcode pairwise correlation
#'
#' `plotBarcodePairCorrelation()` visualizes the correlation of each pair of
#'  Barcodes in the `barbieQ` object using a dot plot.
#'  Correlated Barcodes can be identified as showing high correlation in
#'  Barcode proportions across samples using the [clusterCorrelatingBarcodes]
#'  function.
#'  Visualizing the pair wise correlation assists in determining the threshold
#'  for tagging the highly correlated Barcode clusters, which are considered
#'  multiple Barcodes existing in the same clone.
#'
#' @param barbieQ A `barbieQ` object created by the [createBarbieQ] function.
#' @param method A string specifying the correlation method to use.
#'  Defaults to 'pearson'. Options include: 'pearson', 'kendall', 'spearman'.
#' @param yScaleMetric A string indicating what to present against
#'  correlation in the dot plot. Defaults to `mean`, representing the mean CPM
#'  for each pair. The alternative option is 'max'.
#' @param corThresh A numeric value that sets the threshold for high correlation
#'  Defaults to 0.95
#' @param cpmThresh A numeric value that sets the minimum level of
#'  Barcode pair's mean CPM for a Barcode pair to be considered as
#'  highly correlated co-existing Barcodes. Defaults to 2^10
#' @param preDefinedCluster preDefinedCluster A `list` of known groups containing different Barcodes,
#'  or a `vector`/`array` indicating Barcode groups;
#'  or an equivalent `matrix`, `data.frame`, or `DataFrame` with a single column.
#'  Defaults to NULL.
#'
#' @return A `ggplot` S3 object displaying a dot plot of the correlation in CPM
#'  between each pair of Barcodes, plotted against the mean or max of their CPM.
#'
#' @export
#'
#' @import ggplot2
#' @importFrom magrittr %>%
#' @import data.table
#' @importFrom dplyr mutate
#' @importFrom dplyr case_when
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment assays
#' @importFrom SummarizedExperiment rowData
#' @importFrom S4Vectors metadata
#'
#' @examples
#' nbarcodes <- 50
#' nsamples <- 12
#' count <- abs(matrix(rnorm(nbarcodes * nsamples), nbarcodes, nsamples))
#' rownames(count) <- paste0('Barcode', seq_len(nbarcodes))
#' barbieQ <- createBarbieQ(count)
#' plotBarcodePairCorrelation(barbieQ, preDefinedCluster = c(rep(seq_len(10), 5)))
#' plotBarcodePairCorrelation(barbieQ, preDefinedCluster = list(
#'  group1 = c('Barcode1', 'Barcode2', 'Barcode3'), group2 = c('Barcode4', 'Barcode5')))
plotBarcodePairCorrelation <- function(barbieQ, method = "pearson", yScaleMetric = "mean",
    corThresh = 0.95, cpmThresh = 2^10, preDefinedCluster = NULL) {
    ## check yScaleMetric
    yScaleMetric <- match.arg(yScaleMetric, c("mean", "max"))

    ## dispatch preprocessing function extract barcode pairwise correlation and
    ## pre-defined pairs
    barbieQ <- extractBarcodePairs(barbieQ, method = method, preDefinedCluster = preDefinedCluster)

    ## processed info is saved in this DFrame
    corDF <- SummarizedExperiment::rowData(barbieQ)$barcodeCorrelation
    preDefinedDf <- S4Vectors::metadata(corDF)$preDefinedBarcodePair
    meanMatCPM <- S4Vectors::metadata(corDF)$meanMatCPM
    corMat <- as.matrix(corDF)

    ## compute mean CPM for each barcode
    meanCPM <- rowMeans(SummarizedExperiment::assays(barbieQ)$CPM)
    ## N x N symmetric matrix taking the max between each pair on their row means
    maxMatCPM <- outer(meanCPM, meanCPM, FUN = function(x, y) pmax(x, y))

    ## get row and column indices of upper triangle of the square
    upperIdx <- which(upper.tri(corMat), arr.ind = TRUE)
    ## determine correlated pairs by passing correlation coefficient and CPM
    ## thresholds.
    highCorMat <- meanMatCPM >= cpmThresh & corMat >= corThresh
    ## convert the it into a long data
    identifiedPair <- highCorMat[upperIdx]
    ## create a vector indicating correlation groups
    correlationGroup <- vector(mode = "character", length = length(identifiedPair))
    correlationGroup[identifiedPair] <- "Identified"

    ## create a data.frame for identified pairs
    corResults <- data.frame(name1 = rownames(corMat)[upperIdx[, 1]], name2 = colnames(corMat)[upperIdx[,
        2]], mean = meanMatCPM[upperIdx], max = maxMatCPM[upperIdx], coefficient = corMat[upperIdx],
        correlationGroup)

    ## process preDefined pairs
    if (nrow(preDefinedDf) > 0L) {
        preDefinedPairs <- paste(preDefinedDf[, 1], preDefinedDf[, 2], sep = "_")
    } else {
        preDefinedPairs <- ""
    }
    ## tag preDefined pairs in the result data; preserve 'Identified'; 'tag' the rest
    corResults <- corResults %>%
        dplyr::mutate(correlationGroup = dplyr::case_when((paste(name1, name2, sep = "_") %in%
            preDefinedPairs) | (paste(name2, name1, sep = "_") %in% preDefinedPairs) ~ "pre-Defined",
            correlationGroup == "Identified" ~ "Identified", TRUE ~ "non-Corr"))

    ## choose what to present on y axis
    yAxis <- log2(corResults[, yScaleMetric] + 1)
    yTitle <- paste0("log2 (", yScaleMetric, " CPM+1) between each Barcode pair")

    ## plotting correlations
      p <- ggplot(corResults, aes(x = coefficient)) + geom_histogram(aes(y = (after_stat(count))/max(after_stat(count)) *
        max(yAxis)), binwidth = 0.05, alpha = 0.3, fill = "grey") + geom_point(aes(y = yAxis,
        color = correlationGroup)) + stat_ecdf(geom = "step", aes(y = ..y.. * max(yAxis),
        color = correlationGroup), alpha = 0.5) + scale_color_manual(values = c(`pre-Defined` = "#00BFC4",
        Identified = "#F8766D", `non-Corr` = "grey")) + theme_classic() + theme(aspect.ratio = 1) +
        labs(x = paste0(method, " correlation coefficient"), y = yTitle) + scale_x_continuous(expand = c(0.05,
        0), limits = c(-1, 1)) + scale_y_continuous(sec.axis = sec_axis(~./max(yAxis), name = "cummulative freq.")) +
        geom_hline(yintercept = log2(cpmThresh + 1), linetype = "dashed", color = "#7CAE00",
            alpha = 0.8) + geom_vline(xintercept = corThresh, linetype = "dashed", color = "#C77CFF",
        alpha = 0.8) + annotate("text", x = corThresh - 0.1, y = max(yAxis) * 1.03, color = "#C77CFF",
        alpha = 1, size = 3, label = paste0("Cor=", corThresh)) + annotate("text", x = corThresh -
        0.4, y = log2(cpmThresh + 1) - max(yAxis) * 0.02, color = "#7CAE00", alpha = 1,
        size = 3, label = paste0("log2(", yScaleMetric, "+1)=", cpmThresh))


    return(p)
}


#' Cluster correlated Barcodes based on pairwise correlation
#'
#' `clusterCorrelatingBarcodes()` groups Barcode into clusters by
#'  identifying correlated Barcode pairs based on two criteria:
#'   * The proportions of Barcodes across samples exhibit high correlation,
#'    exceeding the threshold specified by `corThresh`.
#'   * The mean CPM of the Barcode pair exceeds the threshold specified
#'    by `cpmThresh`.
#'  These two parameters can be optimized up to the users using the visualization
#'  function [plotBarcodePairCorrelation].
#'
#' @param barbieQ A `barbieQ` object created by the [createBarbieQ] function.
#' @param method A string specifying the correlation method to use.
#'  Defaults to 'pearson'. Options include: 'pearson', 'kendall', 'spearman'.
#' @param corThresh A numeric value that sets the threshold for high correlation
#'  Defaults to 0.95
#' @param cpmThresh A numeric value that sets the minimum level of
#'  Barcode pair's mean CPM for a Barcode pair to be considered as
#'  highly correlated co-existing Barcodes. Defaults to 2^10.
#' @param preDefinedCluster A `list` of known groups containing different Barcodes,
#'  or a `vector`/`array` indicating Barcode groups;
#'  or an equivalent `matrix`, `data.frame`, or `DataFrame` with a single column.
#'  Defaults to NULL.
#'
#' @return A `barbieQ` object updated with `barcodeCorrelatedCluster`
#'  as a `DataFrame` saved to the `rowData` of `barbieQ`,
#'  containing a single column of identified Barcode cluster groups,
#'  where the `metadata` saves a plot of the identified Barcode clusters.
#'
#' @export
#'
#' @importFrom dplyr filter
#' @importFrom igraph graph_from_edgelist
#' @importFrom igraph clusters
#' @importFrom magrittr %>%
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment assays
#' @importFrom SummarizedExperiment rowData
#' @importFrom S4Vectors DataFrame
#' @importFrom S4Vectors metadata
#'
#' @examples
#' nbarcodes <- 50
#' nsamples <- 12
#' count <- abs(matrix(rnorm(nbarcodes * nsamples), nbarcodes, nsamples))
#' rownames(count) <- paste0('Barcode', seq_len(nbarcodes))
#' barbieQ <- createBarbieQ(count)
#' clusterCorrelatingBarcodes(barbieQ, preDefinedCluster = c(rep(seq_len(10), 5)))
#' clusterCorrelatingBarcodes(barbieQ, preDefinedCluster = list(
#'  group1 = c('Barcode1', 'Barcode2', 'Barcode3'), group2 = c('Barcode4', 'Barcode5')))
clusterCorrelatingBarcodes <- function(barbieQ, method = "pearson", corThresh = 0.95, cpmThresh = 2^10,
    preDefinedCluster = NULL) {
    ## dispatch preprocessing function extract barcode pairwise correlation and
    ## pre-defined pairs
    barbieQ <- extractBarcodePairs(barbieQ, method = method, preDefinedCluster = preDefinedCluster)

    ## processed info is saved in this DFrame
    corDF <- SummarizedExperiment::rowData(barbieQ)$barcodeCorrelation
    preDefinedDf <- S4Vectors::metadata(corDF)$preDefinedBarcodePair
    meanMatCPM <- S4Vectors::metadata(corDF)$meanMatCPM
    corMat <- as.matrix(corDF)

    ## determine correlated pairs by passing correlation coefficient and CPM
    ## thresholds.
    highCorMat <- meanMatCPM >= cpmThresh & corMat >= corThresh
    ## arrow to extract the upper triangle, while also defined as 'correlated pairs'
    identifiedMat <- upper.tri(corMat) & highCorMat
    ## get row and column indices of identified pairs
    upperIdx <- which(identifiedMat, arr.ind = TRUE)

    ## create a data.frame for identified pairs
    highCorDf <- data.frame(name1 = rownames(corMat)[upperIdx[, 1]], name2 = colnames(corMat)[upperIdx[,
        2]], meanCPM = meanMatCPM[identifiedMat], coefficient = corMat[identifiedMat])

    ## combine pre-defined pairs and identified pairs
    totalPairs <- rbind(preDefinedDf, highCorDf[, c("name1", "name2")])

    ## plotting highly correlating Barcode pairs
    graphObject <- igraph::graph_from_edgelist(totalPairs %>%
        as.matrix(), directed = FALSE)
    clustersPlot <- plot(graphObject, vertex.label.cex = 0.4)
    ## cluster Barcode groups based on correlating pairs
    groups <- igraph::components(graphObject)$membership
    identifiedList <- list()
    if (length(groups) > 0L)
        identifiedList <- split(names(groups), groups)
    ## create an array indicating which group Barcodes belong to
    BarcodeGroupArray <- setNames(numeric(nrow(barbieQ)), rownames(barbieQ))
    BarcodeGroupArray[names(groups)] <- groups
    ## assign single Barcodes by unique group names
    tagNoGroup <- which(BarcodeGroupArray == 0)
    BarcodeGroupArray[tagNoGroup] <- -seq_along(tagNoGroup)

    ## save identified Barcode clusters to barbieQ object as a rowData
    barcodeCorrelatedCluster <- S4Vectors::DataFrame(cluster = BarcodeGroupArray)
    SummarizedExperiment::rowData(barbieQ)$barcodeCorrelatedCluster <- barcodeCorrelatedCluster
    ## save the graph to the metadata of `barcodeCorrelatedCluster`
    S4Vectors::metadata(barcodeCorrelatedCluster)$clustersPlot <- clustersPlot

    ## message discovered clusters
    message("identified ", length(identifiedList), " clusters, including ", length(groups),
        " Barcodes.")

    return(barbieQ)
}

#' Generate pairwise Barcode correlation and standardize
#'  known Barcode cluster information
#'
#' @param barbieQ A `barbieQ` object created by the [createBarbieQ] function.
#' @param method A string specifying the correlation method to use.
#'  Defaults to 'pearson'. Options include: 'pearson', 'kendall', 'spearman'.
#' @param preDefinedCluster A `list` of known groups containing different Barcodes,
#'  e.g. `list(group1 = c('Barcode1', 'Barcode2', 'Barcode3'), group2 = c('Barcode4', 'Barcode5'))`
#'  or a `vector`/`array` indicating Barcode groups;
#'  or an equivalent `matrix`, `data.frame`, or `DataFrame` with a single column.
#'  Defaults to NULL.
#'
#' @return `barbieQ` object, updated with a `rowData` called `barcodeCorrelation`: 
#'  a symmetric `DataFrame` with all Barcodes in both rows and columns,
#'  representing pairwise Barcode correlations; where `metadata` 
#'  called `preDefinedBarcodePair` is a `data.frame` of two columns,
#'  where each row saves the names of a pair of pre-defined correlated Barcodes.
#'
#' @importFrom utils combn
#' @importFrom dplyr setequal
#' @importFrom magrittr %>%
#' @importFrom stats cor
#' @importFrom stats p.adjust
#' @import data.table
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment assays
#' @importFrom SummarizedExperiment rowData
#' @importFrom S4Vectors DataFrame
#' @importFrom S4Vectors metadata
#'
#' @noRd
#'
#' @examples \donttest{
#' nbarcodes <- 50
#' nsamples <- 12
#' count <- abs(matrix(rnorm(nbarcodes * nsamples), nbarcodes, nsamples))
#' rownames(count) <- paste0('Barcode', seq_len(nbarcodes))
#' barbieQ <- createBarbieQ(count)
#' barbieQ:::extractBarcodePairs(
#'   barbieQ,
#'   preDefinedCluster = c(rep(seq_len(10), 5))
#' )
#' }
extractBarcodePairs <- function(barbieQ, method = "pearson", preDefinedCluster = NULL) {
    ## check method
    method <- match.arg(method, c("pearson", "kendall", "spearman"))
    ## extract data
    mat <- SummarizedExperiment::assays(barbieQ)$CPM
    ## confirm Barcode IDs
    if (is.null(rownames(mat)))
        rownames(mat) <- rownames(barbieQ)

    ## case when `preDefinedCluster` unspecified, assign each barcode to a unique
    ## cluster number
    if (is.null(preDefinedCluster)) {
        preDefinedCluster <- S4Vectors::DataFrame(cluster = seq_along(rownames(barbieQ)))
        rownames(preDefinedCluster) <- rownames(barbieQ)
    }
    ## check `preDefinedCluster` structure
    if (inherits(preDefinedCluster, "matrix") || inherits(preDefinedCluster, "data.frame") ||
        is(preDefinedCluster, "DataFrame")) {
        ## case when `preDefinedCluster` is a d.f., DF, or matrix of single column
        if (ncol(preDefinedCluster) == 1 && nrow(preDefinedCluster) == nrow(barbieQ)) {
            ## take it when dimension is right
            preDefinedCluster <- S4Vectors::DataFrame(preDefinedCluster)
            colnames(preDefinedCluster) <- "cluster"
        } else {
            stop("`preDefinedCluster` must be a single column with nrow matching `barbiQ`.")
        }
    } else if ((is.vector(preDefinedCluster) || is.factor(preDefinedCluster)) && !is.list(preDefinedCluster)) {
        ## case when `preDefinedCluster` is a vector or array
        if (length(preDefinedCluster) == nrow(barbieQ)) {
            preDefinedCluster <- S4Vectors::DataFrame(cluster = preDefinedCluster)
        } else {
            stop("`preDefinedCluster` length must match number of Barcodes in `barbiQ`.")
        }
    } else if (is.list(preDefinedCluster)) {
        ## case when it's a list
        preDefinedClusterList <- preDefinedCluster
    } else {
        stop("wrong format for `preDefinedCluster`!")
    }

    ## now `preDefinedCluster` should be either a DFrame or a list
    if (!is.list(preDefinedCluster)) {
        ## for DFrames, check rownames
        if (is.null(rownames(preDefinedCluster))) {
            ## if barcode names unspecified
            rownames(preDefinedCluster) <- rownames(barbieQ)
            message("Barcode names of `preDefinedCluster` set up in line with `barbieQ`.")
        } else if (setequal(rownames(preDefinedCluster), rownames(barbieQ))) {
            ## if barcode names overlap but in different order
            preDefinedCluster <- preDefinedCluster[match(rownames(barbieQ), rownames(preDefinedCluster)),
                , drop = FALSE]
        } else {
            rownames(preDefinedCluster) <- rownames(barbieQ)
            warning("Barcode names don't match! Assigning Barcode names from `barbieQ` to `preDefinedCluster.`")
        }
        ## convert the DFrame into a list of barcode groups
        preDefinedClusterList <- split(rownames(preDefinedCluster), preDefinedCluster$cluster)
    }

    ## convert cluster list into data.frame of each pair
    preDefinedClusterList <- lapply(preDefinedClusterList, function(xList) {
        if (length(xList) >= 2)
            utils::combn(xList, 2) %>%
                t()
    })
    preDefinedDf <- do.call(rbind, preDefinedClusterList)
    # if `preDefinedDf` is null, make it an empty dataframe of two columns
    if (is.null(preDefinedDf))
        preDefinedDf <- data.frame()

    ## transpose barcode into columns for calculating barcode pairwise correlation
    corDF <- stats::cor(t(mat), method = method) |>
        S4Vectors::DataFrame()
    ## save the results of a nrow x nrow matrix into DFrame save the correlation
    ## method in the metadata of DFrame
    S4Vectors::metadata(corDF)$method <- method
    message("processing Barcode pairwise ", method, " correlation.")

    ## compute mean CPM for each barcode
    meanCPM <- rowMeans(mat)
    ## N x N symmetric matrix with average between each pair on their row means
    meanMat <- outer(meanCPM, meanCPM, FUN = function(x, y) (x + y)/2)
    ## N x N symmetric matrix with LFC between each pair on their CPM
    # logMat <- log2(mat + 0.5)
    # LFCMat <- outer(logMat, logMat, FUN = function(x, y) (x - y))
    ## save meanMat to DFrame as a metadata
    S4Vectors::metadata(corDF)$meanMatCPM <- meanMat
    # S4Vectors::metadata(corDF)$LFCMatCPM <- LFCMat

    ## save preDefinedDf to corDF metadata
    S4Vectors::metadata(corDF)$preDefinedBarcodePair <- preDefinedDf
    ## save corDF to rowData of barbieQ
    SummarizedExperiment::rowData(barbieQ)$barcodeCorrelation <- corDF

    return(barbieQ)
}
