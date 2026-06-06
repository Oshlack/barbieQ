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
#' @param transformation A string specifying the transformation method for proportion 
#'  to be used for computing barcode pair-wise correlation, and mean proportion.
#'  Options include: 'asin-sqrt', 'logit', and 'none'.
#'  Defaults to 'asin-sqrt'.
#' @param showRawProportion A logical value, deciding whether to label y axis 
#'  by raw proportion when correlation is calculated based on transformed data specified by
#'  `transformation`. Default to FALSE.
#' @param yScaleMetric A string indicating what to present against
#'  correlation in the dot plot. Defaults to `mean`, representing the mean proportion
#'  for each pair. The alternative option is 'max'.
#' @param corThresh A numeric value that sets the threshold for high correlation
#'  Defaults to 0.95
#' @param propThresh A numeric value that sets the minimum level of
#'  Barcode pair's mean proportion (raw) for a Barcode pair to be considered as
#'  highly correlated co-existing Barcodes. Defaults to 0.001.
#' @param preDefinedCluster preDefinedCluster A `list` of known groups containing different Barcodes,
#'  or a `vector`/`array` indicating Barcode groups;
#'  or an equivalent `matrix`, `data.frame`, or `DataFrame` with a single column.
#'  Defaults to NULL.
#'
#' @return A `ggplot` S3 object displaying a dot plot of the correlation in proportion
#'  between each pair of Barcodes, plotted against the mean or max of their proportion.
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
plotBarcodePairCorrelation <- function(barbieQ, method = "pearson", transformation = "asin-sqrt", showRawProportion = FALSE,
    yScaleMetric = "mean", corThresh = 0.95, propThresh = 0.001, preDefinedCluster = NULL) {
    ## check yScaleMetric
    yScaleMetric <- match.arg(yScaleMetric, c("mean", "max"))

    ## dispatch preprocessing function extract barcode pairwise correlation and
    ## pre-defined pairs
    barbieQ <- extractBarcodePairs(barbieQ, method = method, transformation = transformation, preDefinedCluster = preDefinedCluster)

    ## processed info is saved in this DFrame
    corDF <- SummarizedExperiment::rowData(barbieQ)$barcodeCorrelation
    preDefinedDf <- S4Vectors::metadata(corDF)$preDefinedBarcodePair
    meanMatProp <- S4Vectors::metadata(corDF)$meanMatProp
    corMat <- as.matrix(corDF)
    S4Vectors::metadata(corDF)$transformation

    ## compute mean proportion for each barcode
    meanProp <- rowMeans(SummarizedExperiment::assays(barbieQ)$proportion)
    ## N x N symmetric matrix taking the max between each pair on their row means
    maxMatProp <- outer(meanProp, meanProp, FUN = function(x, y) pmax(x, y))

    ## get row and column indices of upper triangle of the square
    upperIdx <- which(upper.tri(corMat), arr.ind = TRUE)

    ## determine correlated pairs by passing correlation coefficient and proportion
    ## thresholds.
    ## scale the propThresh 
    scaled_propThresh <- ifelse(
      transformation == "none", propThresh,
      ifelse(transformation == "asin-sqrt", asin(sqrt(propThresh)), log(propThresh/(1-propThresh))))
    highCorMat <- meanMatProp >= scaled_propThresh & corMat >= corThresh
    ## convert the it into a long data
    identifiedPair <- highCorMat[upperIdx]
    ## create a vector indicating correlation groups
    correlationGroup <- vector(mode = "character", length = length(identifiedPair))
    correlationGroup[identifiedPair] <- "Identified"

    ## create a data.frame for identified pairs
    corResults <- data.frame(name1 = rownames(corMat)[upperIdx[, 1]], name2 = colnames(corMat)[upperIdx[,
        2]], mean = meanMatProp[upperIdx], max = maxMatProp[upperIdx], coefficient = corMat[upperIdx],
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
    yAxis <- corResults[, yScaleMetric]
    yTitle <- paste0("Barcode pair ", yScaleMetric, " proportion")
    if(!showRawProportion && transformation != "none") {
      yTitle <- paste0(yScaleMetric, " ", transformation, " proportion")
    } else if(showRawProportion && transformation != "none") {
      yTitle <- paste0(yScaleMetric, " ", transformation, " proportion\n(inverse-labelled in raw scale)")
    }
    
    ## scale the propThresh for visualisation
    scaled_propThresh <- ifelse(
      transformation == "none", propThresh,
      ifelse(transformation == "asin-sqrt", asin(sqrt(propThresh)), log(propThresh/(1-propThresh))))

    ## plotting correlations
      p <- ggplot(corResults, aes(x = coefficient)) + geom_histogram(aes(y = (after_stat(count))/max(after_stat(count)) *
        max(yAxis)), binwidth = 0.05, alpha = 0.3, fill = "grey") + geom_point(aes(y = yAxis,
        color = correlationGroup)) + stat_ecdf(geom = "step", aes(y = after_stat(y) * max(yAxis),
        color = correlationGroup), alpha = 0.5) + scale_color_manual(values = c(`pre-Defined` = "#00BFC4",
        Identified = "#F8766D", `non-Corr` = "grey")) + theme_classic() +
        labs(x = paste0(method, " correlation coefficient"), y = yTitle, color = "barcode pair") + scale_x_continuous(expand = c(0.05,
        0), limits = c(-1, 1)) + scale_y_continuous(sec.axis = sec_axis(~./max(yAxis), name = "cummulative freq.")) +
        geom_hline(yintercept = scaled_propThresh, linetype = "dashed", color = "#7CAE00",
            alpha = 0.8) + geom_vline(xintercept = corThresh, linetype = "dashed", color = "#C77CFF",
        alpha = 0.8) + annotate("text", x = corThresh - 0.1, y = max(yAxis) * 1.03, color = "#C77CFF",
        alpha = 1, size = 3, label = paste0("Cor=", corThresh)) + annotate("text", x = corThresh -
        0.4, y = scaled_propThresh - max(yAxis) * 0.02, color = "#7CAE00", alpha = 1,
        size = 3, label = paste0("Raw prop=", propThresh))
      
    ## ---- when transformation == the transformed, and showRawProp is true ----
    if((transformation == "asin-sqrt" || transformation == "logit") && showRawProportion){
      pb <- ggplot_build(p)
      # Extract y scale information
      breaks_trans <- pb$layout$panel_params[[1]]$y$breaks
      # Convert breaks back to original scale using inverse: sin^2(x)
      if(transformation == "asin-sqrt") {
        breaks_orig <- sin(breaks_trans)^2
      } else {
        ## inverse: e^x/(1+e^x) ---- but this is reversed to the proportion based on countPlus
        breaks_orig <- exp(breaks_trans)/(1+ exp(breaks_trans))
      }
      ## re-apply reversed y-axis labels
      p <- p +
        scale_y_continuous(
          breaks = breaks_trans,
          labels = format(round(breaks_orig, 4), scientific = FALSE),
          sec.axis = sec_axis(~./max(yAxis), name = "cumulative freq.")
        )
      message(paste0("showing raw proportion in y, but plotted by transfromation"))
    }
      

    return(p)
}


#' Cluster correlated Barcodes based on pairwise correlation
#'
#' `clusterCorrelatingBarcodes()` groups Barcode into clusters by
#'  identifying correlated Barcode pairs based on two criteria:
#'   * The proportions of Barcodes across samples exhibit high correlation,
#'    exceeding the threshold specified by `corThresh`.
#'   * The mean proportion of the Barcode pair exceeds the threshold specified
#'    by `propThresh`.
#'  These two parameters can be optimized up to the users using the visualization
#'  function [plotBarcodePairCorrelation].
#'
#' @param barbieQ A `barbieQ` object created by the [createBarbieQ] function.
#' @param method A string specifying the correlation method to use.
#'  Defaults to 'pearson'. Options include: 'pearson', 'kendall', 'spearman'.
#' @param transformation A string specifying the transformation method for proportion 
#'  to be used for computing barcode pair-wise correlation, and mean proportion.
#'  Options include: 'asin-sqrt', 'logit', and 'none'.
#'  Defaults to 'asin-sqrt'.
#' @param corThresh A numeric value that sets the threshold for high correlation
#'  Defaults to 0.95
#' @param propThresh A numeric value that sets the minimum level of
#'  Barcode pair's mean proportion for a Barcode pair to be considered as
#'  highly correlated co-existing Barcodes. Defaults to 0.001
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
clusterCorrelatingBarcodes <- function(barbieQ, method = "pearson", transformation = "asin-sqrt", 
    corThresh = 0.95, propThresh = 0.001, preDefinedCluster = NULL) {
    ## dispatch preprocessing function extract barcode pairwise correlation and
    ## pre-defined pairs
    barbieQ <- extractBarcodePairs(barbieQ, method = method, transformation = transformation, preDefinedCluster = preDefinedCluster)

    ## processed info is saved in this DFrame
    corDF <- SummarizedExperiment::rowData(barbieQ)$barcodeCorrelation
    preDefinedDf <- S4Vectors::metadata(corDF)$preDefinedBarcodePair
    meanMatProp <- S4Vectors::metadata(corDF)$meanMatProp
    corMat <- as.matrix(corDF)

    ## determine correlated pairs by passing correlation coefficient and proportion
    ## thresholds.
    ## scale the propThresh 
    scaled_propThresh <- ifelse(
      transformation == "none", propThresh,
      ifelse(transformation == "asin-sqrt", asin(sqrt(propThresh)), log(propThresh/(1-propThresh))))
    highCorMat <- meanMatProp >= scaled_propThresh & corMat >= corThresh
    ## arrow to extract the upper triangle, while also defined as 'correlated pairs'
    identifiedMat <- upper.tri(corMat) & highCorMat
    ## get row and column indices of identified pairs
    upperIdx <- which(identifiedMat, arr.ind = TRUE)

    ## create a data.frame for identified pairs
    highCorDf <- data.frame(name1 = rownames(corMat)[upperIdx[, 1]], name2 = colnames(corMat)[upperIdx[,
        2]], meanProp = meanMatProp[identifiedMat], coefficient = corMat[identifiedMat])

    ## combine pre-defined pairs and identified pairs
    if(ncol(preDefinedDf) == 2) {
      colnames(preDefinedDf) <- c("name1", "name2")
    }
    totalPairs <- rbind(preDefinedDf, highCorDf[, c("name1", "name2")])

    ## plotting highly correlating Barcode pairs
    clusterStructure <- igraph::graph_from_edgelist(totalPairs %>%
        as.matrix(), directed = FALSE)
    ## cluster Barcode groups based on correlating pairs
    groups <- igraph::components(clusterStructure)$membership
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
    ## save the graph to the metadata of `barcodeCorrelatedCluster`
    S4Vectors::metadata(barcodeCorrelatedCluster)$clusterStructure <- clusterStructure
    S4Vectors::metadata(barcodeCorrelatedCluster)$transformation <- transformation
    ## save object
    SummarizedExperiment::rowData(barbieQ)$barcodeCorrelatedCluster <- barcodeCorrelatedCluster

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
#' @param transformation A string specifying the transformation method for proportion 
#'  to be used for computing barcode pair-wise correlation, and mean proportion.
#'  Options include: 'asin-sqrt', 'logit', and 'none'.
#'  Defaults to 'asin-sqrt'.
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
extractBarcodePairs <- function(barbieQ, method = "pearson", transformation = "asin-sqrt", preDefinedCluster = NULL) {
    ## check method
    method <- match.arg(method, c("pearson", "kendall", "spearman"))
    ## check transformation
    transformation <- match.arg(transformation, c("asin-sqrt", "logit", "none"))
    ## extract data
    if(transformation == "none"){
      mat <- SummarizedExperiment::assays(barbieQ)$proportion %>% as.matrix()
    } else if(transformation == "asin-sqrt") {
      mat <- SummarizedExperiment::assays(barbieQ)$proportion %>% as.matrix()
      mat <- asin(sqrt(mat))
    } else if(transformation == "logit") {
      countPlus <- SummarizedExperiment::assay(barbieQ) + 0.5
      proportionPlus <- countPlus / colSums(countPlus)
      mat <- log((proportionPlus)/(1 - proportionPlus)) %>% as.matrix()
    }
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
    S4Vectors::metadata(corDF)$transformation <- transformation
    message("processing Barcode pairwise ", method, " correlation on propotion (", transformation, " transformation).")

    ## compute mean proportion for each barcode
    meanProp <- rowMeans(mat)
    ## N x N symmetric matrix with average between each pair on their row means
    meanMat <- outer(meanProp, meanProp, FUN = function(x, y) (x + y)/2)
    ## N x N symmetric matrix with LFC between each pair on their proportion
    # logMat <- log2(mat + 0.5)
    # LFCMat <- outer(logMat, logMat, FUN = function(x, y) (x - y))
    ## save meanMat to DFrame as a metadata
    S4Vectors::metadata(corDF)$meanMatProp <- meanMat
    # S4Vectors::metadata(corDF)$LFCMatProp <- LFCMat

    ## save preDefinedDf to corDF metadata
    S4Vectors::metadata(corDF)$preDefinedBarcodePair <- preDefinedDf
    ## save corDF to rowData of barbieQ
    SummarizedExperiment::rowData(barbieQ)$barcodeCorrelation <- corDF

    return(barbieQ)
}

#' Visualising the barcode clusters, the cluster sizes, and barcode proportion 
#'  within each cluster, after running the `clusterCorrelatingBarcodes` function.
#'
#' `inspectCorrelatingBarcodes()` visualize the cluster details generated by 
#'  calling `clusterCorrelatingBarcodes()`.
#' 
#' @param barbieQ A `barbieQ` object with information of clusters 
#'  of correlating barcodes, which should be predicted by the 
#'  [clusterCorrelatingBarcodes] function, after initialized the barbieQ object
#'  by the [createBarbieQ] function.
#'
#' @returns A list containing three ggplot objects, each for visualizing 
#'  the cluster structures, number of barcodes per cluster, and barcode proportion per cluster.
#'  
#' @export
#' 
#' @import igraph
#' @import scales 
#' @import ggplot2 
#' @import ggbreak
#' @import ggraph
#' @importFrom stats setNames
#'
#' @examples
#' nbarcodes <- 50
#' nsamples <- 12
#' count <- abs(matrix(rnorm(nbarcodes * nsamples), nbarcodes, nsamples))
#' rownames(count) <- paste0('Barcode', seq_len(nbarcodes))
#' barbieQ <- createBarbieQ(count)
#' barbieQ <- clusterCorrelatingBarcodes(barbieQ, preDefinedCluster = list(
#'  group1 = c('Barcode1', 'Barcode2', 'Barcode3'), group2 = c('Barcode4', 'Barcode5')))
#' inspectCorrelatingBarcodes(barbieQ)
inspectCorrelatingBarcodes <- function(barbieQ) {
  ## extract the clusterStructure from barbieQ
  clusterStructure <- barbieQ@elementMetadata$barcodeCorrelatedCluster@metadata$clusterStructure
  transformation <- barbieQ@elementMetadata$barcodeCorrelatedCluster@metadata$transformation
  transformation <- "none"
  ## check if the clusterStructure exist
  if(is.null(clusterStructure)) {
    stop("clusterStructure does not exist. Please run `clusterCorrelatingBarcodes` first.")
  }
  
  ## get cluster membership as factor, preserving level order
  memberships <- igraph::components(clusterStructure)$membership
  groups <- factor(memberships)
  ## get ggplot colors
  n_groups <- length(levels(groups))
  gg_colors <- scales::hue_pal()(n_groups)
  ## map colors to factor levels
  color_map_cluster <- stats::setNames(gg_colors, levels(groups))
  ## assign colors to vertices
  igraph::V(clusterStructure)$color <- color_map_cluster[as.character(groups)]
  
  ## alternative plotting as a basic R plot
  # ## create layout coordinates first
  # layout_coords <- igraph::layout_with_kk(clusterStructure)
  # ## extract cluster membership
  # memberships <- igraph::components(clusterStructure)$membership
  # ## plot igraph cluster centers are too far away so we plot barcodes name here
  # igraph::plot.igraph(
  #   clusterStructure,
  #   vertex.color = igraph::V(clusterStructure)$color,
  #   edge.color = "grey60",
  #   vertex.size = 12,         # ↑ bigger dots (nodes)
  #   vertex.frame.color = "grey60",
  #   vertex.label.cex = 0.2,   # ↓ smaller text
  #   vertex.label.dist = 0.5,  # adjust label distance if needed
  #   vertex.label.color = "grey30",
  #   layout = layout_coords
  # )
  
  ## alternative plotting as a ggplot object
  V(clusterStructure)$cluster <- groups
  ## Fruchterman-Reingold layout
  p_cluster <- ggraph::ggraph(clusterStructure, layout = "fr") +
    ggraph::geom_edge_link(color = "grey80") +
    ggraph::geom_node_point(aes(color = cluster), size = 5) +
    ggraph::geom_node_text(aes(label = name), repel = TRUE, size = 2) +
    theme_void()
  
  ## extract proportion based on transformation
  if(transformation == "none"){
    mat <- SummarizedExperiment::assays(barbieQ)$proportion %>% as.matrix()
  } else if(transformation == "asin-sqrt") {
    mat <- SummarizedExperiment::assays(barbieQ)$proportion %>% as.matrix()
    mat <- asin(sqrt(mat))
  } else if(transformation == "logit") {
    countPlus <- SummarizedExperiment::assay(barbieQ) + 0.5
    proportionPlus <- countPlus / colSums(countPlus)
    mat <- log((proportionPlus)/(1 - proportionPlus)) %>% as.matrix()
  }
  
  ## compute barcode feature
  BC_feature <- data.frame(
    Barcode = rownames(barbieQ),
    Cluster = barbieQ@elementMetadata$barcodeCorrelatedCluster$cluster,
    proportion = rowMeans(mat)
  )
  
  ## rename the clusters and unclustered groups
  BC_feature <- BC_feature %>% 
    mutate(Group = ifelse(
      Cluster < 0, "Unclustered \n barcodes", paste0("Cluster ", Cluster))) %>% 
    mutate(Group = factor(
      Group, levels = c(
        "Unclustered \n barcodes", 
        paste0("Cluster ", sort(unique(Cluster[Cluster>0]))))))
  
  ## set up color scale
  color_map_group <- setNames(
    c("grey", color_map_cluster), 
    c("Unclustered \n barcodes", paste0("Cluster ", names(color_map_cluster))))
  
  # Step 1: compute counts per Group
  group_counts <- BC_feature %>%
    count(Group, name = "count")
  # Step 2: find second largest
  second_largest <- sort(unique(group_counts$count), decreasing = TRUE)[2]
  # Step 3: total max
  max_count <- max(group_counts$count)
  
  p_cluster_size <- ggplot2::ggplot(data = BC_feature) +
    geom_bar(aes(x = Group, fill = Group), color = "white") +
    scale_fill_manual(values = color_map_group) +
    theme_classic() +
    labs(y = "Num. of barcodes") +
    theme(
      # aspect.ratio = 1, 
      legend.position = "none", 
      axis.text.x = element_text(angle = 45, hjust =1, vjust = 1), 
      axis.title.x = element_blank()) + 
    ggbreak::scale_y_break(
      c(second_largest+1, round((max_count)*0.95)), 
      scales = 0.1, 
      ticklabels = c(round((max_count)*0.95)+1, max_count))
  
  p_cluster_prop <- ggplot2::ggplot(data = BC_feature) +
    geom_point(aes(x = Group, y = proportion), color = "grey", shape = 1, size = 4) +
    geom_point(
      data = filter(BC_feature, Cluster > 0), 
      aes(x = Group, y = proportion, color = Group), shape = 1, size = 4, stroke = 1) +
    scale_color_manual(values = color_map_group) +
    theme_classic() +
    theme(
      # aspect.ratio = 1, 
      legend.position = "none", 
      axis.text.x = element_text(angle = 45, hjust =1, vjust = 1), 
      axis.title.x = element_blank()) + 
    labs(y = paste0(
      "Mean ", ifelse(transformation=="none", "", transformation), " proportion"))
  
  if((transformation == "asin-sqrt" || transformation == "logit")){
    pb <- ggplot_build(p_cluster_prop)
    # Extract y scale information
    breaks_trans <- pb$layout$panel_params[[1]]$y$breaks
    # Convert breaks back to original scale using inverse: sin^2(x)
    if(transformation == "asin-sqrt") {
      breaks_orig <- sin(breaks_trans)^2
    } else {
      ## inverse: e^x/(1+e^x) ---- but this is reversed to the proportion based on countPlus
      breaks_orig <- exp(breaks_trans)/(1+ exp(breaks_trans))
    }
    ## re-apply reversed y-axis labels
    p_cluster_prop <- p_cluster_prop +
      scale_y_continuous(
        breaks = breaks_trans,
        labels = format(round(breaks_orig, 4), scientific = FALSE)
      ) +
      labs(y = paste0(
        "Mean ", ifelse(transformation=="none", "", transformation), " proportion",
        "\n(inverse-labelled in raw scale)"))
    message(paste0("showing raw proportion in y, but plotted by transfromation"))
  }
  
  
  ## return a list of three plots
  p_list <- list(
    p_cluster = p_cluster, 
    p_cluster_size = p_cluster_size, 
    p_cluster_prop = p_cluster_prop)
  
  return(p_list)
}

#' Visualising the proportion trajectories of barcodes within each cluster,
#'  after running the `clusterCorrelatingBarcodes` function.
#'
#' `inspectCorrelatingBarcodesDetailed()` visualizes per-barcode proportion
#'  trajectories across samples, faceted by cluster, for clustered barcodes
#'  generated by calling `clusterCorrelatingBarcodes()`.
#'
#' @param barbieQ A `barbieQ` object with information of clusters
#'  of correlating barcodes, which should be predicted by the
#'  [clusterCorrelatingBarcodes] function, after initialized the barbieQ object
#'  by the [createBarbieQ] function.
#' @param transformation A character string specifying the transformation
#'  applied to barcode proportions before plotting. One of `"none"` (default),
#'  `"asin-sqrt"`, or `"logit"`. Should match the transformation used in
#'  [clusterCorrelatingBarcodes] for consistency.
#' @param ncol An integer specifying the number of columns in the facet layout.
#'  Defaults to `NULL`, which automatically arranges facets in an approximately
#'  square grid based on the number of clusters.
#'
#' @returns A `ggplot` object with one facet per cluster, showing per-barcode
#'  proportion trajectories across samples. Samples are ordered by mean
#'  proportion within each cluster. Unclustered barcodes (cluster <= 0) are
#'  excluded.
#'
#' @export
#'
#' @import ggplot2
#' @import scales
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr group_by mutate ungroup distinct arrange pull filter
#' @importFrom ggh4x facet_wrap2 facetted_pos_scales
#' @importFrom stats setNames
#' @importFrom SummarizedExperiment assays assay
#'
#' @examples
#' nbarcodes <- 50
#' nsamples <- 12
#' count <- abs(matrix(rnorm(nbarcodes * nsamples), nbarcodes, nsamples))
#' rownames(count) <- paste0('Barcode', seq_len(nbarcodes))
#' barbieQ <- createBarbieQ(count)
#' barbieQ <- clusterCorrelatingBarcodes(barbieQ, preDefinedCluster = list(
#'  group1 = c('Barcode1', 'Barcode2', 'Barcode3'), group2 = c('Barcode4', 'Barcode5')))
#' inspectCorrelatingBarcodesDetailed(barbieQ)
inspectCorrelatingBarcodesDetailed <- function(barbieQ, transformation = "none", ncol = NULL) {
  
  ## check if cluster exists
  clusterInfo <- barbieQ@elementMetadata$barcodeCorrelatedCluster
  if (is.null(clusterInfo)) {
    stop("Cluster information does not exist. Please run `clusterCorrelatingBarcodes` first.")
  }
  
  ## extract cluster membership
  cluster_vec <- clusterInfo$cluster
  
  ## filter to clustered barcodes only (cluster > 0)
  flag_clustered <- cluster_vec > 0
  if (sum(flag_clustered) == 0) {
    stop("No clustered barcodes found (all cluster values <= 0).")
  }
  cluster_BC <- barbieQ[flag_clustered, ]
  cluster_group <- cluster_vec[flag_clustered]
  
  ## extract proportion matrix based on transformation
  if (transformation == "none") {
    mat <- SummarizedExperiment::assays(cluster_BC)$proportion %>% as.matrix()
  } else if (transformation == "asin-sqrt") {
    mat <- SummarizedExperiment::assays(cluster_BC)$proportion %>% as.matrix()
    mat <- asin(sqrt(mat))
  } else if (transformation == "logit") {
    countPlus <- SummarizedExperiment::assay(cluster_BC) + 0.5
    proportionPlus <- countPlus / colSums(countPlus)
    mat <- log(proportionPlus / (1 - proportionPlus)) %>% as.matrix()
  } else {
    stop("transformation must be one of: 'none', 'asin-sqrt', 'logit'")
  }
  
  ## pivot to long format
  df_long <- as.data.frame(mat) %>%
    tibble::rownames_to_column("barcode") %>%
    tidyr::pivot_longer(
      cols = -barcode,
      names_to = "sample",
      values_to = "prop"
    )
  
  ## add cluster group label
  cluster_group_named <- stats::setNames(cluster_group, rownames(cluster_BC))
  df_long$cluster_group <- paste0("Cluster ", cluster_group_named[df_long$barcode])

  ## reorder samples by mean proportion within each cluster
  df_long <- df_long %>%
    dplyr::group_by(cluster_group, sample) %>%
    dplyr::mutate(sample_mean = mean(prop)) %>%
    dplyr::ungroup() 
    
  ## get consistent cluster colors (matching inspectCorrelatingBarcodes)
  n_groups <- length(unique(cluster_group))
  gg_colors <- scales::hue_pal()(n_groups)
  color_map <- stats::setNames(gg_colors, paste0("Cluster ", sort(unique(cluster_group))))
  
  ## build y-axis label
  y_label <- paste0(
    "Barcode proportion",
    ifelse(transformation == "none", "", paste0("\n(", transformation, " transformed)"))
  )
  
  ## default ncol to approximate square layout
  n_clusters <- length(unique(cluster_group))
  if (is.null(ncol)) {
    ncol <- ceiling(sqrt(n_clusters))
  }
  
  ## build per-facet x scale list
  scale_list <- lapply(unique(df_long$cluster_group), function(cg) {
    samples_ordered <- df_long %>%
      dplyr::filter(cluster_group == cg) %>%
      dplyr::distinct(sample, sample_mean) %>%
      dplyr::arrange(sample_mean) %>%
      dplyr::pull(sample) %>%
      as.character()
    ggplot2::scale_x_discrete(limits = samples_ordered)
  })
  
  
  ## plot
  p <- ggplot2::ggplot(df_long, aes(x = sample, y = prop, group = barcode)) +
    ggplot2::geom_line(alpha = 0.7, aes(colour = cluster_group)) +
    ggplot2::geom_point(size = 0.8, color = "grey30") +
    ggh4x::facet_wrap2(~ cluster_group, scales = "free", ncol = ncol) +
    ggh4x::facetted_pos_scales(x = scale_list) +
    ggplot2::scale_color_manual(values = color_map) +
    ggplot2::labs(x = "Sample", y = y_label) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(p)
}

#' Merging correlating barcodes, by keeping a representative barcode
#'  and removing the rest for each cluster. The representative barcode is 
#'  the one with the highest mean proportion across samples.
#'
#' @param barbieQ_clustered A `barbieQ` object with information of clusters 
#'  of correlating barcodes, which should be predicted by the 
#'  [clusterCorrelatingBarcodes] function, after initialized the barbieQ object
#'  by the [createBarbieQ] function.
#' @param barbieQ_toMerge A `barbieQ` object, whose barcodes will be merged 
#'  based on clusters specified by the `barbieQ_clustered` parameter. 
#'  If not specified, defaults to the same `barbieQ` object as `barbieQ_clustered`.
#'
#' @returns A `barbieQ` object updated by merging the correlating barcodes into one 
#'  by selecting a representative from each cluster and deleting the rest. 
#'  The processing information saved in `metadata$predicted_barcode_clusters` and 
#'  `metadata$removed_barcodes`.
#' 
#' @import dplyr
#' @importFrom SummarizedExperiment assay
#' @importFrom S4Vectors metadata
#' 
#' @export
#'
#' @examples
#' nbarcodes <- 50
#' nsamples <- 12
#' count <- abs(matrix(rnorm(nbarcodes * nsamples), nbarcodes, nsamples))
#' rownames(count) <- paste0('Barcode', seq_len(nbarcodes))
#' barbieQ <- createBarbieQ(count)
#' barbieQ <- clusterCorrelatingBarcodes(barbieQ, preDefinedCluster = list(
#'  group1 = c('Barcode1', 'Barcode2', 'Barcode3'), group2 = c('Barcode4', 'Barcode5')))
#' mergeCorrelatingBarcodes(barbieQ, barbieQ)
mergeCorrelatingBarcodes <- function(barbieQ_clustered, barbieQ_toMerge=NULL) {
  ## extract the predictedClusters from barbieQ_clustered
  predictedClusters <- barbieQ_clustered@elementMetadata$barcodeCorrelatedCluster
  ## check if the predictedClusters exist
  if(is.null(predictedClusters)) {
    stop("predictedClusters does not exist. Please run `clusterCorrelatingBarcodes` first.")
  }
  ## if `barbieQ_toMerge` is unspecified, default to merging the `barbieQ_clustered`
  if(is.null(barbieQ_toMerge)) {
    barbieQ_toMerge <- barbieQ_clustered
  }
  ## extract cluster information
  clusterArray <- predictedClusters$cluster
  ## fetch the clustered barcodes to be merged
  flag_clustered <- clusterArray > 0
  barbieQ_clustered_trimmed <- barbieQ_clustered[flag_clustered,]
  ## features of them
  BC_feature <- data.frame(
    Barcode = rownames(barbieQ_clustered_trimmed),
    meanCPM = rowMeans(barbieQ_clustered_trimmed@assays@data$CPM),
    cluster = clusterArray[flag_clustered])
  ## locate the barcode with maximum meanCPM
  BC_feature <- BC_feature %>% 
    dplyr::group_by(cluster) %>% 
    ##this select one barcode per cluster even with ties
    dplyr::mutate(selectMax = row_number() == which.max(meanCPM)) %>%
    dplyr::ungroup() %>% 
    as.data.frame()
  rownames(BC_feature) <- BC_feature$Barcode
  
  ## print barcodes to be merged
  removed_barcodes <- BC_feature$Barcode[!BC_feature$selectMax]
  print(paste0("printing removed barcodes: ", removed_barcodes))
  
  ## fetch barcodes in the dataset to be merged
  BC_pre_merge <- data.frame(
    Barcode = rownames(barbieQ_toMerge),
    toKeep = TRUE)
  rownames(BC_pre_merge) <- rownames(barbieQ_toMerge)
  
  ## find overlaps between the barbieQ_clustered and barbieQ_toMerge
  common_rows <- intersect(BC_feature$Barcode, BC_pre_merge$Barcode)
  BC_pre_merge[common_rows, "toKeep"] <- BC_feature[common_rows,"selectMax"]
  
  ## merging by selecting the representative barcodes
  ## !!! re-calculating the proportion by calling `createBarbieQ`
  barbieQ_merged <- createBarbieQ(
    object = SummarizedExperiment::assay(barbieQ_toMerge)[BC_pre_merge$toKeep,], 
    sampleMetadata = barbieQ_toMerge$sampleMetadata, 
    factorColors =  S4Vectors::metadata(barbieQ_toMerge)$factorColors)
  message("!! re-computing barcode proportion, CPM, rank... from the selected barcodes.")
  S4Vectors::metadata(barbieQ_merged) <- S4Vectors::metadata(barbieQ_toMerge)
  SummarizedExperiment::colData(barbieQ_merged) <- SummarizedExperiment::colData(barbieQ_toMerge) 
  
  ## saving clustering information
  barbieQ_merged@metadata$predicted_barcode_clusters <- predictedClusters@metadata$clusterStructure
  barbieQ_merged@metadata$removed_barcodes <- removed_barcodes
  
  
  return(barbieQ_merged)
}