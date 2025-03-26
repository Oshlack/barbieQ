#' Tag each Barcode as being part of the major contributors or not
#'
#' `tagTopBarcodes()` tags \emph{top} Barcodes that collectively contribute to 
#'  the majority of counts across the dataset.
#'  It is designed for subsequently filtering out the background noise,
#'  i.e., filtering out Barcodes that consistently have low contributions across samples.
#'  In each sample, Barcodes are tagged as major contributing Barcodes
#'  (\emph{top} Barcodes) or not, based on whether their
#'  combined proportion passes a defined threshold in the sample.
#'  Across the entire dataset, a Barcode is considered \emph{top}
#'  if it is tagged as \emph{top} in a number of samples,
#'  meeting a specified appearance threshold across all selected samples.
#'
#' @param barbieQ A `SummarizedExperiment` object created by the [createBarbieQ] function.
#' @param activeSamples A logical vector indicating individual samples (columns)
#'  to consider or avoid when determining the \emph{top} Barcodes across the entire
#'  dataset. Default to considering all samples in the `barbieQ` object.
#' @param proportionThreshold A numeric value ranging from 0 to 1,
#'  used as a threshold for determining \emph{top} Barcodes in each sample
#'  based on their combined proportion in theat sample. Default to 0.99.
#' @param nSampleThreshold An integer specifying the minimum number of times
#'  a Barcode must be tagged as \emph{top} across the selected samples
#'  (specified by `activeSamples`) in order to be considered \emph{top}
#'  for the entire dataset. Default to 1.
#'  
#' @importClassesFrom S4Vectors DataFrame
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment assays
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment rowData
#'
#' @return A `barbieQ` object with slots `isTopAssay` and `isTopBarcode` added or updated, 
#' while inheriting other components from the `barbieQ` object passed into the function.
#' See `barbieQ` structure (as a `SummarizedExperiment` object) in [createBarbieQ].
#'  * `isTopAssay`: a logical matrix stored in `assays(barbieQ)` 
#'    tagging Barcodes as \emph{top} in each sample.
#'  * `isTopBarcode`: a `DataFrame` stored in `rowData(barbieQ)` 
#'    with a single logical column tagging Barcodes as \emph{top} or not across the entire dataset.
#'
#' @export
#'
#' @examples
#' ## create a `barbieQ` object
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
#' barcodeCount <- abs(matrix(rnorm(nbarcodes * nsamples), nbarcodes, nsamples))
#' rownames(barcodeCount) <- paste0('Barcode', seq_len(nbarcodes))
#' ## create a `barbieQ` object
#' myBarbieQ <- createBarbieQ(barcodeCount, sampleConditions, conditionColor)
#' ## tag top Barcodes
#' tagTopBarcodes(myBarbieQ)
tagTopBarcodes <- function(barbieQ, activeSamples = NULL, proportionThreshold = 0.99, nSampleThreshold = 1) {
    ## extract raw count
    mat <- SummarizedExperiment::assay(barbieQ)
    ## dispatch 'returnNumMat' function to ensure the object is a numeric matrix apart
    ## from NAs.
    if (inherits(mat, "data.frame") || inherits(mat, "matrix") || is.vector(mat)) {
        mat <- returnNumMat(mat)
    } else {
        stop("`assay(barbieQ)` should be a data.frame, matrix or vector of Barcode counts.")
    }
    ## if `activeSamples` is not specified, assume all columns are active.
    if (is.null(activeSamples))
        activeSamples <- rep(TRUE, ncol(mat))
    ## dispatch `tagTopEachColumn` to determine 'top Barcodes' in each column
    ## individually.
    topsInMat <- apply(mat, 2, function(col) tagTopEachColumn(col, proportionThreshold))
    ## select the active columns used to determine 'top Barcodes' out of the samples.
    subTopMat <- topsInMat[, activeSamples, drop = FALSE]
    ## determine the 'top Barcodes' based on the numbers of being true out of all
    ## 'active' samples.
    topOverall <- rowSums(subTopMat, na.rm = TRUE) >= nSampleThreshold
    ## store `topOverall` and `topsInMat`.
    rownames(topsInMat) <- rownames(barbieQ)
    colnames(topsInMat) <- colnames(barbieQ)
    ## assay slots in `assays` must have proper rownames and colnames.
    SummarizedExperiment::assays(barbieQ)$isTopAssay <- topsInMat
    SummarizedExperiment::rowData(barbieQ)$isTopBarcode <- S4Vectors::DataFrame(isTop = topOverall)

    return(barbieQ)
}

#' Tag \emph{top} Barcodes in each sample.
#'
#' @param tempCol A numeric vector of Barcode counts in a sample,
#'  usually being a column in a count matrix
#' @param proportionThreshold A numeric value ranging from 0 to 1,
#'  used as a threshold for determining \emph{top} Barcodes in each sample
#'  based on their combined proportion in theat sample. Default to 0.99.
#'
#' @return a logical vector tagging Barcodes as \emph{top} in the sample.
#'
#' @noRd
#'
#' @importFrom stats na.omit
#'
#' @examples \donttest{
#' myCount <- c(1, 2, 3, 98, NA, NA)
#' barbieQ:::tagTopEachColumn(myCount, 0.99)
#' }
tagTopEachColumn <- function(tempCol, proportionThreshold = 0.99) {
    ## set up default threshold at 0.99
    if (is.null(proportionThreshold))
        proportionThreshold <- 0.99
    ## NAs are retained as NA in 'barbieQ$rank', but we want NAs to be placed last in
    ## the ranking here.  setting 'na.last=TRUE' will place NAs last in the ranking.
    ## setting 'ties.method='first'' will handle equal values by assigning the first
    ## the lowest rank.  setting 'minus x' will rank 'x' by decreasing order.
    rankCol <- rank(-tempCol, ties.method = "first", na.last = TRUE)
    ## sort Barcodes in the column by decreasing order.  setting 'na.last=TRUE' will
    ## place NAs last in the order.
    sortedValue <- sort(tempCol, decreasing = TRUE, na.last = TRUE)  # sort the column by decreasing order

    ## find the position where the cumulative sum is > threshold and the previous
    ## value not exceeds it step1: calculate cumulative values from the first to the
    ## current value in the sorted values, ordered in decreasing order.  NAs will be
    ## retained as NA in the cumulative values-cumSum.
    cumSum <- vapply(rankCol, function(x) sum(sortedValue[seq_len(x)]), numeric(1))
    ## step2: calculate cumulative values from the first to the (current ranking + 1)
    ## value in the sorted values.
    cumSumMinus1 <- vapply(rankCol, function(x) {
        if (x > 1) {
            sum(sortedValue[seq_len(x - 1)])
        } else if (x == 1)
            sum(sortedValue[1])
    }, numeric(1))
    ## step3: calculate sum of all values from the first to the last, omitting NAs
    sums <- sum(stats::na.omit(tempCol))
    ## avoid dividing by zero
    if (sums == 0)
        sums <- Inf
    ## step4: find the position where the cumulative sum is >= proportionThreshold and
    ## the previous value not exceeds it
    position <- which(replace(cumSum >= proportionThreshold * sums, is.na(cumSum), FALSE) &
        replace(cumSumMinus1 < proportionThreshold * sums, is.na(cumSumMinus1), FALSE))
    ## step5: determine all values ranking higher than the position is 'top'.
    if (length(position) == 0L) {
        isTopInColumn <- rep(FALSE, length(rankCol))
    } else {
        isTopInColumn <- rankCol <= rankCol[position]
    }

    return(isTopInColumn)
}
