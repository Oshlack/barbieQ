#' Tag each Barcode as being part of the major contributors or not
#'
#' `tagTopBarcodes()` is designed for better filtering out the background noise,
#'  i.e., Barcodes that consistently have low contributions across samples.
#'  Each Barcode is tagged as one of the major contributing Barcodes
#'  (\emph{top} Barcodes) or not in each sample, based on whether their
#'  combined proportion passes a defined threshold in the sample.
#'  Across the entire dataset, a Barcode is considered \emph{top}
#'  if it is tagged as \emph{top} in a number of samples,
#'  meeting a specified appearance threshold across all selected samples.
#'
#' @param barbieQ An object created by the [createBarbieQ] function.
#' @param activeSamples a logical vector indicating which samples (columns)
#'  to consider when determining the \emph{top} Barcodes across the entire
#'  dataset. Default to considering all samples in the `barbieQ` object.
#' @param proportionThreshold A numeric value ranging from 0 to 1,
#'  used as a threshold for determining \emph{top} Barcodes in each sample
#'  based on their combined proportion in theat sample. Default to 0.99.
#' @param minTopFrequency An integer specifying the minimum number of times
#'  a Barcode must be tagged as \emph{top} across the selected samples
#'  (specified by `activeSamples`) in order to be considered \emph{top}
#'  for the entire dataset. Default to 1.
#'
#' @return A `barbieQ` object updating `isTop` while inheriting other components
#' from the `barbieQ` object passed into the function.
#' See `barbieQ` structure in [createBarbieQ].
#' `isTop` is a list including:
#'  * `vec`: a logical vector tagging Barcodes as \emph{top}
#'    across the entire dataset.
#'  * `mat`: a logical matrix tagging Barcodes as \emph{top} in each sample.
#'
#' @export
#'
#' @examples
#' ## create a `barbieQ` object
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
#' barcodeCount <- abs(matrix(rnorm(nbarcodes * nsamples), nbarcodes, nsamples))
#' rownames(barcodeCount) <- paste0("Barcode", seq_len(nbarcodes))
#' ## create a `barbieQ` object
#' myBarbieQ <- createBarbieQ(barcodeCount, sampleConditions, conditionColor)
#' ## tag top Barcodes
#' tagTopBarcodes(myBarbieQ)
tagTopBarcodes <- function(barbieQ, activeSamples = NULL,
                           proportionThreshold = 0.99, minTopFrequency = 1) {
  mat <- barbieQ$assay
  ## dispatch 'returnNumMat' function to ensure the object is a numeric matrix apart from NAs.
  if (inherits(mat, "data.frame") || inherits(mat, "matrix") || is.vector(mat)) {
    mat <- returnNumMat(mat)
  } else {
    stop("`barbieQ$assay` should be a data.frame, matrix or vector of Barcode counts.")
  }
  ## if 'activeSamples' is not specified, assume all columns are active.
  if (is.null(activeSamples)) activeSamples <- rep(TRUE, ncol(mat))
  ## dispatch 'tagTopEachColumn' to individully determine 'top Barcodes' in each column.
  topsInMat <- apply(
    mat, 2, function(col) tagTopEachColumn(col, proportionThreshold)
  )
  ## select the active columns used to determine 'top Barcodes' out of the samples.
  subTopMat <- topsInMat[, activeSamples, drop = FALSE]
  ## determine the 'top Barcodes' based on the numbers of being true out of all 'active' samples.
  topOverall <- rowSums(subTopMat, na.rm = TRUE) >= minTopFrequency
  ## store
  updatedObject <- list(
    ## retain other components
    assay = barbieQ$assay,
    metadata = barbieQ$metadata,
    proportion = barbieQ$proportion,
    CPM = barbieQ$CPM,
    occurrence = barbieQ$occurrence,
    rank = barbieQ$rank,
    ## update `isTop`
    isTop = list(
      vec = topOverall,
      mat = topsInMat
    ),
    ## retain other components
    clusters = barbieQ$clusters,
    factorColors = barbieQ$factorColors
  )

  # include other existing elements in 'barbieQ' that are not mentioned above.
  for (elementName in names(barbieQ)) {
    if (!(elementName %in% names(updatedObject))) {
      updatedObject[[elementName]] <- barbieQ[[elementName]]
    }
  }

  return(updatedObject)
}

#' Tag \emph{top} Barcodes in each sample.
#'
#' @param tempCol a numeric vector of Barcode counts in a sample,
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
  if (is.null(proportionThreshold)) proportionThreshold <- 0.99
  ## NAs are retained as NA in 'barbieQ$rank', but we want NAs to be placed last in the ranking here.
  ## setting 'na.last=TRUE' will place NAs last in the ranking.
  ## setting 'ties.method="first"' will handle equal values by assigning the first the lowest rank.
  ## setting 'minus x' will rank 'x' by decreasing order.
  rankCol <- rank(-tempCol, ties.method = "first", na.last = TRUE)
  ## sort Barcodes in the column by decreasing order.
  ## setting 'na.last=TRUE' will place NAs last in the order.
  sortedValue <- sort(tempCol, decreasing = TRUE, na.last = TRUE) # sort the column by decreasing order

  ## find the position where the cumulative sum is > threshold and the previous value not exceeds it
  ## step1: calculate cumulative values from the first to the current value in the sorted values, ordered in decreasing order.
  ## NAs will be retained as NA in the cumulative values-cumSum.
  cumSum <- vapply(rankCol, function(x) sum(sortedValue[seq_len(x)]), numeric(1))
  ## step2: calculate cumulative values from the first to the (current ranking + 1) value in the sorted values.
  cumSumMinus1 <- vapply(rankCol, function(x) {
    if (x > 1) {
      sum(sortedValue[seq_len(x - 1)])
    } else if (x == 1) sum(sortedValue[1])
  }, numeric(1))
  ## step3: calculate sum of all values from the first to the last, omitting NAs
  sums <- sum(stats::na.omit(tempCol))
  ## avoid dividing by zero
  if (sums == 0) sums <- Inf
  ## step4: find the position where the cumulative sum is >= proportionThreshold and the previous value not exceeds it
  position <- which(
    replace(cumSum >= proportionThreshold * sums, is.na(cumSum), FALSE) &
      replace(
        cumSumMinus1 < proportionThreshold * sums, is.na(cumSumMinus1), FALSE
      )
  )
  ## step5: determine all values ranking higher than the position is 'top'.
  if (length(position) == 0L) {
    isTopInColumn <- rep(FALSE, length(rankCol))
  } else {
    isTopInColumn <- rankCol <= rankCol[position]
  }

  return(isTopInColumn)
}
