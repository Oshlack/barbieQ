#' Tag each Barcode as being part of the major contributors or not
#'
#' `tagTopBarcodes` is designed for better filtering out the background noise,
#'  i.e., Barcodes that consistently have low contributions across samples.
#'  Each Barcode is tagged as one of the major contributing Barcodes or not
#'  in each sample, referred to *top* Barcodes. In the entire dataset,
#'  a Barcode is considered *top* if it is tagged as *top*
#'  in a number of samples passing a defined threshold
#'  across all selected samples.
#'
#' @param Barbie an object created by the [createBarbie] function
#' @param activeColumns a logical vector indicating which column (sample)
#'  to be considered when determining the *top* Barcodes
#' @param threshold a number ranging from 0 to 1 used for thresholding the Barcodes' contribution
#' @param minTopAppearance An integer indicating the minimum times a barcode must be 'top' in each column to be considered a 'top barcode' across samples
#'
#' @return A 'Barbie' object, including components:
#'  * `isTop` (updated): a list containing a vector `vec` and a matrix `mat`,
#'    tagging each Barcode as being part of the major contributors or not.
#'  * Other components inherited from the `Barbie` object passed into the
#'    function. See `Barbie` structure in [createBarbie].
#'
#' @export
#'
#' @examples
#' ## create a `Barbie` object
#' ## sample conditions and color palettes
#' sampleConditions <- data.frame(
#'   Treat=factor(rep(c("ctrl", "drug"), each=6)),
#'   Time=rep(rep(1:2, each=3), 2))
#' conditionColor <- list(
#'   Treat=c(ctrl="#999999", drug="#112233"),
#'   Time=c("1"="#778899", "2"="#998877"))
#' ## Barcode count data
#' nbarcodes <- 50
#' nsamples <- 12
#' barcodeCount <- abs(matrix(rnorm(nbarcodes*nsamples), nbarcodes, nsamples))
#' rownames(barcodeCount) <- paste0("Barcode", 1:nbarcodes)
#' ## create a `Barbie` object
#' myBarbie <- createBarbie(barcodeCount, sampleConditions, conditionColor)
#' ## tag top Barcodes
#' tagTopBarcodes(myBarbie)
tagTopBarcodes <- function(Barbie, activeColumns=NULL,
                           threshold=0.99, minTopAppearance=1){
  mat <- Barbie$assay
  ## dispatch 'returnNumMat' function to ensure the object is a numeric matrix apart from NAs.
  if(inherits(mat, "data.frame") || inherits(mat, "matrix") || is.vector(mat))
    mat <- returnNumMat(mat)
  else stop("`Barbie$assay` should be a data.frame, matrix or vector of Barcode counts.")
  ## if 'activeColumns' is not specified, assume all columns are active.
  if(is.null(activeColumns)) activeColumns <- rep(TRUE, ncol(mat))
  ## dispatch 'tagTopEachColumn' to individully determine 'top Barcodes' in each column.
  topsInMat <- apply(mat, 2, function(col) tagTopEachColumn(col, threshold))
  ## select the active columns used to determine 'top Barcodes' out of the samples.
  subTopMat <- topsInMat[, activeColumns, drop = FALSE]
  ## determine the 'top Barcodes' based on the numbers of being true out of all 'active' samples.
  topOverall <- rowSums(subTopMat, na.rm = TRUE) >= minTopAppearance
  ## store
  updatedObject <- list(
    ## retain other components
    assay = Barbie$assay,
    metadata = Barbie$metadata,
    proportion = Barbie$proportion,
    CPM = Barbie$CPM,
    occurrence = Barbie$occurrence,
    rank = Barbie$rank,
    ## update `isTop`
    isTop = list(vec = topOverall,
                 mat = topsInMat),
    ## retain other components
    clusters = Barbie$clusters,
    factorColors = Barbie$factorColors)

  # include other existing elements in 'Barbie' that are not mentioned above.
  for (elementName in names(Barbie)) {
    if (!(elementName %in% names(updatedObject)))
      updatedObject[[elementName]] <- Barbie[[elementName]]
  }

  return(updatedObject)
}

#' Tag top Barcodes in each column (sample).
#'
#' @param tempCol a numeric vector of Barcode counts in a sample, usually being a column in a count matrix
#' @param threshold a number ranging from 0 to 1 used for thresholding the Barcodes' contribution
#'
#' @return a logical vector indicating whether each Barcode is considered as 'top' in this vector
#'
#' @noMd
#'
#' @import stats
#' @examples
#' myCount <- c(1, 2, 3, 98, NA, NA)
#' tagTopEachColumn(myCount, 0.99)
tagTopEachColumn <- function(tempCol, threshold=0.99){
  ## need update this function: with a more robust threshold instead of empirical 0.99

  ## set up default threshold at 0.99
  if(is.null(threshold)) threshold <- 0.99
  ## NAs are retained as NA in 'Barbie$rank', but we want NAs to be placed last in the ranking here.
  ## setting 'na.last=TRUE' will place NAs last in the ranking.
  ## setting 'ties.method="first"' will handle equal values by assigning the first the lowest rank.
  ## setting 'minus x' will rank 'x' by decreasing order.
  rankCol <- rank(-tempCol, ties.method = "first", na.last = TRUE)
  ## sort Barcodes in the column by decreasing order.
  ## setting 'na.last=TRUE' will place NAs last in the order.
  sortedValue <- sort(tempCol, decreasing = TRUE, na.last = TRUE) #sort the column by decreasing order

  ## find the position where the cumulative sum is > threshold and the previous value not exceeds it
  ## step1: calculate cumulative values from the first to the current value in the sorted values, ordered in decreasing order.
  ## NAs will be retained as NA in the cumulative values-cumSum.
  cumSum <- vapply(rankCol, function(x) sum(sortedValue[1:x]), numeric(1))
  ## step2: calculate cumulative values from the first to the (current ranking + 1) value in the sorted values.
  cumSumMinus1 <- vapply(rankCol, function(x) {
    if(x > 1) sum(sortedValue[1:(x-1)])
    else if (x == 1) sum(sortedValue[1])
    }, numeric(1))
  ## step3: calculate sum of all values from the first to the last, omitting NAs
  sums <- sum(na.omit(tempCol))
  ## step4: find the position where the cumulative sum is >= threshold and the previous value not exceeds it
  position <- which(
    replace(cumSum >= threshold*sums, is.na(cumSum), FALSE) &
      replace(cumSumMinus1 < threshold*sums, is.na(cumSumMinus1), FALSE))
  ## step5: determine all values ranking higher than the position is 'top'.
  isTopInColumn <- rankCol <= rankCol[position]

  return(isTopInColumn)
}

