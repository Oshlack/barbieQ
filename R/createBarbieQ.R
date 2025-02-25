#' Extract Barcode count data and build into a `barbieQ` object for
#' subsequent analysis
#'
#' `createBarbieQ()` constructs an object-centred list called `barbieQ` object,
#' which is designed to process Barcode count data gained from
#' cell clonal tracking experiments.
#'
#' @param object A numeric `matrix` of Barcode counts, with Barcodes in rows
#'  and samples in columns;
#'  Alternatively, you can pass an existing `barbieQ` object, from which
#'  barcode counts, sample conditions, and the mapping colors will be inherited.
#' @param target A `matrix` or `data.frame` where rows represent samples 
#'  and columns represent factors, each containing the specific condition assigned to a sample.
#'  If not specified, all samples are assigned the same condition.
#'  If a `barbieQ` object is passed to `object`, the sample conditions
#'  will be updated by specifying `target`.
#' @param factorColors A `list` mapping factors in `target` to color palettes, 
#'  where each element is a named `vector` assigning colors to specific sample conditions.
#'  If not specified, defaults to `NULL`.
#'  If a `barbieQ` object is passed to `object`, the mapping colors
#'  will be updated by specifying `factorColors`.
#'
#' @return A `barbieQ` object, including several components:
#'  * `assay`: a `data.frame` of Barcode counts, with Barcodes in rows and samples in columns.
#'  * `metadata`: a `data.frame` of sample conditions, with samples in rows and
#'    experimental factors in columns.
#'  * `factorColors`: A `list` mapping factors in `target` to color palettes,
#'    where each element is a named `vector` assigning colors to specific sample conditions.
#'  * Matrices of `proportion`, `CPM`, `occurrence`, and `rank` representing
#'    Barcode proportion, Counts Per Million (CPM), presence/absence, and ranking
#'    in each sample, based on the Barcode counts in `assay`.
#'  * `isTop`: a `list` containing a `vector` (`vec`) and a `matrix` (`mat`),
#'    indicating whether each Barcode is a major contributor, determined by calling [tagTopBarcodes].
#'  * Additional processed information for further analysis.
#'
#' @export
#'
#' @importFrom utils menu
#'
#' @examples
#' ## Sample conditions and color palettes
#' sampleConditions <- data.frame(
#'   Treat = factor(rep(c("ctrl", "drug"), each = 6)),
#'   Time = rep(rep(seq_len(2), each = 3), 2)
#' )
#' conditionColor <- list(
#'   Treat = c(ctrl = "#999999", drug = "#112233"),
#'   Time = c("1" = "#778899", "2" = "#998877")
#' )
#'
#' ## Barcode count data
#' nbarcodes <- 50
#' nsamples <- 12
#' barcodeCount <- abs(matrix(rnorm(nbarcodes * nsamples), nbarcodes, nsamples))
#' rownames(barcodeCount) <- paste0("Barcode", seq_len(nbarcodes))
#'
#' ## Passing `object`, `target` and `factorColors`
#' createBarbieQ(object = barcodeCount)
#' createBarbieQ(object = barcodeCount, target = sampleConditions)
#' createBarbieQ(
#'   object = barcodeCount, target = sampleConditions, factorColors = conditionColor
#' )
#'
#' ## Updating a `barbieQ` object by passing new `target` and `factorColors`
#' myBarbieQ <- createBarbieQ(barcodeCount, sampleConditions, conditionColor)
#' createBarbieQ(
#'   object = myBarbieQ, target = data.frame(Mouse = rep(seq_len(4), each = 3))
#' )
#' createBarbieQ(
#'   object = myBarbieQ, target = data.frame(Mouse = rep(seq_len(4), each = 3)),
#'   factorColors = list(
#'     Mouse = c("1" = "#111199", "2" = "#112200", "3" = "#441111", "4" = "#000000")
#'   )
#' )
createBarbieQ <- function(object, target = NULL, factorColors = NULL) {
  ## define the barbieQ object structure
  barbieQ <- list(
    ## assay stores Barcode counts.
    assay = NULL,
    ## metadata stores experimental design of each sample, as the equivalent to 'target file' in single cell data.
    metadata = NULL,
    proportion = matrix(),
    CPM = matrix(),
    occurrence = matrix(),
    rank = matrix(),
    isTop = list(),
    clusters = factor(levels = character(), labels = character()),
    ## factorColors stores color palettes corresponding to the factors in 'metadata'.
    factorColors = NULL
  )

  ## extract Barcode counts from `object`
  ## dispatch `returnNumMat` function to ensure the object is a numeric matrix apart from NAs.
  if (inherits(object, "data.frame") ||
    inherits(object, "matrix") ||
    is.vector(object) && !(is.list(object))) {
    barbieQ$assay <- returnNumMat(object)
  } else if (is.list(object)) {
    ## inherit components if a `barbieQ` object is passed to `object`
    checkBarbieQDimensions(object)
    barbieQ <- object
    ## extract barcode counts from other objects. Need a function for it.
    # if Seurat, extract assay and metadata.
    # if BARtab object, extract ...
  }
  if (!nrow(barbieQ$assay)) {
    stop("Barcode count matrix extracted from 'object' has zero rows.")
  }
  if (anyNA(barbieQ$assay)) {
    warning("Barcode count matrix extracted from 'object' includes NAs.")
  }
  ## extract metadata from `target`
  ## if `target` not specified, inherit metadata from the passed `barbieQ` object
  if (is.null(target)) target <- barbieQ$metadata
  ## if 'barbieQ$metadata' is NULL either, default a homogeneous setting in the metadata.
  if (is.null(target)) {
    target <- as.data.frame(matrix(1, ncol(barbieQ$assay), 1))
    warning("no valid 'target' provided. Assigning a homogeneous condition to all samples.")
  } else {
    ## with valid `target` provided by either @param `target` or `object$metadata`, check the structure
    if ((!(is.list(target)) && is.vector(target)) || is.factor(target)) {
      target <- data.frame(target)
    } else if (!(inherits(target, "data.frame") || inherits(target, "matrix"))) {
      stop("'target' must be a vector, matrix or data.frame describing sample conditions.")
    }
  }

  ## check if sample size match between 'barbieQ$metadata' and 'barbieQ$assay'
  if (identical(ncol(barbieQ$assay), nrow(target))) {
    barbieQ$metadata <- as.data.frame(target)
  } else {
    ## if sample size not matching, set a homogeneous setting in metadata.
    barbieQ$metadata <- as.data.frame(matrix("1", ncol(barbieQ$assay), 1))
    warning("sample size mismatch between conditions and barcode counts! Assigning a homogeneous condition to all samples.")
  }
  if (anyNA(barbieQ$metadata)) warning("sample conditions includes NAs.")

  ## extract factor color list from @param 'factorColors'
  ## if @param factorColors is not provided, inherit factorColors from the 'barbieQ' object created.
  if (is.null(factorColors)) factorColors <- barbieQ$factorColors
  ## if 'barbieQ$factorColors' is NULL either, remind users and continue.
  if (is.null(factorColors)) {
    message("continuing with missing factorColors.")
  } else {
    barbieQ$factorColors <- factorColors
  }

  ## now 'barbieQ$assay' should be a matrix already, check if it's numeric apart from NAs.
  ## find NAs in the 'barbieQ$assay'
  ColumnHasNA <- vapply(
    as.data.frame(barbieQ$assay), 
    function(x) any(is.na(x)), 
    FUN.VALUE = logical(1L)
    )
  WhichColumnHasNA <- which(ColumnHasNA)
  RowHasNA <- vapply(
    seq_len(nrow(barbieQ$assay)),
    function(i) any(is.na(barbieQ$assay[i, ])),
    FUN.VALUE = logical(1L)
  )
  WhichRowHasNA <- which(RowHasNA)
  if (any(ColumnHasNA)) {
    warning("Barcode count matrix has NA, across ", length(WhichRowHasNA), " Barcodes (rows), and ", length(WhichColumnHasNA), " samples (columns).")
    ## display a menu for the user to choose how to handle NAs.
    choice <- utils::menu(c("yes, continue with zeros", "no, keep NA and continue", "stop"), title = "Do you want to convert all NA into zero and continue?")
    if (choice == 1) {
      ## replace NAs by zeros.
      barbieQ$assay[is.na(barbieQ$assay)] <- 0
      message("all NAs in Barcode count table are converted into zero.")
    } else if (choice == 2) {
      ## proceed without preprocessing calculations
      message("NAs are retained. barbieQ object is created without auto-preprocessing. Note that NAs may cause issues in subsequent analyses.")
    } else {
      stop("please address the NAs in the barcode count matrix first.")
    }
  }

  ## dispatch 'returnNumMat' function to ensure 'barbieQ$assay' is a numeric matrix apart from NAs.
  barbieQ$assay <- as.data.frame(returnNumMat(barbieQ$assay))

  ## preprocessing.
  ## calculate Barcode proportion, CPM, occurrence, rank ...
  ## set na.rm=TRUE in colSums to handle NA values.
  libSize <- colSums(barbieQ$assay, na.rm = TRUE)
  ## avoid division by zero by setting zero columns to NA.
  libSize[libSize == 0] <- NA

  barbieQ$proportion <- (t(barbieQ$assay) / libSize) |> t()
  barbieQ$CPM <- barbieQ$proportion * 1e6
  ## determine Barcode occurrence by default threshold.
  barbieQ$occurrence <- barbieQ$CPM >= 1
  ## calculate the rank of Barcodes in each column
  ## NAs will be placed last in the order.
  ## setting na.last="keep" will retain NAs as NA instead of assigning a rank.
  ## setting ties.method="first" will handle equal values by assigning the first the lowest rank.
  ## setting 'minus x' will rank 'x' by decreasing value.
  barbieQ$rank <- apply(barbieQ$assay, 2, function(x) {
    rank(-x, ties.method = "first", na.last = "keep")
  })

  ## dispatch function to decide if each Barcode is "top Barcode" across the samples.
  ## barbieQ$isTop <- getTopBar(barbieQ = barbieQ)

  return(barbieQ)
}
