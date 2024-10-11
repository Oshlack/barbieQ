#' Extract Barcode count data and build into a `Barbie` object for
#' subsequent analysis
#'
#' `createbarbie()` constructs an object-centred list called `Barbie` object,
#' which is designed to process Barcode count data gained from
#' cell clonal tracking experiments.
#'
#' @param object A numeric matrix of Barcode counts, with Barcodes in rows
#'  and samples in columns;
#'  Alternatively, you can pass an existing `Barbie` object, from which
#'  barcode counts, sample conditions, and color palettes will be inherited.
#' @param target A `matrix` or `data.frame` of sample conditions,
#'  with each factor in a separate column.
#'  If a `Barbie` object is passed to `object`, the sample conditions
#'  in `Barbie` will be updated by `target`.
#'  If not specified, all samples are assigned the same condition.
#' @param factorColors A `list` of colour palettes corresponding to
#'  the factors in `target`.
#'  If a `Barbie` object is passed to `object`, the color palettes
#'  in `Barbie`will be updated by `factorColors`.
#'  If not specified, defaults to NULL.
#'
#' @return A `Barbie` object, including several components:
#'  * `assay`: a `data.frame` containing Barcode counts across samples.
#'  * `metadata`: a `data.frame` representing sample conditions,
#'    organized by various experimental factors.
#'  * `factorColors`: a `list` of color palettes corresponding to
#'    the sample conditions.
#'  * Matrices of `proportion`, `CPM`, `occurrence`, and `rank` representing
#'    the proportion, Counts Per Million (CPM), presence or absence, and ranking
#'    of Barcodes in each sample, based on the Barcode counts in `assay`.
#'  * `isTop`: a list containing a vector `vec` and a matrix `mat`, tagging
#'    each Barcode as being part of the major contributors or not,
#'    by calling the function [tagTopBarcodes].
#'  * Additional processed information for further analysis.
#'
#' @export
#'
#' @importFrom utils menu
#'
#' @examples
#' ## Sample conditions and color palettes
#' sampleConditions <- data.frame(
#'   Treat=factor(rep(c("ctrl", "drug"), each=6)),
#'   Time=rep(rep(1:2, each=3), 2))
#' conditionColor <- list(
#'   Treat=c(ctrl="#999999", drug="#112233"),
#'   Time=c("1"="#778899", "2"="#998877"))
#'
#' ## Barcode count data
#' nbarcodes <- 50
#' nsamples <- 12
#' barcodeCount <- abs(matrix(rnorm(nbarcodes*nsamples), nbarcodes, nsamples))
#' rownames(barcodeCount) <- paste0("Barcode", 1:nbarcodes)
#'
#' ## Passing `object`, `target` and `factorColors`
#' createBarbie(object=barcodeCount)
#' createBarbie(object=barcodeCount, target=sampleConditions)
#' createBarbie(
#'   object=barcodeCount, target=sampleConditions, factorColors=conditionColor)
#'
#' ## Updating a `Barbie` object by passing new `target` and `factorColors`
#' myBarbie <- createBarbie(barcodeCount, sampleConditions, conditionColor)
#' createBarbie(object=myBarbie, target=data.frame(Mouse=rep(1:4, each=3)))
#' createBarbie(
#'   object=myBarbie, target=data.frame(Mouse=rep(1:4, each=3)),
#'   factorColors=list(
#'     Mouse=c("1"="#111199", "2"="#112200", "3"="#441111", "4"="#000000")))
createBarbie <- function(object, target=NULL, factorColors=NULL) {
  ## define the Barbie object structure
  Barbie <- list(
    ## assay stores Barcode counts.
    assay=NULL,
    ## metadata stores experimental design of each sample, as the equivalent to 'target file' in single cell data.
    metadata=NULL,
    proportion=matrix(),
    CPM=matrix(),
    occurrence=matrix(),
    rank=matrix(),
    isTop=list(),
    clusters=factor(levels=character(), labels=character()),
    ## factorColors stores color palettes corresponding to the factors in 'metadata'.
    factorColors=NULL
  )

  ## extract Barcode counts from `object`
  ## dispatch `returnNumMat` function to ensure the object is a numeric matrix apart from NAs.
  if (inherits(object, "data.frame") ||
      inherits(object, "matrix") ||
      is.vector(object) && !(is.list(object))) {
    Barbie$assay <- returnNumMat(object)
    } else if (is.list(object)) {
      ## inherit components if a `Barbie` object is passed to `object`
      checkBarbieDimensions(object)
      Barbie <- object
    ## extract barcode counts from other objects. Need a function for it.
    # if Seurat, extract assay and metadata.
    # if BARtab object, extract ...
  }
  if (!nrow(Barbie$assay)) {
    stop("Barcode count matrix extracted from 'object' has zero rows.")
  }
  if (anyNA(Barbie$assay)) {
    warning("Barcode count matrix extracted from 'object' includes NAs.")
  }
  ## extract metadata from `target`
  ## if `target` not specified, inherit metadata from the passed `Barbie` object
  if(is.null(target)) target <- Barbie$metadata
  ## if 'Barbie$metadata' is NULL either, assume a homogeneous setting in the metadata.
  if(is.null(target)) {
    target <- as.data.frame(matrix(1, ncol(Barbie$assay), 1))
    message("no 'target' provided.
now creating a pseudo uni-group with a homogeneous setting in 'Barbie$metadata'.")
  } else {
    ## now target is provided by either @param target or Barbie$metadata, check target format
    if(is.vector(target) || is.factor(target)) target <- matrix(target, ncol = 1)
    else if(!(inherits(target, "data.frame") || inherits(target, "matrix")))
      stop("'target' should be a vector, matrix or data.frame describing sample conditions.")
  }

  ## check if the sample sizes match between 'Barbie$metadata' and 'Barbie$assay'
  if(identical(ncol(Barbie$assay), nrow(target)))
    Barbie$metadata <- as.data.frame(target)
  else {
    ## if sample sizes don't match, create a homogeneous setting in metadata.
    Barbie$metadata <- as.data.frame(matrix("1", ncol(Barbie$assay), 1))
    warning("sample size of 'target' (or 'Barbie$metadata') doesn't match sample size of Barcode count extracted in 'object'!
now creating a pseudo uni-group with a homogeneous setting in 'Barbie$metadata'.")
  }
  if(anyNA(Barbie$metadata)) warning("Barbie$metadata includes NAs.")

  ## extract factor color list from @param 'factorColors'
  ## if @param factorColors is not provided, inherit factorColors from the 'Barbie' object created.
  if(is.null(factorColors)) factorColors <- Barbie$factorColors
  ## if 'Barbie$factorColors' is NULL either, remind users and continue.
  if(is.null(factorColors)) message("continuing with missing factorColors.")
  else {
    Barbie$factorColors <- factorColors
  }

  ## now 'Barbie$assay' should be a matrix already, check if it's numeric apart from NAs.
  ## find NAs in the 'Barbie$assay'
  ColumnHasNA <- vapply(as.data.frame(Barbie$assay), function(x) any(is.na(x)), FUN.VALUE = logical(1L))
  WhichColumnHasNA <- which(ColumnHasNA)
  RowHasNA <- vapply(seq_len(nrow(Barbie$assay)),
                     function(i) any(is.na(Barbie$assay[i, ])),
                     FUN.VALUE = logical(1L))
  WhichRowHasNA <- which(RowHasNA)
  if(any(ColumnHasNA)) {
    warning("In Barcode count matrix (Barbie$assay), ", length(WhichRowHasNA), " Barcodes (rows) show NAs across ", length(WhichColumnHasNA), " samples (columns).")
    ## display a menu for the user to choose how to handle NAs.
    choice <- utils::menu(c("yes, continue with zeros", "no, keep NA and continue", "stop"), title = "Do you want to convert all NA into zero and continue?")
    if(choice == 1) {
      ## replace NAs by zeros.
      Barbie$assay[is.na(Barbie$assay)] <- 0
      message("all NAs in Barcode count table are converted into zero.")
    } else if(choice == 2) {
      ## proceed without preprocessing calculations
      message("NAs are retained. Barbie object is created without preprocessing.
note that NAs may cause issues in subsequent analyses.")
    } else {
      stop("please address the NAs in the 'object' parameter first.")
    }
  }

  ## dispatch 'returnNumMat' function to ensure 'Barbie$assay' is a numeric matrix apart from NAs.
  Barbie$assay <- as.data.frame(returnNumMat(Barbie$assay))

  ## preprocessing.
  ## calculate Barcode proportion, CPM, occurrence, rank ...
  ## set na.rm=TRUE in colSums to handle NA values.
  libSize <- colSums(Barbie$assay, na.rm=TRUE)
  ## avoid division by zero by setting zero columns to NA.
  libSize[libSize == 0] <- NA

  Barbie$proportion <- (t(Barbie$assay) / libSize) |> t()
  Barbie$CPM <- Barbie$proportion * 1e6
  ## determine Barcode occurrence by default threshold.
  Barbie$occurrence <- Barbie$CPM >= 1
  ## calculate the rank of Barcodes in each column
  ## NAs will be placed last in the order.
  ## setting na.last="keep" will retain NAs as NA instead of assigning a rank.
  ## setting ties.method="first" will handle equal values by assigning the first the lowest rank.
  ## setting 'minus x' will rank 'x' by decreasing value.
  Barbie$rank <- apply(Barbie$assay, 2, function(x) {
    rank(-x, ties.method="first", na.last="keep")})

  ## dispatch function to decide if each Barcode is "top Barcode" across the samples.
  ## Barbie$isTop <- getTopBar(Barbie = Barbie)

  return(Barbie)
}
