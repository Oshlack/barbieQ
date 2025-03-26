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
#' @param sampleMetadata A `matrix` or `data.frame` where rows represent samples 
#'  and columns represent factors, each containing the specific condition assigned to a sample.
#'  If not specified, all samples are assigned the same condition.
#'  If a `barbieQ` object is passed to `object`, the sample conditions
#'  will be updated by specifying `sampleMetadata`.
#' @param factorColors A `list` mapping factors in `sampleMetadata` to color palettes, 
#'  where each element is a named `vector` assigning colors to specific sample conditions.
#'  If not specified, defaults to `NULL`.
#'  If a `barbieQ` object is passed to `object`, the mapping colors
#'  will be updated by specifying `factorColors`.
#'
#' @return A `barbieQ` is a `SummarizedExperimnt` object, including several components:
#'  * `assay`:A `data.frame` in `assays` storing Barcode counts, with Barcodes in rows and samples in columns.
#'  * `sampleMetadata`: A `DataFrame` in `colData` storing sample conditions, with samples in rows and
#'    experimental factors in columns.
#'  * `factorColors`: A `list` mapping factors in `sampleMetadata` to color palettes,
#'    where each element is a named `vector` assigning colors to specific sample conditions.
#'  * Other matrices in `assays`: `proportion`, `CPM`, `occurrence`, and `rank` representing
#'    Barcode proportion, Counts Per Million (CPM), presence/absence, and ranking
#'    in each sample, based on the Barcode counts in `assay`.
#'  * `isTopBarcode`: A `DataFrame` in `elementMetadata` (`rowData`) containing a single logical column
#'    indicating whether each Barcode is a major contributor. It is initially set to `TRUE` for all Barcodes, 
#'    but should be refined by calling [tagTopBarcodes].
#'  * Additional processed information during further analysis.
#'
#' @export
#'
#' @importFrom utils menu
#' @importClassesFrom S4Vectors DataFrame
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment assays
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment rowData
#' @importFrom S4Vectors metadata
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
#' ## Passing `object`, `sampleMetadata` and `factorColors`
#' createBarbieQ(object = barcodeCount)
#' createBarbieQ(object = barcodeCount, sampleMetadata = sampleConditions)
#' createBarbieQ(
#'   object = barcodeCount, sampleMetadata = sampleConditions, factorColors = conditionColor
#' )
#'
#' ## Updating a `barbieQ` object by passing new `sampleMetadata` and `factorColors`
#' myBarbieQ <- createBarbieQ(barcodeCount, sampleConditions, conditionColor)
#' createBarbieQ(
#'   object = myBarbieQ, sampleMetadata = data.frame(Mouse = rep(seq_len(4), each = 3))
#' )
#' createBarbieQ(
#'   object = myBarbieQ, sampleMetadata = data.frame(Mouse = rep(seq_len(4), each = 3)),
#'   factorColors = list(
#'     Mouse = c("1" = "#111199", "2" = "#112200", "3" = "#441111", "4" = "#000000")
#'   )
#' )
createBarbieQ <- function(object, sampleMetadata = NULL, factorColors = NULL) {
  ## barbieQ object relies on S4 class SummarizedExperiment object
  ## define the components to be included
    ## `assays`:
      ## `assay` for raw counts
      ## `proportion`, `CPM`, `occurrence`, `rank`, `isTopArray`
    ## `colData`:
      ## `sampleMetadata`
    ## `elementMetadata` (accessor `rowData`):
      ## `isTopBarcode`
    ## `metadata`:
      ## `factorColorts` storing color palettes corresponding to the factors in `sampleMetadata`

  ## extract Barcode counts from `object`
  ## dispatch `returnNumMat` function to ensure the object is a numeric matrix apart from NAs.
  ## taking `assays(SE)$assay` or `assay(SE)` for `object` is considered as `matrix` or `data.frame`
  if (inherits(object, "data.frame") ||
    inherits(object, "matrix") ||
    is.vector(object) && !(is.list(object))) {
    ## initiating barbieQ as a SE S4 class object
    ## basic components for now: `assays$assay`, `colData` as NULL, `metadata`as list(), `elementMetadata` as `DataFrame` with 0 columns, and `NAMES` saving barcode names
    if(is.null(colnames(object))) {
      colnames(object) <- paste0("sample", seq_len(ncol(object)))
    }
    if(is.null(rownames(object))) {
      rownames(object) <- paste0("barcode", seq_len(nrow(object)))
    }
    barbieQ <- SummarizedExperiment::SummarizedExperiment(
      ## the matrix has to have colnames and rownames - otherwise not allowed in SE
      assays = list(assay = returnNumMat(object))
    )
  } else if (inherits(object, "SummarizedExperiment")) {
    ## if a SE object is passed to `object` - no need to check dimensions
    ## inherit colData, rowData, metadata 
    ## for assays in the object, only pass `assay` - the raw count
    ## other slots in assays should be recalculated, such as CPM, etc...
    barbieQ <- SummarizedExperiment::SummarizedExperiment(
      assays = list(assay = SummarizedExperiment::assay(object)),
      colData = SummarizedExperiment::colData(object),
      rowData = SummarizedExperiment::rowData(object),
      metadata = S4Vectors::metadata(object)
    )
    message("inheriting `assay`, `colData`, `roweData` and `metadata`. Other slots including `proportion`, `CPM`, `occurrence` and `rank` will be recalculated based on `assay`.")
  }
  ## no need to check row dimension - even passing an empty matrix, SE initiates a matrix of [1,1] with NA
  ## check NAs in `assay`
  if (anyNA(SummarizedExperiment::assay(barbieQ))) {
    warning("Barcode count matrix extracted from 'object' includes NAs.")
  }
  ## extract metadata from `sampleMetadata`
  ## if `sampleMetadata` not specified, inherit metadata from the passed `barbieQ` object
  ## `colData(barbieQ)$sampleMetadata` returns NULL if there's no such thing in `barbieQ` object
  if (is.null(sampleMetadata)) {
    sampleMetadata <- SummarizedExperiment::colData(barbieQ)$sampleMetadata
  }
  ## if `sampleMetadata` not provided by `barbieQ` either, default a homogeneous setting for `sampleMetadata`.
  if (is.null(sampleMetadata)) {
    sampleMetadata <- S4Vectors::DataFrame(matrix(1, ncol(barbieQ), 1))
    warning("no valid 'sampleMetadata' provided. Assigning a homogeneous condition to all samples.")
  } else {
    ## for `sampleMetadata` specified either way, check the format
    if ((!(is.list(sampleMetadata)) && is.vector(sampleMetadata)) 
        || is.factor(sampleMetadata)) {
      sampleMetadata <- S4Vectors::DataFrame(sampleMetadata)
    } else if (!(inherits(sampleMetadata, "data.frame") || 
                 inherits(sampleMetadata, "matrix") || 
                 is(sampleMetadata, "DataFrame"))) {
      stop("'sampleMetadata' must be a vector, matrix, data.frame or DataFrame describing sample conditions.")
    }
  }
  
  ## check if sample size match between `sampleMetadata` and `assay`
  if (identical(ncol(barbieQ), nrow(sampleMetadata))) {
    ## pass `sampleMetadata` to `barbieQ` if sample names are in line
    ## check if sample names in line
    if (is.null(rownames(SummarizedExperiment::colData(barbieQ)))) {
      if (is.null(rownames(S4Vectors::DataFrame(sampleMetadata)))) {
        SummarizedExperiment::colData(barbieQ)$sampleMetadata <- S4Vectors::DataFrame(sampleMetadata)
        warning("sample names unspecified in either `object` or `sampleMetadata`.")
      } else {
        SummarizedExperiment::colData(barbieQ)$sampleMetadata <- S4Vectors::DataFrame(sampleMetadata)
        message("sample names set up in line with `sampleMetadata`")
      }
    } else {
      if (is.null(rownames(S4Vectors::DataFrame(sampleMetadata)))) {
        SummarizedExperiment::colData(barbieQ)$sampleMetadata <- S4Vectors::DataFrame(sampleMetadata)
        message("sample names set up in line with `object`")
      } else if (identical(rownames(colData(barbieQ)), 
                           rownames(S4Vectors::DataFrame(sampleMetadata)))) {
        ## if sample names align between `assay` and `sampleMetadata`
        ## assign `sampleMetadata` to `colData$sampleMetadata`
        SummarizedExperiment::colData(barbieQ)$sampleMetadata <- S4Vectors::DataFrame(sampleMetadata)
      } else if (setequal(rownames(colData(barbieQ)), 
                          rownames(S4Vectors::DataFrame(sampleMetadata)))) {
        ## if sample names overlap but in a different order
        SummarizedExperiment::colData(barbieQ)$sampleMetadata <- S4Vectors::DataFrame(sampleMetadata)[rownames(colData(barbieQ)),]
        message("sample order following `object`; `sampleMetadata` reordered.")
      } else {
        ## sample names don't overlap, assign `assay` sample names to `sampleMetadata`
        rownames(S4Vectors::DataFrame(sampleMetadata)) <- rownames(colData())
        SummarizedExperiment::colData(barbieQ)$sampleMetadata <- S4Vectors::DataFrame(sampleMetadata)
        warning("sample names don't match! Assigning sample names from `object` to `sampleMetadata`")
      }
    }
  } else {
    ## if sample size not matching, set a homogeneous setting in metadata.
    SummarizedExperiment::colData(barbieQ)$sampleMetadata <- S4Vectors::DataFrame(matrix("1", ncol(barbieQ), 1))
    warning("sample size don't match! Assigning a homogeneous condition to all samples.")
  }
  ## check if any NAs exist in sampleMetadata
  if (anyNA(SummarizedExperiment::colData(barbieQ)$sampleMetadata)) {
    warning("sample conditions includes NAs.")
  }
  
  ## extract factor color list from @param `factorColors`
  ## if `factorColors` is not provided, inherit factorColors from the `barbieQ` object
  if (is.null(factorColors)) {
    factorColors <- S4Vectors::metadata(barbieQ)$factorColors
  }
  ## if not provided by `object`, remind users and continue.
  if (is.null(factorColors)) {
    message("continuing with missing `factorColors`.")
  } else {
    S4Vectors::metadata(barbieQ)$factorColors <- factorColors
  }

  ## now `assay(barbieQ)` should be a matrix already, check if it's numeric apart from NAs
  ## find NAs in the 'assay(barbieQ)'
  ColumnHasNA <- vapply(
    as.data.frame(SummarizedExperiment::assay(barbieQ)), 
    function(x) any(is.na(x)), 
    FUN.VALUE = logical(1L)
    )
  WhichColumnHasNA <- which(ColumnHasNA)
  RowHasNA <- vapply(
    seq_len(nrow(barbieQ)),
    function(i) any(is.na(SummarizedExperiment::assay(barbieQ)[i, ])),
    FUN.VALUE = logical(1L)
  )
  WhichRowHasNA <- which(RowHasNA)
  if (any(ColumnHasNA)) {
    warning("Barcode count matrix has NA, across ", length(WhichRowHasNA), " Barcodes (rows), and ", length(WhichColumnHasNA), " samples (columns).")
    ## display a menu for the user to choose how to handle NAs.
    choice <- utils::menu(c("yes, continue with zeros", "no, keep NA and continue", "stop"), title = "Do you want to convert all NA into zero and continue?")
    if (choice == 1) {
      ## replace NAs by zeros.
      SummarizedExperiment::assay(barbieQ)[is.na(SummarizedExperiment::assay(barbieQ))] <- 0
      message("all NAs in Barcode count table are converted into zero.")
    } else if (choice == 2) {
      ## proceed without preprocessing calculations
      message("NAs are retained. barbieQ object is created without auto-preprocessing. Note that NAs may cause issues in subsequent analyses.")
    } else {
      stop("please address the NAs in the barcode count matrix first.")
    }
  }

  ## dispatch 'returnNumMat' function to ensure 'assay(barbieQ)' is a numeric matrix apart from NAs.
  SummarizedExperiment::assay(barbieQ) <- as.data.frame(returnNumMat(SummarizedExperiment::assay(barbieQ)))

  ## data transformation: processing other slots to be included into assays.
  ## calculate Barcode proportion, CPM, occurrence, rank ...
  ## set na.rm=TRUE in colSums to handle NA values.
  libSize <- colSums(SummarizedExperiment::assay(barbieQ), na.rm = TRUE)
  ## avoid division by zero by setting zero columns to NA.
  libSize[libSize == 0] <- NA
  
  ## store transformed slots into `assays`
  SummarizedExperiment::assays(barbieQ)$proportion <- (t(SummarizedExperiment::assay(barbieQ)) / libSize) |> t()
  SummarizedExperiment::assays(barbieQ)$CPM <- SummarizedExperiment::assays(barbieQ)$proportion * 1e6
  ## determine Barcode occurrence by default threshold.
  SummarizedExperiment::assays(barbieQ)$occurrence <- SummarizedExperiment::assays(barbieQ)$CPM >= 1
  ## calculate the rank of Barcodes in each column
  ## NAs will be placed last in the order.
  ## setting na.last="keep" will retain NAs as NA instead of assigning a rank.
  ## setting ties.method="first" will handle equal values by assigning the first the lowest rank.
  ## setting 'minus x' will rank 'x' by decreasing value.
  SummarizedExperiment::assays(barbieQ)$rank <- apply(SummarizedExperiment::assay(barbieQ), 2, function(x) {
    rank(-x, ties.method = "first", na.last = "keep")
  })

  ## initiating `isTopBarcode` by seting all barcodes as "top"
  SummarizedExperiment::rowData(barbieQ)$isTopBarcode <- S4Vectors::DataFrame(matrix(rep(TRUE, nrow(barbieQ)), nrow(barbieQ), 1))
  ## users should dispatch proper function to update `isTopAssay` and `isTopBarcode`

  return(barbieQ)
}
