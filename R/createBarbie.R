#' create an object-centred structure-a 'Barbie' object-based on the input Barcode count and other information.
#'
#'\code{createBarbie}
#' @param object a numeric matrix or equivalent; or a Barbie object or equivalent
#' @param target a matrix of sample conditions
#' @param factorColors a list containing colour palettes corresponding to the factors in the target parameter
#'
#' @return a Barbie object
#' @export
#'
#' @examples
#' createBarbie(mymat)
#' createBarbie(mymat, mytarget)
#' createBarbie(mymat, mytarget, mycolors)
#'
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
  ## extract Barcode counts from @param object
  ## dispatch 'returnNumMat' function to ensure the object is a numeric matrix apart from NAs.
  if(inherits(object, "data.frame") || inherits(object, "matrix") || is.vector(object))
    Barbie$assay <- returnNumMat(object)
  else{
    ###$
    # extract barcode counts from objects. Need a function for it.
    # if Seurat, extract assay and metadata.
    # if BARtab object, extract ...
    # if Matrix object, ...
    # if Barbie itself, ...[get target and factorColors]
  }
  if(!nrow(Barbie$assay)) stop("Barcode count matrix extracted from 'object' has zero rows.")
  if(anyNA(Barbie$assay)) warning("Barcode count matrix extracted from 'object' includes NAs.")

  ## extract metadata from @param target
  ## if @param target is not provided, inherit metadata from the 'Barbie' object created.
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
    choice <- menu(c("yes, continue with zeros", "no, keep NA and continue", "stop"), title = "Do you want to convert all NA into zero and continue?")
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
