#' extracting a numeric matrix out of the input object
#'
#' @param object a matrix, a data.frame, or a vector
#'
#' @return a numeric matrix with NAs retained
#'
#' @examples
#' mat <- data.frame(id = letters[1:5], matrix(1:25,5,5))
#' returnNumMat(mat)
returnNumMat <- function(object) {
  ## if the object is a vector, treating it as a one-column matrix
  if(is.vector(object)){
    object <- as.matrix(object, ncol=1L)
    message("object is a vecotr, converted it into a one-column matrix.")
  }
  ## check if the object is a matrix or data.frame.
  if(inherits(object, "data.frame") || inherits(object, "matrix")) {
    ## check if all columns are numeric, excluding NAs
    ColumnIsNumeric <- vapply(
      as.data.frame(object),
      function(x) is.numeric(na.omit(x)), FUN.VALUE=logical(1L))
    if(all(ColumnIsNumeric))
      ## update object as a numeric matrix
      objectUpdated <- as.matrix(object)
    else {
      ## find the non-numeric column
      WhichNotNumeric <- which(!ColumnIsNumeric)
      ## if the first column is the only non-numeric column, treating it as row IDs
      if(identical(sum(WhichNotNumeric), 1L) && length(ColumnIsNumeric) > 1L) {
        ## only update rownames when it's NULL
        if(is.null(rownames(object)) || all(rownames(object) == seq(nrow(object)))) {
          message("attempting to convert first column as Barcode IDs.")
          ## extract first non-numeric column
          firstCol <- object[,1,drop = T]
          ## if NAs exist in the first column, name the NAs
          if(any(is.na(firstCol))) {
            na_indices <- which(is.na(firstCol))
            firstCol[na_indices] <- paste0("Add_Barcode_", seq_along(na_indices))
            message("NAs exist in Barcode IDs, replaced by new names.")
          }
          ## if first column values (Barcode IDs) are unique, make them unique
          if(any(duplicated(firstCol))) {
            firstCol <- make.unique(as.character(firstCol))
            message("duplicated Barcode IDs exist, making them unique.")
          }
          ## update object
          objectUpdated <- as.matrix(object[,-1,drop=FALSE])
          rownames(objectUpdated) <- firstCol
          message("treating first column as barcode IDs; column number -1")
        } else {
          stop("cannot rename Barcode IDs by first column, as object rownames already exist.")
        }
      } else {
        stop(
          "object should be numeric, instead it is a data.frame (or matrix) with ",
          length(WhichNotNumeric), " non-numeric colunms.")
      }
    }
  }
  return(objectUpdated)
}

