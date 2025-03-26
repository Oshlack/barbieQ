#' extracting a numeric matrix out of the input object
#'
#' @param object a matrix, a data.frame, or a vector
#'
#' @return a numeric matrix with NAs retained
#'
#' @noRd
#'
#' @examples \donttest{
#' mat <- data.frame(id = letters[1:5], matrix(seq_len(25), 5, 5))
#' barbieQ:::returnNumMat(mat)
#' }
returnNumMat <- function(object) {
    ## if the object is a vector, treating it as a one-column matrix
    if (is.vector(object)) {
        object <- as.matrix(object, ncol = 1L)
        message("object is a vecotr, converted it into a one-column matrix.")
    }
    ## check if the object is a matrix or data.frame.
    if (inherits(object, "data.frame") || inherits(object, "matrix")) {
        ## check if all columns are numeric, excluding NAs
        ColumnIsNumeric <- vapply(as.data.frame(object), function(x) is.numeric(na.omit(x)),
            FUN.VALUE = logical(1L))
        if (all(ColumnIsNumeric)) {
            ## update object as a numeric matrix
            objectUpdated <- as.matrix(object)
        } else {
            ## find the non-numeric column
            WhichNotNumeric <- which(!ColumnIsNumeric)
            ## if the first column is the only non-numeric column, treating it as row
            ## IDs
            if (identical(sum(WhichNotNumeric), 1L) && length(ColumnIsNumeric) > 1L) {
                ## only update rownames when it's NULL
                if (is.null(rownames(object)) || all(rownames(object) == seq(nrow(object)))) {
                  message("attempting to convert first column as Barcode IDs.")
                  ## extract first non-numeric column
                  firstCol <- object[, 1, drop = TRUE]
                  ## if NAs exist in the first column, name the NAs
                  if (any(is.na(firstCol))) {
                    na_indices <- which(is.na(firstCol))
                    firstCol[na_indices] <- paste0("Add_Barcode_", seq_along(na_indices))
                    message("NAs exist in Barcode IDs, replaced by new names.")
                  }
                  ## if first column values (Barcode IDs) are unique, make them unique
                  if (any(duplicated(firstCol))) {
                    firstCol <- make.unique(as.character(firstCol))
                    message("duplicated Barcode IDs exist, making them unique.")
                  }
                  ## update object
                  objectUpdated <- as.matrix(object[, -1, drop = FALSE])
                  rownames(objectUpdated) <- firstCol
                  message("treating first column as barcode IDs; column number -1")
                } else {
                  stop("cannot rename Barcode IDs by first column, as object rownames already exist.")
                }
            } else {
                stop("object should be numeric, instead it is a data.frame (or matrix) with ",
                  length(WhichNotNumeric), " non-numeric colunms.")
            }
        }
    }
    return(objectUpdated)
}

#' Check if `barbieQ` object is in right format
#'
#' `checkBarbieQDimensions()` ensures that the `proportion`, `CPM`, `occurrence`,
#'  and `rank` components of `barbieQ` have the same number of samples and
#'  Barcodes as `assay`.
#'
#' @param barbieQ A `barbieQ` object created by the [createBarbieQ] function.
#'
#' @return a logical value
#'  * Returns TRUE if components in the `barbieQ` object are correctly formatted.
#'  * Otherwise, the function throws an error with specific a reason.
#'
#' @noRd
#'
#' @examples \donttest{
#' Treat <- factor(rep(seq_len(2), each = 6))
#' Time <- rep(rep(seq_len(2), each = 3), 2)
#' nbarcodes <- 50
#' nsamples <- 12
#' count <- abs(matrix(rnorm(nbarcodes * nsamples), nbarcodes, nsamples))
#' rownames(count) <- paste0('Barcode', seq_len(nbarcodes))
#' barbieQ <- createBarbieQ(count, data.frame(Treat = Treat, Time = Time))
#' barbieQ:::checkBarbieQDimensions(barbieQ)
#' }
checkBarbieQDimensions <- function(barbieQ) {
    ## check barbieQ$assay's format
    if (is.data.frame(barbieQ$assay) || is.matrix(barbieQ$assay)) {
        NumDim <- dim(barbieQ$assay)
    } else {
        stop("barbieQ$assay wrong format.
call createBarbieQ() to generate proper `barbieQ` object - don't modify by hand..")
    }
    ## check barbieQ$metadata's format and dimension
    if (is.data.frame(barbieQ$metadata) || is.matrix(barbieQ$metadata)) {
        if (nrow(barbieQ$metadata) != NumDim[2]) {
            stop("row dimension of barbieQ$metadata doesn't match column dimention of barbieQ$assay.
call createBarbieQ() to generate proper `barbieQ` object - don't modify by hand..")
        }
    } else {
        stop("BarbieQ$assay wrong format.
call createBarbieQ() to generate proper `barbieQ` object - don't modify by hand..")
    }
    ## check other components' format dimensions
    elements <- c("proportion", "CPM", "occurrence", "rank")
    for (elem in elements) {
        if (is.data.frame(barbieQ[[elem]]) || is.matrix(barbieQ[[elem]])) {
            elemDim <- dim(barbieQ[[elem]])
            if (any(elemDim != NumDim)) {
                stop("dimensions of barbieQ component-", elem, "don't match dimentions of barbieQ$assay.
call createBarbieQ() to generate proper `barbieQ` object - don't modify by hand..")
            }
        } else {
            stop("barbieQ component-", elem, " isn't in right format.
call createBarbieQ() to generate proper `barbieQ` object - don't modify by hand..")
        }
    }
    ## check barbieQ$isTop format and dimensions
    if (!is.null(barbieQ$isTop$vec)) {
        if (is.vector(barbieQ$isTop$vec)) {
            if (length(barbieQ$isTop$vec) != NumDim[1]) {
                stop("length of barbieQ$isTop$vec doesn't match row dimention of barbieQ$assay.
call createBarbieQ() to generate proper `barbieQ` object - don't modify by hand..")
            }
        } else {
            stop("BarbieQ$isTop$vec isn't in right format.
call createBarbieQ() to generate proper `barbieQ` object - don't modify by hand..")
        }
    }
    if (!is.null(barbieQ$isTop$mat)) {
        if (is.matrix(barbieQ$isTop$mat) || is.data.frame(barbieQ$isTop$mat)) {
            if (nrow(barbieQ$isTop$mat) != NumDim[1]) {
                stop("row dimension barbieQ$isTop$mat doesn't match row dimention of barbieQ$assay.
call createBarbieQ() to generate proper `barbieQ` object - don't modify by hand.")
            }
        } else {
            stop("barbieQ$isTop$mat isn't in right format.
call createBarbieQ() to generate proper `barbieQ` object - don't modify by hand..")
        }
    }

    return(TRUE)
}

#' Extract sampleMetadata and primary factor
#'
#' `extractSampleMetadataAndPrimaryFactor()` extracts sampleMetadata `data.frame` from the
#' specified `sampleMetadata` (prioritised) and `sampleGroup`, or inherits from
#' `barbieQ`, and identifies the primary factor based on the `sampleGroup`.
#'
#' @param barbieQ A `barbieQ` object based on `SummarizedExperiment` created by the [createBarbieQ] function.
#' @param sampleMetadata A `matrix`, `data.frame` or `DataFrame` of sample conditions,
#'  where each factor is represented in a separate column. Defaults to NULL,
#'  in which case sample conditions are inherited from `colData(barbieQ)$sampleMetadata`.
#' @param sampleGroup A string representing the name of a factor from the
#'  sample conditions passed by `barbieQ` or `sampleMetadata`, or a `vector` of
#'  sample conditions, indicating the primary factor to be evaluated.
#'  Defaults to the first factor in the sample conditions, or a `list` with a vector like this,
#'  or a `data.frame` or `DataFrame` with a single column of sample conditions.
#'
#' @return A `DataFrame` object:
#'  * `listData` with samples in rows and conditions in columns
#'  * `medadata$primaryFactor` storing the name of the primary factor
#'
#' @noRd
#' 
#' @importClassesFrom S4Vectors DataFrame
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment colData
#' @importFrom S4Vectors metadata
#'
#' @examples \donttest{
#' Treat <- factor(rep(seq_len(2), each = 6))
#' Time <- rep(rep(seq_len(2), each = 3), 2)
#' nbarcodes <- 50
#' nsamples <- 12
#' count <- abs(matrix(rnorm(nbarcodes * nsamples), nbarcodes, nsamples))
#' rownames(count) <- paste0('Barcode', seq_len(nbarcodes))
#' barbieQ <- createBarbieQ(count, data.frame(Treat = Treat, Time = Time))
#' barbieQ:::extarctSampleMetadataAndPrimaryFactor(
#'   barbieQ = barbieQ, sampleGroup = 'Treat'
#' )
#' }
extractSampleMetadataAndPrimaryFactor <- function(barbieQ, sampleMetadata = NULL, sampleGroup = NULL) {
    ## check sampleMetadata: if 'sampleMetadata' is not specified, inherit from
    ## barbieQ (could still be NULL)
    if (is.null(sampleMetadata)) {
        sampleMetadata <- SummarizedExperiment::colData(barbieQ)$sampleMetadata
    }
    ## case when sampleMetadata is specified or provided by barbieQ in right format
    if (is.vector(sampleMetadata) || is.factor(sampleMetadata)) {
        sampleMetadata <- S4Vectors::DataFrame(V1 = sampleMetadata)
        ## check dimensions - sample size
        if (nrow(sampleMetadata) != ncol(barbieQ)) {
            stop("sample size don't match between `sampleMetadata` and `barbieQ`.")
        }
    } else if (is.matrix(sampleMetadata) || is.data.frame(sampleMetadata) || is(sampleMetadata,
        "DataFrame")) {
        sampleMetadata <- S4Vectors::DataFrame(sampleMetadata)
        ## check dimensions - sample size
        if (nrow(sampleMetadata) != ncol(barbieQ)) {
            stop("sample size don't match between `sampleMetadata` and `barbieQ`.")
        }
    } else {
        ## case when sampleMetadata is not specified by `object` or `sampleMetadata`
        ## or not any of the above formats when group is a vector, or factor with
        ## length of sample size assign it to sampleMetadata
        if (((is.vector(sampleGroup) && !is.list(sampleGroup)) || is.factor(sampleGroup)) &&
            length(sampleGroup) > 1L) {
            ## check sampleGroup length
            if (length(sampleGroup) != ncol(barbieQ)) {
                stop("sample size don't match between `sampleGroup` and `barbieQ`.")
            } else {
                sampleMetadata <- S4Vectors::DataFrame(sampleGroup = sampleGroup)
                ## as `sampleGroup` will be the only columnn in `sampleMetadata` set
                ## it up as primary factor and store in metadata of `sampleMetadata`
                S4Vectors::metadata(sampleMetadata)$primaryFactor <- "sampleGroup"
                message("assigning `sampleGroup` to `sampleMetadata`, set up as primary factor.")
            }
        } else if ((is.list(sampleGroup) && length(sampleGroup) == 1L) || (is(sampleGroup,
            "DataFrame") && ncol(sampleGroup) == 1L)) {
            ## case when `sampleGroup` is a list with a named vector of proper length
            ## or data.frame or DataFrame of a single column of proper rows
            if (length(sampleGroup[[1]]) != ncol(barbieQ)) {
                stop("sample size don't match between `sampleGroup` and `barbieQ`.")
            } else {
                ## assign `sampleGroup` to `sampleMetadata`
                sampleMetadata <- S4Vectors::DataFrame(sampleGroup)
                S4Vectors::metadata(sampleMetadata)$primaryFactor <- names(sampleGroup)[1]
                message("assigning ", names(sampleGroup)[1], " to `sampleMetadata`, set up as primary factor.")
            }
        } else {
            ## case when `sampleMetadata` is still NULL or not properly specified
            stop("no valid `sampleMetadata` or `sampleGroup` specified.")
        }
    }

    ## now sampleMetadata should be a matrix, data.frame or DataFrame already this
    ## section if to set up the primary factor and save it to sampleGroup if
    ## unspecified, default sampleGroup to the first factor in sampleMetadata
    if (is.null(sampleGroup)) {
        S4Vectors::metadata(sampleMetadata)$primaryFactor <- colnames(sampleMetadata)[1]
        message("setting ", colnames(sampleMetadata)[1], " as the primary factor in `sampleMetadata`.")
    }
    ## if `sampleGroup` is a specified column name of the `sampleMetadata`, extract
    ## the entire vector
    if (is.character(sampleGroup) && length(sampleGroup) == 1L) {
        ## case when `sampleGroup` is the name of a factor in `sampleMetadata`
        if (sampleGroup %in% colnames(sampleMetadata)) {
            S4Vectors::metadata(sampleMetadata)$primaryFactor <- sampleGroup
            message("setting ", sampleGroup, " as the primary factor in `sampleMetadata`.")
        } else {
            stop("cannot find `sampleGroup` in factor names in `sampleMetadata`.")
        }
    } else if (((is.vector(sampleGroup) && !is.list(sampleGroup)) || is.factor(sampleGroup)) &&
        length(sampleGroup) > 1L) {
        ## case when `sampleGroup` is specified with the existence of proper
        ## `sampleMetadata`
        if (length(sampleGroup) != ncol(barbieQ)) {
            stop("sample size don't match between `sampleGroup` and `barbieQ`.")
        } else {
            ## cbind `sampleGroup` to `sampleMetadata`
            sampleMetadata <- cbind(sampleGroup, sampleMetadata)
            S4Vectors::metadata(sampleMetadata)$primaryFactor <- "sampleGroup"
            message("binding `sampleGroup` to `sampleMetadata`, set up as primary factor in `sampleMetadata`.")
        }
    } else if ((is.list(sampleGroup) && length(sampleGroup) == 1L) || (is(sampleGroup, "DataFrame") &&
        ncol(sampleGroup) == 1L)) {
        ## case when `sampleGroup` is a list with a named vector of proper length or
        ## data.frame or DataFrame of a single column of proper rows
        if (length(sampleGroup[[1]]) != ncol(barbieQ)) {
            stop("sample size don't match between `sampleGroup` and `barbieQ`.")
        } else {
            ## cbind `sampleGroup` to `sampleMetadata`
            sampleMetadata <- S4Vectors::DataFrame(cbind(sampleGroup, sampleMetadata))
            S4Vectors::metadata(sampleMetadata)$primaryFactor <- names(sampleGroup)[1]
            message("binding ", names(sampleGroup)[1], " to `sampleMetadata`, set up as primary factor in `sampleMetadata`.")
        }
    } else if (!is.null(sampleGroup)) {
        ## case when `sampleGroup` is not properly specified
        S4Vectors::metadata(sampleMetadata)$primaryFactor <- colnames(sampleMetadata)[1]
        message("setting ", colnames(sampleMetadata)[1], " as the primary factor in `sampleMetadata`.")
    }

    return(sampleMetadata)
}
