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
                  ## if first column values (Barcode IDs) are unique, make them
                  ## unique
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

#' Extract targets and primary factor
#'
#' `extractTargetsAndPrimaryFactor()` extracts targets `data.frame` from the
#' specified `targets` (prioritised) and `sampleGroups`, or inherits from
#' `barbieQ`, and identifies the primary factor based on the `sampleGroups`.
#'
#' @param barbieQ A `barbieQ` object created by the [createBarbieQ] function.
#' @param targets A `matrix` or `data.frame` of sample conditions,
#'  where each factor is represented in a separate column. Defaults to NULL,
#'  in which case sample conditions are inherited from `barbieQ$metadata`.
#' @param sampleGroups A string representing the name of a factor from the
#'  sample conditions passed by `barbieQ` or `targets`, or a vector of
#'  sample conditions, indicating the primary factor to be evaluated.
#'  Defaults to the first factor in the sample conditions.
#'
#' @return A list including:
#'  * A `data.frame` of targets (sample conditions) with each factor
#'    in a separate column.
#'  * A vector indicating which column in the targets represents the
#'    primary factor.
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
#' barbieQ:::extarctTargetsAndPrimaryFactor(
#'   barbieQ = barbieQ, sampleGroups = 'Treat'
#' )
#' }
extractTargetsAndPrimaryFactor <- function(barbieQ, targets = NULL, sampleGroups = NULL) {
    ## check targets: if 'targets' is not specified, assign barbieQ$metadata (could
    ## still be NULL)
    if (is.null(targets))
        targets <- barbieQ$metadata
    ## case when targets is specified or provided by barbieQ$metadata in right format
    if (is.vector(targets) || is.factor(targets)) {
        targets <- data.frame(V1 = targets)
    } else if (is.matrix(targets) || is.data.frame(targets)) {
        if (nrow(targets) != ncol(barbieQ$assay)) {
            stop("the row dimension of 'targets' doesn't match the column dimension (sample size) of 'barbieQ$assay'.")
        }
    } else {
        ## case when targets is still NULL or not in right format if group is a
        ## vector or factor of correct length, add it to targets
        if ((is.vector(sampleGroups) || is.factor(sampleGroups)) && length(sampleGroups) >
            1L) {
            ## check sampleGroups length
            if (length(sampleGroups) != ncol(barbieQ$assay)) {
                stop("the length of 'sampleGroups' doesn't match the column dimention (sample size) of 'barbieQ$assay'.")
            } else {
                mytargets <- data.frame(sampleGroups = sampleGroups)
                pointer <- which(colnames(mytargets) == "sampleGroups")
                message("adding sampleGroups to targets.")
            }
        } else {
            stop("target not properly specified; barbieQ$metadata not provided; sampleGroups not properly specified.
           at least one of them is needed in right format.")
        }
    }

    ## now targets should be a matrix or data.frame already default sampleGroups to
    ## the first factor in targets
    if (is.null(sampleGroups)) {
        sampleGroups <- colnames(targets)[1]
    }
    ## check sampleGroups: if 'sampleGroups' is a specified effector name, extract
    ## the entire vector
    if (is.character(sampleGroups) && length(sampleGroups) == 1L) {
        if (sampleGroups %in% colnames(targets)) {
            pointer <- which(colnames(targets) == sampleGroups)
            sampleGroups <- targets[, sampleGroups]
            mytargets <- targets
            message("setting ", colnames(targets)[pointer], " as the primary effector of sample conditions.")
        } else {
            stop("sampleGroups not correspond to an effector name of sample conditionss. please ensure it is spelled correctly.")
        }
    } else if (is.vector(sampleGroups) || is.factor(sampleGroups)) {
        if (length(sampleGroups) != ncol(barbieQ$assay)) {
            stop("the length of 'sampleGroups' doesn't match the column dimention (sample size) of 'targets' or'barbieQ$assay'.")
        } else {
            mytargets <- data.frame(sampleGroups = sampleGroups, targets)
            pointer <- which(colnames(mytargets) == "sampleGroups")
            message("binding 'sampleGroups' to 'targets'.")
        }
    } else {
        sampleGroups <- rep(1, ncol(barbieQ$assay))
        mytargets <- data.frame(sampleGroups = sampleGroups, targets)
        pointer <- which(colnames(mytargets) == "sampleGroups")
        message("no properly specified 'sampleGroups'. setting samples by homogenenous group.")
    }

    return(list(mytargets = mytargets, pointer = pointer))
}
