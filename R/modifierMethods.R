#' Trim samples of a `Barbie` object based on specified specified
#'  sample conditions
#'
#' In the `matrix` and `data.frame` components of a `Barbie` object,
#'  samples are organized in columns and Barcodes in rows. The `metadata`
#'  component in `Barbie` is a `data.frame` representing sample conditions,
#'  organized by various experimental factors.
#'  `subsetSamplesByMetadata()` trims the columns based on specified
#'  `specifiedConditions` in the specified `factor` in `Barbie$metadata`.
#'
#' @param Barbie A `Barbie` object created by the [createBarbie] function.
#' @param factor A string representing the name of a factor in
#'  `Barbie$metadata`. Defaults to the first factor in `Barbie$metadata`.
#' @param specifiedConditions A string vector specifying the conditions
#'  in the specified factor. Defaults to all conditions in the specified factor.
#' @param keep A logical value: TRUE to retain or FALSE to exclude
#'  the columns (samples) of the specified conditions. Defaults to TRUE.
#'
#' @return A `Barbie` object with:
#'  * Updated components reflecting the subsetted samples, including `assay`,
#'    `metadata`, `proportion`, `CPM`, `occurrence`, `rank`, and `isTop$mat`.
#'  * Other components are inherited from the original `Barbie` object.
#'
#' @export
#'
#' @examples
#' myBarbie <- createBarbie(
#'   object = matrix(seq_len(20), nrow = 5, ncol = 4),
#'   target = data.frame(fac = c(1, 1, 2, 2))
#' )
#' subsetSamplesByMetadata(
#'   Barbie = myBarbie, factor = "fac", specifiedConditions = 1, keep = TRUE
#' )
subsetSamplesByMetadata <- function(
    Barbie, factor = NULL, specifiedConditions = NULL, keep = TRUE) {
  ## check `Barbie` format and default factor to all samples
  checkBarbieDimensions(Barbie)
  if (is.null(factor)) {
    factor <- colnames(Barbie$metadata)[1]
  }
  ## stop if the specified factor and conditions are not in the metadata of the object.
  if (!factor %in% colnames(Barbie$metadata)) {
    stop("specified factor in 'Barbie$metadata' is not found.")
  }
  if (is.null(specifiedConditions)) {
    specifiedConditions <- unique(Barbie$metadata[[factor]])
  }
  if (!all(specifiedConditions %in% Barbie$metadata[[factor]])) {
    stop("specified conditions in 'Barbie$metadata$", factor, "' is not found.")
  }
  ## choose columns specified by the factor and conditions, retaining or excluding them based on 'keep'.
  specifiedColumns <- Barbie$metadata[[factor]] %in% specifiedConditions
  retainedColumns <- if (keep) {
    specifiedColumns
  } else {
    !specifiedColumns
  }
  ## dispatch 'subsetColumns()' to update all components in the object based on the retained columns (samples).
  subsetObject <- subsetSamples(Barbie = Barbie, retainedColumns = retainedColumns)

  return(subsetObject)
}

#' Trim samples of a `Barbie` object based on specified columns
#'
#' In the `matrix` and `data.frame` components of a `Barbie` object,
#'  samples are organized in columns and Barcodes in rows. `subsetSamples()`
#'  trims the columns based on specified `retainedColumns`.
#'
#' @param Barbie A `Barbie` object created by the [createBarbie] function.
#' @param retainedColumns A logical vector indicating each column
#'  to be retained or excluded. Default to retaining all samples.
#'
#' @return A `Barbie` object with:
#'  * Updated components reflecting the subsetted samples, including `assay`,
#'    `metadata`, `proportion`, `CPM`, `occurrence`, `rank`, and `isTop$mat`.
#'  * Other components are inherited from the original `Barbie` object.
#'
#' @export
#'
#' @examples
#' myBarbie <- createBarbie(
#'   object = matrix(seq_len(20), nrow = 5, ncol = 4),
#'   target = data.frame(fac = c(1, 1, 2, 2))
#' )
#' subsetSamples(Barbie = myBarbie, retainedColumns = c(TRUE, TRUE, FALSE, FALSE))
subsetSamples <- function(Barbie, retainedColumns = NULL) {
  ## check `Barbie` format and default retainedColumns to all samples
  if (checkBarbieDimensions(Barbie)) {
    if (is.null(retainedColumns)) {
      retainedColumns <- rep(TRUE, ncol(Barbie$assay))
    }
  }
  ## update all components in the object based on the retained columns.
  subsetObject <- list(
    assay = Barbie$assay[, retainedColumns, drop = FALSE],
    metadata = Barbie$metadata[retainedColumns, , drop = FALSE],
    proportion = Barbie$proportion[, retainedColumns, drop = FALSE],
    CPM = Barbie$CPM[, retainedColumns, drop = FALSE],
    occurrence = Barbie$occurrence[, retainedColumns, drop = FALSE],
    rank = Barbie$rank[, retainedColumns, drop = FALSE],
    isTop = list(
      vec = Barbie$isTop$vec,
      mat = Barbie$isTop$mat[, retainedColumns, drop = FALSE]
    ),
    clusters = Barbie$clusters,
    factorColors = Barbie$factorColors
  )

  ## include other existing elements in 'Barbie' that are not mentioned above.
  for (elementName in names(Barbie)) {
    if (!(elementName %in% names(subsetObject))) {
      subsetObject[[elementName]] <- Barbie[[elementName]]
    }
  }

  return(subsetObject)
}

#' Merge several `Barbie` objects into one object organized by samples
#'
#' `combineBarbies` merges multiple `Barbie` objects, aligning their components
#'  based on samples. Components such as `assay`, `metadata`, `proportion`,
#'  `CPM`, `occurrence`, `rank`, and `isTop$mat` are combined by binding their
#'  columns (samples). The `isTop$vec` component is merged using a logical
#'  OR (`||`) operator, while other inherited components are passed unchanged.
#'
#' @param ... A list of `Barbie` objects created by the [createBarbie] function.
#'
#' @return A `Barbie` object including all the samples passed `Barbie` objects.
#'
#' @note All input `Barbie` objects must have the same number of Barcodes
#'  (rows) in the same order.
#'
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
#' Barbie1 <- createBarbie(
#'   object = matrix(seq_len(20), nrow = 5, ncol = 4),
#'   target = data.frame(
#'     a = c("a", "a", "b", "b"),
#'     b = c("x", "x", "x", "x")
#'   ),
#'   factorColors = list(
#'     a = c("a" = "#FFFFFF", "b" = "#CCCCFF"),
#'     b = c("x" = "#000099")
#'   )
#' )
#' Barbie2 <- createBarbie(
#'   object = matrix(seq(21, 20), nrow = 5, ncol = 2),
#'   target = data.frame(
#'     a = c("a", "c"),
#'     c = c("u", "u")
#'   ),
#'   factorColors = list(
#'     a = c("a" = "#FFFFFF", "c" = "#000000"),
#'     c = c("u" = "#CCCCCC")
#'   )
#' )
#' Barbie3 <- createBarbie(
#'   object = matrix(seq(31, 35), nrow = 5, ncol = 1),
#'   target = data.frame(b = c("x")),
#'   factorColors = list(b = c("x" = "#999999"))
#' )
#' combineBarbies(Barbie1, Barbie2, Barbie3)
combineBarbies <- function(...) {
  ## need fix unmatched rownames (Barcode ID)
  ## collect all Barbie objects into a list
  BarbieList <- list(...)
  ## start with the first Barbie object in the list
  if (checkBarbieDimensions(BarbieList[[1]])) {
    combinedObject <- BarbieList[[1]]
  }
  ## loop through the remaining Barbie objects in the list
  for (i in 2:length(BarbieList)) {
    ## check Barbie dimensions
    checkBarbieDimensions(BarbieList[[i]])
    Barbie <- BarbieList[[i]]
    ## check number of Barcodes
    if (nrow(Barbie$assay) != nrow(combinedObject$assay)) {
      stop("all `Barbie` object must contain same number of Barcodes.")
    }
    ## combine the elements of the current Barbie with the combinedObject
    combinedObject <- list(
      assay = cbind(combinedObject$assay, Barbie$assay) %>% as.data.frame(),
      ## dispatch `mergeMetadata()`
      metadata = mergeMetadata(combinedObject$metadata, Barbie$metadata) %>%
        as.data.frame(),
      proportion = cbind(combinedObject$proportion, Barbie$proportion),
      CPM = cbind(combinedObject$CPM, Barbie$CPM),
      occurrence = cbind(combinedObject$occurrence, Barbie$occurrence),
      rank = cbind(combinedObject$rank, Barbie$rank),
      isTop = list(
        vec = combinedObject$isTop$vec | Barbie$isTop$vec,
        mat = cbind(combinedObject$isTop$mat, Barbie$isTop$mat)
      ),
      clusters = paste(combinedObject$clusters, Barbie$clusters, sep = "_"),
      ## dispatch 'mergeLists()' to merge factor color lists from the objects
      factorColors = mergeLists(
        combinedObject$factorColors, Barbie$factorColors
      )
    )

    ## check and include elements present in both combinedObject and the currentBarbie but not mentioned above
    allElementNames <- unique(c(names(combinedObject), names(Barbie)))
    for (elementName in allElementNames) {
      if (!elementName %in% names(combinedObject)) {
        if (elementName %in% names(combinedObject) &&
          elementName %in% names(Barbie)) {
          ## if the element exists in both combinedObject and the currentBarbie
          ## bind columns if the element is a matrix or data.frame, otherwise treating the elements as lists
          if (is.matrix(combinedObject[[elementName]]) ||
            is.data.frame(combinedObject[[elementName]])) {
            combinedObject[[elementName]] <- cbind(
              combinedObject[[elementName]], Barbie[[elementName]]
            )
          } else {
            mergeLists(combinedObject[[elementName]], Barbie[[elementName]])
          }
        } else if (elementName %in% names(combinedObject)) {
          ## include element from combinedObject if not present in the currentBarbie
          combinedObject[[elementName]] <- combinedObject[[elementName]]
        } else if (elementName %in% names(Barbie)) {
          ## include element from Barbie if not present in combinedObject
          combinedObject[[elementName]] <- Barbie[[elementName]]
        }
      }
    }
  }

  return(combinedObject)
}

#' Merge multiple Barbie$metadata
#'
#' @param ... A list of Barbie$metadata
#'
#' @return A `data.frame`
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr bind_cols
#'
#' @noRd
#'
#' @examples \donttest{
#' Barbie:::mergeMetadata(Barbie1$metadata, Barbie2$metadata)
#' }
mergeMetadata <- function(...) {
  dfs <- list(...)
  sampleSizes <- vapply(dfs, nrow, integer(1))
  totalSample <- sum(sampleSizes)
  extendedDfs <- vector("list", length = length(dfs))
  for (i in seq(sampleSizes)) {
    exDf <- as.data.frame(matrix(NA, nrow = totalSample, ncol = ncol(dfs[[i]])))
    if (i == 1) {
      exDf[seq_len(sum(sampleSizes[1])), ] <- dfs[[i]]
    } else {
      exDf[(sum(sampleSizes[seq_len(i) - 1]) + 1):
      sum(sampleSizes[seq_len(i)]), ] <- dfs[[i]]
    }
    colnames(exDf) <- colnames(dfs[[i]])
    extendedDfs[[i]] <- exDf
  }
  factors <- unlist(lapply(dfs, colnames))

  boundDf <- dplyr::bind_cols(extendedDfs)
  boundDf[is.na(boundDf)] <- ""
  mergedDf <- tapply(names(boundDf), factors, function(cols) {
    subsetDf <- boundDf[, cols, drop = FALSE]
    apply(subsetDf, 1, paste0, collapse = "")
  }) %>%
    dplyr::bind_cols() %>%
    as.data.frame()

  rownames(mergedDf) <- lapply(dfs, rownames) %>%
    unlist() %>%
    make.unique()

  return(mergedDf)
}


#' Merge multiple `list`s into one and include all elements in the old `list`s
#'
#' @param ... A list of `Barbie` list components
#'
#' @return A `list`
#'
#' @noRd
#'
#' @examples \donttest{
#' list1 <- list(
#'   a = c("a" = "#FFFFFF", "b" = "#CCCCFF"),
#'   b = c("x" = "#000099")
#' )
#' list2 <- list(
#'   a = c("a" = "#FFFFFF", "c" = "#000000"),
#'   c = c("u" = "#CCCCCC")
#' )
#' list3 <- list(b = c("x" = "#999999"))
#' Barbie:::mergeLists(list1, list2, list3)
#' }
mergeLists <- function(...) {
  ## combine all lists into a single list of lists
  lists <- list(...)
  ## start with the first list in the argument list
  combinedList <- lists[[1]]
  ## loop through the remaining lists
  for (i in 2:length(lists)) {
    listAdding <- lists[[i]]
    ## add elements from listAdding to combinedList
    for (name in names(listAdding)) {
      if (name %in% names(combinedList)) {
        ## if the element exists in both, combine and keep unique values
        combined <- c(combinedList[[name]], listAdding[[name]])
        uniqueCombined <- combined[!duplicated(combined)]
        combinedList[[name]] <- uniqueCombined
      } else {
        ## if the element doesn't exist, simply add it
        combinedList[[name]] <- listAdding[[name]]
      }
    }
  }
  return(combinedList)
}

#' Trim Barcodes of a `Barbie` object based on specified rows
#'
#' In the `matrix` and `data.frame` components of a `Barbie` object,
#'  samples are organized in columns and Barcodes in rows. `subsetBarcodes()`
#'  trims the rows based on specified `retainedRows`.
#'
#' @param Barbie A `Barbie` object created by the [createBarbie] function.
#' @param retainedRows A logical vector indicating each row
#'  to be retained or excluded. Default to retaining all Barcodes.
#'
#' @return A `Barbie` object with:
#'  * Updated components reflecting the subsetted Barcodes, including `assay`,
#'    `proportion`, `CPM`, `occurrence`, `rank`, `isTop$mat`, `isTop$vec`,
#'    and `clusters`.
#'  * Other components are inherited from the original `Barbie` object.
#'
#' @export
#'
#' @examples
#' myBarbie <- createBarbie(
#'   object = matrix(seq_len(20), nrow = 5, ncol = 4),
#'   target = data.frame(fac = c(1, 1, 2, 2))
#' )
#' subsetBarcodes(
#'   Barbie = myBarbie, retainedRows = c(TRUE, TRUE, FALSE, FALSE, FALSE)
#' )
subsetBarcodes <- function(Barbie, retainedRows = NULL) {
  ## check `Barbie` format and default retainedRows to all Barcodes
  checkBarbieDimensions(Barbie)
  if (is.null(retainedRows)) {
    retainedRows <- rep(TRUE, nrow(Barbie$assay))
  }
  ## update all components in the object based on the retained rows.
  subsetObject <- list(
    assay = Barbie$assay[retainedRows, , drop = FALSE],
    metadata = Barbie$metadata,
    proportion = Barbie$proportion[retainedRows, , drop = FALSE],
    CPM = Barbie$CPM[retainedRows, , drop = FALSE],
    occurrence = Barbie$occurrence[retainedRows, , drop = FALSE],
    rank = Barbie$rank[retainedRows, , drop = FALSE],
    isTop = list(
      vec = Barbie$isTop$vec[retainedRows],
      mat = Barbie$isTop$mat[retainedRows, , drop = FALSE]
    ),
    clusters = Barbie$clusters[retainedRows],
    factorColors = Barbie$factorColors
  )

  ## include other existing elements in 'Barbie' that are not mentioned above.
  ## subset these elements based on 'retainedRows'
  for (elementName in names(Barbie)) {
    if (!(elementName %in% names(subsetObject))) {
      if (is.data.frame(Barbie[[elementName]]) ||
        is.matrix(Barbie[[elementName]])) {
        subsetObject[[elementName]] <- Barbie[[elementName]][retainedRows, ]
      } else if (is.vector(Barbie[[elementName]])) {
        subsetObject[[elementName]] <- Barbie[[elementName]][retainedRows]
      }
    }
  }

  return(subsetObject)
}

#' Collapse samples of a `Barbie` object based on sample grouping
#'
#' In the `matrix` and `data.frame` components of a `Barbie` object,
#'  samples are organized in columns and Barcodes in rows. `collapseSamples()`
#'  merges the columns based on specified `groupArray`. Samples within the same
#'  group are collapsed into a single column using specified aggregation method.
#'
#' @param Barbie A `Barbie` object created by the [createBarbie] function.
#' @param groupArray A vector indicating the group to which each sample belongs.
#'  Defaults to each sample being treated as a unique group.
#' @param method The name of an aggregation function used to collapse samples
#'  within the same group. Defaults to `mean`.
#'
#' @return A `Barbie` object with updated components reflecting the collapsed
#'  sample groups, including:
#'  * `assay`, `proportion`, and `CPM` aggregated by the specified `method`
#'  * `metadata` merged for all sample conditions within the group
#'  * `occurrence` aggregated by `any` within the group
#'  * `rank` aggregated by `min` within the group
#'  * `isTop$mat` aggregated by `any` within the group
#'  * Other components inherited from the original `Barbie` object.
#'
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
#' myBarbie <- createBarbie(
#'   object = matrix(seq_len(20), nrow = 5, ncol = 4),
#'   target = data.frame(fac = c(1, 1, 2, 2))
#' )
#' collapseSamples(Barbie = myBarbie, groupArray = c(1, 1, 2, 2))
collapseSamples <- function(Barbie, groupArray = NULL, method = mean) {
  ## check method
  if (!(is.function(method))) {
    stop("'method' must be a valid function.")
  }
  testVector <- c(1, 2)
  testResult <- tryCatch(method(testVector), error = function(e) NULL)
  if (is.null(testResult) || length(testResult) != 1) {
    stop("'method' must be an aggregation function that returns a single value.")
  }
  ## check object format
  checkBarbieDimensions(Barbie)
  if (is.null(groupArray)) {
    groupArray <- seq_len(ncol(Barbie$assay))
  }
  collapsedObject <- list(
    assay = Barbie$assay %>%
      collapseColumnsByArray(groupArray = groupArray, method = method) %>%
      as.data.frame(),
    metadata = Barbie$metadata %>% t() %>%
      collapseColumnsByArray(groupArray = groupArray, method = function(x) {
        checkDup <- duplicated(x)
        ifelse(
          length(checkDup) - 1 == sum(checkDup), x[[1]],
          paste(x, collapse = ".")
        )
      }) %>% t() %>% as.data.frame(),
    proportion = Barbie$proportion %>%
      collapseColumnsByArray(groupArray = groupArray, method = method),
    CPM = Barbie$CPM %>%
      collapseColumnsByArray(groupArray = groupArray, method = method),
    occurrence = Barbie$occurrence %>%
      collapseColumnsByArray(groupArray = groupArray, method = any),
    rank = Barbie$rank %>%
      collapseColumnsByArray(groupArray = groupArray, method = min),
    isTop = list(
      vec = Barbie$isTop$vec,
      mat = Barbie$isTop$mat %>%
        collapseColumnsByArray(groupArray = groupArray, method = any)
    ),
    clusters = Barbie$clusters,
    factorColors = Barbie$factorColors
  )

  ## include other existing elements in 'Barbie' that are not mentioned above.
  for (elementName in names(Barbie)) {
    if (!(elementName %in% names(collapsedObject))) {
      collapsedObject[[elementName]] <- Barbie[[elementName]]
    }
  }

  return(collapsedObject)
}

#' Collapse a matrix by column grouping
#'
#' @param mat A matrix, data.frame, or vector.
#' @param groupArray A vector indication the group to which each column belongs.
#' @param method The name of an aggregation function used to collapse columns
#'  within the same group. Defaults to `mean`.
#'
#' @return A matrix, data.frame, or vector aggregated by the specified `method`.
#'
#' @noRd
#'
#' @importFrom magrittr %>%
#'
#' @examples \donttest{
#' mymatrix <- matrix(seq_len(25), nrow = 5, ncol = 5)
#' Barbie:::collapseColumnsByArray(
#'   mymatrix,
#'   groupArray = c(1, 1, 1, 2, 2), method = mean
#' )
#' }
collapseColumnsByArray <- function(mat, groupArray, method = mean) {
  collapsedMat <- NULL
  if (is.matrix(mat) || is.data.frame(mat)) {
    subMat <- tapply(seq_len(ncol(mat)), groupArray, function(indices) {
      apply(mat[, indices, drop = FALSE], 1, method)
    })
    if (nrow(mat) > 1) {
      collapsedMat <- do.call(cbind, subMat)
    } else {
      collapsedMat <- unlist(subMat)
    }
  } else if (is.vector(mat) || is.factor(mat) || is.array(mat)) {
    subMat <- tapply(seq_along(mat), groupArray, function(indices) {
      mat[indices] %>% method()
    })
    collapsedMat <- unlist(subMat)
  }

  return(collapsedMat)
}

#' Collapse Barcodes of a `Barbie` object based on Barcode grouping
#'
#' In the `matrix` and `data.frame` components of a `Barbie` object,
#'  samples are organized in columns and Barcodes in rows. `collapseBarcodes()`
#'  merges the rows based on specified `groupArray`. Barcodes within the same
#'  group are collapsed into a single row using specified aggregation method.
#'
#' @param Barbie A `Barbie` object created by the [createBarbie] function.
#' @param groupArray A vector indicating the group to which each Barcode
#'  belongs. Defaults to each Barcode being treated as a unique group.
#' @param method The name of an aggregation function used to collapse barcodes
#'  within the same group. Defaults to `max`.
#'
#' @return A `Barbie` object with updated components reflecting the collapsed
#'  Barcode groups, including:
#'  * `assay`, `proportion`, and `CPM` aggregated by the specified `method`
#'  * `occurrence` aggregated by `any` within the group
#'  * `rank` aggregated by `min` within the group
#'  * `isTop$mat` and `isTop$vec` aggregated by `any` within the group
#'  * `clusters` aggregated by the most abundant cluster name within the group
#'  * Other components inherited from the original `Barbie` object.
#'
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
#' myBarbie <- createBarbie(
#'   object = matrix(seq_len(20), nrow = 5, ncol = 4),
#'   target = data.frame(fac = c(1, 1, 2, 2))
#' )
#' collapseBarcodes(Barbie = myBarbie, groupArray = c(1, 1, 2, 2, 2))
collapseBarcodes <- function(Barbie, groupArray = NULL, method = max) {
  ## check method
  if (!(is.function(method))) {
    stop("'method' must be a valid function.")
  }
  testVector <- c(1, 2)
  testResult <- tryCatch(method(testVector), error = function(e) NULL)
  if (is.null(testResult) || length(testResult) != 1) {
    stop("'method' must be an aggregation function that returns a single value.")
  }
  ## check object format
  checkBarbieDimensions(Barbie)
  if (is.null(groupArray)) {
    groupArray <- seq_len(nrow(Barbie$assay))
  }
  collapsedObject <- list(
    assay = Barbie$assay %>%
      collapseRowsByArray(groupArray = groupArray, method = method) %>% as.data.frame(),
    metadata = Barbie$metadata,
    proportion = Barbie$proportion %>%
      collapseRowsByArray(groupArray = groupArray, method = method),
    CPM = Barbie$CPM %>%
      collapseRowsByArray(groupArray = groupArray, method = method),
    occurrence = Barbie$occurrence %>%
      collapseRowsByArray(groupArray = groupArray, method = any),
    rank = Barbie$rank %>%
      collapseRowsByArray(groupArray = groupArray, method = min),
    isTop = list(
      vec = Barbie$isTop$vec %>%
        collapseRowsByArray(groupArray = groupArray, method = any),
      mat = Barbie$isTop$mat %>%
        collapseRowsByArray(groupArray = groupArray, method = any)
    ),
    clusters = Barbie$clusters %>%
      collapseRowsByArray(
        groupArray = groupArray,
        method = function(x) names(which.max(table(x)))
      ),
    factorColors = Barbie$factorColors
  )

  ## include existing elements in the collapsedObject
  for (elementName in names(Barbie)) {
    if (!(elementName %in% names(collapsedObject))) {
      collapsedObject[[elementName]] <- Barbie[[elementName]] %>%
        collapseRowsByArray(
          groupArray = groupArray,
          method = function(x) {
            x[[1]]
          }
        )
    }
  }

  return(collapsedObject)
}

#' collapse a matrix by row grouping
#'
#' @param mat A matrix, data.frame, or vector.
#' @param groupArray A vector indication the group to which each row belongs.
#' @param method The name of an aggregation function used to collapse rows
#'  within the same group. Defaults to `max`.
#'
#' @return A matrix, data.frame, or vector aggregated by the specified `method`.
#'
#' @noRd
#'
#' @importFrom magrittr %>%
#'
#' @examples \donttest{
#' mymatrix <- matrix(seq_len(25), nrow = 5, ncol = 5)
#' Barbie:::collapseRowsByArray(
#'   mymatrix,
#'   groupArray = c(1, 1, 1, 2, 2), method = mean
#' )
#' }
collapseRowsByArray <- function(mat, groupArray, method = max) {
  collapsedMat <- NULL
  if (is.matrix(mat) || is.data.frame(mat)) {
    subMat <- tapply(seq_len(nrow(mat)), groupArray, function(indices) {
      apply(mat[indices, , drop = FALSE], 2, method)
    })
    if (ncol(mat) > 1) {
      collapsedMat <- do.call(rbind, subMat)
    } else {
      collapsedMat <- unlist(subMat)
    }
  } else if (is.vector(mat)) {
    subMat <- tapply(seq_along(mat), groupArray, function(indices) {
      mat[indices] %>% method()
    })
    collapsedMat <- subMat
  }

  return(collapsedMat)
}
