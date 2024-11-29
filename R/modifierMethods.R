#' Trim samples of a `barbieQ` object based on specified specified
#'  sample conditions
#'
#' In the `matrix` and `data.frame` components of a `barbieQ` object,
#'  samples are organized in columns and Barcodes in rows. The `metadata`
#'  component in `barbieQ` is a `data.frame` representing sample conditions,
#'  organized by various experimental factors.
#'  `subsetSamplesByMetadata()` trims the columns based on specified
#'  `specifiedConditions` in the specified `factor` in `barbieQ$metadata`.
#'
#' @param barbieQ A `barbieQ` object created by the [createBarbieQ] function.
#' @param factor A string representing the name of a factor in
#'  `barbieQ$metadata`. Defaults to the first factor in `barbieQ$metadata`.
#' @param specifiedConditions A string vector specifying the conditions
#'  in the specified factor. Defaults to all conditions in the specified factor.
#' @param keep A logical value: TRUE to retain or FALSE to exclude
#'  the columns (samples) of the specified conditions. Defaults to TRUE.
#'
#' @return A `barbieQ` object with:
#'  * Updated components reflecting the subsetted samples, including `assay`,
#'    `metadata`, `proportion`, `CPM`, `occurrence`, `rank`, and `isTop$mat`.
#'  * Other components are inherited from the original `barbieQ` object.
#'
#' @export
#'
#' @examples
#' myBarbieQ <- createBarbieQ(
#'   object = matrix(seq_len(20), nrow = 5, ncol = 4),
#'   target = data.frame(fac = c(1, 1, 2, 2))
#' )
#' subsetSamplesByMetadata(
#'   barbieQ = myBarbieQ, factor = "fac", specifiedConditions = 1, keep = TRUE
#' )
subsetSamplesByMetadata <- function(
    barbieQ, factor = NULL, specifiedConditions = NULL, keep = TRUE) {
  ## check `barbieQ` format and default factor to all samples
  checkBarbieQDimensions(barbieQ)
  if (is.null(factor)) {
    factor <- colnames(barbieQ$metadata)[1]
  }
  ## stop if the specified factor and conditions are not in the metadata of the object.
  if (!factor %in% colnames(barbieQ$metadata)) {
    stop("specified factor in 'barbieQ$metadata' is not found.")
  }
  if (is.null(specifiedConditions)) {
    specifiedConditions <- unique(barbieQ$metadata[[factor]])
  }
  if (!all(specifiedConditions %in% barbieQ$metadata[[factor]])) {
    stop("specified conditions in 'barbieQ$metadata$", factor, "' is not found.")
  }
  ## choose columns specified by the factor and conditions, retaining or excluding them based on 'keep'.
  specifiedColumns <- barbieQ$metadata[[factor]] %in% specifiedConditions
  retainedColumns <- if (keep) {
    specifiedColumns
  } else {
    !specifiedColumns
  }
  ## dispatch 'subsetColumns()' to update all components in the object based on the retained columns (samples).
  subsetObject <- subsetSamples(barbieQ = barbieQ, retainedColumns = retainedColumns)

  return(subsetObject)
}

#' Trim samples of a `barbieQ` object based on specified columns
#'
#' In the `matrix` and `data.frame` components of a `barbieQ` object,
#'  samples are organized in columns and Barcodes in rows. `subsetSamples()`
#'  trims the columns based on specified `retainedColumns`.
#'
#' @param barbieQ A `barbieQ` object created by the [createBarbieQ] function.
#' @param retainedColumns A logical vector indicating each column
#'  to be retained or excluded. Default to retaining all samples.
#'
#' @return A `barbieQ` object with:
#'  * Updated components reflecting the subsetted samples, including `assay`,
#'    `metadata`, `proportion`, `CPM`, `occurrence`, `rank`, and `isTop$mat`.
#'  * Other components are inherited from the original `barbieQ` object.
#'
#' @export
#'
#' @examples
#' myBarbieQ <- createBarbieQ(
#'   object = matrix(seq_len(20), nrow = 5, ncol = 4),
#'   target = data.frame(fac = c(1, 1, 2, 2))
#' )
#' subsetSamples(barbieQ = myBarbieQ, retainedColumns = c(TRUE, TRUE, FALSE, FALSE))
subsetSamples <- function(barbieQ, retainedColumns = NULL) {
  ## check `barbieQ` format and default retainedColumns to all samples
  if (checkBarbieQDimensions(barbieQ)) {
    if (is.null(retainedColumns)) {
      retainedColumns <- rep(TRUE, ncol(barbieQ$assay))
    }
  }
  ## update all components in the object based on the retained columns.
  subsetObject <- list(
    assay = barbieQ$assay[, retainedColumns, drop = FALSE],
    metadata = barbieQ$metadata[retainedColumns, , drop = FALSE],
    proportion = barbieQ$proportion[, retainedColumns, drop = FALSE],
    CPM = barbieQ$CPM[, retainedColumns, drop = FALSE],
    occurrence = barbieQ$occurrence[, retainedColumns, drop = FALSE],
    rank = barbieQ$rank[, retainedColumns, drop = FALSE],
    isTop = list(
      vec = barbieQ$isTop$vec,
      mat = barbieQ$isTop$mat[, retainedColumns, drop = FALSE]
    ),
    clusters = barbieQ$clusters,
    factorColors = barbieQ$factorColors
  )

  ## include other existing elements in 'barbieQ' that are not mentioned above.
  for (elementName in names(barbieQ)) {
    if (!(elementName %in% names(subsetObject))) {
      subsetObject[[elementName]] <- barbieQ[[elementName]]
    }
  }

  return(subsetObject)
}

#' Merge several `barbieQ` objects into one object organized by samples
#'
#' `combineBarbieQs` merges multiple `barbieQ` objects, aligning their components
#'  based on samples. Components such as `assay`, `metadata`, `proportion`,
#'  `CPM`, `occurrence`, `rank`, and `isTop$mat` are combined by binding their
#'  columns (samples). The `isTop$vec` component is merged using a logical
#'  OR (`||`) operator, while other inherited components are passed unchanged.
#'
#' @param ... A list of `barbieQ` objects created by the [createBarbieQ] function.
#'
#' @return A `barbieQ` object including all the samples passed `barbieQ` objects.
#'
#' @note All input `barbieQ` objects must have the same number of Barcodes
#'  (rows) in the same order.
#'
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
#' barbieQ1 <- createBarbieQ(
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
#' barbieQ2 <- createBarbieQ(
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
#' barbieQ3 <- createBarbieQ(
#'   object = matrix(seq(31, 35), nrow = 5, ncol = 1),
#'   target = data.frame(b = c("x")),
#'   factorColors = list(b = c("x" = "#999999"))
#' )
#' combineBarbieQs(barbieQ1, barbieQ2, barbieQ3)
combineBarbieQs <- function(...) {
  ## need fix unmatched rownames (Barcode ID)
  ## collect all barbieQ objects into a list
  BarbieQList <- list(...)
  ## start with the first barbieQ object in the list
  if (checkBarbieQDimensions(BarbieQList[[1]])) {
    combinedObject <- BarbieQList[[1]]
  }
  ## loop through the remaining barbieQ objects in the list
  for (i in 2:length(BarbieQList)) {
    ## check barbieQ dimensions
    checkBarbieQDimensions(BarbieQList[[i]])
    barbieQ <- BarbieQList[[i]]
    ## check number of Barcodes
    if (nrow(barbieQ$assay) != nrow(combinedObject$assay)) {
      stop("all `barbieQ` object must contain same number of Barcodes.")
    }
    ## combine the elements of the current barbieQ with the combinedObject
    combinedObject <- list(
      assay = cbind(combinedObject$assay, barbieQ$assay) %>% as.data.frame(),
      ## dispatch `mergeMetadata()`
      metadata = mergeMetadata(combinedObject$metadata, barbieQ$metadata) %>%
        as.data.frame(),
      proportion = cbind(combinedObject$proportion, barbieQ$proportion),
      CPM = cbind(combinedObject$CPM, barbieQ$CPM),
      occurrence = cbind(combinedObject$occurrence, barbieQ$occurrence),
      rank = cbind(combinedObject$rank, barbieQ$rank),
      isTop = list(
        vec = combinedObject$isTop$vec | barbieQ$isTop$vec,
        mat = cbind(combinedObject$isTop$mat, barbieQ$isTop$mat)
      ),
      clusters = paste(combinedObject$clusters, barbieQ$clusters, sep = "_"),
      ## dispatch 'mergeLists()' to merge factor color lists from the objects
      factorColors = mergeLists(
        combinedObject$factorColors, barbieQ$factorColors
      )
    )

    ## check and include elements present in both combinedObject and the currentBarbieQ but not mentioned above
    allElementNames <- unique(c(names(combinedObject), names(barbieQ)))
    for (elementName in allElementNames) {
      if (!elementName %in% names(combinedObject)) {
        if (elementName %in% names(combinedObject) &&
          elementName %in% names(barbieQ)) {
          ## if the element exists in both combinedObject and the currentBarbieQ
          ## bind columns if the element is a matrix or data.frame, otherwise treating the elements as lists
          if (is.matrix(combinedObject[[elementName]]) ||
            is.data.frame(combinedObject[[elementName]])) {
            combinedObject[[elementName]] <- cbind(
              combinedObject[[elementName]], barbieQ[[elementName]]
            )
          } else {
            mergeLists(combinedObject[[elementName]], barbieQ[[elementName]])
          }
        } else if (elementName %in% names(combinedObject)) {
          ## include element from combinedObject if not present in the currentBarbieQ
          combinedObject[[elementName]] <- combinedObject[[elementName]]
        } else if (elementName %in% names(barbieQ)) {
          ## include element from barbieQ if not present in combinedObject
          combinedObject[[elementName]] <- barbieQ[[elementName]]
        }
      }
    }
  }

  return(combinedObject)
}

#' Merge multiple barbieQ$metadata
#'
#' @param ... A list of barbieQ$metadata
#'
#' @return A `data.frame`
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr bind_cols
#'
#' @noRd
#'
#' @examples \donttest{
#' barbieQ:::mergeMetadata(barbieQ1$metadata, barbieQ2$metadata)
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
#' @param ... A list of `barbieQ` list components
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
#' barbieQ:::mergeLists(list1, list2, list3)
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

#' Trim Barcodes of a `barbieQ` object based on specified rows
#'
#' In the `matrix` and `data.frame` components of a `barbieQ` object,
#'  samples are organized in columns and Barcodes in rows. `subsetBarcodes()`
#'  trims the rows based on specified `retainedRows`.
#'
#' @param barbieQ A `barbieQ` object created by the [createBarbieQ] function.
#' @param retainedRows A logical vector indicating each row
#'  to be retained or excluded. Default to retaining all Barcodes.
#'
#' @return A `barbieQ` object with:
#'  * Updated components reflecting the subsetted Barcodes, including `assay`,
#'    `proportion`, `CPM`, `occurrence`, `rank`, `isTop$mat`, `isTop$vec`,
#'    and `clusters`.
#'  * Other components are inherited from the original `barbieQ` object.
#'
#' @export
#'
#' @examples
#' myBarbieQ <- createBarbieQ(
#'   object = matrix(seq_len(20), nrow = 5, ncol = 4),
#'   target = data.frame(fac = c(1, 1, 2, 2))
#' )
#' subsetBarcodes(
#'   barbieQ = myBarbieQ, retainedRows = c(TRUE, TRUE, FALSE, FALSE, FALSE)
#' )
subsetBarcodes <- function(barbieQ, retainedRows = NULL) {
  ## check `barbieQ` format and default retainedRows to all Barcodes
  checkBarbieQDimensions(barbieQ)
  if (is.null(retainedRows)) {
    retainedRows <- rep(TRUE, nrow(barbieQ$assay))
  }
  ## update all components in the object based on the retained rows.
  subsetObject <- list(
    assay = barbieQ$assay[retainedRows, , drop = FALSE],
    metadata = barbieQ$metadata,
    proportion = barbieQ$proportion[retainedRows, , drop = FALSE],
    CPM = barbieQ$CPM[retainedRows, , drop = FALSE],
    occurrence = barbieQ$occurrence[retainedRows, , drop = FALSE],
    rank = barbieQ$rank[retainedRows, , drop = FALSE],
    isTop = list(
      vec = barbieQ$isTop$vec[retainedRows],
      mat = barbieQ$isTop$mat[retainedRows, , drop = FALSE]
    ),
    clusters = barbieQ$clusters[retainedRows],
    factorColors = barbieQ$factorColors
  )

  ## include other existing elements in 'barbieQ' that are not mentioned above.
  ## subset these elements based on 'retainedRows'
  for (elementName in names(barbieQ)) {
    if (!(elementName %in% names(subsetObject))) {
      if (is.data.frame(barbieQ[[elementName]]) ||
        is.matrix(barbieQ[[elementName]])) {
        subsetObject[[elementName]] <- barbieQ[[elementName]][retainedRows, ]
      } else if (is.vector(barbieQ[[elementName]])) {
        subsetObject[[elementName]] <- barbieQ[[elementName]][retainedRows]
      }
    }
  }

  return(subsetObject)
}

#' Collapse samples of a `barbieQ` object based on sample grouping
#'
#' In the `matrix` and `data.frame` components of a `barbieQ` object,
#'  samples are organized in columns and Barcodes in rows. `collapseSamples()`
#'  merges the columns based on specified `groupArray`. Samples within the same
#'  group are collapsed into a single column using specified aggregation method.
#'
#' @param barbieQ A `barbieQ` object created by the [createBarbieQ] function.
#' @param groupArray A vector indicating the group to which each sample belongs.
#'  Defaults to each sample being treated as a unique group.
#' @param method The name of an aggregation function used to collapse samples
#'  within the same group. Defaults to `mean`.
#'
#' @return A `barbieQ` object with updated components reflecting the collapsed
#'  sample groups, including:
#'  * `assay`, `proportion`, and `CPM` aggregated by the specified `method`
#'  * `metadata` merged for all sample conditions within the group
#'  * `occurrence` aggregated by `any` within the group
#'  * `rank` aggregated by `min` within the group
#'  * `isTop$mat` aggregated by `any` within the group
#'  * Other components inherited from the original `barbieQ` object.
#'
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
#' myBarbieQ <- createBarbieQ(
#'   object = matrix(seq_len(20), nrow = 5, ncol = 4),
#'   target = data.frame(fac = c(1, 1, 2, 2))
#' )
#' collapseSamples(barbieQ = myBarbieQ, groupArray = c(1, 1, 2, 2))
collapseSamples <- function(barbieQ, groupArray = NULL, method = mean) {
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
  checkBarbieQDimensions(barbieQ)
  if (is.null(groupArray)) {
    groupArray <- seq_len(ncol(barbieQ$assay))
  }
  collapsedObject <- list(
    assay = barbieQ$assay %>%
      collapseColumnsByArray(groupArray = groupArray, method = method) %>%
      as.data.frame(),
    metadata = barbieQ$metadata %>% t() %>%
      collapseColumnsByArray(groupArray = groupArray, method = function(x) {
        checkDup <- duplicated(x)
        ifelse(
          length(checkDup) - 1 == sum(checkDup), x[[1]],
          paste(x, collapse = ".")
        )
      }) %>% t() %>% as.data.frame(),
    proportion = barbieQ$proportion %>%
      collapseColumnsByArray(groupArray = groupArray, method = method),
    CPM = barbieQ$CPM %>%
      collapseColumnsByArray(groupArray = groupArray, method = method),
    occurrence = barbieQ$occurrence %>%
      collapseColumnsByArray(groupArray = groupArray, method = any),
    rank = barbieQ$rank %>%
      collapseColumnsByArray(groupArray = groupArray, method = min),
    isTop = list(
      vec = barbieQ$isTop$vec,
      mat = barbieQ$isTop$mat %>%
        collapseColumnsByArray(groupArray = groupArray, method = any)
    ),
    clusters = barbieQ$clusters,
    factorColors = barbieQ$factorColors
  )

  ## include other existing elements in 'barbieQ' that are not mentioned above.
  for (elementName in names(barbieQ)) {
    if (!(elementName %in% names(collapsedObject))) {
      collapsedObject[[elementName]] <- barbieQ[[elementName]]
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
#' barbieQ:::collapseColumnsByArray(
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

#' Collapse Barcodes of a `barbieQ` object based on Barcode grouping
#'
#' In the `matrix` and `data.frame` components of a `barbieQ` object,
#'  samples are organized in columns and Barcodes in rows. `collapseBarcodes()`
#'  merges the rows based on specified `groupArray`. Barcodes within the same
#'  group are collapsed into a single row using specified aggregation method.
#'
#' @param barbieQ A `barbieQ` object created by the [createBarbieQ] function.
#' @param groupArray A vector indicating the group to which each Barcode
#'  belongs. Defaults to each Barcode being treated as a unique group.
#' @param method The name of an aggregation function used to collapse barcodes
#'  within the same group. Defaults to `max`.
#'
#' @return A `barbieQ` object with updated components reflecting the collapsed
#'  Barcode groups, including:
#'  * `assay`, `proportion`, and `CPM` aggregated by the specified `method`
#'  * `occurrence` aggregated by `any` within the group
#'  * `rank` aggregated by `min` within the group
#'  * `isTop$mat` and `isTop$vec` aggregated by `any` within the group
#'  * `clusters` aggregated by the most abundant cluster name within the group
#'  * Other components inherited from the original `barbieQ` object.
#'
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
#' myBarbieQ <- createBarbieQ(
#'   object = matrix(seq_len(20), nrow = 5, ncol = 4),
#'   target = data.frame(fac = c(1, 1, 2, 2))
#' )
#' collapseBarcodes(barbieQ = myBarbieQ, groupArray = c(1, 1, 2, 2, 2))
collapseBarcodes <- function(barbieQ, groupArray = NULL, method = max) {
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
  checkBarbieQDimensions(barbieQ)
  if (is.null(groupArray)) {
    groupArray <- seq_len(nrow(barbieQ$assay))
  }
  collapsedObject <- list(
    assay = barbieQ$assay %>%
      collapseRowsByArray(groupArray = groupArray, method = method) %>% as.data.frame(),
    metadata = barbieQ$metadata,
    proportion = barbieQ$proportion %>%
      collapseRowsByArray(groupArray = groupArray, method = method),
    CPM = barbieQ$CPM %>%
      collapseRowsByArray(groupArray = groupArray, method = method),
    occurrence = barbieQ$occurrence %>%
      collapseRowsByArray(groupArray = groupArray, method = any),
    rank = barbieQ$rank %>%
      collapseRowsByArray(groupArray = groupArray, method = min),
    isTop = list(
      vec = barbieQ$isTop$vec %>%
        collapseRowsByArray(groupArray = groupArray, method = any),
      mat = barbieQ$isTop$mat %>%
        collapseRowsByArray(groupArray = groupArray, method = any)
    ),
    clusters = barbieQ$clusters %>%
      collapseRowsByArray(
        groupArray = groupArray,
        method = function(x) names(which.max(table(x)))
      ),
    factorColors = barbieQ$factorColors
  )

  ## include existing elements in the collapsedObject
  for (elementName in names(barbieQ)) {
    if (!(elementName %in% names(collapsedObject))) {
      collapsedObject[[elementName]] <- barbieQ[[elementName]] %>%
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
#' barbieQ:::collapseRowsByArray(
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
