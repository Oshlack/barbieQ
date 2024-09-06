#' Title make a subset of columns from all components of a 'Barbie' object based on specified conditions in the metadata
#'
#' @param Barbie an object created by the 'createBarbie()' function in the 'Barbie' package
#' @param factor a string specifying the factor name in 'Barbie$metadata'
#' @param specifiedConditions a string or a string vector specifying the conditions in the specified factor
#' @param keep a logical value indicating whether to retain or exclude the columns (samples) of the specified conditions
#'
#' @return an updated Barbie object with subsets of its components
#' @export
#'
#' @examples
#' myBarbie <- createBarbie(
#'   object=matrix(c(1:20), nrow=5, ncol=4),
#'   target=data.frame(fac=c(1, 1, 2, 2)))
#' subsetColumnsByMetadata(Barbie=myBarbie, factor="fac", specifiedConditions=1, keep=TRUE)
subsetSamplesByMetadata <- function(Barbie, factor, specifiedConditions, keep = TRUE) {
  ## stop if the specified factor and conditions are not in the metadata of the object.
  if (!factor %in% colnames(Barbie$metadata))
    stop("specified factor in 'Barbie$metadata' is not found.")
  if (!all(specifiedConditions %in% Barbie$metadata[[factor]]))
    stop("specified conditions in 'Barbie$metadata$", factor, "' is not found.")
  ## choose columns specified by the factor and conditions, retaining or excluding them based on 'keep'.
  specifiedColumns <- Barbie$metadata[[factor]] %in% specifiedConditions
  retainedColumns <- if(keep) specifiedColumns
    else !specifiedColumns
  ## dispatch 'subsetColumns()' to update all components in the object based on the retained columns (samples).
  subsetObject <- subsetColumns(Barbie=Barbie, retainedColumns=retainedColumns)

  return(subsetObject)
}

#' Title make a subset of columns from all components of a 'Barbie' object based on specified columns
#'
#' @param Barbie an object created by the 'createBarbie()' function in the 'Barbie' package
#' @param retainedColumns a logical vector indicating each column to be retained or excluded
#'
#' @return an updated Barbie object with subsets of its components
#' @export
#'
#' @examples
#' myBarbie <- createBarbie(
#'   object=matrix(c(1:20), nrow=5, ncol=4),
#'   target=data.frame(fac=c(1, 1, 2, 2)))
#' subsetColumns(Barbie=myBarbie, retainedColumns=c(TRUE, TRUE, FALSE, FALSE))
subsetSamples <- function(Barbie, retainedColumns) {
  ## update all components in the object based on the retained columns (samples).
  subsetObject <- list(
    assay = Barbie$assay[, retainedColumns, drop=FALSE],
    metadata = Barbie$metadata[retainedColumns, , drop=FALSE],
    proportion = Barbie$proportion[, retainedColumns, drop=FALSE],
    CPM = Barbie$CPM[, retainedColumns, drop=FALSE],
    occurrence = Barbie$occurrence[, retainedColumns, drop=FALSE],
    rank = Barbie$rank[, retainedColumns,drop=FALSE],
    isTop = list(vec = Barbie$isTop$vec,
                 mat = Barbie$isTop$mat[, retainedColumns, drop=FALSE]),
    clusters = Barbie$clusters,
    factorColors = Barbie$factorColors)

  ## include other existing elements in 'Barbie' that are not mentioned above.
  for (elementName in names(Barbie)) {
    if (!(elementName %in% names(subsetObject)))
      subsetObject[[elementName]] <- Barbie[[elementName]]
  }

  return(subsetObject)
}

#' Title combining several 'Barbie' object by binding matrix components by samples and merging list components
#'
#' @param ... a list of 'Barbie' object
#'
#' @return a 'Barbie' object
#' @export
#'
#' @import dplyr
#' @import magrittr
#'
#' @examples
#' combineBarbies(Barbie1, Barbie2, Barbie3)
combineBarbies <- function(...) {
  ## need fix unmatched rownames (Barcode ID)
  ## collect all Barbie objects into a list
  BarbieList <- list(...)
  ## start with the first Barbie object in the list
  combinedObject <- BarbieList[[1]]
  ## loop through the remaining Barbie objects in the list
  for (i in 2:length(BarbieList)) {
    Barbie <- BarbieList[[i]]
    ## combine the elements of the current Barbie with the combinedObject
    combinedObject <- list(
      assay = cbind(combinedObject$assay, Barbie$assay) %>% as.data.frame(),
      ## dispatch 'mergeMetadata()'
      metadata = mergeMetadata(combinedObject$metadata, Barbie$metadata) %>% as.data.frame(),
      proportion = cbind(combinedObject$proportion, Barbie$proportion),
      CPM = cbind(combinedObject$CPM, Barbie$CPM),
      occurrence = cbind(combinedObject$occurrence, Barbie$occurrence),
      rank = cbind(combinedObject$rank, Barbie$rank),
      isTop = list(vec=combinedObject$isTop$vec | Barbie$isTop$vec,
                   mat=cbind(combinedObject$isTop$mat, Barbie$isTop$mat)),
      clusters = paste(combinedObject$clusters, Barbie$clusters, sep = "_"),
      ## dispatch 'mergeLists()' to merge factor color lists from the objects
      factorColors = mergeLists(combinedObject$factorColors, Barbie$factorColors))

    ## check and include elements present in both combinedObject and the currentBarbie but not mentioned above
    allElementNames <- unique(c(names(combinedObject), names(Barbie)))
    for (elementName in allElementNames) {
      if (!elementName %in% names(combinedObject))
        if (elementName %in% names(combinedObject) && elementName %in% names(Barbie)) {
        ## if the element exists in both combinedObject and the currentBarbie
        ## bind columns if the element is a matrix or data.frame, otherwise treating the elements as lists
        if (is.matrix(combinedObject[[elementName]]) || is.data.frame(combinedObject[[elementName]]))
          combinedObject[[elementName]] <- cbind(combinedObject[[elementName]], Barbie[[elementName]])
        else mergeLists(combinedObject[[elementName]], Barbie[[elementName]])
      } else if (elementName %in% names(combinedObject)) {
        ## include element from combinedObject if not present in the currentBarbie
        combinedObject[[elementName]] <- combinedObject[[elementName]]
      } else if (elementName %in% names(Barbie)) {
        ## include element from Barbie if not present in combinedObject
        combinedObject[[elementName]] <- Barbie[[elementName]]
      }
    }
  }

  return(combinedObject)
}

#' Title merging several Barbie$metadata
#'
#' @param ... a list of metadata
#'
#' @return a data.frame
#'
#' @import magrittr
#'
#' @examples
#' mergeMetadata(Barbie1$metadata, Barbie2$metadata)
mergeMetadata <- function(...) {
  dfs <- list(...)
  sampleSizes <- vapply(dfs, nrow, integer(1))
  totalSample <- sum(sampleSizes)
  extendedDfs <- vector("list", length=length(dfs))
  for(i in seq(sampleSizes)) {
    exDf <- as.data.frame(matrix(NA, nrow=totalSample, ncol=ncol(dfs[[i]])))
    if(i == 1) exDf[1:sum(sampleSizes[1]),] <- dfs[[i]]
    else exDf[(sum(sampleSizes[1:i-1])+1):sum(sampleSizes[1:i]),] <- dfs[[i]]
    colnames(exDf) <- colnames(dfs[[i]])
    extendedDfs[[i]] <- exDf
  }
  factors <- unlist(lapply(dfs, colnames))

  boundDf <- bind_cols(extendedDfs)
  boundDf[is.na(boundDf)] <- ""
  mergedDf <- tapply(names(boundDf), factors, function(cols) {
    subsetDf <- boundDf[, cols, drop=FALSE]
    apply(subsetDf, 1, paste0, collapse = "")
  }) %>% bind_cols() %>%
    as.data.frame()

  rownames(mergedDf) <- lapply(dfs, rownames) %>%
    unlist() %>%
    make.unique()

  return(mergedDf)
}


#' Title Merging several lists into one list and include all elements in the old lists
#'
#' @param ... a list of lists
#'
#' @return a list
#'
#' @examples
#' list1 <- list(a = c("a"="#FFFFFF", "b" = "#CCCCFF"),
#'               b = c("x"="#000099"))
#' list2 <- list(a = c("a" = "#FFFFFF", "c"="#000000"),
#'               c = c("u" = "#CCCCCC"))
#' list3 <- list(b = c("x"="#999999"))
#' mergeLists(list1, list2, list3)

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

#' Title make a subset of rows from all components of a 'Barbie' object based on specified rows
#'
#' @param Barbiean object created by the 'createBarbie()' function in the 'Barbie' package
#' @param retainedRows a logical vector indicating each row to be retained or excluded
#'
#' @return an updated Barbie object with subsets of its components
#' @export
#'
#' @examples
#' myBarbie <- createBarbie(
#'   object=matrix(c(1:20), nrow=5, ncol=4),
#'   target=data.frame(fac=c(1, 1, 2, 2)))
#' subsetBarcodes(Barbie=myBarbie, retainedRows=c(TRUE, TRUE, FALSE, FALSE, FALSE))
subsetBarcodes <- function(Barbie, retainedRows) {
  subsetObject <- list(
    assay = Barbie$assay[retainedRows, ,drop=FALSE],
    metadata = Barbie$metadata,
    proportion = Barbie$proportion[retainedRows, ,drop=FALSE],
    CPM = Barbie$CPM[retainedRows, ,drop=FALSE],
    occurrence = Barbie$occurrence[retainedRows, ,drop=FALSE],
    rank = Barbie$rank[retainedRows, ,drop=FALSE],
    isTop = list(vec = Barbie$isTop$vec[retainedRows],
                 mat = Barbie$isTop$mat[retainedRows, , drop=FALSE]),
    clusters = Barbie$clusters[retainedRows],
    factorColors = Barbie$factorColors)

  ## include other existing elements in 'Barbie' that are not mentioned above.
  ## subset these elements based on 'retainedRows'
  for (elementName in names(Barbie)) {
    if (!(elementName %in% names(subsetObject))) {
      if(is.data.frame(Barbie[[elementName]]) || is.matrix(Barbie[[elementName]]))
        subsetObject[[elementName]] <- Barbie[[elementName]][retainedRows,]
      else if(is.vector(Barbie[[elementName]]))
        subsetObject[[elementName]] <- Barbie[[elementName]][retainedRows]
    }
  }

  return(subsetObject)
}

#' Title collapse a 'Barbie' object by sample groups
#'
#' @param Barbie an object created by the 'createBarbie()' function in the 'Barbie' package
#' @param groupArray a vector indication the group of each sample
#'
#' @return a 'Barbie' object with samples collapsed by 'groupArray'
#' @export
#'
#' @import magrittr
#'
#' @examples
#' myBarbie <- createBarbie(
#'   object=matrix(c(1:20), nrow=5, ncol=4),
#'   target=data.frame(fac=c(1, 1, 2, 2)))
#' collapseSamples(Barbie=myBarbie, groupArray=c(1,1,2,2))
collapseSamples <- function(Barbie, groupArray) {
  collapsedObject <- list(
    assay = Barbie$assay %>%
      collapseColumnsByArray(groupArray=groupArray) %>% as.data.frame(),
    metadata = Barbie$metadata %>% t() %>%
      collapseColumnsByArray(groupArray=groupArray, method=function(x) {
        checkDup <- duplicated(x)
        ifelse(length(checkDup) -1 == sum(checkDup), x[[1]], paste(x, collapse = "."))
      }) %>% t() %>% as.data.frame(),
    proportion = Barbie$proportion %>%
      collapseColumnsByArray(groupArray=groupArray),
    CPM = Barbie$CPM %>%
      collapseColumnsByArray(groupArray=groupArray),
    occurrence = Barbie$occurrence %>%
      collapseColumnsByArray(groupArray=groupArray, method=max),
    rank = Barbie$rank %>%
      collapseColumnsByArray(groupArray=groupArray, method=min),
    isTop = list(vec = Barbie$isTop$vec,
                 mat = Barbie$isTop$mat %>%
                   collapseColumnsByArray(groupArray=groupArray, method=any)),
    clusters = Barbie$clusters,
    factorColors = Barbie$factorColors)

  ## include other existing elements in 'Barbie' that are not mentioned above.
  for (elementName in names(Barbie)) {
    if (!(elementName %in% names(collapsedObject)))
      collapsedObject[[elementName]] <- Barbie[[elementName]]
  }

  return(collapsedObject)
}

#' Title collapse a matrix by column groups
#'
#' @param mat a matrix, data.frame, or vector
#' @param groupArray a vector indication the groups of each column in 'mat'
#' @param method a function indicating the method for collapsing each Barcodes in the subset matrix of each sample group
#'
#' @return a matrix, data.frame, or vector collapsed by column groups
#'
#' @import magrittr
#'
#' @examples
#' mymatrix <- matrix(1:25, nrow=5, ncol=5)
#' collapseColumnsByArray(mymatrix, groupArray=c(1,1,1,2,2), method=mean)
collapseColumnsByArray <- function(mat, groupArray, method=mean) {
  if (is.matrix(mat) || is.data.frame(mat)) {
    subMat <- tapply(seq_len(ncol(mat)), groupArray, function(indices)
      apply(mat[, indices, drop=F], 1, method))
    if(nrow(mat) > 1) {
      collapsedMat <- do.call(cbind, subMat)
      } else collapsedMat <- unlist(subMat)
  } else if(is.vector(mat) || is.factor(mat) || is.array(mat)) {
    subMat <- tapply(seq_along(mat), groupArray, function(indices)
      mat[indices] %>% method)
    collapsedMat <- unlist(subMat)
  }

  return(collapsedMat)
}

collapseBarcodes <- function(Barbie, groupArray) {
  collapsedObject <- list(
    assay = Barbie$assay %>%
      collapseRowsByArray(groupArray=groupArray) %>% as.data.frame(),
    metadata = Barbie$metadata,
    proportion = Barbie$proportion %>%
      collapseRowsByArray(groupArray=groupArray),
    CPM = Barbie$CPM %>%
      collapseRowsByArray(groupArray=groupArray),
    occurrence = Barbie$occurrence %>%
      collapseRowsByArray(groupArray=groupArray),
    rank = Barbie$rank %>%
      collapseRowsByArray(groupArray=groupArray, methods=min),
    isTop = list(vec = Barbie$isTop$vec,
                 mat = Barbie$isTop$mat %>%
                   collapseRowsByArray(groupArray=groupArray, method=any)),
    clusters = Barbie$clusters,
    factorColors = Barbie$factorColors)

  ## include existing elements in the collapsedObject
  for (elementName in names(Barbie)) {
    if (!(elementName %in% names(collapsedObject)))
      collapsedObject[[elementName]] <- Barbie[[elementName]] %>%
        collapseRowsByArray(groupArray = groupArray, methods = function(x) x[[1]])
  }

  return(collapsedObject)
}

collapseRowsByArray <- function(mat, groupArray, method=max) {
  if (is.matrix(mat) || is.data.frame(mat)) {
    subMat <- tapply(seq_len(nrow(mat)), groupArray, function(indices) {
      apply(mat[indices, , drop = F], 2, method)
    })
    if(ncol(mat) > 1) {
      collapsedMat <- do.call(rbind, subMat)
    } else collapsedMat <- unlist(subMat)
  } else if (is.vector(mat)) {
    subMat <- tapply(seq_along(mat), groupArray, function(indices)
      mat[indices] %>% method)
    collapsedMat <- subMat
  }

  return(collapsedMat)
}

