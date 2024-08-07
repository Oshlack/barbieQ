#' Building a central data structure-based object similar to Seurat object.

#' Create a Barbie object
#'
#' @param object A data frame or matrix containing barcode counts.
#' @param target Optional parameter for experimental information.
#' @param color_panel Optional color panel for the Barbie object.
#' @param ... Additional arguments.
#' @return A list with the Barbie object structure.

createBarbie <- function(object, target = NULL, color_panel = NULL, ...) {
# define the object structure
  Barbie <- list(
    assay = NULL, #matrix. barcodes in rows, samples in columns
    metadata = NULL, #d.f. samples in rows, factors in columns. like targets file.
    proportion = matrix(),
    CPM = matrix(),
    occurrence = matrix(),
    rank = matrix(), #d.f. barcodes rank in each sample
    is_top = logical(), #logical vector indicating top or bottom contributing barcodes
    clusters = factor(levels = character(), labels = character()),
    color_panel = NULL #list of factors, each contains colors for each condition.
  )

# Create a constructor function
# extract barcode counts from object
  if(inherits(object, "data.frame") || inherits(object, "matrix")) {
#   check if all columns are numeric, excluding NAs from the check
    ColumnIsNumeric <- vapply(object, function(x) is.numeric(na.omit(x)), FUN.VALUE = logical(1L))
    if(all(ColumnIsNumeric)) Barbie$assay <- as.matrix(object)
    else {
      WhichNotNumeric <- which(!ColumnIsNumeric)
#     if the first column is the only non-numeric column, treating it as barcode IDs
      if(identical(sum(WhichNotNumeric), 1L) && length(ColumnIsNumeric) > 1L) {
#       only update row(barcode) names when it's NULL
        if(is.null(rownames(object))) {
          message("Missing Barcode IDs and detecting first column non-numeric, so converting first column as Barcode IDs.")
#         extract first non-numeric column
          firstCol <- object[,1,drop = T]

#         check if NAs exist in the first column
          if(any(is.na(firstCol))) {
            na_indices <- which(is.na(firstCol))
            firstCol[na_indices] <- paste0("Add_Barcode_", seq_along(na_indices))
            message("NAs exist in Barcode IDs, replaced by new names.")
          }
#         check if first column values (Barcode IDs) are unique
          if(any(duplicated(firstCol))) {
            firstCol <- make.unique(as.character(firstCol))
            message("duplicated Barcode IDs exist, making them unique.")
          }
#         update object
          Barbie$assay <- as.matrix(object[,-1,drop = FALSE])
          rownames(Barbie$assay) <- firstCol
          message("treating first column as barcode IDs; column number -1")

        } else {
          stop("cannot rename Barcodes by first column, as object rownames already exist.")
        }
      } else {
        stop("object should be numeric, instead it is a data.frame with ", length(WhichNotNumeric), " non-numeric colunms.")
      }
    }
  } else if(is.vector(object))
    stop("please convert object into a matrix, make sure Barcodes in rows and samples in columns.")
  else{
###$
    # extract barcode counts from objects. Need a function for it.
    # if Seurat, extract assay and metadata.
    # if BARtab object, extract ...
    # if Matrix object, ...
    # if Barbie itself, ...[get metadata!!!]
  }
  if(!nrow(Barbie$assay)) stop("Barcode count matrix has zero rows.")
  if(anyNA(Barbie$assay)) warning("Barcode count includes NAs.")

# extract metadata from target
# check if @param target is not provided - then inherit metadata from @param object such as Barbie objects.
  if(is.null(target)) target <- Barbie$metadata
# check if @param object does not provide metadata, then make up homogeneous settings.
  if(is.null(target)) {
    target <- as.data.frame(matrix(1, ncol(Barbie$assay), 1))
    message("no target files provided. making a homogeneous setting.")
  } else {
#   check @param target format
    if(is.vector(target)) target <- matrix(target, ncol = 1)
    else if(!(inherits(target, "data.frame") || inherits(target, "matrix")))
      stop("target should be a vector, matrix or data.frame of sample conditions.")
  }

# check if dimensions of assay and metadata match
  if(identical(ncol(Barbie$assay), nrow(target)))
    Barbie$metadata <- as.data.frame (target)
  else {
    stop("row dimension of target doesn;t match column dimension of Barcode count.")
  }
  if(anyNA(Barbie$assay)) warning("Target file (metadata) includes NAs.")

# now assay (Barcode count) should be a matrix, check if it's numeric apart from possible NAs.
# if any NAs, return Barbie object without further processing and with warning.
  ColumnHasNA <- vapply(as.data.frame(Barbie$assay), function(x) any(is.na(x)), FUN.VALUE = logical(1L))
  WhichColumnHasNA <- which(ColumnHasNA)
  RowHasNA <- vapply(seq_len(nrow(Barbie$assay)),
                     function(i) any(is.na(Barbie$assay[i, ])),
                     FUN.VALUE = logical(1L))
  WhichRowHasNA <- which(RowHasNA)
  if(any(ColumnHasNA)) {
    warning("In Barcode count matrix, ", length(WhichRowHasNA), " Barcodes (rows) show NAs across ", length(WhichColumnHasNA), " samples (columns).")
#   Display a menu for the user to choose
    choice <- menu(c("yes, continue with zeros", "no, keep NA and continue", "stop"), title = "Do you want to convert all NA into zero and continue?")
    if(choice == 1) {
#     replace NAs by zeros.
      Barbie$assay[is.na(Barbie$assay)] <- 0
      message("all NAs in Barcode count table are converted into zero.")
    } else if(choice == 2) {
      message("NAs may cause issues in the following analysis.")
    } else {
      stop("please deal with the NAs in the Barcode count table first.")
    }
  }

# confirm Barcode count table is numeric apart from NAs.
  ColumnIsNumeric <- vapply(as.data.frame(Barbie$assay), function(x) is.numeric(na.omit(x)), FUN.VALUE = logical(1L))
  if(!all(ColumnIsNumeric))
    stop("The assay, Barcode count table, should be a numeric matrix.")

# extract colors from @param color_panel or existing objects.
  if(is.null(color_panel)) color_panel <- Barbie$color_panel
  if(is.null(color_panel)) {
###$   dispatch a color picker function
    message("@param color_panel not provided. generating a new list of colors corresponding to metadata.")
  } else {
#   check if all factors in metadata has a corresponding color panel.
    factors <- colnames(Barbie$metadata)
    WhichFactorNoColor <- which(!(factors %in% names(colors)))
    if(length(WhichFactorNoColor) > 0L) {
###$     dispatch a color picker function
      message(length(WhichFactorNoColor), " factors in metadata doesn't have equivalent color panel, generating new colors for them.")
    }
  }
# it's OK if a color panel doesn't cover all conditions within the corresponding factor.

# calculate proportion, CPM, occurrence, rank ...
# use na.rm = TRUE in colSums to handle NA values
  col_sums <- colSums(Barbie$assay, na.rm = TRUE)
# avoid division by zero by setting zero columns to NA
  col_sums[col_sums == 0] <- NA

  Barbie$proportion <- (t(Barbie$assay) / col_sums) |> t()
  Barbie$CPM <- Barbie$proportion * 1e6
# calculate default Barcode occurrence
  Barbie$occurrence <- (Barbie$CPM >= 1) +1 -1
# calculate rank of Barcodes in each column
# NA will be the last ones in the order. letting na.last = "keep" will keep NAs as NA instead of a rank.
  Barbie$rank <- apply(Barbie$assay, 2, function(x) {rank(-x, ties.method = "first", na.last = "keep")})

# calculate is.top matrix and is.top vector
  Barbie$is_top <- getTopBar(Barbie = Barbie)

###$ are we updating clusters here? barcode clusaters? sample clusters?

###$ are we predicting occurrence against background noise? based on the bottom Barcode counts.


  return(Barbie)
}

# define a function for flagging top rows in a column
topsInVec <- function(temp_col, cut_off=NULL){

###$ think of a more robust cutoff instead of am empirical .99
# set up default cutoff at 0.99
  if(is.null(cut_off)){
    cut_off = 0.99
  }

  #   rank NAs as in the bottom of the vector
  temp_rank <- rank(-temp_col, ties.method = "first", na.last = TRUE) #rank the column with equal elements separated, by deceasing order
  #   sort NAs in the last
  temp_sort <- sort(temp_col, decreasing = TRUE, na.last = TRUE) #sort the column by decreasing order
  #   NAs will get NA in the cumulative percentile
  temp_cumsum <- sapply(temp_rank, function(x) sum(temp_sort[1:x])) #calculating the cumulative contribution
  #   NAs will stay NA after assessment
  temp_top <- temp_cumsum <= cut_off * sum(na.omit(temp_col)) #flag the top 95% of total counts

  return(temp_top)
}

## filterTopBar
getTopBar <- function(Barbie, column = NULL, cut_off = 0.99, num_appearance = 1){
  x <- Barbie$assay

###$ check x is a numeric matrix.


# column: flag the testing columns for selecting top rows.
  if(is.null(column)){
    column <- rep(TRUE, ncol(x))
  }

  topsInMat <- apply(x, 2, function(j) topsInVec(j, cut_off)) # a matrix of flagging top rows in each column

  subTopMat <- topsInMat[,column,drop = FALSE] # selecting columns

  top_flag <- rowSums(subTopMat, na.rm = TRUE) >= num_appearance # flag the rows defined as top at least once in the testing columns

  updatedObject <- list(
    assay = Barbie$assay,
    metadata = Barbie$metadata,
    proportion = Barbie$proportion,
    CPM = Barbie$CPM,
    occurrence = Barbie$occurrence,
    rank = Barbie$rank,
#   update!
    is_top = top_flag,
    topsInMat = topsInMat,
#
    clusters = Barbie$clusters,
    color_panel = Barbie$color_panel
    # Add other components as needed
  )

  # Include existing elements in the updatedObject
  for (element_name in names(Barbie)) {
    if (!(element_name %in% names(updatedObject))) {
      updatedObject[[element_name]] <- Barbie[[element_name]]
    }
  }

  return(updatedObject)
}

## trim data
trimObjectByMetadata <- function(Barbie, condition, specified, keep = TRUE) {
  if (!condition %in% colnames(Barbie$metadata)) {
    stop("Specified metadata column not found.")
  }

  if (!all(specified %in% Barbie$metadata[[condition]])) {
    stop("Specified value in metadata condition not found.")
  }

  keep_columns <- Barbie$metadata[[condition]] %in% specified
  if (!keep) {
    keep_columns <- !keep_columns
  }

  trimmedObject <- list(
    assay = Barbie$assay[,keep_columns,drop=FALSE],
    metadata = Barbie$metadata[keep_columns,,drop=FALSE],
    prop = Barbie$proportion[,keep_columns,drop=FALSE],
    CPM = Barbie$CPM[,keep_columns,drop=FALSE],
    occurrence = Barbie$occurrence[,keep_columns,drop=FALSE],
    rank = Barbie$rank[,keep_columns,drop=FALSE],
    is_top = Barbie$is_top,
    clusters = Barbie$clusters[keep_columns,drop=FALSE],
    color_panel = Barbie$color_panel
    # Add other components as needed
  )

  # Include existing elements in the trimmedObject
  for (element_name in names(Barbie)) {
    if (!(element_name %in% names(trimmedObject))) {
      trimmedObject[[element_name]] <- Barbie[[element_name]]
    }
  }

  return(trimmedObject)
}

trimColumn <- function(Barbie, keep_columns) {

  trimmedObject <- list(
    assay = Barbie$assay[,keep_columns,drop=FALSE],
    metadata = Barbie$metadata[keep_columns,,drop=FALSE],
    proportion = Barbie$proportion[,keep_columns,drop=FALSE],
    CPM = Barbie$CPM[,keep_columns,drop=FALSE],
    occurrence = Barbie$occurrence[,keep_columns,drop=FALSE],
    rank = Barbie$rank[,keep_columns,drop=FALSE],
    is_top = Barbie$is_top,
    clusters = Barbie$clusters[keep_columns,drop=FALSE],
    color_panel = Barbie$color_panel
    # Add other components as needed
  )

  # Include existing elements in the trimmedObject
  for (element_name in names(Barbie)) {
    if (!(element_name %in% names(trimmedObject))) {
      trimmedObject[[element_name]] <- Barbie[[element_name]]
    }
  }

  return(trimmedObject)
}

reorderColumn <- function(Barbie, order_array) {
  orderedObject <- list(
    assay = Barbie$assay[,order_array],
    metadata = Barbie$metadata[order_array,],
    proportion = Barbie$proportion[,order_array],
    CPM = Barbie$CPM[,order_array],
    occurrence = Barbie$occurrence[,order_array],
    rank = Barbie$rank[,order_array],
    is_top = Barbie$is_top,
    clusters = Barbie$clusters[order_array],
    color_panel = Barbie$color_panel
    # Add other components as needed
  )

  # Include existing elements in the orderedObject
  for (element_name in names(Barbie)) {
    if (!(element_name %in% names(orderedObject))) {
      orderedObject[[element_name]] <- Barbie[[element_name]]
    }
  }

  return(orderedObject)
}

combineColumn <- function(Barbie_1, Barbie_2) {
  combinedObject <- list(
    assay = cbind(Barbie_1$assay, Barbie_2$assay) %>% as.data.frame(),
    metadata = rbind(Barbie_1$metadata, Barbie_2$metadata) %>% as.data.frame(),
    proportion = cbind(Barbie_1$proportion, Barbie_2$proportion),
    CPM = cbind(Barbie_1$CPM, Barbie_2$CPM),
    occurrence = cbind(Barbie_1$occurrence, Barbie_2$occurrence),
    rank = cbind(Barbie_1$rank, Barbie_2$rank),
    is_top = Barbie_1$is_top | Barbie_2$is_top,
    clusters = cbind(Barbie_1$clusters, Barbie_2$clusters),
    color_panel = Barbie$color_panel
    # Add other components as needed
  )

  # Include existing elements in the combinedObject
  for (element_name in names(Barbie_1)) {
    if (!(element_name %in% names(combinedObject))) {
      combinedObject[[element_name]] <- Barbie_1[[element_name]]
    }
  }
  for (element_name in names(Barbie_2)) {
    if (!(element_name %in% names(combinedObject))) {
      combinedObject[[element_name]] <- Barbie_2[[element_name]]
    }
  }

  return(combinedObject)
}

trimRow <- function(Barbie, keep_rows) {

  trimmedObject <- list(
    assay = Barbie$assay[keep_rows,,drop=FALSE],
    metadata = Barbie$metadata,
    proportion = Barbie$proportion[keep_rows,,drop=FALSE],
    CPM = Barbie$CPM[keep_rows,,drop=FALSE],
    occurrence = Barbie$occurrence[keep_rows,,drop=FALSE],
    rank = Barbie$rank[keep_rows,,drop=FALSE],
    is_top = Barbie$is_top[keep_rows,drop=FALSE],s
    clusters = Barbie$clusters,
    color_panel = Barbie$color_panel
    # Add other components as needed
  )

  # Include existing elements in the trimmedObject
  for (element_name in names(Barbie)) {

    if (!(element_name %in% names(trimmedObject))) {

      if(is.data.frame(Barbie[[element_name]])) {
        trimmedObject[[element_name]] <- Barbie[[element_name]][keep_rows,]
      }

      if(is.vector(Barbie[[element_name]])) {
        trimmedObject[[element_name]] <- Barbie[[element_name]][keep_rows]
      }
    }
  }

  return(trimmedObject)
}

# Collapsing Rows by Array

collapse_rows_by_array <- function(mat, group_array, methods = max) {
  # get the max of each column in each submatrix (clone) decided by group array
  if (is.matrix(mat) | is.data.frame(mat)) { # mat is a matrix
    c_mat <- tapply(seq_len(nrow(mat)), group_array, function(indices) {
      apply(mat[indices, , drop = F], 2, methods)
    })

    # bind max count of each clone
    collapsed_mat <- do.call(rbind, c_mat)
  } else if (is.vector(mat)) { # mat is a vector
    c_mat <- tapply(seq_along(mat), group_array, function(indices) {
      mat[indices] %>% methods
    })

    collapsed_mat <- c_mat
  }

  return(collapsed_mat)
}

CollapseRow <- function(Barbie, group_array) {

  collapsedObject <- list(
    assay = Barbie$assay %>% collapse_rows_by_array(group_array = group_array) %>% as.data.frame(),
    metadata = Barbie$metadata,
    proportion = Barbie$proportion %>% collapse_rows_by_array(group_array = group_array),
    CPM = Barbie$CPM %>% collapse_rows_by_array(group_array = group_array),
    occurrence = Barbie$occurrence %>% collapse_rows_by_array(group_array = group_array),
    rank = Barbie$rank %>% collapse_rows_by_array(group_array = group_array, methods = min),
    is_top = Barbie$is_top %>% collapse_rows_by_array(group_array = group_array, methods = any),
    clusters = Barbie$clusters,
    color_panel = Barbie$color_panel
    # Add other components as needed
  )

  # Include existing elements in the collapsedObject
  for (element_name in names(Barbie)) {
    if (!(element_name %in% names(collapsedObject))) {
      collapsedObject[[element_name]] <- Barbie[[element_name]] %>%
        collapse_rows_by_array(group_array = group_array, methods = function(x) x[[1]])
    }
  }

  return(collapsedObject)
}

# Collapsing Columns by Array
collapse_cols_by_array <- function(mat, group_array, methods = mean) {
  if (is.matrix(mat) | is.data.frame(mat)) { # mat is a matrix
    c_mat <- tapply(seq_len(ncol(mat)), group_array, function(indices) {
      apply(mat[, indices, drop = F], 1, methods)
    })

    collapsed_mat <- do.call(cbind, c_mat)
  } else if(is.vector(mat) | is.factor(mat) | is.array(mat)) { # mat is a vector
    c_mat <- tapply(seq_along(mat), group_array, function(indices) {
      mat[indices] %>% methods
    })

    collapsed_mat <- c_mat
  }

  return(collapsed_mat)
}

CollapseColumn <- function(Barbie, group_array) {
  collapsedObject <- list(
    assay = Barbie$assay %>% collapse_cols_by_array(group_array = group_array) %>% as.data.frame(),
    metadata = Barbie$metadata %>% t() %>%
      collapse_cols_by_array(group_array = group_array, methods = function(x) {
        check_dup <- duplicated(x)
        ifelse(length(check_dup) -1 == sum(check_dup),
               x[[1]],
               paste(x, collapse = "."))
      }) %>% t() %>%
      as.data.frame(),
    proportion = Barbie$proportion %>% collapse_cols_by_array(group_array = group_array),
    CPM = Barbie$CPM %>% collapse_cols_by_array(group_array = group_array),
    occurrence = Barbie$occurrence %>% collapse_cols_by_array(group_array = group_array, methods = max),
    rank = Barbie$rank %>% collapse_cols_by_array(group_array = group_array, methods = min),
    is_top = Barbie$is_top,
    clusters = Barbie$clusters %>% collapse_cols_by_array(group_array = group_array, methods = function(x) x[[1]]),
    color_panel = Barbie$color_panel
    # Add other components as needed
  )

  for (element_name in names(Barbie)) {
    if (!(element_name %in% names(collapsedObject))) {
      collapsedObject[[element_name]] <- Barbie[[element_name]]
    }
  }

  return(collapsedObject)
}

# add elements in Barbie
addToBarbie <- function(Barbie, new_element_name, new_values) {
  # Check if the new element already exists in Barbie
  if (new_element_name %in% names(Barbie)) {
    stop("Element already exists. Choose a different name.")
  }

  # Add the new element to Barbie
  Barbie[[new_element_name]] <- new_values
  return(Barbie)
} #example usage: Barbie <- addToBarbie(Barbie, "new_element", "new_values")

# add color gradient panel

GetColorGradient <- function(color.ls) {
  togradient <- function(maxcolors, mincolor = "#FFFFFF"){ #maxcolors is a vector containing elements
    #mincolor = "#FFFFFF"
    #maxcolors <- color.ls$lineage
    grad_ls <- list()
    for (i in 1:length(maxcolors)) {
      grad_ls[[i]] <- colorRamp2(c(0, 1), c(mincolor, maxcolors[i]))
    }
    names(grad_ls) <- names(maxcolors)
    return(grad_ls)
  }

  color.gr.ls <- lapply(color.ls, togradient)

  return(color.gr.ls)
}

#UpdateColor
UpdateColor <- function(Barbie, condition, new_color_ls){
  for(i in seq_along(condition)){
    Barbie$color_panel[[condition[i]]] <- new_color_ls[[condition[i]]]
  }

  return(Barbie)
}

# Step4: update methods

# Step5: create accessor functions

