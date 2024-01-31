# Building a central data structure-based object similar to Seurat object.

# Step1: define the object structure
Barbie <- list(
  assay = NULL, #barcodes in rows, samples in columns
  metadata = data.frame(), #samples in rows, factors in columns. like targets file.
  prop = NULL,
  CPM = NULL,
  presence = NULL,
  rank = matrix(),
  is_top = logical(),
  clusters = factor(levels = character(), labels = character()),
  color_panel = list()
)

# Step2: create a constructor function
createBarbie <- function(counts, metadata, color_panel, ...) {
  prop <- (t(counts) / colSums(counts)) |> t()
  CPM <- prop * 1000000
  presence <- CPM >= 1 +1 -1

  Barbie <- list(
    assay = counts,
    metadata = metadata,
    prop = prop,
    CPM = CPM,
    presence = presence,
    rank = apply(counts, 2, function(x) {cbind(rank(-x, ties.method = "first"))}),
    is_top = logical(length = nrow(counts)),
    clusters = factor(colnames(counts)),
    color_panel = color_panel,
    ...
  )
  return(Barbie)
}

# Step3: add methods
## filterTopBar
getTopBar <- function(Barbie, column = NULL){
  x <- Barbie$assay

  # define a function for flagging top rows in each column
  tops <- function(temp_col, cutoff_top = 0.99){

    if(is.null(cutoff_top)){
      cutoff_top = 0.99 #set up default cutoff at 0.99
    }

    temp_rank <- rank(-temp_col, ties.method = "first") #rank the column with equal elements separated, by deceasing order

    temp_sort <- sort(temp_col, decreasing = TRUE) #sort the column by decreasing order

    temp_cumsum <- sapply(temp_rank, function(x) sum(temp_sort[1:x])) #calculating the cumulative contribution

    temp_top <- temp_cumsum <= cutoff_top * sum(temp_col) #flag the top 95% of total counts

    return(temp_top)
  }

  # column: flag the testing columns for selecting top rows.
  if(is.null(column)){
    column <- rep(TRUE, ncol(x))
  }

  top_mx <- apply(x, 2, tops) # a matrix of flagging top rows in each column

  top_flag <- rowSums(top_mx[, column]) > 0 # flag the rows defined as top at least once in the testing columns

  updatedObject <- list(
    assay = Barbie$assay,
    metadata = Barbie$metadata,
    prop = Barbie$prop,
    CPM = Barbie$CPM,
    presence = Barbie$presence,
    rank = Barbie$rank,
    is_top = top_flag,
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
    assay = Barbie$assay[,keep_columns],
    metadata = Barbie$metadata[keep_columns,],
    prop = Barbie$prop[,keep_columns],
    CPM = Barbie$CPM[,keep_columns],
    presence = Barbie$presence[,keep_columns],
    rank = Barbie$rank[,keep_columns],
    is_top = Barbie$is_top,
    clusters = Barbie$clusters[keep_columns],
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
    assay = Barbie$assay[,keep_columns],
    metadata = Barbie$metadata[keep_columns,],
    prop = Barbie$prop[,keep_columns],
    CPM = Barbie$CPM[,keep_columns],
    presence = Barbie$presence[,keep_columns],
    rank = Barbie$rank[,keep_columns],
    is_top = Barbie$is_top,
    clusters = Barbie$clusters[keep_columns],
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

trimRow <- function(Barbie, keep_rows) {

  trimmedObject <- list(
    assay = Barbie$assay[keep_rows,],
    metadata = Barbie$metadata,
    prop = Barbie$prop[keep_rows,],
    CPM = Barbie$CPM[keep_rows,],
    presence = Barbie$presence[keep_rows,],
    rank = Barbie$rank[keep_rows,],
    is_top = Barbie$is_top[keep_rows],
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
    prop = Barbie$prop %>% collapse_rows_by_array(group_array = group_array),
    CPM = Barbie$CPM %>% collapse_rows_by_array(group_array = group_array),
    presence = Barbie$presence %>% collapse_rows_by_array(group_array = group_array),
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
  } else if(is.vector(mat) | is.factor(mat)) { # mat is a vector
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
        if(length(check_dup) -1 == sum(check_dup))
          paste(x, collapse = ".")
      })
         %>% t() %>%
      as.data.frame(),
    prop = Barbie$prop %>% collapse_cols_by_array(group_array = group_array),
    CPM = Barbie$CPM %>% collapse_cols_by_array(group_array = group_array),
    presence = Barbie$presence %>% collapse_cols_by_array(group_array = group_array, methods = max),
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

