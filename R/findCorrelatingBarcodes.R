# get mean total counts of each pair
cross.mean <- function(vec) {
  size <- length(vec)
  mean_mx <- matrix(nrow = size, ncol = size)
  rownames(mean_mx) <- names(vec)
  colnames(mean_mx) <- names(vec)

  for (i in 1:size) {
    for (j in 1:size) {
      mean_mx[i,j] <- mean(vec[i], vec[j])
    }
  }

  return(mean_mx)
}

cross.name <- function(vec) {
  size <- length(vec)
  name_mx <- matrix(nrow = size, ncol = size)
  rownames(name_mx) <- names(vec)
  colnames(name_mx) <- names(vec)

  for (i in 1:size) {
    for (j in 1:size) {
      name_mx[i,j] <- paste(vec[i], vec[j], sep = ".")
    }
  }

  return(name_mx)
}

cross.cor.test <- function(mx) {
  size <- ncol(mx)
  test_mx <- matrix(nrow = size, ncol = size)
  rownames(test_mx) <- names(mx)
  colnames(test_mx) <- names(mx)

  for (i in 1:size) {
    for (j in 1:size) {
      test_mx[i,j] <- cor.test(mx[,i], mx[,j], method = "pearson")$p.value
    }
  }

  return(test_mx)
}

# separate correlation within / between clones
# flag the inner clone correlations by group vector
defined.diag <- function(vec) {
  size <- length(vec)
  flag_mx <- matrix(nrow = size, ncol = size)
  rownames(flag_mx) <- vec
  colnames(flag_mx) <- vec

  for (i in 1:size) {
    for (j in 1:size) {
      ifelse(vec[i] == vec[j],
             flag_mx[i,j] <- TRUE,
             flag_mx[i,j] <- FALSE)
    }
  }

  return(flag_mx)
}

getLongFromMatrix <- function(mx, group_mx){
  inner <- mx[defined.diag(group_mx) & lower.tri(mx, diag = F)]
  inter <- mx[!defined.diag(group_mx) & lower.tri(mx, diag = F)]

  inner.df <- data.frame(Value = inner, group = rep("within-clone barcode pair", length = length(inner)))
  inter.df <- data.frame(Value = inter, group = rep("undefined barcode pair", length = length(inter)))

  total.df <- rbind(inner.df, inter.df)

  return(total.df)
}


plotBarcodePairCor_Count <- function(count.mx, group.vec = NULL, plot = "mean") {
  #read barcode count, make sure barcodes in columns

  if(is.null(group.vec)){
    group.vec <- seq(nrow(count.mx))
  }

  count.mx <- t(count.mx)

  count.mx <- count.mx[,order(group.vec)]
  group.vec <- group.vec[order(group.vec)]

  bar_cor <- cor(count.mx) # correlation matrix

  bar_cor.test <- cross.cor.test(count.mx) # correlation test matrix
  bar_cor.test <- bar_cor.test < 0.05 # reject null hypo (0 cor), accept correlated

  bar_Cmean <- colMeans(count.mx) # get average barcode count across samples
  names(bar_Cmean) <- colnames(count.mx)

  bar_Cmax <- apply(count.mx, 2, max) # get the max of barcode count across sample
  names(bar_Cmax) <- colnames(bar_Cmax)

  bar_mean <- cross.mean(bar_Cmean) # mean cross matrix of each pair

  bar_max <- cross.mean(bar_Cmax) # max cross matrix of each pair

  bar_barID <- cross.name(names(bar_Cmean)) # cross matrix of barcode ID

  bar_cloneID <- cross.name(group.vec) # cross matrix of clone ID

  #convert all matrix into long vector
  all_bar_cor <- getLongFromMatrix(mx = bar_cor, group_mx = group.vec)

  all_bar_mean <- getLongFromMatrix(mx = bar_mean, group_mx = group.vec)

  all_bar_max <- getLongFromMatrix(mx = bar_max, group_mx = group.vec)

  all_bar_barID <- getLongFromMatrix(mx = bar_barID, group_mx = group.vec)

  all_bar_cloneID <- getLongFromMatrix(mx = bar_cloneID, group_mx = group.vec)

  all_bar_cor_sig <- getLongFromMatrix(mx = bar_cor.test, group_mx = group.vec)

  all_bar <- data.frame(group = all_bar_barID$group,
                        barID = all_bar_barID$Value,
                        cloneID = all_bar_cloneID$Value,
                        Cor = all_bar_cor$Value,
                        Mean = all_bar_mean$Value,
                        Max = all_bar_max$Value,
                        sig. = all_bar_cor_sig$Value)

  # plot mean vs Cor.
  meanCor <- ggplot(all_bar, aes(x = Cor, y = log2(Mean+1), color = group, text = paste(cloneID, barID, sep = " : "))) +
    geom_point(shape = 1)

  # plot max vs Cor.
  maxCor <- ggplot(all_bar, aes(x = Cor, y = log2(Max+1), color = group, text = paste(cloneID, barID, sep = " : "))) +
    geom_point(shape = 1)

  # plot correlation hgistogram
  hisCor <- ggplot(all_bar, aes(x = Cor, fill = group)) +
    geom_histogram()

  return(
    if(plot == "mean") {
      meanCor
    } else if(plot == "max") {
      maxCor
    } else if(plot == "histogram") {
      hisCor
    } else if(plot == "table") {
      all_bar
    } else {
      warning("Invalid 'plot' argument. Returning NULL.")
      NULL
    }
  )
}

# Create new barcode-clone list
createBarcodeCloneRef <- function(new_pairlist) {
  # old_reflist is the old barcode-clone reference list
  # new_pairlist is the new determined pairwise correlated barcode list

  total_list <- new_pairlist

  g <- igraph::graph_from_edgelist(do.call(rbind, total_list), directed = F)

  plot(g, vertex.label.cex=0.4)

  groups <- igraph::clusters(g)$membership

  groups_list <- split(names(groups), groups)

  return(groups_list)
}


# make new or update barcode-clone reference list
updateBarcodeCloneRef <- function(old_reflist, new_pairlist) {
  # old_reflist is the old barcode-clone reference list
  # new_pairlist is the new determined pairwise correlated barcode list

    # remove single points in ovbar
    ovbar2 <- old_reflist[lengths(old_reflist) >= 2]

    # convert groups into pairs
    ovbar_df <- lapply(names(ovbar2), function(x) {
      data.frame(Node = ovbar2[[x]],
                 Group = rep(x, length(ovbar2[[x]])))
    })

    ovbar_df2 <- do.call(rbind, ovbar_df)

    #function to create pairwise relationships
    create_pairwise <- function(group_data) {
      edges <- combn(group_data, 2)
      split_list <- lapply(seq_len(ncol(edges)), function(i) edges[,i])
    }

    pairwise_list <- tapply(ovbar_df2$Node, ovbar_df2$Group, create_pairwise)

    pairwise_list2 <- do.call(c, pairwise_list) # bind all the sublists

    total_list <- c(pairwise_list2, new_pairlist) # bind the list with pool

  g <- igraph::graph_from_edgelist(do.call(rbind, total_list), directed = F)

  plot(g, vertex.label.cex=0.4)

  groups <- igraph::clusters(g)$membership

  groups_list <- split(names(groups), groups)

  return(groups_list)
}

# for one clone in the newlist, find the overlapping clone in the old list
find_ref_cloneID <- function(newlisti, reflist) {

  flag <- sapply(reflist, function(x) any(newlisti %in% x))

  if(sum(flag) == 1) {
    find_name <- names(reflist)[flag]
  } else if(sum(flag == 0)) {
    find_name <- paste0("Clone_", -new_clone_num)
    new_clone_num <<- new_clone_num + 1
  } else {
    find_name <- paste(names(reflist)[flag], collapse = ".")
  }

  return(find_name)

}

# update cloneID for the newlist
updateCloneIDRef <- function(newlist, reflist) {

  names(newlist) <- sapply(newlist, function(x) {find_ref_cloneID(x, reflist)})

  print(paste0("next new_clone_num: ", new_clone_num))

  flag_left <- ! names(reflist) %in% names(newlist)

  update_list <- c(newlist, reflist[flag_left])

  return(update_list)
}

# for each barcode, find cloneID in the refernce list
match_cloneID <- function(barIDi, reflist) {
  flag <- sapply(reflist, function(x) barIDi %in% x)

  if(sum(flag) == 1) {
    find_name <- names(reflist)[flag]
  } else if(sum(flag == 0)) {
    find_name <- paste0("Clone_", -new_clone_num)
    new_clone_num <<- new_clone_num + 1
  } else {
    find_name <- paste(names(reflist)[flag], collapse = ".")
  }

  return(find_name)
}

# update cloneID in Barbie according to the ref list
updateCloneIDBarbie <- function(Barbie, reflist) {
  barID <- rownames(Barbie$assay)

  group_cloneID <- sapply(barID, function(x) {match_cloneID(x, reflist)})
  print(paste0("next new_clone_num: ", new_clone_num))

  return(group_cloneID)
}
