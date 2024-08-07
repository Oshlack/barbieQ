# Barbie <- week4_IVIF
# Factor <- "lineage"
# Vector_customized <- Barbie$metadata$lineage
# Vector_customized[Vector_customized %in% c("Bcell", "Tcell")] <- "Lympoid"

GetPotential <- function(Barbie, Factor = "lineage", Vector_customized = NULL, method = "hgeom"){
  #Vector_customized is a categorical vector with length of columns (sample number).

  mx <- Barbie$presence
  facs <- Barbie$metadata[[Factor]]

  if(!is.null(Vector_customized)) {
    facs <- Vector_customized
  }

  cdt <- table(facs) # extract all conditions in analysis, names(cdt) is the condition names
  m_n <- sum(cdt) # total condition number

  #in hypergeometric distribution, k is the total outcome
  k <- rowSums(mx)
  #get total non-detected numbers
  k_non <- rowSums(!mx)

  #k_x, k_y, ... , is the outcome of each condition
  ki <- lapply(names(cdt), function(x){rowSums(mx[,facs == x])})
  names(ki) <- names(cdt)
  #non_detected outcomes
  ki_non <- lapply(names(cdt), function(x){rowSums(!mx[,facs == x])})
  names(ki_non) <- names(cdt)

  #hypergeometric distribution
  bias.ls <- lapply(names(cdt), function(x){
    ifelse(k>0,
           phyper(q = ki[[x]],
                  m = cdt[[x]],
                  n = m_n - cdt[[x]],
                  k = k),
           0)
  })
  names(bias.ls) <- names(cdt)

  #simple probability of outcome / condition
  prob.ls <- lapply(names(cdt), function(x){
    ifelse(k>0,
           ki[[x]] / cdt[[x]],
           0)
    })
  names(prob.ls) <- names(cdt)

  #convert list into matrix
  bias.mx <- do.call(rbind, bias.ls) |> t() |> as.data.frame()
  prob.mx <- do.call(rbind, prob.ls) |> t() |> as.data.frame()

  if(method == "hgeom"){
    return(bias.mx)
  }

  if(method == "freq"){
    return(prob.mx)
  }

}

# get contingency table from the presence data

GetContingencyTable <- function(Barbie, Factor = "lineage", Vector_customized = NULL){
  #Vector_customized is a categorical vector with length of columns (sample number).

  mx <- Barbie$presence
  facs <- Barbie$metadata[[Factor]]

  if(!is.null(Vector_customized)) {
    facs <- Vector_customized
  }

  cdt <- table(facs) # extract all conditions in analysis, names(cdt) is the condition names
  m_n <- sum(cdt) # total condition number

  #in hypergeometric distribution, k is the total outcome
  k <- rowSums(mx)
  #get total non-detected numbers
  k_non <- rowSums(!mx)

  #k_x, k_y, ... , is the outcome of each condition
  ki <- lapply(names(cdt), function(x){rowSums(mx[,facs == x])})
  names(ki) <- names(cdt)
  #non_detected outcomes
  ki_non <- lapply(names(cdt), function(x){rowSums(!mx[,facs == x])})
  names(ki_non) <- names(cdt)

  #check if ki_non[[1]] + ki[[1]] == cdt[1]
  #check if ki_non[[2]] + ki[[2]] == cdt[2]

  #check if sum(ki) == k
  # all(ki[[1]] + ki[[2]] == k)
  # all(ki_non[[1]] + ki_non[[2]] == k_non)

  # contigincy table
  c_tables <- lapply(rownames(mx), function(x) {
    c_table <- matrix(nrow = length(cdt)+1, ncol = 3,
                      dimnames = list(c("num_presence", "num_absence", "total"),
                                      c(names(cdt), "total")),
                      byrow = TRUE)
    # Populate the matrix with data
    for (i in seq_along(cdt)) {
      c_table[, i] <- c(ki[[i]][x], ki_non[[i]][x], cdt[i])
    }
    c_table[1, length(cdt) + 1] <- k[x]
    c_table[2, length(cdt) + 1] <- k_non[x]
    c_table[3, length(cdt) + 1] <- sum(cdt)

    return(c_table)
  })

  names(c_tables) <- rownames(mx)

  return(c_tables)
}

# example usage:
# Vector_customized <- week4_IVIF$metadata$lineage
# Vector_customized[Vector_customized %in% c("Bcell", "Tcell")] <- "Lympoid"
# c_tables_week4 <- GetContingencyTable(week4_IVIF, Vector_customized = Vector_customized)

#Plot potential by heatmap
PotentialHP <- function(Barbie, potential, dataFormat = "CPM") {
  design.mx <- Barbie$metadata
  color.ls <- Barbie$color_panel
  color.gr.ls <- Barbie$color_panel_gradient

  column_ha = HeatmapAnnotation(Treat = design.mx$treat,
                                Mouse = design.mx$mouse,
                                Tissue = design.mx$tissue,
                                Lineage = design.mx$lineage,
                                Cell = design.mx$celltype,
                                annotation_name_side = "left",
                                annotation_name_gp = gpar(fontsize = 10),
                                col = list(Lineage = color.ls$lineage,
                                           Mouse = color.ls$mouse,
                                           Tissue = color.ls$tissue,
                                           Treat = color.ls$treat,
                                           Cell = color.ls$celltype
                                )
  )

  tpcol <- lapply(colnames(potential) , FUN = function(x){color.gr.ls$lineage[[x]]})
  names(tpcol) <- colnames(potential)

  row_ha = rowAnnotation(df = potential,
                         show_annotation_name = T,
                         col = tpcol
  )

  # heatmap
  if(dataFormat == "CPM") {
    bar_heat <- log2(as.matrix(Barbie$CPM) + 1)
    ht <- Heatmap(bar_heat, name = "log2(CPM+1)", width = unit(5,"cm"), height = unit(5,"cm"),
                  column_title = "Samples",
                  row_title = paste0(nrow(bar_heat), " Clones"),
                  column_title_gp = gpar(fontsize = 10), row_title_gp = gpar(fontzise = 10),
                  show_row_names = FALSE, show_column_names = FALSE,
                  row_labels = rownames(bar_heat), column_labels = colnames(bar_heat),
                  top_annotation = column_ha,
                  left_annotation = row_ha,
                  show_heatmap_legend = T
    )
  }

  col_fun = structure(c(2,4), names = c("1","0")) # black, red, green, blue

  if(dataFormat == "presence") {
    bar_heat <- Barbie$presence + 1 - 1
    ht <- Heatmap(bar_heat, name = "Presence", width = unit(5,"cm"), height = unit(5,"cm"),
                  column_title = "Samples",
                  row_title = paste0(nrow(bar_heat), " Clones"),
                  column_title_gp = gpar(fontsize = 10), row_title_gp = gpar(fontzise = 10),
                  show_row_names = FALSE, show_column_names = FALSE,
                  row_labels = rownames(bar_heat), column_labels = colnames(bar_heat),
                  top_annotation = column_ha,
                  left_annotation = row_ha,
                  show_heatmap_legend = T,
                  col = col_fun
    )
  }

  return(ht)
}

# Get Fisher test result
GetOccBiasGroup <- function(Barbie, contingency_table_ls) {

  fisher_table <- lapply(contingency_table_ls,  function(x) x[1:2, 1:2]) #remove the "total" column and "total" row.
  fisher_results <- lapply(fisher_table, function(x) fisher.test(x)) #default: two.sided

  # Extract p-values and log-odds ratios
  fisher_Pvalue <- sapply(fisher_results, function(x) x$p.value)
  fisher_Odds_Ratio <- sapply(fisher_results, function(x) x$estimate)

  fisher_Pvalue.gr <- lapply(fisher_table, function(x) fisher.test(x, alternative = "greater")$p.value) |> unlist()#default: two.sided

  fisher_Pvalue.le <- lapply(fisher_table, function(x) fisher.test(x, alternative = "less")$p.value) |> unlist()#default: two.sided

  #figure out bias
  # LymGrMye <- fisher_Pvalue.gr < 0.05
  # MyeGrLym <- fisher_Pvalue.le < 0.05
  #
  # Bias_Fisher <- data.frame(
  #   Lymphoid = LymGrMye,
  #   Myeloid = MyeGrLym,
  #   Unbiased = LymGrMye == FALSE & MyeGrLym == FALSE)

  Bias_Fisher <- data.frame(
    Lymphoid = (fisher_Pvalue < 0.05) & (fisher_Odds_Ratio > 1),
    Myeloid = (fisher_Pvalue < 0.05) & (fisher_Odds_Ratio < 1),
    Unbiased = fisher_Pvalue >= 0.05
  )

  con_name <- colnames(fisher_table[[1]])[1:2]
  colnames(Bias_Fisher)[1:2] <- con_name

  #check
  all(rowSums(Bias_Fisher) == 1)

  Bias_Fisher_group <- numeric(length = nrow(Bias_Fisher))
  for (i in 1:nrow(Bias_Fisher)) {
    tp <- as.matrix(Bias_Fisher)
    Bias_Fisher_group[i] <- colnames(tp)[tp[i,]]
  }

  Barbie$Bias_Occ <- data.frame(pvalue = fisher_Pvalue,
                            pvalue.LymGrMye = fisher_Pvalue.gr,
                            pvalue.MyeGrLym = fisher_Pvalue.le,
                            group = Bias_Fisher_group)

  if(is.null(Barbie$color_panel$bias_group)) {
    Barbie$color_panel$bias_group <- c("C1" = "#33AAFF", "C2" = "#FF5959", "Unbiased" = "#FFC000")
    names(Barbie$color_panel$bias_group)[1:2] <- con_name
  }

  return(Barbie)

}

GetFisherBiasGroup <- GetOccBiasGroup

# Plot clone as vs average Rank colored by bias group.

PlotBiasVsRank <- function(Barbie, passing_data = "rank", bias_group = NULL, bias_pvalue = NULL) {
  # Define a custom color/shape palette for groups
  custom_shape <- c("C1" = 21, "C2" = 24, "Unbiased" = 23)
  names(custom_shape) <- names(Barbie$color_panel$bias_group)

  output <- Barbie$presence |> rowSums()
  avgRank <- Barbie$rank |> rowMeans()

  # check out consistency of ranks across samples.
  rank_reorder <- apply(Barbie$rank, 2, rank)

  #transformation
  avgRank_reorder <- rowMeans(rank_reorder)
  varRank_reorder <- apply(rank_reorder, 1, var)

  df <- data.frame(# bias = -log10(Barbie$Bias_Fisher$pvalue),
                   output = output,
                   rank = avgRank,
                   rank_fixed = avgRank_reorder,
                   rank_fixed_var = varRank_reorder,
                   # group = Barbie$Bias_Occ$group,
                   clone = rownames(Barbie$CPM),
                   meanprop = log2(rowMeans((Barbie$CPM + 1)))
                   )

  # p_oneside_min <- apply(Barbie$Bias_Occ[, 2:3], 1, min)

  if(is.null(bias_group)) {
    df$group <- Barbie$Bias_Occ$group # default
    df$bias <- -log10(Barbie$Bias_Occ$pvalue) # p_oneside_min
  } else {
    df$group <- bias_group
    df$bias <- -log10(bias_pvalue)
  }

  if(passing_data == "output") {
    p <- ggplot(df, aes(x = output, y = bias, text = clone)) +
      geom_point(aes(color = group, shape = group, fill = group), size = 4, stroke = 1) +
      theme_classic() +
      theme(aspect.ratio = 1) +
      labs(title = "Lym vs Mye",
           y = "Bias = -log10(p.value)",
           x = "Output = number of samples in which clone is detected") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
      annotate("text", x = max(df$output)*0.9, y = -log10(0.05), label = "p.value = 0.05", vjust = 1.5, hjust = 1) +
      scale_color_manual(values = Barbie$color_panel$bias_group) +
      scale_shape_manual(values = custom_shape) +
      scale_fill_manual(values = alpha(Barbie$color_panel$bias_group, 0.2))
  } else if(passing_data == "rank") {
    p <- ggplot(df, aes(x = nrow(df) - rank_fixed, y = bias, text = clone)) +
      geom_point(aes(color = group, shape = group, fill = group), size = 4, stroke = 1) +
      theme_classic() +
      theme(aspect.ratio = 1) +
      labs(title = "Lym vs Mye",
           y = "Bias = -log10(p.value)",
           x = paste0("Genral Contribution = ", nrow(df), " - average rank of contributions across samples")) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
      annotate("text", x = max(nrow(df))*0.9, y = -log10(0.05), label = "p.value = 0.05", vjust = 1.5, hjust = 1) +
      scale_color_manual(values = Barbie$color_panel$bias_group) +
      scale_shape_manual(values = custom_shape) +
      scale_fill_manual(values = alpha(Barbie$color_panel$bias_group, 0.2))
  } else if(passing_data == "rank_var"){
    p <- ggplot(df, aes(x = nrow(df) - rank_fixed, y = rank_fixed_var, text = clone)) +
      geom_point(aes(color = group, shape = group, fill = group), size = 4, stroke = 1) +
      theme_classic() +
      theme(aspect.ratio = 1) +
      labs(title = "Mean-Variance of clone ranks",
           y = "Variance of rank",
           x = paste0(nrow(df), " - Mean of rank")) +
      scale_color_manual(values = Barbie$color_panel$bias_group) +
      scale_shape_manual(values = custom_shape) +
      scale_fill_manual(values = alpha(Barbie$color_panel$bias_group, 0.2))
  } else if(passing_data == "contribution") {
    p <- ggplot(df, aes(x = meanprop, y = bias, text = clone)) +
      geom_point(aes(color = group, shape = group, fill = group), size = 4, stroke = 1) +
      theme_classic() +
      theme(aspect.ratio = 1) +
      labs(title = "Lym vs Mye",
           y = "Bias = -log10(p.value)",
           x = "log Average Clone Contribution across Samples") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
      annotate("text", x = max(df$meanprop)*0.9, y = -log10(0.05), label = "p.value = 0.05", vjust = 1.5, hjust = 1) +
      scale_color_manual(values = Barbie$color_panel$bias_group) +
      scale_shape_manual(values = custom_shape) +
      scale_fill_manual(values = alpha(Barbie$color_panel$bias_group, 0.2))
  }

  return(p)

}
