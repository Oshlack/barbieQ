# plot heatmap of CPM or pre
plotBarbieHeatmap <- function(Barbie, showTest = FALSE, rowHA = NULL, splitBy = NULL,
                        groupBy = NULL, biasDirection = NULL) {

  column_ha = HeatmapAnnotation(Sample_Group = groupBy,
                                annotation_name_side = "left",
                                annotation_name_gp = gpar(fontsize = 10),
                                col = list(Sample_Group = Barbie$color_panel$biasDirection)
  )

  mat <- log2(Barbie$CPM +1) %>% as.matrix()
  mat_name <- "logCPM"
  col_fun <- circlize::colorRamp2(c(min(mat), median(mat), max(mat)), c("blue", "white", "red"))

  if(is.null(biasDirection)){
    biasDirection <- Barbie$Bias_Occ$group
  }

  if(!showTest) {
    rowHA <- NULL
  } else if (is.null(rowHA)) {
    rowHA <- rowAnnotation(Bias = biasDirection,
                            annotation_name_side = "top",
                            annotation_name_gp = gpar(fontsize = 10),
                            col = list(Bias = Barbie$color_panel$biasDirection),
                            show_legend = TRUE,
                            show_annotation_name = TRUE)
  }


  hp <- Heatmap(mat,
                name = mat_name, width = unit(6,"cm"), height = unit(6,"cm"),
                cluster_rows = TRUE, cluster_columns = TRUE,
                show_row_names = FALSE, show_column_names = FALSE,
                column_title = "Samples", row_title = paste(nrow(mat), " Top Clones"),
                # col = col_fun,
                right_annotation = rowHA,
                top_annotation = column_ha,
                column_split = splitBy,
                cluster_column_slices = FALSE
  )

  return(hp)
}

# plot presence
PlotPreHP_0 <- function(Barbie, showTest = FALSE, rowHA = NULL, splitBy = NULL,
                        groupBy = NULL, biasDirection = NULL) {

  column_ha = HeatmapAnnotation(Sample_Group = groupBy,
                                annotation_name_side = "left",
                                annotation_name_gp = gpar(fontsize = 10),
                                col = list(Sample_Group = Barbie$color_panel$biasDirection)
  )


  mat <- (Barbie$presence +1 -1) %>% as.matrix()
  mat_name <- "Presence"
  col_fun <- structure(c(2,4), names = c("1","0"))

  if(is.null(biasDirection)){
    biasDirection <- Barbie$Bias_Occ$group
  }

  if(!showTest) {
    rowHA <- NULL
  } else if (is.null(rowHA)) {
    rowHA <- rowAnnotation(Bias = biasDirection,
                            annotation_name_side = "top",
                            annotation_name_gp = gpar(fontsize = 10),
                            col = list(Bias = Barbie$color_panel$biasDirection),
                            show_legend = TRUE,
                            show_annotation_name = TRUE)
  }

  hp <- Heatmap(mat,
                name = mat_name, width = unit(6,"cm"), height = unit(6,"cm"),
                cluster_rows = TRUE, cluster_columns = TRUE,
                show_row_names = FALSE, show_column_names = FALSE,
                column_title = "Samples", row_title = paste(nrow(mat), " Top Clones"),
                col = col_fun,
                right_annotation = rowHA,
                column_split = splitBy,
                top_annotation = column_ha,
                cluster_column_slices = FALSE
  )

  return(hp)
}


