# plot heatmap of CPM or pre
PlotCpmHP_0 <- function(Barbie, show_bias = FALSE, row_ha = NULL, split_by = NULL,
                        Vector_customized = NULL, bias_group = NULL) {

  column_ha = HeatmapAnnotation(Sample_Group = Vector_customized,
                                annotation_name_side = "left",
                                annotation_name_gp = gpar(fontsize = 10),
                                col = list(Sample_Group = Barbie$color_panel$bias_group)
  )

  mat <- log2(Barbie$CPM +1) %>% as.matrix()
  mat_name <- "logCPM"
  col_fun <- circlize::colorRamp2(c(min(mat), median(mat), max(mat)), c("blue", "white", "red"))

  if(is.null(bias_group)){
    bias_group <- Barbie$Bias_Occ$group
  }

  if(!show_bias) {
    row_ha <- NULL
  } else if (is.null(row_ha)) {
    row_ha <- rowAnnotation(Bias = bias_group,
                            annotation_name_side = "top",
                            annotation_name_gp = gpar(fontsize = 10),
                            col = list(Bias = Barbie$color_panel$bias_group),
                            show_legend = TRUE,
                            show_annotation_name = TRUE)
  }


  hp <- Heatmap(mat,
                name = mat_name, width = unit(6,"cm"), height = unit(6,"cm"),
                cluster_rows = TRUE, cluster_columns = TRUE,
                show_row_names = FALSE, show_column_names = FALSE,
                column_title = "Samples", row_title = paste(nrow(mat), " Top Clones"),
                # col = col_fun,
                right_annotation = row_ha,
                top_annotation = column_ha,
                column_split = split_by,
                cluster_column_slices = FALSE
  )

  return(hp)
}

# plot presence
PlotPreHP_0 <- function(Barbie, show_bias = FALSE, row_ha = NULL, split_by = NULL,
                        Vector_customized = NULL, bias_group = NULL) {

  column_ha = HeatmapAnnotation(Sample_Group = Vector_customized,
                                annotation_name_side = "left",
                                annotation_name_gp = gpar(fontsize = 10),
                                col = list(Sample_Group = Barbie$color_panel$bias_group)
  )


  mat <- (Barbie$presence +1 -1) %>% as.matrix()
  mat_name <- "Presence"
  col_fun <- structure(c(2,4), names = c("1","0"))

  if(is.null(bias_group)){
    bias_group <- Barbie$Bias_Occ$group
  }

  if(!show_bias) {
    row_ha <- NULL
  } else if (is.null(row_ha)) {
    row_ha <- rowAnnotation(Bias = bias_group,
                            annotation_name_side = "top",
                            annotation_name_gp = gpar(fontsize = 10),
                            col = list(Bias = Barbie$color_panel$bias_group),
                            show_legend = TRUE,
                            show_annotation_name = TRUE)
  }

  hp <- Heatmap(mat,
                name = mat_name, width = unit(6,"cm"), height = unit(6,"cm"),
                cluster_rows = TRUE, cluster_columns = TRUE,
                show_row_names = FALSE, show_column_names = FALSE,
                column_title = "Samples", row_title = paste(nrow(mat), " Top Clones"),
                col = col_fun,
                right_annotation = row_ha,
                column_split = split_by,
                top_annotation = column_ha,
                cluster_column_slices = FALSE
  )

  return(hp)
}


