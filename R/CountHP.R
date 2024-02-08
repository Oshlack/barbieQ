PlotCpmHP_0 <- function(Barbie, show_bias = FALSE, row_ha = NULL, split_by = NULL,
                        Vector_customized = NULL) {

  column_ha = HeatmapAnnotation(Sample_Group = Vector_customized,
                                annotation_name_side = "left",
                                annotation_name_gp = gpar(fontsize = 10),
                                col = list(Sample_Group = Barbie$color_panel$bias_group)
  )

  mat <- log2(Barbie$assay +1) %>% as.matrix()
  mat_name <- "logCPM"
  col_fun <- circlize::colorRamp2(c(min(mat), median(mat), max(mat)), c("blue", "white", "red"))


  if(!show_bias) {
    row_ha <- NULL
  } else if (is.null(row_ha)) {
    row_ha <- rowAnnotation(Bias = Barbie$Bias$group,
                            annotation_name_side = "bottom",
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
                column_split = split_by
  )

  return(hp)
}


# plot heatmap of CPM or pre

PlotCpmHP <- function(Barbie, show_bias = FALSE, row_ha = NULL, split_by = NULL) {

  if (is.null(Barbie$color_panel)) {
    Barbie$color_panel <- list(lineage = c(Unknown = "black"), mouse = c(Unknown = "black"),
                               tissue = c(Unknown = "black"), treat = c(Unknown = "black"),
                               celltype = c(Unknown = "black"), time = c(Unknown = "black"))
  }


  if (is.null(Barbie$metadata)) {
    Barbie$metadata <- data.frame(time = rep("Unknown", ncol(Barbie$assay)), treat = rep("Unknown", ncol(Barbie$assay)),
                                  lineage = rep("Unknown", ncol(Barbie$assay)), tissue = rep("Unknown", ncol(Barbie$assay)),
                                  celltype = rep("Unknown", ncol(Barbie$assay)), mouse = rep("Unknown", ncol(Barbie$assay)))
  }

  column_ha = HeatmapAnnotation(Time = factor(Barbie$metadata$time, levels = c("week4", "week8")),
                                Treat = Barbie$metadata$treat,
                                Lineage = Barbie$metadata$lineage,
                                Tissue = Barbie$metadata$tissue,
                                Cell = Barbie$metadata$celltype,
                                Mouse = Barbie$metadata$mouse,
                                annotation_name_side = "left",
                                annotation_name_gp = gpar(fontsize = 10),
                                col = list(Lineage = Barbie$color_panel$lineage,
                                           Mouse = Barbie$color_panel$mouse,
                                           Tissue = Barbie$color_panel$tissue,
                                           Treat = Barbie$color_panel$treat,
                                           Cell = Barbie$color_panel$celltype,
                                           Time = Barbie$color_panel$time
                                )
  )


  mat <- log2(Barbie$assay +1) %>% as.matrix()
  mat_name <- "logCPM"
  col_fun <- circlize::colorRamp2(c(min(mat), median(mat), max(mat)), c("blue", "white", "red"))


  if(!show_bias) {
    row_ha <- NULL
    } else if (is.null(row_ha)) {
    row_ha <- rowAnnotation(Bias = Barbie$Bias$group,
                              annotation_name_side = "bottom",
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
                column_split = split_by,
                top_annotation = column_ha
  )

  return(hp)
}


# ploy presence
PlotPreHP_0 <- function(Barbie, show_bias = FALSE, row_ha = NULL, split_by = NULL,
                        Vector_customized = NULL) {

  column_ha = HeatmapAnnotation(Sample_Group = Vector_customized,
                                annotation_name_side = "left",
                                annotation_name_gp = gpar(fontsize = 10),
                                col = list(Sample_Group = Barbie$color_panel$bias_group)
  )


  mat <- (Barbie$presence +1 -1) %>% as.matrix()
  mat_name <- "Presence"
  col_fun <- structure(c(2,4), names = c("1","0"))

  if(!show_bias) {
    row_ha <- NULL
  } else if (is.null(row_ha)) {
    row_ha <- rowAnnotation(Bias = Barbie$Bias$group,
                            annotation_name_side = "bottom",
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
                top_annotation = column_ha
  )

  return(hp)
}

# ploy presence
PlotPreHP <- function(Barbie, show_bias = FALSE, row_ha = NULL, split_by = NULL) {

  column_ha = HeatmapAnnotation(Time = Barbie$metadata$time,
                                Treat = Barbie$metadata$treat,
                                Lineage = Barbie$metadata$lineage,
                                Tissue = Barbie$metadata$tissue,
                                Cell = Barbie$metadata$celltype,
                                Mouse = Barbie$metadata$mouse,
                                annotation_name_side = "left",
                                annotation_name_gp = gpar(fontsize = 10),
                                col = list(Lineage = Barbie$color_panel$lineage,
                                           Mouse = Barbie$color_panel$mouse,
                                           Tissue = Barbie$color_panel$tissue,
                                           Treat = Barbie$color_panel$treat,
                                           Cell = Barbie$color_panel$celltype,
                                           Time = Barbie$color_panel$time
                                )
  )


  mat <- (Barbie$presence +1 -1) %>% as.matrix()
  mat_name <- "Presence"
  col_fun <- structure(c(2,4), names = c("1","0"))

  if(!show_bias) {
    row_ha <- NULL
  } else if (is.null(row_ha)) {
    row_ha <- rowAnnotation(Bias = Barbie$Bias$group,
                            annotation_name_side = "bottom",
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
                  top_annotation = column_ha
    )

  return(hp)
}




