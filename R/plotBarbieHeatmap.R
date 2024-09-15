
# @import ComplexHeatmap

# plot heatmap of CPM or pre
plotBarbieHeatmap <- function(Barbie, value="proportion", showTest=FALSE, splitBy=NULL, groupBy=NULL,
                              rowAnnotation=NULL, columnAnnotation=NULL,
                              biasDirection=NULL) {
  ## check Barbie dimensions
  if(!checkBarbieDimensions(Barbie))
    stop("Barbie components are not in right format or dimensions don't match.
         use createBarbie() and other functions in the 'Barbie' package to modify the object and don't do it manually.")

  ## check which value to visualize
  value <- match.arg(value, c("proportion", "occurrence"))

  ## check groupBy: if 'groupBy' is a specified effector name, extract the entire vector
  if(is.character(groupBy)) {
    if(groupBy %in% colnames(targets)) {
      groupBy <- targets[,groupBy]
      mytargets <- targets
      pointer <- which(colnames(mytargets) == groupBy)
      message("found", groupBy, "as an effector in targets or Barbie$metadata.")
    } else {stop("the groupBy specified is a charactor value,
                 but it's not an effector name found in targets or Barbie$metadata.
                 make sure you spell it correctly.")}
  } else if(is.vector(groupBy) || is.factor(groupBy)) {
    if(length(groupBy) != ncol(Barbie$assay)) {
      stop("the length of 'groupBy' doesn't match the column dimention (sample size) of 'targets' or'Barbie$assay'.")
    } else {
      mytargets <- data.frame(groupBy=groupBy, targets)
      pointer <- which(colnames(mytargets) == "groupBy")
      message("binding 'groupBy' to 'targets'.")
    }
  }

  columnHA = HeatmapAnnotation(
    Sample_Group = groupBy,
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 10),
    col = list(
      Sample_Group = Barbie$colorFactors$biasDirection)
  )

  mat <- log2(Barbie$CPM +1) %>% as.matrix()
  mat_name <- "logCPM"
  col_fun <- circlize::colorRamp2(c(min(mat), median(mat), max(mat)), c("blue", "white", "red"))

  if(is.null(biasDirection)){
    biasDirection <- Barbie$Bias_Occ$group
  }

  if(!showTest) {
    rowAnnotation <- NULL
  } else if (is.null(rowAnnotation)) {
    rowAnnotation <- rowAnnotation(Bias = biasDirection,
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
                right_annotation = rowAnnotation,
                top_annotation = columnHA,
                column_split = splitBy,
                cluster_column_slices = FALSE
  )

  return(hp)
}

# plot presence
PlotPreHP_0 <- function(Barbie, showTest = FALSE, rowAnnotation = NULL, splitBy = NULL,
                        groupBy = NULL, biasDirection = NULL) {

  columnHA = HeatmapAnnotation(Sample_Group = groupBy,
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
    rowAnnotation <- NULL
  } else if (is.null(rowAnnotation)) {
    rowAnnotation <- rowAnnotation(Bias = biasDirection,
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
                right_annotation = rowAnnotation,
                column_split = splitBy,
                top_annotation = columnHA,
                cluster_column_slices = FALSE
  )

  return(hp)
}


