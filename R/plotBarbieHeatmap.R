#' plotting Barcode CPM (or occurrence) in a Heatmap with sample annotations
#'
#' @param Barbie a Barbie object created by createBarbie()
#' @param value a string value choosing to present "CPM" or "occurrence"
#' @param splitSamples a logical value deciding whether to split samples into slices
#' @param targets a data.frame containing sample conditions - each effector is a column
#' @param sampleGroups a vector or a string value indicating the primary effector in targets
#' @param barcodeAnnotation a row Annotation object created by ComplexHeatmap::rowAnnotation()
#' @param sampleAnnotation a column Annotation object created by ComplexHeatmap::HeatmapAnnotation()
#'
#' @return a "Heatmap" S4 object
#' @export
#'
#' @import ComplexHeatmap
#' @importFrom circlize colorRamp2
#' @import grid
#'
#' @examples
#' HSC <- Barbie::HSC
#' plotBarbieHeatmap(Barbie = HSC)
plotBarbieHeatmap <- function(Barbie, value="CPM", splitSamples=FALSE,
                              targets=NULL, sampleGroups=NULL,
                              barcodeAnnotation=NULL, sampleAnnotation=NULL) {
  ## check Barbie dimensions
  if(!checkBarbieDimensions(Barbie))
    stop("Barbie components are not in right format or dimensions don't match.
please start with Barbie::createBarbie() and use proper functions to modify the object - don't do it manually.")

  ## check which value to visualize
  value <- match.arg(value, c("CPM", "occurrence"))

  ## extract targets and primary effector based on arguments
  targetsInfo <- extractTargetsAndPrimaryFactor(Barbie=Barbie, targets=targets, sampleGroups=sampleGroups)
  mytargets <- targetsInfo$mytargets
  pointer <- targetsInfo$pointer
  ## set the primary effector as sample splitter displaying at bottom
  bottomTargets <- mytargets[, pointer, drop=FALSE]
  ## the rest of effectors displayed at top
  topTargets <- mytargets[, setdiff(seq_along(colnames(mytargets)), pointer), drop=FALSE]

  sampleAnnotation = HeatmapAnnotation(
    df = topTargets,
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 10),
    col = Barbie$factorColors
  )

  if(splitSamples) {
      groupAnnotation = HeatmapAnnotation(
        df = bottomTargets,
        annotation_name_side = "left",
        annotation_name_gp = gpar(fontsize = 10),
        col = Barbie$factorColors
      )
    splitBy <- bottomTargets
  } else {
    groupAnnotation <- NULL
    splitBy <- NULL}

  ## choose values to be visualised
  if(value == "CPM") {
    mat <- log2(Barbie$CPM +1) %>% as.matrix()
    matTitle <- "log2 CPM+1"
    colorFun <- circlize::colorRamp2(
      c(min(mat), mean(mat), max(mat)),
      c("blue", "white", "red"))
  } else {
    mat <- (Barbie$occurrence +1 -1) %>% as.matrix()
    matTitle <- "occurrence"
    colorFun <- structure(c(2,4), names = c("1","0"))
  }

  hp <- Heatmap(mat,
                name = matTitle, width = unit(6,"cm"), height = unit(6,"cm"),
                cluster_rows = TRUE, cluster_columns = TRUE,
                show_row_names = FALSE, show_column_names = FALSE,
                column_title = paste0(ncol(mat), " Samples"),
                row_title = paste0(nrow(mat), " Barcodes"),
                col = colorFun,
                right_annotation = barcodeAnnotation,
                top_annotation = sampleAnnotation,
                bottom_annotation = groupAnnotation,
                column_split = splitBy,
                cluster_column_slices = FALSE)

  return(hp)
}

