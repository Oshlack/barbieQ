#' Plot sample pairwise correlation in a Heatmap
#'
#' `plotSamplePairCorrelation()` visualizes the pairwise correlation of
#'  CPM values across samples in a `Barbie` object,
#'  using a heatmap with a checkerboard-like pattern.
#'
#' @param Barbie A `Barbie` object created by the [createBarbie] function.
#' @param sampleOrder A character vector of names of the factors in the
#'  sample conditions provided by `Barbie` or `targets`,
#'  specifying the order in which samples should be arranged.
#'  Defaults to the original order of factors in the sample conditions.
#' @param targets A `matrix` or `data.frame` of sample conditions,
#'  where each factor is represented in a separate column. Defaults to NULL,
#'  in which case sample conditions are inherited from `Barbie$metadata`.
#' @param sampleGroups A string representing the name of a factor from the
#'  sample conditions passed by `Barbie` or `targets`, or a vector of
#'  sample conditions, indicating the primary factor to split sample slices.
#' @param method A string specifying the correlation method to use.
#'  Defaults to "pearson". Options include: "pearson", "spearman".
#'
#' @return A "Heatmap" S4 object displaying the pairwise correlation between
#'  samples in a Heatmap.
#'
#' @export
#'
#' @import ComplexHeatmap
#' @importFrom circlize colorRamp2
#' @importFrom grid gpar
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr pull
#' @importFrom dplyr arrange
#' @importFrom dplyr across
#' @importFrom dplyr row_number
#' @importFrom dplyr all_of
#' @importFrom stats cor
#' @import data.table
#'
#' @examples
#' ## sample conditions and color palettes
#' sampleConditions <- data.frame(
#'   Treat=factor(rep(c("ctrl", "drug"), each=6)),
#'   Time=rep(rep(1:2, each=3), 2))
#' conditionColor <- list(
#'   Treat=c(ctrl="#999999", drug="#112233"),
#'   Time=c("1"="#778899", "2"="#998877"))
#' ## Barcode count data
#' nbarcodes <- 50
#' nsamples <- 12
#' barcodeCount <- abs(matrix(10, nbarcodes, nsamples))
#' barcodeCount[21:50,] <- 0.0001
#' rownames(barcodeCount) <- paste0("Barcode", 1:nbarcodes)
#' ## create a `Barbie` object
#' myBarbie <- createBarbie(barcodeCount, sampleConditions, conditionColor)
#' plotSamplePairCorrelation(myBarbie)
plotSamplePairCorrelation <- function(
    Barbie, sampleOrder=NULL, targets=NULL, sampleGroups=NULL,
    method="pearson") {
  ## check Barbie dimensions
  checkBarbieDimensions(Barbie)
  ## extract targets and primary effector based on arguments
  targetsInfo <- extractTargetsAndPrimaryFactor(
    Barbie=Barbie, targets=targets, sampleGroups=sampleGroups)
  mytargets <- targetsInfo$mytargets
  pointer <- targetsInfo$pointer
  ## set the primary effector as sample splitter displaying at bottom
  bottomTargets <- mytargets[, pointer, drop=FALSE]
  ## the rest of effectors displayed at top
  topTargets <- mytargets[, dplyr::setdiff(
    seq_along(colnames(mytargets)), pointer), drop=FALSE]
  ## check sampleOrder
  if(is.null(sampleOrder)) {
    sampleOrder <- c(colnames(bottomTargets), colnames(topTargets))
  }
  ## check method
  method <- match.arg(method, c("pearson", "spearman"))

  ## reorder the columns of the data.frame based on the specified order
  mytargets <- mytargets[ , sampleOrder, drop = FALSE]
  ## extract sample order
  rowOrder <- mytargets %>%
    dplyr::mutate(rowNumber = dplyr::row_number()) %>%
    dplyr::arrange(dplyr::across(dplyr::all_of(sampleOrder))) %>%
    dplyr::pull(rowNumber)
  ## if sampleGroups column is automatically added, remove it
  if(all(mytargets[,pointer] == 1) &&
     colnames(mytargets)[pointer] == "sampleGroups") {
    mytargets <- mytargets[ , -pointer, drop = FALSE]
  }
  ## compute annotation obejct
  sampleAnnotationColumn = HeatmapAnnotation(
    df = mytargets,
    annotation_name_side = "left",
    annotation_name_gp = grid::gpar(fontsize = 10),
    col = Barbie$factorColors
  )
  sampleAnnotationRow = rowAnnotation(
    df = mytargets,
    annotation_name_side = "bottom",
    annotation_name_gp = grid::gpar(fontsize = 10),
    col = Barbie$factorColors ,
    show_legend = FALSE,
    show_annotation_name = F
  )
  ## calculate correlation
  corMat <- stats::cor(log2(Barbie$CPM +1), method = method)
  message("displaying", method, "correlation of Barcode log2 CPM+1.")

  hp <- Heatmap(
    corMat,
    name = "correlation", width = unit(6,"cm"), height = unit(6,"cm"),
    col = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
    heatmap_legend_param = list(at = c(-1, 0, 1)),
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = FALSE, show_column_names = FALSE,
    top_annotation = sampleAnnotationColumn,
    left_annotation = sampleAnnotationRow,
    column_order = rowOrder, row_order = rowOrder)

  return(hp)
}
