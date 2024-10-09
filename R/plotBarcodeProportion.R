#' plotting Barcode contributions (percentage of average Barcode proportion across samples)
#'
#' @param Barbie  a Barbie object created by createBarbie()
#' @param coordFixed a logical value to decide whether to coordinate x and y scale
#' @param colorGradient a logical value to choose bar colors being consistant or gradient
#'
#' @return a "ggplot" S3 class object
#' @export
#'
#' @import ggplot2
#' @importFrom circlize colorRamp2
#' @import dplyr
#'
#' @examples
#' HSC <- Barbie::HSC
#' plotBarcodeProportion(HSC)
#' plotBarcodeProportion(HSC, colorGradient=TRUE)
plotBarcodeProportion <- function(Barbie, coordFixed=FALSE, colorGradient=FALSE) {
  ## check Barbie dimensions
  if(!checkBarbieDimensions(Barbie))
    stop("Barbie components are not in right format or dimensions don't match.
please start with Barbie::createBarbie() and use proper functions to modify the object - don't do it manually.")
  ## compute mean Barcode proportion across samples as contribution
  contribution <- rowSums(Barbie$proportion) / ncol(Barbie$proportion)
  relativeContribution <- contribution / sum(contribution) *100
  ## rank the relative contribution of each Barcode
  rank <- rank(-relativeContribution, ties.method = "first")
  ## extraxt plotting data
  data <- data.frame(
    individual = rownames(Barbie$proportion),
    percentage = relativeContribution,
    rank = as.numeric(rank)
  )
  ## define color function and assign bar colors
  colorFun <- colorRamp2(
    c(min(relativeContribution), max(relativeContribution)),
    c("#FFC1E0", "#FF3399"))
  if(colorGradient) {
    barColor <- colorFun(relativeContribution)
  } else {
    barColor <- "#FF3399"
  }

  p <- ggplot(data, aes(x=rank, y=percentage)) +
    geom_bar(stat="identity", alpha= 1, color = alpha(barColor, 0.2),
             fill = barColor, linewidth = 1) +
    labs(title = "Barcode Contribution",
         x = paste0(nrow(data), " Barcodes"),
         y = "Relative Average Barcode Proportion (%)") +
    theme_minimal() +
    theme(axis.ticks = element_blank(),
          panel.background = element_blank(),  # Remove plot background
          panel.grid = element_blank(),
          axis.text.x = element_blank() # Remove x-axis tick labels
    )

  if(coordFixed){
    p <- p + coord_fixed()
  } else {
    p <- p + theme(aspect.ratio = 1)
  }

  return(p)
}
