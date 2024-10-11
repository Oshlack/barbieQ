#' Plot total contributions of Barcodes grouped by \emph{top} vs. \emph{bottom}
#' in a stacked bar plot.
#'
#' After the [tagTopBarcodes] function tags Barcodes as either
#'  \emph{top} or \emph{bottom}, `plotBarcodeSankey()` visualizes their
#'  relative frequency and total contribution across all samples using
#'  a stacked bar plot, resembling a Sankey plot.
#'
#' @param Barbie A `Barbie` object created by the [createBarbie] function.
#'
#' @return A `ggplot` S3 class object displaying the Sankey-like stacked bar
#'  plot, where Barcodes are categorized as either \emph{top} or \emph{bottom}.
#'
#' @export
#'
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
#' @import ggplot2
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
#' myBarbie <- tagTopBarcodes(myBarbie)
#' plotBarcodeSankey(myBarbie)
plotBarcodeSankey <- function(Barbie) {
  flag <- Barbie$isTop$vec
  contribution <- Barbie$CPM
  ## sankey contribution
  ntop <- sum(flag) / length(flag) * 100
  nbottom <- sum(!flag) / length(flag) * 100
  totalSum <- rowSums(contribution) / sum(rowSums(contribution)) *100
  sumtop <- totalSum[flag] |> sum()
  sumbottom <- totalSum[!flag] |> sum()
  ## define sample data
  ppdata <- data.frame(
    category = c("num of Barcodes", "relative total proportion") %>%
      factor(levels = c("num of Barcodes", "relative total proportion")),
    top = c(ntop, sumtop),
    bottom = c(nbottom, sumbottom)
  )
  ## reshape data for stacked bar plot
  dataLong <- tidyr::gather(ppdata, key = "variable", value = "value", -category)
  dataLong$variable <- factor(dataLong$variable, levels = c("top", "bottom"))

  suppressWarnings({
    p <- ggplot(dataLong, aes(x = category, y = value, fill = variable)) +
      geom_bar(stat = "identity", width = 0.5) +
      facet_grid(~category, scales = "free_x", space = "free_x") +
      geom_text(aes(label = paste0(round(value), "%")), position = position_stack(vjust = 0.5)) +
      labs(x = " ", y = "relative total proportion (%)") +
      labs(fill = "Barcodes") +
      scale_fill_manual(
        values = c("top" = "#FF3399", "bottom" = "#0066FF"),
        labels = c(paste0("Top ", sum(flag)), "Others")
      ) +
      theme(axis.ticks = element_blank(),
            panel.background = element_blank(),
            panel.grid = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
  })


  return(p)
}
