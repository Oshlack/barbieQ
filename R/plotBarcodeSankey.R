#' plotting relative Barcode total contribution in Sankey plot - top vs bottom Barcodes
#'
#' @param Barbie
#'
#' @return a "ggplot" S3 class object
#' @export
#'
#' @import tidyr
#'
#' @examples
#' HSC <- Barbie::HSC
#' plotBarcodeSankey(HSC)
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

  return(p)
}
