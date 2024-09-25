#' plotting Barcode proportion by a circular bar plot
#'
#' @param Barbie a Barbie object created by createBarbie()
#'
#' @return a "ggplot" S3 class object
#' @export
#'
#' @import dplyr
#'
#' @examples
#' HSC <- Barbie::HSC
#' plotBarcodePareto(HSC)
plotBarcodePareto <- function(Barbie) {
  flag <- Barbie$isTop$vec
  contribution <- Barbie$CPM
  ## create dataset
  ntop <- sum(flag)
  nbottom <- sum(!flag)
  topTag <- factor(flag, levels = c(TRUE, FALSE))
  data <- data.frame(
    individual=rownames(contribution),
    group=topTag,
    value=rowSums(contribution) / sum(rowSums(contribution)) *100
  )
  ## set a number of 'empty bar' to add at the end of each group
  emptyBar <- 5
  toAdd <- data.frame( matrix(NA, emptyBar*nlevels(data$group), ncol(data)))
  colnames(toAdd) <- colnames(data)
  toAdd$group <- rep(levels(data$group), each=emptyBar)
  data <- rbind(data, toAdd)
  ## order data
  data = data %>% arrange(group, value)
  data$id <- nrow(data) |> seq()
  ## add label
  labelData <- data
  numberBar <- nrow(labelData)
  angle <- 90 - 360*(labelData$id-0.5) / numberBar
  labelData$hjust <- ifelse(angle < -90, 1, 0)
  labelData$angle <- ifelse(angle < -90, angle+180, angle)

  ## prepare a data frame for base lines
  baseData <- data %>%
    dplyr::group_by(group) %>%
    dplyr::summarize(start = min(id), end = max(id) - emptyBar) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(title = mean(c(start, end)))
  ## calculate the number of clones for each group
  baseData$labs <- ifelse(baseData$group == "TRUE", ntop, nbottom)
  ## prepare a data frame for grid (scales)
  gridData <- baseData
  gridData$end <- gridData$end[ c(nrow(gridData), 1:nrow(gridData)-1)] + 1
  gridData$start <- gridData$start - 1
  gridData <- gridData[-1,]

  p <- ggplot(data, aes(x=as.factor(id), y=value)) +
    ## add the bars with a blue color
    geom_bar(stat="identity", alpha=1, color="#999966") +
    ## add a val=100/75/50/25 lines - compute it first to make sure barplots are OVER it.
    geom_segment(data=gridData, aes(x = end, y = 20, xend = start, yend = 20),
                 colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=gridData, aes(x = end, y = 15, xend = start, yend = 15),
                 colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=gridData, aes(x = end, y = 10, xend = start, yend = 10),
                 colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=gridData, aes(x = end, y = 5, xend = start, yend = 5),
                 colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    ## add text showing the value of each 100/75/50/25 lines
    annotate("text", x = rep(max(data$id),5), y = c(5, 10, 15, 20, 25),
             label = c("5%", "10%", "15%", "20%", "relative mean Barcode Proportion") ,
             color="#999966", size=3 , angle=0, fontface="bold", hjust=1) +
    ## the negative value controls the size of the central circle, the positive to add size over each bar
    ylim(-15,25) +
    ## custom the theme: no axis title and no cartesian grid
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(0,4), "cm") # This remove unnecessary margin around plot
    ) +
    ## make the coordinate polar instead of cartesian.
    coord_polar(start = 0) +
    #  geom_text(data = labelData, aes(x = id, y = value +30, label=individual, hjust=hjust),
    #            color = "black", fontface = "bold", alpha=0.5, size=1, angle= labelData$angle, inherit.aes = F) +
    ## add base line information
    geom_segment(data=baseData, aes(x = start, y = -2, xend = end, yend = -2, color=group),
                 alpha=0.8, size=2 , inherit.aes = FALSE )  +
    geom_text(data=baseData, aes(x = title+15, y = c(-5,-10), label=labs, color=group),
              hjust=c(1,1), alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE, show.legend = F) +
    labs(color = "Num of clones") +
    scale_color_manual(
      values = c("TRUE" = "#FF3399", "FALSE" = "#0066FF"),
      labels = c("Top Barcodes", "Others")
    )

  return(p)
}
