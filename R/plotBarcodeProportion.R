# @import colorRamp2

# plot a histogram contribution

DesignColorPalette <- function(Barbie, light = "#FFC1E0", dark = "#FF3399") {
  contribution <- Barbie$CPM
  avg_contribution <- rowSums(contribution) / sum(rowSums(contribution)) *100

  rank <- rank(-avg_contribution, ties.method = "first")

  color_palette <- colorRamp2(c(min(avg_contribution), max(avg_contribution)),
                              c(light, dark)) #FF99CC

  color_group <- color_palette(avg_contribution)

  names(color_group) <- rownames(contribution)

  return(color_group)
}

PlotBarContribution <- function(Barbie, coord = FALSE, color_defaull = TRUE, color_group) {
  contribution <- Barbie$CPM
  avg_contribution <- rowSums(contribution) / sum(rowSums(contribution)) *100

  rank <- rank(-avg_contribution, ties.method = "first")

  data <- data.frame(
    individual = rownames(contribution),
    value = avg_contribution,
    ranking = as.numeric(rank)
  )

  # #order data
  # data = data %>% arrange(rank)
  # data$id <- nrow(data) |> seq()

  if(color_defaull) {
    color_group <- "#FF3399"
    # color_board <- "#FFC1E0"
  } else {
    color_group <- color_group
    # color_board <- color_group[which.max(rank)]
  }

  #make the plot
  p <- ggplot(data, aes(x=ranking, y=value)) +
    geom_bar(stat="identity", alpha= 1, color = alpha(color_group, 0.2), fill = color_group, size = 1) +

    labs(title = "Contribution of Top Clones",
         x = paste0(nrow(data), " Clones"),
         y = "% Contribution") +
    theme_minimal() +
    theme(axis.ticks = element_blank(),
          panel.background = element_blank(),  # Remove plot background
          panel.grid = element_blank(),
          axis.text.x = element_blank() # Remove x-axis tick labels
    )

  p

  if(coord){
    p <- p + coord_fixed()
  }

  if(!coord){
    p <- p + theme(aspect.ratio = 1)
  }

  return(p)
}


# example usage:
# PlotBarContribution(week4_top, coord = FALSE)
# PlotBarContribution(week8_top)



# plot two sets of contribution
PlotBarContribution_Double <- function(BarbieA, BarbieB, coord = FALSE) {
  contributionA <- BarbieA$CPM
  contributionB <- BarbieB$CPM

  valueA <- rowSums(contributionA) / sum(rowSums(contributionA)) *100
  valueB <- rowSums(contributionB) / sum(rowSums(contributionB)) *100

  dfA <- data.frame(rownames = rownames(contributionA), valueA = valueA, stringsAsFactors = FALSE)
  dfB <- data.frame(rownames = rownames(contributionB), valueB = valueB, stringsAsFactors = FALSE)

  contribution <- merge(dfA, dfB, by = "rownames", all = TRUE)
  rownames(contribution) <- contribution$rownames
  contribution <- contribution[,-1]
  contribution[is.na(contribution)] <- 0

  contribution$mean <- rowMeans(contribution)

  rank <- rank(-contribution$valueA, ties.method = "first")

  contribution$ranking <- rank

  # contribution = contribution %>% arrange(-contribution$mean)
  # contribution$id <- nrow(contribution) |> seq()

  color_palette <- colorRamp2(c(min(contribution$valueA), max(contribution$valueA)),
                              c("#FFC1E0", "#FF3399")) #FF99CC

  color_group <- color_palette(contribution$valueA)

  names(color_group) <- rownames(contribution)

  color_group[! rownames(contribution) %in% rownames(contributionB)] <- "#3CDA2D" #A only

  color_group[! rownames(contribution) %in% rownames(contributionA)] <- "#7833FF" #B only

  data <- contribution

  p <- ggplot(data) +
    geom_bar(aes(x=ranking, y=valueA),
             stat="identity", alpha=1, color = alpha(color_group, 0.2), fill = color_group, size = 1) +
    geom_bar(aes(x=ranking, y=-valueB),
             stat="identity", alpha=1, color = alpha(color_group, 0.2), fill = color_group, size = 1) +
    labs(title = "Contribution of Top Clones",
         x = paste0(nrow(data), " Clones"),
         y = "% Contribution") +
    theme_minimal() +
    theme(axis.ticks = element_blank(),
          panel.background = element_blank(),  # Remove plot background
          panel.grid = element_blank(),
          axis.text.x = element_blank() # Remove x-axis tick labels
    ) +
    geom_hline(yintercept = 0, linetype = "solid", color = "white", size = 0.1)  # Add horizontal line at x = 0

  if(coord){
    p <- p + coord_fixed()
  }

  if(!coord){
    p <- p + theme(aspect.ratio = 1)
  }

  return(p)

}
