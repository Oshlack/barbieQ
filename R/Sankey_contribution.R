# plot total contribution: top v.s. bottom clones

PlotTotalContribution <- function(Barbie) {
  flag <- Barbie$is_top
  contribution <- Barbie$CPM

  # sankey contribution
  ntop <- sum(flag) / length(flag) * 100
  nbottom <- sum(!flag) / length(flag) * 100

  total_sum <- rowSums(contribution) / sum(rowSums(contribution)) *100
  sumtop <- total_sum[flag] |> sum()
  sumbottom <- total_sum[!flag] |> sum()

  # Sample data
  ppdata <- data.frame(
    category = c("Num of Clones", "Contribution in Progeny Samples") %>% factor(levels = c("Num of Clones", "Contribution in Progeny Samples")),
    top = c(ntop, sumtop),
    bottom = c(nbottom, sumbottom)
  )

  # Reshape data for stacked bar plot
  data_long <- tidyr::gather(ppdata, key = "variable", value = "value", -category)
  data_long$variable <- factor(data_long$variable, levels = c("top", "bottom"))

  # Create the stacked bar plot
  p <- ggplot(data_long, aes(x = category, y = value, fill = variable)) +
    geom_bar(stat = "identity", width = 0.5) +
    facet_grid(~category, scales = "free_x", space = "free_x") +
    geom_text(aes(label = paste0(round(value), "%")), position = position_stack(vjust = 0.5)) +
    labs(x = " ", y = "Percentage", title = "Clone Contributions") +

    labs(fill = "Clones") +
    #set color
    scale_fill_manual(
      values = c("top" = "#FF3399", "bottom" = "#0066FF"),
      labels = c(paste0("Top ", sum(flag)), "Others")
    ) +

    theme(axis.ticks = element_blank(),
          panel.background = element_blank(),  # Remove plot background
          panel.grid = element_blank(),
          axis.title.y = element_blank(),  # Remove y-axis label
          axis.text.y = element_blank(),   # Remove y-axis tick labels
          axis.ticks.y = element_blank())

  return(p)
}

# example usage:
# plot(PlotTotalContribution(week4_ba))
# plot(PlotTotalContribution(week8_ba))
