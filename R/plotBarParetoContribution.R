# Plot circular bar plot of clone contribution.

PlotCircularContribution <- function(Barbie) {

  flag <- Barbie$is_top
  contribution <- Barbie$CPM
  # Create dataset
  ntop <- sum(flag)
  nbottom <- sum(!flag)
  top_factor <- factor(flag, levels = c(TRUE, FALSE))

  data <- data.frame(
    individual=rownames(contribution),
    group=top_factor,
    value=rowSums(contribution) / sum(rowSums(contribution)) *100
  )

  # Set a number of 'empty bar' to add at the end of each group
  empty_bar <- 5
  to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)))
  colnames(to_add) <- colnames(data)
  to_add$group <- rep(levels(data$group), each=empty_bar)
  data <- rbind(data, to_add)

  #order data
  data = data %>% arrange(group, value)
  data$id <- nrow(data) |> seq()

  #add label
  label_data <- data
  number_of_bar <- nrow(label_data)
  angle <- 90 - 360*(label_data$id-0.5) / number_of_bar
  label_data$hjust <- ifelse(angle < -90, 1, 0)
  label_data$angle <- ifelse(angle < -90, angle+180, angle)

  # prepare a data frame for base lines
  base_data <- data %>%
    dplyr::group_by(group) %>%
    dplyr::summarize(start = min(id), end = max(id) - empty_bar) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(title = mean(c(start, end)))

  # Calculate the number of clones for each group
  base_data$labs <- ifelse(base_data$group == "TRUE", ntop, nbottom)

  # prepare a data frame for grid (scales)
  grid_data <- base_data
  grid_data$end <- grid_data$end[ c(nrow(grid_data), 1:nrow(grid_data)-1)] + 1
  grid_data$start <- grid_data$start - 1
  grid_data <- grid_data[-1,]


  # Make the plot
  p <- ggplot(data, aes(x=as.factor(id), y=value)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar

    # This add the bars with a blue color
    geom_bar(stat="identity", alpha=1, color="#999966") +

    # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
    geom_segment(data=grid_data, aes(x = end, y = 20, xend = start, yend = 20), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 15, xend = start, yend = 15), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 10, xend = start, yend = 10), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 5, xend = start, yend = 5), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +

    # Add text showing the value of each 100/75/50/25 lines
    annotate("text", x = rep(max(data$id),5), y = c(5, 10, 15, 20, 25), label = c("5%", "10%", "15%", "20%", "Clone Contribution") , color="#999966", size=3 , angle=0, fontface="bold", hjust=1) +

    # Limits of the plot = very important. The negative value controls the size of the inner circle, the positive one is useful to add size over each bar
    ylim(-15,25) +

    # Custom the theme: no axis title and no cartesian grid
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(0,4), "cm") # This remove unnecessary margin around plot
    ) +

    # This makes the coordinate polar instead of cartesian.
    coord_polar(start = 0) +
    #  geom_text(data = label_data, aes(x = id, y = value +30, label=individual, hjust=hjust),
    #            color = "black", fontface = "bold", alpha=0.5, size=1, angle= label_data$angle, inherit.aes = F) +

    # Add base line information
    geom_segment(data=base_data, aes(x = start, y = -2, xend = end, yend = -2, color=group), alpha=0.8, size=2 , inherit.aes = FALSE )  +
    geom_text(data=base_data, aes(x = title+15, y = c(-5,-10), label=labs, color=group), hjust=c(1,1), alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE, show.legend = F) +

    labs(color = "Num of clones") +
    scale_color_manual(
      values = c("TRUE" = "#FF3399", "FALSE" = "#0066FF"),
      labels = c("Top clones", "Others")
    )

  return(p)
}

# example usage:
# plot(PlotCircularContribution(Barbie = week4_ba))


