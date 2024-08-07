# fit logistic model
fit_logistic_model <- function(Barbie, row) {
  dat <- data.frame(
    Group = as.factor(Barbie$metadata$group),
    Barcode_i = factor(row)
  )

  logiM <- glm(
    Barcode_i ~ Group,
    data = dat,
    family = "binomial"
  )

  plot_model(logiM,
             type = "pred",
             terms = "Group"
  ) +
    labs(y = "Prob(Presence)")

  return(summary(logiM))
}

# fit firth model
firth_model <- function(Barbie, row) {
  dat <- data.frame(
    Group = as.factor(Barbie$metadata$group),
    Barcode_i = factor(row)
  )

  firthM <- logistf(Barcode_i ~ Group, data = dat)

  return(invisible(firthM))

}

# Apply the function to each row of outcome and store the results
GetLogitOcc <- function(Barbie, penalized = F, group_names) {
  # get occurence data
  outcome <- Barbie$presence +1 -1
  # get group names for comparison
  if(is.null(group_names)){
    group_names <- unique(Barbie$metadata$group)
  }

  if(!penalized){ # fit classic logistic model

    results <- lapply(
      1:nrow(outcome),
      function(i)
        fit_logistic_model(Barbie, outcome[i, ])
      )

    # View the summary of the first model as an example
    # results[[1]]$coefficients[2,4]

    p.values.logi <- lapply(results, function(x)
      x$coefficients[2,4]) %>% unlist()

    beta.logi <- lapply(results, function(x)
      x$coefficients[2,1]) %>% unlist()

    signif.logi <- ifelse(p.values.logi >= 0.05,
                          "Unbiased",
                          ifelse(beta.logi > 0,
                                 group_names[2], group_names[1]))

    # get stat output
    LR_Occ <- data.frame(
      p.values = p.values.logi,
      group = signif.logi
    )

    Barbie$LR_Occ_classic <- LR_Occ

  } else { # fit penalized logistic model

    # Apply the function to each row of outcome and store the results
    results_firth <- lapply(
      1:nrow(outcome),
      function(i)
        invisible(firth_model(Barbie, outcome[i, ]))
      )

    p.values.firth <- lapply(results_firth, function(x) x$prob[2]) %>%
      unlist()

    beta.firth <- lapply(
      results_firth, function(x)
        x$coefficients[2]) %>% unlist()

    signif.firth <- ifelse(p.values.firth >= 0.05,
                           "Unbiased",
                           ifelse(beta.firth > 0,
                                  group_names[2], group_names[1]))

    # get stat output
    LR_Occ <- data.frame(
      p.values = p.values.firth,
      group = signif.firth
    )

    Barbie$LR_Occ_penalized <- LR_Occ

  }

  #return
  if(is.null(Barbie$color_panel$bias_group)) {
    Barbie$color_panel$bias_group <- c("C1" = "#33AAFF", "C2" = "#FF5959", "Unbiased" = "#FFC000")
    names(Barbie$color_panel$bias_group)[1:2] <- group_names
  }

  return(Barbie)
}





