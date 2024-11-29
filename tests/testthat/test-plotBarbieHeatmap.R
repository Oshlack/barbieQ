test_that("plotting barcode Heatmap works", {
  ## Sample conditions and color palettes
  sampleConditions <- data.frame(
    Treat = factor(rep(c("ctrl", "drug"), each = 6)),
    Time = rep(rep(seq_len(2), each = 3), 2)
  )
  conditionColor <- list(
    Treat = c(ctrl = "#999999", drug = "#112233"),
    Time = c("1" = "#778899", "2" = "#998877")
  )

  ## Barcode count data
  nbarcodes <- 50
  nsamples <- 12
  barcodeCount <- abs(matrix(rnorm(nbarcodes * nsamples), nbarcodes, nsamples))
  rownames(barcodeCount) <- paste0("Barcode", seq_len(nbarcodes))

  object1 <- createBarbieQ(
    object = barcodeCount, target = sampleConditions, factorColors = conditionColor
  )

  hp <- plotBarcodeHeatmap(barbieQ = object1)
  expect_s4_class(hp, "Heatmap")
  expect_equal(hp@name, "log2 CPM+1")
  hp <- plotBarcodeHeatmap(barbieQ = object1, splitSamples = TRUE)
  expect_equal(
    hp@top_annotation@anno_list %>% names(),
    c("Time")
  )
  hp <- plotBarcodeHeatmap(
    barbieQ = object1, splitSamples = TRUE, sampleGroups = "Treat"
  )
  expect_equal(hp@bottom_annotation@anno_list %>% names(), "Treat")

  object1$metadata$Treat <- factor(
    object1$metadata$Treat,
    levels = c("drug", "ctrl")
  )
  hp <- plotBarcodeHeatmap(
    barbieQ = object1, splitSamples = TRUE, sampleGroups = "Treat"
  )
  expect_equal(
    hp@matrix_param[["column_split"]][[1]] %>% levels(),
    c("drug", "ctrl")
  )

  hp <- plotBarcodeHeatmap(barbieQ = object1, value = "occurrence")
  expect_equal(hp@name, "occurrence")
})
