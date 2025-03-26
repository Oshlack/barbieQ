test_that("plotting sample pair-wise correlation in Heatmap works", {
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
    object = barcodeCount, sampleMetadata = sampleConditions, factorColors = conditionColor
  )

  p <- plotSamplePairCorrelation(barbieQ = object1)
  expect_s4_class(p, "Heatmap")

  SummarizedExperiment::colData(object1)$sampleMetadata$Treat <- factor(
    SummarizedExperiment::colData(object1)$sampleMetadata$Treat,
    levels = c("drug", "ctrl")
  )
  SummarizedExperiment::colData(object1)$sampleMetadata$Time <- factor(
    SummarizedExperiment::colData(object1)$sampleMetadata$Time,
    levels = c(2, 1)
  )
  sampleOrder <- c("Time", "Treat")
  p <- plotSamplePairCorrelation(barbieQ = object1, sampleOrder = sampleOrder)
  expect_equal(
    SummarizedExperiment::colData(object1)$sampleMetadata$Time[p@column_order] %>% table() %>% names(),
    c("2", "1")
  )
})
