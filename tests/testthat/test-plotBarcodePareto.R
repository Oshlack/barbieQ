test_that("plotting circular bar plot of Barcode proportion works.", {
  ## sample conditions and color palettes
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
  barcodeCount <- abs(matrix(10, nbarcodes, nsamples))
  barcodeCount[seq(21, 50), ] <- 0.0001
  rownames(barcodeCount) <- paste0("Barcode", seq_len(nbarcodes))
  ## create a `barbieQ` object
  object1 <- createBarbieQ(barcodeCount, sampleConditions, conditionColor)
  object1 <- tagTopBarcodes(object1)
  p <- plotBarcodePareto(object1)
  expect_s3_class(p, "ggplot")
  ## check correct *top* and *bottom* numbers, adding 5 spaces
  expect_equal(table(p$data$group), c(25, 35), ignore_attr = TRUE)
})
