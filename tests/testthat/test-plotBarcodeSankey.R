test_that("plotting total Barcode proportion in Sankey plot works", {
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
  p <- plotBarcodeSankey(object1)
  expect_s3_class(p, "ggplot")
  ## check correct frequency of *bottom* Barcodes
  expect_equal(p$data[3, 3], 60)
  ## check correct total proportion of *bottom* Barcodes
  expect_equal(
    p$data[4, 3], 0.0001 * 30 / (0.0001 * 30 + 10 * 20) * 100,
    tolerance = 1e-10
  )
})
