test_that("plotting Barcode average proportion works", {
  ## sample conditions and color palettes
  sampleConditions <- data.frame(
    Treat=factor(rep(c("ctrl", "drug"), each=6)),
    Time=rep(rep(1:2, each=3), 2))
  conditionColor <- list(
    Treat=c(ctrl="#999999", drug="#112233"),
    Time=c("1"="#778899", "2"="#998877"))
  ## Barcode count data
  nbarcodes <- 50
  nsamples <- 12
  barcodeCount <- abs(matrix(10, nbarcodes, nsamples))
  barcodeCount[21:50,] <- 0.0001
  rownames(barcodeCount) <- paste0("Barcode", 1:nbarcodes)
  ## create a `Barbie` object
  object1 <- createBarbie(barcodeCount, sampleConditions, conditionColor)
  p <- plotBarcodeProportion(object1)
  ## check correct ggplot object
  expect_s3_class(p, "ggplot")
  expect_equal(sum(p$data$percentage), 100)
  expect_equal(mean(p$data$rank), mean(1:nrow(barcodeCount)))
})
