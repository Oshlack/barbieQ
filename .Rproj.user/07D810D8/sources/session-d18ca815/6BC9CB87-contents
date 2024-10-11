test_that("plotting sample pair-wise correlation in Heatmap works", {
  ## Sample conditions and color palettes
  sampleConditions <- data.frame(
    Treat=factor(rep(c("ctrl", "drug"), each=6)),
    Time=rep(rep(1:2, each=3), 2))
  conditionColor <- list(
    Treat=c(ctrl="#999999", drug="#112233"),
    Time=c("1"="#778899", "2"="#998877"))

  ## Barcode count data
  nbarcodes <- 50
  nsamples <- 12
  barcodeCount <- abs(matrix(rnorm(nbarcodes*nsamples), nbarcodes, nsamples))
  rownames(barcodeCount) <- paste0("Barcode", 1:nbarcodes)

  object1 <- createBarbie(
    object=barcodeCount, target=sampleConditions, factorColors=conditionColor)

  p <- plotSamplePairCorrelation(Barbie = object1)
  expect_s4_class(p, "Heatmap")

  object1$metadata$Treat <- factor(
    object1$metadata$Treat, levels = c("drug", "ctrl"))
  object1$metadata$Time <- factor(
    object1$metadata$Time, levels = 2:1)
  sampleOrder <- c("Time", "Treat")
  p <- plotSamplePairCorrelation(Barbie = object1, sampleOrder = sampleOrder)
  expect_equal(object1$metadata$Time[p@column_order] %>% table() %>% names(),
               c("2", "1"))
})
