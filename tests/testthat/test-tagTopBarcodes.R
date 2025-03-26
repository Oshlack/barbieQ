test_that("tagging top Barcodes across samples works", {
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
  barcodeCount[50, ] <- 0
  rownames(barcodeCount) <- paste0("Barcode", seq_len(nbarcodes))
  ## create a `barbieQ` object
  object1 <- createBarbieQ(barcodeCount, sampleConditions, conditionColor)
  object2 <- tagTopBarcodes(barbieQ = object1)
  ## check `vec` and `mat` in correct format
  expect_equal(SummarizedExperiment::rowData(object2)$isTopBarcode$isTop, 
               c(rep(TRUE, 49), FALSE), ignore_attr = TRUE)
  expect_equal(SummarizedExperiment::assays(object2)$isTopAssay[50, ], rep(FALSE, 12), ignore_attr = TRUE)
  ## check `activeSamples` correctly choose samples to consider
  barcodeCount[, seq_len(6)] <- 0
  object1 <- createBarbieQ(barcodeCount, sampleConditions, conditionColor)
  object3 <- tagTopBarcodes(
    barbieQ = object1, activeSamples = rep(c(TRUE, FALSE), each = 6)
  )
  expect_equal(SummarizedExperiment::rowData(object3)$isTopBarcode$isTop, rep(FALSE, 50))
  ## cehck `proportionThreshold` setting correct cutoff for *top* Barcodes
  barcodeCount[, 1] <- 10
  object1 <- createBarbieQ(barcodeCount, sampleConditions, conditionColor)
  object4 <- tagTopBarcodes(
    barbieQ = object1, activeSamples = rep(c(TRUE, FALSE), each = 6),
    proportionThreshold = 0.5
  )
  expect_equal(SummarizedExperiment::rowData(object4)$isTopBarcode$isTop, rep(c(TRUE, FALSE), each = 25))
  ## check `nSampleThreshold` setting up correct minimum frequency of *top*
  object5 <- tagTopBarcodes(
    barbieQ = object1, activeSamples = rep(c(TRUE, FALSE), each = 6),
    proportionThreshold = 0.5, nSampleThreshold = 2
  )
  expect_equal(SummarizedExperiment::rowData(object5)$isTopBarcode$isTop, rep(FALSE, 50))
})

test_that("tagging top Barcodes in one sample works", {
  vec <- c(1, 2, 3, 4, NA, NA)
  top <- c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE)
  ## check NAs are not tagged as *top*
  expect_equal(tagTopEachColumn(vec), top)
  vecDump <- c(1, 1.5, 2, 96, NA, NA)
  topDump <- c(FALSE, TRUE, TRUE, TRUE, FALSE, FALSE)
  ## check low proportions are not tagged as *top*
  expect_equal(tagTopEachColumn(vecDump), topDump)
  ## case when all zeros in a vector - return all FALSE
  zeros <- rep(0, 10)
  expect_equal(tagTopEachColumn(zeros), rep(FALSE, 10))
})
