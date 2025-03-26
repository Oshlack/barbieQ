test_that("plot Barcode test results - scatter plot - works", {
  Block <- c(1, 1, 2, 3, 3, 4, 1, 1, 2, 3, 3, 4)
  Treat <- factor(rep(c("ctrl", "drug"), each = 6))
  Time <- rep(rep(seq_len(2), each = 3), 2)
  nbarcodes <- 50
  nsamples <- 12
  count <- matrix(rnorm(nbarcodes * nsamples), nbarcodes, nsamples) |> abs()
  rownames(count) <- paste0("Barcode", seq_len(nbarcodes))
  barbieQ <- createBarbieQ(count, data.frame(Treat = Treat, Time = Time))

  testBB <- testBarcodeBias(barbieQ, sampleGroup = "Treat")
  p <- plotBarcodePValue(barbieQ = testBB)

  expect_s3_class(p, "ggplot")
  expect_equal(length(p$layers), 3L)
  expect_equal(nrow(p$data), nrow(testBB))
  plotPoint <- ggplot_build(p)$data[[1]]
  expect_equal(plotPoint$y, -log10(p$data$adj.P.Val))
  expect_equal(plotPoint$x, -p$data$avgRank)
})

test_that("plot Barcode test results - MA plot - works", {
  Block <- c(1, 1, 2, 3, 3, 4, 1, 1, 2, 3, 3, 4)
  Treat <- factor(rep(c("ctrl", "drug"), each = 6))
  Time <- rep(rep(seq_len(2), each = 3), 2)
  nbarcodes <- 50
  nsamples <- 12
  count <- matrix(rnorm(nbarcodes * nsamples), nbarcodes, nsamples) |> abs()
  rownames(count) <- paste0("Barcode", seq_len(nbarcodes))
  barbieQ <- createBarbieQ(count, data.frame(Treat = Treat, Time = Time))
  
  testBB <- testBarcodeBias(barbieQ, sampleGroup = "Treat")
  p <- plotBarcodeMA(barbieQ = testBB)
  
  expect_s3_class(p, "ggplot")
  expect_equal(length(p$layers), 1L)
  expect_equal(nrow(p$data), nrow(testBB))
  plotPoint <- ggplot_build(p)$data[[1]]
  expect_equal(plotPoint$y, p$data$meanDiff)
  expect_equal(plotPoint$x, p$data$Amean)
  
  testBB <- testBarcodeBias(barbieQ, sampleGroup = "Treat", method = "diffOcc")
  p <- plotBarcodeMA(barbieQ = testBB)
  
  expect_s3_class(p, "ggplot")
  expect_equal(length(p$layers), 1L)
  expect_equal(nrow(p$data), nrow(testBB))
  plotPoint <- ggplot_build(p)$data[[1]]
  expect_equal(plotPoint$y, p$data$logOR)
  expect_equal(plotPoint$x, p$data$totalOcc)
})

test_that("plot Barcode test results - heatmap works", {
  Block <- c(1, 1, 2, 3, 3, 4, 1, 1, 2, 3, 3, 4)
  Treat <- factor(rep(c("ctrl", "drug"), each = 6))
  Time <- rep(rep(seq_len(2), each = 3), 2)
  nbarcodes <- 50
  nsamples <- 12
  count <- matrix(rnorm(nbarcodes * nsamples), nbarcodes, nsamples) |> abs()
  rownames(count) <- paste0("Barcode", seq_len(nbarcodes))
  barbieQ <- createBarbieQ(count, data.frame(Treat = Treat, Time = Time))

  testBB <- testBarcodeBias(barbieQ, sampleGroup = "Treat")

  hp <- plotSignifBarcodeHeatmap(testBB)
  expect_s4_class(hp, "Heatmap")
  expect_equal(
    levels(as.factor(hp@matrix_param$column_split$testingBarcode)),
    c("Treatctrl", "Treatdrug")
  )

})
