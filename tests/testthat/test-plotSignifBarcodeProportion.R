test_that("multiplication works", {
  Block <- c(1, 1, 2, 3, 3, 4, 1, 1, 2, 3, 3, 4)
  Treat <- factor(rep(c("ctrl", "drug"), each = 6))
  Time <- rep(rep(seq_len(2), each = 3), 2)
  nbarcodes <- 50
  nsamples <- 12
  count <- matrix(rnorm(nbarcodes * nsamples), nbarcodes, nsamples) |> abs()
  
  set.seed(2025)
  for (i in 1:10) {
    count[i, Treat == "ctrl"] <- 0
    count[i, Treat == "drug"] <- 1
  }
  for (i in 11:20) {
    count[i, Treat == "ctrl"] <- 1
    count[i, Treat == "drug"] <- 0
  }
  for (i in 21:30) {
    count[i, Treat == "ctrl"] <- 1
    count[i, Treat == "drug"] <- 1
  }
  for (i in 31:50) {
    count[i, Treat == "ctrl"] <- 0
    count[i, Treat == "drug"] <- 0
  }
  rownames(count) <- paste0("Barcode", seq_len(nbarcodes))
  barbieQ <- createBarbieQ(count, data.frame(Treat = Treat, Time = Time))
  
  testBB <- testBarcodeSignif(barbieQ, sampleGroup = "Treat")
  p <- plotSignifBarcodeProportion(barbieQ = testBB)
  
  expect_s3_class(p, "ggplot")
  expect_equal(length(p$layers), 3L)
  expect_equal(nrow(p$data), 3*ncol(testBB))
  plotPoint <- ggplot_build(p)$data[[1]]
  expect_true(all((unique(plotPoint$fill) %in% c("#33AAFF","#FF5959","#FFC000"))))
})
