test_that("plot Barcode test results - scatter plot - works", {
  Block <- c(1,1,2,3,3,4,1,1,2,3,3,4)
  Treat <- factor(rep(c("ctrl", "drug"), each=6))
  Time <- rep(rep(1:2, each=3), 2)
  nbarcodes <- 50
  nsamples <- 12
  count <- matrix(rnorm(nbarcodes*nsamples), nbarcodes, nsamples) %>% abs()
  rownames(count) <- paste0("Barcode", 1:nbarcodes)
  Barbie <- createBarbie(count, data.frame(Treat=Treat, Time=Time))

  testBB <- testBarcodeBias(Barbie, sampleGroups = "Treat")
  p <- plotBarcodeBiasScatterPlot(
    Barbie = testBB, elementName = "diffProp_Treat")

  expect_s3_class(p, "ggplot")
  expect_equal(length(p$layers), 3L)
  expect_equal(nrow(p$data), nrow(testBB$proportion))
  plotPoint <- ggplot_build(p)$data[[1]]
  expect_equal(plotPoint$y, p$data$minusLogP)
  expect_equal(plotPoint$x, -p$data$avgRank)
})

test_that("plot Barcode test results - heatmap works", {
  Block <- c(1,1,2,3,3,4,1,1,2,3,3,4)
  Treat <- factor(rep(c("ctrl", "drug"), each=6))
  Time <- rep(rep(1:2, each=3), 2)
  nbarcodes <- 50
  nsamples <- 12
  count <- matrix(rnorm(nbarcodes*nsamples), nbarcodes, nsamples) %>% abs()
  rownames(count) <- paste0("Barcode", 1:nbarcodes)
  Barbie <- createBarbie(count, data.frame(Treat=Treat, Time=Time))

  testBB <- testBarcodeBias(Barbie, sampleGroups = "Treat")

  hp <- plotBarcodeBiasHeatmap(testBB)
  expect_s4_class(hp, "Heatmap")
  expect_equal(levels(hp@matrix_param$column_split$Treat),
               c("ctrl", "drug"))

  testBB <- testBarcodeBias(
    Barbie, sampleGroups = "Treat", contrastLevels = c("drug", "ctrl"))
  hp <- plotBarcodeBiasHeatmap(testBB)
  expect_equal(levels(hp@matrix_param$column_split$Treat),
               c("drug", "ctrl"))
})

