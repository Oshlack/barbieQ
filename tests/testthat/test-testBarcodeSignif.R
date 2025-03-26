test_that("testing differential proportions works", {
  Block <- c(1, 1, 2, 3, 3, 4, 1, 1, 2, 3, 3, 4)
  Treat <- factor(rep(c("ctrl", "drug"), each = 6))
  Time <- rep(rep(1:2, each = 3), 2)
  nbarcodes <- 50
  nsamples <- 12
  set.seed(2025)
  count <- matrix(rnorm(nbarcodes * nsamples), nbarcodes, nsamples) |> abs()
  for (i in 1:20) {
    count[i, Treat == "ctrl"] <- rnorm(sum(Treat == "ctrl"), mean = 10, sd = 2)  # ctrl around 10
    count[i, Treat == "drug"] <- rnorm(sum(Treat == "drug"), mean = 20, sd = 2)  # drug around 20
  }
  for (i in 21:40) {
    count[i, Treat == "ctrl"] <- rnorm(sum(Treat == "ctrl"), mean = 20, sd = 2)  # ctrl around 10
    count[i, Treat == "drug"] <- rnorm(sum(Treat == "drug"), mean = 10, sd = 2)  # drug around 20
  }
  for (i in 41:50) {
    count[i, Treat == "ctrl"] <- rnorm(sum(Treat == "ctrl"), mean = 10, sd = 2)  # ctrl around 10
    count[i, Treat == "drug"] <- rnorm(sum(Treat == "drug"), mean = 10, sd = 2)  # drug around 20
  }
  
  rownames(count) <- paste0("Barcode", 1:nbarcodes)
  barbieQ <- createBarbieQ(count, data.frame(Treat = Treat, Time = Time))

  resultStat <- testDiffProp(
    proportion = SummarizedExperiment::assays(barbieQ)$proportion,
    mycontrasts = setNames(c(1,-1,0), c("Treatctrl", "Treatdrug", "Time")),
    designMatrix = model.matrix(~ 0 + Treat + Time),
    transformation = "asin-sqrt",
    block = Block
  )
  expect_equal(rownames(resultStat), rownames(barbieQ))
  
  # expect_equal(resultStat$direction, c(rep(-1,20), rep(1, 20), rep(0, 10)), ignore_attr = TRUE)

  resultStat2 <- testDiffProp(
    proportion = SummarizedExperiment::assays(barbieQ)$proportion,
    mycontrasts = c(0, 0, 1),
    designMatrix = model.matrix(~ 0 + Treat + Time),
    transformation = "asin-sqrt",
    block = Block
  )
  expect_equal(rownames(resultStat2), rownames(barbieQ))
  
  expect_equal(resultStat2$direction, rep(0, 50), ignore_attr = TRUE)
  
})

test_that("testing differential occurrence works", {
  Block <- c(1, 1, 2, 3, 3, 4, 1, 1, 2, 3, 3, 4)
  Treat <- factor(rep(c("ctrl", "drug"), each = 6))
  Time <- rep(rep(1:2, each = 3), 2)
  nbarcodes <- 50
  nsamples <- 12
  count <- matrix(rnorm(nbarcodes * nsamples), nbarcodes, nsamples) |> abs()
  
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
  
  rownames(count) <- paste0("Barcode", 1:nbarcodes)
  barbieQ <- createBarbieQ(count, data.frame(Treat = Treat, Time = Time))

  resultStat <- testDiffOcc(
    occurrence = SummarizedExperiment::assays(barbieQ)$occurrence,
    regularization = "firth",
    mycontrasts = setNames(c(1,-1,0), c("Treatctrl", "Treatdrug", "Time")),
    designMatrix = model.matrix(~ 0 + Treat + Time)
  )
  expect_equal(rownames(resultStat$statsDf), rownames(barbieQ))
  expect_equal(resultStat$statsDf$direction, c(rep(-1, 10), rep(1, 10), rep(0, 30)))
  expect_equal(colnames(resultStat$design), c("(Intercept)", "Treatctrl", "Time"))
  expect_equal(resultStat$variables, list("Treatdrug", "Treatctrl", "Time"), ignore_attr = TRUE)

  ## check one factor
  resultStat <- testDiffOcc(
    occurrence = SummarizedExperiment::assays(barbieQ)$occurrence,
    regularization = "firth",
    mycontrasts = setNames(c(-1,1), c("Treatctrl", "Treatdrug")),
    designMatrix = model.matrix(~ 0 + Treat)
  )
  # expect_equal(resultStat$statsDf$direction, c(rep(1, 10), rep(-1, 10), rep(0, 30)))

  resultStat <- testDiffOcc(
    occurrence = SummarizedExperiment::assays(barbieQ)$occurrence,
    regularization = "firth",
    mycontrasts = setNames(c(1), c("Time")),
    designMatrix = model.matrix(~ 0 + Time)
  )
  expect_equal(resultStat$statsDf$direction, rep(0, 50))

  ## this will cause warings, because of
  resultStat <- testDiffOcc(
    occurrence = SummarizedExperiment::assays(barbieQ)$occurrence,
    regularization = "firth",
    mycontrasts = setNames(c(1,0,0), c("Time", "Treatctrl", "Treatdrug")),
    designMatrix = model.matrix(~ 0 + Time + Treat)
  )
  expect_equal(resultStat$statsDf$direction, rep(0, 50), ignore_attr = TRUE)
})

test_that("barcode test extracting correct arguments, dispatching right function", {
  Block <- c(1, 1, 2, 3, 3, 4, 1, 1, 2, 3, 3, 4)
  Treat <- factor(rep(c("ctrl", "drug"), each = 6))
  Time <- rep(rep(seq_len(2), each = 3), 2)
  nbarcodes <- 50
  nsamples <- 12
  ## Generate the base counts
  base_counts <- rnorm(nbarcodes, mean = 10, sd = 5) |> abs()
  count <- matrix(0, nbarcodes, nsamples)
  count <- matrix(rnorm(nbarcodes * nsamples), nbarcodes, nsamples) |> abs()
  set.seed(2025)
  for (i in 1:20) {
    count[i, Treat == "ctrl"] <- rnorm(sum(Treat == "ctrl"), mean = 10, sd = 2)  # ctrl around 10
    count[i, Treat == "drug"] <- rnorm(sum(Treat == "drug"), mean = 20, sd = 2)  # drug around 20
  }
  for (i in 21:40) {
    count[i, Treat == "ctrl"] <- rnorm(sum(Treat == "ctrl"), mean = 20, sd = 2)  # ctrl around 10
    count[i, Treat == "drug"] <- rnorm(sum(Treat == "drug"), mean = 10, sd = 2)  # drug around 20
  }
  for (i in 41:50) {
    count[i, Treat == "ctrl"] <- rnorm(sum(Treat == "ctrl"), mean = 10, sd = 2)  # ctrl around 10
    count[i, Treat == "drug"] <- rnorm(sum(Treat == "drug"), mean = 10, sd = 2)  # drug around 20
  }
  # ## Decrease counts based on time
  # for (i in seq_len(nbarcodes)) {
  #   count[i, ] <- base_counts[i] * (1 - (Time - 1) * 0.1)
  # }
  rownames(count) <- paste0("Barcode", seq_len(nbarcodes))
  barbieQ <- createBarbieQ(count, data.frame(Treat = Treat, Time = Time))

  testBB1 <- testBarcodeSignif(barbieQ, sampleGroup = "Treat")
  testBB1 <- testBarcodeSignif(barbieQ, sampleMetadata = data.frame(Treat = Treat) ,sampleGroup = "Treat")
  testBB11 <- testBarcodeSignif(
    barbieQ,
    sampleGroup = "Treat",
    contrastFormula = ("Treatdrug - Treatctrl")
  )
  expect_equal(SummarizedExperiment::rowData(testBB11)$testingBarcode$tendencyTo,
               c(rep("Treatdrug", 20), rep("Treatctrl", 20), rep("n.s.", 10)))
  testBB111 <- testBarcodeSignif(
    barbieQ,
    sampleGroup = "Treat",
    contrastFormula = ("Treatctrl - Treatdrug")
  )
  expect_equal(
    SummarizedExperiment::rowData(testBB11)$testingBarcode$tendencyTo,
    SummarizedExperiment::rowData(testBB111)$testingBarcode$tendencyTo
  )

  testBB1111 <- testBarcodeSignif(barbieQ,
    sampleGroup = "Treat", contrastFormula = ("Treatctrl - Treatdrug"),
    designFormula = formula("~0 + Treat + Time")
  )
  testBB11111 <- testBarcodeSignif(barbieQ,
    sampleGroup = "Treat", contrastFormula = ("Treatctrl - Treatdrug"),
    designMatrix = model.matrix(~ 0 + Treat + Time)
  )
  expect_equal(
    SummarizedExperiment::rowData(testBB1111)$testingBarcode,
    SummarizedExperiment::rowData(testBB11111)$testingBarcode
  )

  ## confirm `sampleMetadata` is updated; design matrix and formula are updated accordingly
  testBB12 <- testBarcodeSignif(barbieQ,
    sampleGroup = "Treat", contrastFormula = ("Treatctrl - Treatdrug"),
    sampleMetadata = data.frame(Treat = Treat)
  )
  results <- SummarizedExperiment::rowData(testBB12)$testingBarcode
  expect_equal(colnames(S4Vectors::metadata(results)$design), c("Treatctrl", "Treatdrug"))
  expect_equal(
    S4Vectors::metadata(results)$contrasts,
    matrix(c(1, -1), ncol=1), ignore_attr = TRUE
  )

  ## with specified `designMatrix` used for test, should expect not using default formula
  testBB13 <- testBarcodeSignif(barbieQ,
    sampleGroup = "Treat", contrastFormula = ("Treatctrl - Treatdrug"),
    designMatrix = model.matrix(~ 0 + Treat)
  )
  results <- SummarizedExperiment::rowData(testBB13)$testingBarcode
  expect_equal(
    S4Vectors::metadata(results)$design,
    model.matrix(~ 0 + Treat),
    ignore_attr = TRUE
  )

  ## with specified `designFormula`, should expect updated formula and `designMatrix`
  testBB14 <- testBarcodeSignif(barbieQ,
    sampleGroup = "Treat", contrastFormula = ("Treatctrl - Treatdrug"),
    designFormula = formula("~ 0 + Treat")
  )
  results <- SummarizedExperiment::rowData(testBB14)$testingBarcode
  expect_equal(
    S4Vectors::metadata(results)$design,
    model.matrix(~ 0 + Treat),
    ignore_attr = TRUE
  )

  testBB2 <- testBarcodeSignif(barbieQ, sampleGroup = "Time")
  testBB22 <- testBarcodeSignif(barbieQ, sampleMetadata = data.frame(Time = Time) ,sampleGroup = "Time")

  testBB3 <- testBarcodeSignif(
    barbieQ,
    sampleGroup = "Time", method = "diffOcc",
    designFormula = formula("~ 0 + Time + Treat")
  )

  testBB4 <- testBarcodeSignif(
    barbieQ,
    sampleGroup = "Time", method = "diffOcc",
    designFormula = formula("~ 0 + Time")
  )

  testBB <- testBarcodeSignif(barbieQ, sampleGroup = rep(seq_len(4), each = 3))
  
  ## handling comlex design
  testBB5 <- testBarcodeSignif(
    barbieQ, sampleMetadata = data.frame(Treat = as.factor(rep(seq_len(4), each = 3))),
    sampleGroup = "Treat", contrastLevels = c(1,2),
    designFormula = formula("~ 0 + Treat")
  )
})
