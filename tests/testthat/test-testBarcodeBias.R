test_that("testing differential proportions works", {
  Block <- c(1,1,2,3,3,4,1,1,2,3,3,4)
  Treat <- factor(rep(c("ctrl", "drug"), each=6))
  Time <- rep(rep(1:2, each=3), 2)
  nbarcodes <- 50
  nsamples <- 12
  count <- matrix(rnorm(nbarcodes*nsamples), nbarcodes, nsamples) %>% abs()
  rownames(count) <- paste0("Barcode", 1:nbarcodes)
  Barbie <- Barbie::createBarbie(count, data.frame(Treat=Treat, Time=Time))

  resultStat <- testDiffProp(
    Barbie = Barbie,
    mycontrasts = c(-1, 1, 0),
    contrastLevels = c("ctrl", "drug"),
    designMatrix = model.matrix(~0 + Treat + Time),
    transformation = "asin-sqrt",
    block = Block
    )
  expect_equal(rownames(resultStat), rownames(Barbie$proportion))

  resultStat2 <- testDiffProp(
    Barbie = Barbie,
    mycontrasts = c(0, 0, 1),
    designMatrix = model.matrix(~0 + Treat + Time),
    transformation = "asin-sqrt",
    block = Block
    )
  expect_equal(rownames(resultStat2), rownames(Barbie$proportion))
})

test_that("testing differential proportions works", {
  Block <- c(1,1,2,3,3,4,1,1,2,3,3,4)
  Treat <- factor(rep(c("ctrl", "drug"), each=6))
  Time <- rep(rep(1:2, each=3), 2)
  nbarcodes <- 50
  nsamples <- 12
  count <- matrix(rnorm(nbarcodes*nsamples), nbarcodes, nsamples) %>% abs()
  rownames(count) <- paste0("Barcode", 1:nbarcodes)
  Barbie <- Barbie::createBarbie(count, data.frame(Treat=Treat, Time=Time))

  resultStat <- testDiffOcc(
    Barbie = Barbie,
    regularization = "firth",
    mycontrasts = c(-1, 1, 0),
    contrastLevels = c("ctrl", "drug"),
    designMatrix = model.matrix(~0 + Treat + Time)
  )
  expect_equal(rownames(resultStat), rownames(Barbie$occurrence))

  # resultStat2 <- testDiffOcc(
  #   Barbie = Barbie,
  #   regularization = "firth",
  #   mycontrasts = c(0, 0, 1),
  #   designMatrix = model.matrix(~0 + Treat + Time)
  # )
  # expect_equal(rownames(resultStat2), rownames(Barbie$occurrence))
})

test_that("barcode test extracting correct arguments, dispatching right function", {
  Block <- c(1,1,2,3,3,4,1,1,2,3,3,4)
  Treat <- factor(rep(c("ctrl", "drug"), each=6))
  Time <- rep(rep(1:2, each=3), 2)
  nbarcodes <- 50
  nsamples <- 12
  count <- matrix(rnorm(nbarcodes*nsamples), nbarcodes, nsamples) %>% abs()
  rownames(count) <- paste0("Barcode", 1:nbarcodes)
  Barbie <- Barbie::createBarbie(count, data.frame(Treat=Treat, Time=Time))

  testBB1 <- testBarcodeBias(Barbie, sampleGroups = "Treat")
  testBB11 <- testBarcodeBias(Barbie, sampleGroups = "Treat", contrastLevels = c("ctrl", "drug"))
  testBB111 <- testBarcodeBias(Barbie, sampleGroups = "Treat", contrastLevels = c("drug", "ctrl"))
  expect_equal(testBB11$testBarcodes$diffProp_Treat$results$direction,
               testBB111$testBarcodes$diffProp_Treat$results$direction)

  testBB1111 <- testBarcodeBias(Barbie, sampleGroups = "Treat", contrastLevels = c("drug", "ctrl"),
                                designFormula = formula("~0 + Treat + Time"))
  testBB11111 <- testBarcodeBias(Barbie, sampleGroups = "Treat", contrastLevels = c("drug", "ctrl"),
                               designMatrix = model.matrix(~0 + Treat + Time))
  expect_equal(testBB1111$testBarcodes$diffProp_Treat$result,
               testBB11111$testBarcodes$diffProp_Treat$result)

  testBB2 <- testBarcodeBias(Barbie, sampleGroups = "Time")

  testBB3 <- testBarcodeBias(Barbie, sampleGroups = "Treat", method = "diffOcc")

  testBB <- testBarcodeBias(Barbie, sampleGroups = rep(1:4, each = 3))

  HSC <- Barbie::HSC
  BB <- createBarbie(object = HSC$assay, target = HSC$metadata)
  testBB4 <- testBarcodeBias(BB, sampleGroups = "treat", contrastLevels = c("IV", "IT"))
  testBB5 <- testBarcodeBias(Barbie = BB, sampleGroups = "mouse", contrastLevels = c("M13", "M14"))
})
