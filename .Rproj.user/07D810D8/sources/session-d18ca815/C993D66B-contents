test_that("plot Barcode test results - scatter plot - works", {
  HSC <- Barbie::HSC
  BB <- createBarbie(object = HSC$assay, target = HSC$metadata)
  testBB <- testBarcodeBias(BB, groupBy = "treat", contrastLevels = c("IV", "IT"))
  plotBarcodeBiasScatterPlot(Barbie = testBB)

  testBB <- testBarcodeBias(testBB, groupBy = "lineage", contrastLevels = c("Myeloid", "Tcell"))
  p <- plotBarcodeBiasScatterPlot(Barbie = testBB, elementName = "diffProp_lineage")

  expect_s3_class(p, "ggplot")
  expect_equal(length(p$layers), 3L)
  expect_equal(nrow(p$data), nrow(testBB$proportion))
  plotPoint <- ggplot_build(p)$data[[1]]
  expect_equal(plotPoint$y, p$data$minusLogP)
  expect_equal(plotPoint$x, -p$data$avgRank)
  # testBB <- testBarcodeBias(BB, groupBy = "treat", method = "diffOcc",
  #                           contrastLevels = c("IV", "IT"))
})

test_that("plot Barcode test results - heatmap works", {
  HSC <- Barbie::HSC
  testBB <- testBarcodeBias(HSC , groupBy = "treat", contrastLevels = c("IV", "IT"))

  hp <- plotBarcodeBiasHeatmap(testBB)
  expect_s4_class(hp, "Heatmap")
  expect_equal(levels(hp@matrix_param$column_split$treat), c("IV", "IT", "IF", "T0"))

  testBB <- testBarcodeBias(HSC , groupBy = "treat", contrastLevels = c("IF", "IV"))
  hp <- plotBarcodeBiasHeatmap(testBB)
  expect_equal(levels(hp@matrix_param$column_split$treat), c("IF", "IV", "IT", "T0"))
})

