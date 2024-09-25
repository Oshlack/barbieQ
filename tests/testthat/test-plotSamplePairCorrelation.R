test_that("plotting sample pair-wise correlation in Heatmap works", {
  HSC <- Barbie::HSC
  p <- plotSamplePairCorrelation(Barbie = HSC)
  expect_s4_class(p, "Heatmap")

  HSC$metadata$treat <- factor(HSC$metadata$treat, levels = c("IF", "IT", "IV"))
  HSC$metadata$lineage <- factor(HSC$metadata$lineage, levels = c("Bcell", "Tcell", "immature", "Myeloid"))
  sampleOrder <- c("time", "treat", "lineage", "celltype", "tissue", "mouse")
  p <- plotSamplePairCorrelation(Barbie = HSC, sampleOrder = sampleOrder)
  expect_equal(HSC$metadata$treat[p@column_order] %>% table() %>% names(),
               c("IF", "IT", "IV"))
})
