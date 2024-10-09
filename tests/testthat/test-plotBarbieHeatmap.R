test_that("plotting Barbie Heatmap works", {
  HSC <- Barbie::HSC

  hp <- plotBarbieHeatmap(Barbie = HSC)
  expect_s4_class(hp, "Heatmap")
  expect_equal(hp@name, "log2 CPM+1")
  hp <- plotBarbieHeatmap(Barbie = HSC, splitSamples = TRUE)
  expect_equal(hp@top_annotation@anno_list %>% names(),
               c("treat", "mouse", "tissue", "lineage", "celltype", "PCRrep", "time"))
  hp <- plotBarbieHeatmap(Barbie = HSC, splitSamples = TRUE, sampleGroups = "treat")
  expect_equal(hp@bottom_annotation@anno_list %>% names(), "treat")

  HSC$metadata$treat <- factor(HSC$metadata$treat, levels = c("IT", "IF", "IV"))
  hp <- plotBarbieHeatmap(Barbie = HSC, splitSamples = TRUE, sampleGroups = "treat")
  expect_equal(hp@matrix_param[["column_split"]][[1]] %>% levels(),
               c("IT", "IF", "IV"))

  hp <- plotBarbieHeatmap(Barbie = HSC, value = "occurrence")
  expect_equal(hp@name, "occurrence")
})
