test_that("plotting Barcode average proportion works", {
  HSC <- Barbie::HSC
  p <- plotBarcodeProportion(HSC)
  expect_s3_class(p, "ggplot")
  expect_equal(sum(p$data$percentage), 100)
  expect_equal(mean(p$data$rank), mean(1:nrow(HSC$proportion)))
})
