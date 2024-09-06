test_that("extracting count matrix, target file and factor colors works", {
  dat <- Barbie::HSC
  recreatedDat <- createBarbie(object=dat$assay, target=dat$metadata, factorColors=dat$factorColors)
  expect_equal(dat$assay, recreatedDat$assay)
  expect_equal(dat$metadata, recreatedDat$metadata)
  expect_equal(dat$proportion, recreatedDat$proportion)
  expect_equal(dat$CPM, recreatedDat$CPM)
  expect_equal(dat$occurrence, recreatedDat$occurrence)
  expect_equal(dat$rank, recreatedDat$rank)
  # expect_equal(dat$isTop, recreatedDat$isTop)
  # expect_equal(dat$clusters, recreatedDat$clusters)
  expect_equal(dat$factorColors, recreatedDat$factorColors)
})
