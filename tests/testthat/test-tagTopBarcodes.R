test_that("tagging top Barcodes overall works", {
  dat <- Barbie::HSC
  recreatedDat <- tagTopBarcodes(dat)
  expect_equal(dat$assay, recreatedDat$assay)
  expect_equal(dat$metadata, recreatedDat$metadata)
  expect_equal(dat$proportion, recreatedDat$proportion)
  expect_equal(dat$CPM, recreatedDat$CPM)
  expect_equal(dat$occurrence, recreatedDat$occurrence)
  expect_equal(dat$rank, recreatedDat$rank)
  expect_equal(dat$isTop, recreatedDat$isTop)
  # expect_equal(dat$clusters, recreatedDat$clusters)
  expect_equal(dat$factorColors, recreatedDat$factorColors)
})

test_that("tagging top Barcodes in one sample works", {
  vec <- c(1,2,3,4,NA,NA)
  top <- c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE)
  expect_equal(top, tagTopEachColumn(vec))
  vecDump <- c(1,1.5,2,96,NA,NA)
  topDump <- c(FALSE, TRUE, TRUE, TRUE, FALSE, FALSE)
  expect_equal(topDump, tagTopEachColumn(vecDump))
})
