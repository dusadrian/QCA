local_edition(3)
set.seed(12345)
x <- sample(1:100, size = 15)
gdp <- c(460, 500, 900, 2000, 2100, 2400, 15000, 16000, 20000)

test_that("tests for findTh() have the same output", {
  expect_snapshot(findTh(x))
  expect_snapshot(findTh(x, groups = 3))
  expect_snapshot(findTh(gdp))
  expect_snapshot(findTh(gdp, n = 2))
  expect_snapshot(findTh(gdp, n = 2, hclustm = "complete", distm = "euclidean"))
})
