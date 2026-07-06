local_edition(3)
noflevels <- c(2, 2, 2)
input.combs <- getRow(row.no = c(14, 17), noflevels + 1)

test_that("tests for findSupersets() have the same output", {
  expect_snapshot(findSupersets(input = 14, noflevels = noflevels + 1))
  expect_snapshot(findSupersets(input = input.combs, noflevels + 1))
})
