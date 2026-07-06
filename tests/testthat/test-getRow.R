local_edition(3)
noflevels <- c(2, 2, 2)
rows <- c(2, 4, 5, 7, 8, 10, 11, 13, 14, 16, 17)

mat <- getRow(rows, noflevels + 1)
rownames(mat) <- rows
colnames(mat) <- c("A", "B", "C")

rows2 <- findSubsets(input = 2, noflevels = noflevels + 1)
rows2_20 <- findSubsets(input = 2, noflevels = noflevels + 1, stop = 20)
rows879 <- findSubsets(input = c(8, 79), noflevels = rep(3, 4))
rows14 <- findSupersets(input = 14, noflevels = noflevels + 1)
x <- data.frame(createMatrix(rep(3, 4)), row.names = 1:81)
rows1417 <- findSupersets(input = c(14, 17), noflevels = noflevels + 1)
input.combs <- getRow(row.no = c(14, 17), noflevels + 1)

test_that("tests for getRow() have the same output", {
  expect_snapshot(mat)
  expect_snapshot(getRow(row.no = 2, noflevels = noflevels + 1) - 1)
  expect_snapshot(getRow(row.no = rows2, noflevels = noflevels + 1) - 1)
  expect_snapshot(getRow(row.no = rows2_20, noflevels = noflevels + 1) - 1)
  expect_snapshot(getRow(row.no = rows879, noflevels = rep(3, 4)) - 1)
  expect_snapshot(getRow(row.no = 14, noflevels = noflevels + 1) - 1)
  expect_snapshot(getRow(row.no = rows14, noflevels = noflevels + 1) - 1)
  expect_snapshot(getRow(row.no = c(14, 17), noflevels = noflevels + 1) - 1)
  expect_snapshot(getRow(row.no = rows1417, noflevels = noflevels + 1) - 1)
})
