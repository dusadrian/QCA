local_edition(3)
data(d.represent, package = "QCA")
data(d.jobsecurity, package = "QCA")
data(d.partybans, package = "QCA")
data(CVF, package = "QCA")
chart <- makeChart("a, b, ~c", "abc, a~b~c, a~bc, ~ab~c")
fmin <- unclass(findmin(chart))
fmin.conservative <- unclass(findmin(chart, type = "conservative"))
fmin.hybrid <- unclass(findmin(chart, type = "hybrid"))
fmin.hybrid.solind <- unclass(findmin(chart, type = "hybrid", solind = TRUE))
attributes(fmin) <- NULL
attributes(fmin.conservative) <- NULL
attributes(fmin.hybrid) <- NULL
attributes(fmin.hybrid.solind) <- NULL

represent.chart <- minimize(d.represent, outcome = "WNP", include = "?", details = TRUE)$PIchart
jobsecurity.chart <- minimize(
  d.jobsecurity,
  outcome = "JSR",
  incl.cut = 0.9,
  include = "?",
  details = TRUE
)$PIchart
partybans.chart <- minimize(
  d.partybans,
  outcome = "PB[1]",
  conditions = "C, F, T, V",
  include = "?",
  dir.exp = "1, 1;2, 2, 1",
  details = TRUE
)$PIchart
ttCVF.chart <- minimize(
  truthTable(CVF, outcome = PROTEST, incl.cut = 0.8, show.cases = TRUE, sort.by = "incl, n"),
  include = "?",
  row.dom = TRUE,
  details = TRUE,
  show.cases = TRUE
)$PIchart

representative.charts <- list(
  chart1 = makeChart(
    c("A", "B", "~C", "D"),
    c("A~BCD", "A~BC~D", "~AB~C~D", "~ABCD"),
    snames = c("A", "B", "C", "D")
  ),
  chart3 = makeChart(
    c("AB", "BC", "A~C", "~AC", "~A~B~D", "~B~C~D"),
    c(
      "ABCD", "ABC~D", "AB~CD", "AB~C~D", "A~B~CD", "A~B~C~D",
      "~ABCD", "~ABC~D", "~A~BCD", "~A~BC~D", "~A~B~C~D"
    ),
    snames = c("A", "B", "C", "D")
  ),
  represent = represent.chart,
  jobsecurity = jobsecurity.chart,
  partybans = partybans.chart,
  ttCVF = ttCVF.chart
)

hybrid_cardinality <- function(chart) {
  sum(unclass(findmin(chart, type = "hybrid", solind = TRUE)) > 0)
}

test_that("findmin() works", {
  expect_equal(fmin, 2)
  expect_equal(fmin.conservative, 2)
  expect_equal(fmin.hybrid, 2)
  expect_true(is.numeric(fmin.hybrid.solind))
  expect_equal(sum(fmin.hybrid.solind > 0), 2)
})

test_that("the conservative type replaces the former exact label", {
  expect_error(findmin(chart, type = "exact"), "should be one of")
})

test_that("findmin() defaults to the hybrid backend", {
  .Call("C_resetScpProfile", PACKAGE = "QCA")
  default <- unclass(findmin(chart, solind = TRUE))
  profile <- .Call("C_getScpProfile", PACKAGE = "QCA")

  expect_equal(default, fmin.hybrid.solind)
  expect_gt(profile$original_columns, 0)
})

test_that("hybrid matches conservative on representative PI charts", {
  for (nm in names(representative.charts)) {
    chart <- representative.charts[[nm]]
    conservative_min <- unclass(findmin(chart, type = "conservative"))
    hybrid_min <- unclass(findmin(chart, type = "hybrid"))
    hybrid_solind_min <- hybrid_cardinality(chart)

    expect_equal(
      hybrid_min,
      conservative_min,
      info = sprintf("%s conservative vs hybrid minimum", nm)
    )
    expect_equal(
      hybrid_solind_min,
      conservative_min,
      info = sprintf("%s hybrid solind minimum", nm)
    )
  }
})

test_that("hybrid probes SCP then falls back without coupling conservative", {
  set.seed(1)
  native.chart <- matrix(FALSE, nrow = 80, ncol = 200)
  for (j in seq_len(ncol(native.chart))) {
    native.chart[sample.int(nrow(native.chart), 4), j] <- TRUE
  }
  chart <- t(native.chart)

  .Call("C_resetScpProfile", PACKAGE = "QCA")
  conservative <- unclass(findmin(chart, type = "conservative", gurobi = FALSE, solind = TRUE))
  conservative.profile <- .Call("C_getScpProfile", PACKAGE = "QCA")

  expect_equal(conservative.profile$nodes, 0)
  expect_equal(conservative.profile$original_columns, 0)

  .Call("C_resetScpProfile", PACKAGE = "QCA")
  hybrid <- unclass(findmin(chart, type = "hybrid", solind = TRUE))
  hybrid.profile <- .Call("C_getScpProfile", PACKAGE = "QCA")

  expect_equal(sum(hybrid), sum(conservative))
  expect_true(hybrid.profile$scp_attempted)
  expect_true(hybrid.profile$scp_limited)
  expect_true(hybrid.profile$lpsolve_fallback)
  expect_equal(hybrid.profile$nodes, 5000)
  expect_equal(hybrid.profile$core_columns, nrow(chart))
})

test_that("hybrid removes dominated low-row coverage masks before level 2", {
  set.seed(824)
  nr <- 12L
  base <- matrix(FALSE, nrow = nr, ncol = 40L)
  seen <- character()
  for (j in seq_len(ncol(base))) {
    repeat {
      rows <- sort(sample.int(nr, 6L))
      key <- paste(rows, collapse = ",")
      if (!key %in% seen) break
    }
    seen <- c(seen, key)
    base[rows, j] <- TRUE
  }

  source.column <- rep(seq_len(ncol(base)), each = 4L)
  dominated <- base[, source.column, drop = FALSE]
  for (j in seq_len(ncol(dominated))) {
    covered.rows <- which(dominated[, j])
    dominated[sample(covered.rows, 1L), j] <- FALSE
  }
  native.chart <- cbind(base, dominated)
  chart <- t(native.chart)

  conservative <- unclass(findmin(
    chart, type = "conservative", gurobi = FALSE, solind = TRUE
  ))
  .Call("C_resetScpProfile", PACKAGE = "QCA")
  hybrid <- unclass(findmin(chart, type = "hybrid", solind = TRUE))
  profile <- .Call("C_getScpProfile", PACKAGE = "QCA")

  expect_equal(sum(hybrid), sum(conservative))
  expect_true(all(colSums(native.chart[, hybrid > 0, drop = FALSE]) > 0))
  expect_equal(profile$original_columns, ncol(native.chart))
  expect_equal(profile$presolved_columns, ncol(base))
  expect_false(any(which(conservative > 0) > ncol(base)))
  expect_false(any(which(hybrid > 0) > ncol(base)))
})

test_that("hybrid coverage presolve is not restricted to 20 rows", {
  set.seed(825)
  nr <- 30L
  base <- matrix(FALSE, nrow = nr, ncol = 20L)
  for (j in seq_len(ncol(base))) {
    base[sample.int(nr, 15L), j] <- TRUE
  }
  source.column <- rep(seq_len(ncol(base)), each = 3L)
  dominated <- base[, source.column, drop = FALSE]
  for (j in seq_len(ncol(dominated))) {
    covered.rows <- which(dominated[, j])
    dominated[covered.rows[(j - 1L) %% length(covered.rows) + 1L], j] <- FALSE
  }
  native.chart <- cbind(base, dominated)
  chart <- t(native.chart)

  conservative <- unclass(findmin(
    chart, type = "conservative", gurobi = FALSE, solind = TRUE
  ))
  .Call("C_resetScpProfile", PACKAGE = "QCA")
  hybrid <- unclass(findmin(chart, type = "hybrid", solind = TRUE))
  profile <- .Call("C_getScpProfile", PACKAGE = "QCA")

  expect_equal(sum(hybrid), sum(conservative))
  expect_equal(profile$original_columns, ncol(native.chart))
  expect_equal(profile$presolved_columns, ncol(base))
  expect_false(any(which(conservative > 0) > ncol(base)))
})

test_that("hybrid extends a productive exact-finisher probe", {
  set.seed(2)
  nr <- 36L
  nc <- 800L
  chart <- matrix(FALSE, nrow = nc, ncol = nr)
  for (column in seq_len(nc)) {
    degree <- sample(3L:8L, 1L)
    chart[column, sample.int(nr, degree)] <- TRUE
  }

  conservative <- unclass(findmin(
    chart, type = "conservative", gurobi = FALSE, solind = TRUE
  ))
  .Call("C_resetScpProfile", PACKAGE = "QCA")
  hybrid <- unclass(findmin(chart, type = "hybrid", solind = TRUE))
  profile <- .Call("C_getScpProfile", PACKAGE = "QCA")

  expect_equal(sum(hybrid), sum(conservative))
  expect_true(profile$scp_attempted)
  expect_false(profile$scp_limited)
  expect_false(profile$lpsolve_fallback)
  expect_gte(profile$adaptive_extensions, 1)
  expect_equal(profile$root_branches_completed, profile$root_branches_total)
})

test_that("hybrid stops a stagnant exact-finisher probe", {
  set.seed(4)
  nr <- 48L
  nc <- 600L
  chart <- matrix(FALSE, nrow = nc, ncol = nr)
  for (column in seq_len(nc)) {
    degree <- sample(3L:8L, 1L)
    chart[column, sample.int(nr, degree)] <- TRUE
  }

  conservative <- unclass(findmin(
    chart, type = "conservative", gurobi = FALSE, solind = TRUE
  ))
  .Call("C_resetScpProfile", PACKAGE = "QCA")
  hybrid <- unclass(findmin(chart, type = "hybrid", solind = TRUE))
  profile <- .Call("C_getScpProfile", PACKAGE = "QCA")

  expect_equal(sum(hybrid), sum(conservative))
  expect_true(profile$scp_attempted)
  expect_true(profile$scp_limited)
  expect_true(profile$lpsolve_fallback)
  expect_equal(profile$adaptive_extensions, 0)
  expect_equal(profile$root_branches_completed, 0)
  expect_gt(profile$root_branches_total, 0)
})
