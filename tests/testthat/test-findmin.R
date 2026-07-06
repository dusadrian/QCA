local_edition(3)
data(d.represent, package = "QCA")
data(d.jobsecurity, package = "QCA")
data(d.partybans, package = "QCA")
data(CVF, package = "QCA")
chart <- makeChart("a, b, ~c", "abc, a~b~c, a~bc, ~ab~c")
fmin <- unclass(findmin(chart))
fmin.exact <- unclass(findmin(chart, type = "exact"))
fmin.lagrangian <- unclass(findmin(chart, type = "lagrangian"))
fmin.lagrangian.solind <- unclass(findmin(chart, type = "lagrangian", solind = TRUE))
attributes(fmin) <- NULL
attributes(fmin.exact) <- NULL
attributes(fmin.lagrangian) <- NULL
attributes(fmin.lagrangian.solind) <- NULL

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

lagrangian_cardinality <- function(chart) {
  sum(unclass(findmin(chart, type = "lagrangian", solind = TRUE)) > 0)
}

test_that("findmin() works", {
  expect_equal(fmin, 2)
  expect_equal(fmin.exact, 2)
  expect_equal(fmin.lagrangian, 2)
  expect_true(is.numeric(fmin.lagrangian.solind))
  expect_equal(sum(fmin.lagrangian.solind > 0), 2)
})

test_that("lagrangian matches exact on representative PI charts", {
  for (nm in names(representative.charts)) {
    chart <- representative.charts[[nm]]
    exact_min <- unclass(findmin(chart, type = "exact"))
    lagrangian_min <- unclass(findmin(chart, type = "lagrangian"))
    lagrangian_solind_min <- lagrangian_cardinality(chart)

    expect_equal(
      lagrangian_min,
      exact_min,
      info = sprintf("%s exact vs lagrangian minimum", nm)
    )
    expect_equal(
      lagrangian_solind_min,
      exact_min,
      info = sprintf("%s lagrangian solind minimum", nm)
    )
  }
})
