library(testthat)

test_that("clustering modules define key functions", {
  source(file.path(dirname(getwd()), "src", "R_shared", "io_utils.R"))
  source(file.path(dirname(getwd()), "src", "R_shared", "stats_utils.R"))
  source(file.path(dirname(getwd()), "src", "clustering", "leiden.R"))
  expect_true(is.function(run_leiden))
})
