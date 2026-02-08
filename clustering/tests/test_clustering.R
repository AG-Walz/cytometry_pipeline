library(testthat)

test_that("clustering modules define key functions", {
  source(file.path(dirname(getwd()), "R_shared", "io_utils.R"))
  source(file.path(dirname(getwd()), "R_shared", "stats_utils.R"))
  source(file.path(dirname(getwd()), "clustering", "leiden.R"))
  expect_true(is.function(run_leiden))
})
