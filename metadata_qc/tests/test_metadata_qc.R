library(testthat)

test_that("metadata_qc scripts define key functions", {
  source(file.path(dirname(getwd()), "R_shared", "io_utils.R"))
  source(file.path(dirname(getwd()), "metadata_qc", "check_acquisition.R"))
  expect_true(is.function(get_fcs_keyword))
  expect_true(is.function(check_acquisition_dates))
})
