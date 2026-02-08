library(testthat)

test_that("harmonization scripts define expected functions", {
  source(file.path(dirname(getwd()), "src", "R_shared", "io_utils.R"))
  source(file.path(dirname(getwd()), "src", "harmonization", "panel_homogenization.R"))
  expect_true(is.function(translate_vector))
  expect_true(is.function(homogenize_flowframes))
})
