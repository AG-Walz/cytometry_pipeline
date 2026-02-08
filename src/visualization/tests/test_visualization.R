library(testthat)

test_that("plot_utils sources without error", {
  expect_no_error(source(file.path(dirname(getwd()), "src", "R_shared", "io_utils.R")))
  expect_no_error(source(file.path(dirname(getwd()), "src", "R_shared", "plot_utils.R")))
})

test_that("plot utility functions exist", {
  source(file.path(dirname(getwd()), "src", "R_shared", "plot_utils.R"))
  expect_true(is.function(save_plot_formats))
  expect_true(is.function(create_safe_filename))
})
