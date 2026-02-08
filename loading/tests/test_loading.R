library(testthat)

test_that("loading scripts source without error", {
  # Temporarily set shared_dir
  expect_no_error(source(file.path(dirname(getwd()), "loading", "load_fcs_raw.R")))
})

test_that("key functions exist after sourcing", {
  source(file.path(dirname(getwd()), "R_shared", "io_utils.R"))
  source(file.path(dirname(getwd()), "R_shared", "fcs_utils.R"))
  expect_true(is.function(load_fcs_list))
  expect_true(is.function(load_fcs_folder))
})
