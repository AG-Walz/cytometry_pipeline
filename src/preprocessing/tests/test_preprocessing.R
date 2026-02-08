library(testthat)

test_that("R_shared stats_utils sources without error", {
  expect_no_error(source(file.path(dirname(getwd()), "src", "R_shared", "io_utils.R")))
  expect_no_error(source(file.path(dirname(getwd()), "src", "R_shared", "stats_utils.R")))
})

test_that("key functions exist", {
  source(file.path(dirname(getwd()), "src", "R_shared", "stats_utils.R"))
  expect_true(is.function(run_pca))
  expect_true(is.function(run_umap))
})
