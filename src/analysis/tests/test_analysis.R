library(testthat)

test_that("R_shared cluster_utils sources without error", {
  expect_no_error(source(file.path(dirname(getwd()), "src", "R_shared", "io_utils.R")))
  expect_no_error(source(file.path(dirname(getwd()), "src", "R_shared", "cluster_utils.R")))
})

test_that("cluster utility functions exist", {
  source(file.path(dirname(getwd()), "src", "R_shared", "cluster_utils.R"))
  expect_true(is.function(clean_cluster_name))
  expect_true(is.function(extract_cluster_group))
  expect_true(is.function(convert_binary_to_single_cluster))
})
