library(testthat)

# ---- cytometry src bootstrap (Nextflow-safe); see src/R_shared/paths.R ----
if (!exists(".cyto_paths_loaded", inherits = TRUE)) {
  local({
    root <- Sys.getenv("CYTOMETRY_SRC", "")
    if (!nzchar(root) || !dir.exists(file.path(root, "R_shared"))) {
      of <- NULL
      frames <- sys.frames()
      for (i in rev(seq_along(frames))) {
        f <- frames[[i]]$ofile
        if (!is.null(f) && nzchar(f)) { of <- f; break }
      }
      d <- if (!is.null(of)) dirname(normalizePath(of)) else getwd()
      while (!dir.exists(file.path(d, "R_shared")) && d != dirname(d)) d <- dirname(d)
      root <- d
    }
    source(file.path(root, "R_shared", "paths.R"), local = FALSE)
  })
}

# ---------------------------------------------------------------------------
# Helpers: resolve shared and visualization source directories
# ---------------------------------------------------------------------------
root_dir   <- cyto_repo_root()
shared_dir <- cyto_shared_dir()
viz_dir    <- file.path(cyto_src_root(), "visualization")

# Source shared utilities (required by all visualization scripts)
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "plot_utils.R"))

# =============================================================================
# Tests for create_safe_filename()
# =============================================================================
test_that("create_safe_filename removes special characters", {
  expect_equal(create_safe_filename("hello world!"), "hello_world")
  expect_equal(create_safe_filename("file@#$%^&*name"), "file_name")
  expect_equal(create_safe_filename("a/b\\c:d"), "a_b_c_d")
})

test_that("create_safe_filename collapses multiple underscores", {
  expect_equal(create_safe_filename("a---b___c"), "a-b_c")
  expect_equal(create_safe_filename("lots   of   spaces"), "lots_of_spaces")
})

test_that("create_safe_filename strips leading/trailing underscores", {
  expect_equal(create_safe_filename("__hello__"), "hello")
  expect_equal(create_safe_filename("  spaced  "), "spaced")
})

test_that("create_safe_filename respects max_length", {
  long_name <- paste0(rep("A", 200), collapse = "")
  result <- create_safe_filename(long_name, max_length = 50)
  expect_lte(nchar(result), 50)
  # Default max_length = 100
  result_default <- create_safe_filename(long_name)
  expect_lte(nchar(result_default), 100)
})

test_that("create_safe_filename preserves dots, hyphens, and underscores", {
  expect_equal(create_safe_filename("file.name-v2_final"), "file.name-v2_final")
})

test_that("create_safe_filename handles empty or whitespace input", {
  expect_equal(create_safe_filename(""), "")
  expect_equal(create_safe_filename("   "), "")
})

test_that("create_safe_filename preserves alphanumeric content", {
  expect_equal(create_safe_filename("Sample123_CD4pos"), "Sample123_CD4pos")
})

# =============================================================================
# Tests for create_compact_names()
# =============================================================================
test_that("create_compact_names returns short names unchanged", {
  input <- c("CD4", "CD8", "CD19")
  result <- create_compact_names(input, max_length = 25)
  expect_equal(unname(result), input)
})

test_that("create_compact_names truncates long names with ellipsis", {
  long_name <- "This_is_a_very_long_marker_name_indeed"
  result <- create_compact_names(long_name, max_length = 15)
  expect_equal(nchar(unname(result)), 15)
  expect_true(grepl("\\.\\.\\.$", unname(result)))
})

test_that("create_compact_names handles boundary length exactly", {
  name_25 <- paste0(rep("x", 25), collapse = "")
  result <- create_compact_names(name_25, max_length = 25)
  expect_equal(unname(result), name_25)  # exactly max_length, no truncation
})

test_that("create_compact_names works with vector input", {
  input <- c("short", "a_very_long_name_exceeding_ten_chars", "ok")
  result <- create_compact_names(input, max_length = 10)
  expect_equal(length(result), 3)
  expect_equal(unname(result[1]), "short")
  expect_equal(nchar(unname(result[2])), 10)
  expect_true(endsWith(unname(result[2]), "..."))
  expect_equal(unname(result[3]), "ok")
})

test_that("create_compact_names handles single character max_length gracefully", {
  # Edge case: max_length = 5 (smaller truncation)
  result <- create_compact_names("ABCDEFGHIJ", max_length = 5)
  expect_equal(nchar(unname(result)), 5)
})

# =============================================================================
# Tests for save_plot_formats()
# =============================================================================
test_that("save_plot_formats requires ggplot2", {
  skip_if_not_installed("ggplot2")
  tmp_dir <- tempdir()
  p <- ggplot2::ggplot(data.frame(x = 1:5, y = 1:5), ggplot2::aes(x, y)) +
    ggplot2::geom_point()

  base_path <- file.path(tmp_dir, "test_save_plot")
  result <- save_plot_formats(p, base_path,
                              formats = list(png = TRUE, pdf = FALSE, svg = FALSE))

  expect_true("PNG" %in% result)
  expect_false("PDF" %in% result)
  expect_false("SVG" %in% result)
  expect_true(file.exists(paste0(base_path, ".png")))

  # Cleanup
  unlink(paste0(base_path, ".png"))
})

test_that("save_plot_formats saves multiple formats", {
  skip_if_not_installed("ggplot2")
  tmp_dir <- tempdir()
  p <- ggplot2::ggplot(data.frame(x = 1:5, y = 1:5), ggplot2::aes(x, y)) +
    ggplot2::geom_point()

  base_path <- file.path(tmp_dir, "test_multi_format")
  result <- save_plot_formats(p, base_path,
                              formats = list(png = TRUE, pdf = TRUE, svg = FALSE))

  expect_equal(sort(result), c("PDF", "PNG"))
  expect_true(file.exists(paste0(base_path, ".png")))
  expect_true(file.exists(paste0(base_path, ".pdf")))

  # Cleanup
  unlink(paste0(base_path, ".png"))
  unlink(paste0(base_path, ".pdf"))
})

test_that("save_plot_formats returns empty vector when no formats enabled", {
  skip_if_not_installed("ggplot2")
  p <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
    ggplot2::geom_point()
  result <- save_plot_formats(p, file.path(tempdir(), "no_format"),
                              formats = list(png = FALSE, pdf = FALSE, svg = FALSE))
  expect_length(result, 0)
})

# =============================================================================
# Tests for plot_composition_barplot()
# =============================================================================
test_that("plot_composition_barplot validates missing x_col", {
  skip_if_not_installed("ggplot2")
  source(file.path(viz_dir, "barplots.R"))

  df <- data.frame(sample = c("S1", "S2"), count = c(10, 20))
  expect_error(
    plot_composition_barplot(df, tempdir(), x_col = "cluster"),
    "x_col 'cluster' not found"
  )
})

test_that("plot_composition_barplot validates missing fill_col", {
  skip_if_not_installed("ggplot2")
  source(file.path(viz_dir, "barplots.R"))

  df <- data.frame(cluster = c("A", "B"), count = c(10, 20))
  expect_error(
    plot_composition_barplot(df, tempdir(), fill_col = "sample"),
    "fill_col 'sample' not found"
  )
})

test_that("plot_composition_barplot handles empty data with warning", {
  skip_if_not_installed("ggplot2")
  source(file.path(viz_dir, "barplots.R"))

  df <- data.frame(cluster = character(0), sample = character(0))
  expect_warning(
    result <- plot_composition_barplot(df, tempdir()),
    "Empty data frame"
  )
  expect_null(result)
})

test_that("plot_composition_barplot produces valid ggplot with correct data", {
  skip_if_not_installed("ggplot2")
  source(file.path(viz_dir, "barplots.R"))

  out_dir <- file.path(tempdir(), "barplot_test")
  df <- data.frame(
    cluster = rep(c("A", "B", "C"), each = 10),
    sample  = rep(c("S1", "S2"), 15)
  )
  result <- plot_composition_barplot(df, out_dir)
  expect_true(inherits(result, "gg"))

  # Cleanup
  unlink(out_dir, recursive = TRUE)
})

# =============================================================================
# Tests for plot_contour()
# =============================================================================
test_that("plot_contour validates missing x and y columns", {
  skip_if_not_installed("ggplot2")
  source(file.path(viz_dir, "contour_plots.R"))

  # Missing 'x' and 'y' columns entirely
  coords_bad <- data.frame(a = rnorm(50), b = rnorm(50))
  expect_error(
    plot_contour(coords_bad, output_dir = tempdir()),
    "coords must have 'x' and 'y' columns"
  )
})

test_that("plot_contour validates custom x_col not found", {
  skip_if_not_installed("ggplot2")
  source(file.path(viz_dir, "contour_plots.R"))

  coords <- data.frame(x = rnorm(50), y = rnorm(50))
  expect_error(
    plot_contour(coords, x_col = "UMAP1", output_dir = tempdir()),
    "x_col 'UMAP1' not found"
  )
})

test_that("plot_contour validates custom y_col not found", {
  skip_if_not_installed("ggplot2")
  source(file.path(viz_dir, "contour_plots.R"))

  coords <- data.frame(x = rnorm(50), y = rnorm(50))
  expect_error(
    plot_contour(coords, y_col = "UMAP2", output_dir = tempdir()),
    "y_col 'UMAP2' not found"
  )
})

test_that("plot_contour rejects unknown method", {
  skip_if_not_installed("ggplot2")
  source(file.path(viz_dir, "contour_plots.R"))

  coords <- data.frame(x = rnorm(50), y = rnorm(50))
  expect_error(
    plot_contour(coords, method = "invalid", output_dir = tempdir()),
    "Unknown method"
  )
})

test_that("plot_contour succeeds with valid data and geom_density2d method", {
  skip_if_not_installed("ggplot2")
  source(file.path(viz_dir, "contour_plots.R"))

  out_dir <- file.path(tempdir(), "contour_test")
  coords <- data.frame(x = rnorm(100), y = rnorm(100))
  result <- plot_contour(coords, method = "geom_density2d", output_dir = out_dir,
                         output_formats = list(png = TRUE, pdf = FALSE, svg = FALSE))
  expect_true(inherits(result, "gg"))

  # Cleanup
  unlink(out_dir, recursive = TRUE)
})

test_that("plot_contour succeeds with group_labels", {
  skip_if_not_installed("ggplot2")
  source(file.path(viz_dir, "contour_plots.R"))

  out_dir <- file.path(tempdir(), "contour_grouped_test")
  coords <- data.frame(x = rnorm(100), y = rnorm(100))
  groups <- rep(c("A", "B"), each = 50)
  result <- plot_contour(coords, group_labels = groups,
                         method = "geom_density2d", output_dir = out_dir,
                         output_formats = list(png = TRUE, pdf = FALSE, svg = FALSE))
  expect_true(inherits(result, "gg"))

  # Cleanup
  unlink(out_dir, recursive = TRUE)
})

# =============================================================================
# Tests for plot_cluster_heatmap()
# =============================================================================
test_that("plot_cluster_heatmap validates missing cluster column", {
  skip_if_not_installed("pheatmap")
  source(file.path(viz_dir, "heatmaps.R"))

  df <- data.frame(CD4 = rnorm(20), CD8 = rnorm(20))
  expect_error(
    plot_cluster_heatmap(df, cluster_col = "cluster",
                         marker_cols = c("CD4", "CD8"),
                         output_path = file.path(tempdir(), "hm_test.png")),
    "cluster_col 'cluster' not found"
  )
})

test_that("plot_cluster_heatmap validates missing marker columns", {
  skip_if_not_installed("pheatmap")
  source(file.path(viz_dir, "heatmaps.R"))

  df <- data.frame(cluster = rep(1:3, each = 10), CD4 = rnorm(30))
  expect_error(
    plot_cluster_heatmap(df, cluster_col = "cluster",
                         marker_cols = c("CD4", "CD8", "CD19"),
                         output_path = file.path(tempdir(), "hm_test.png")),
    "marker_cols not found"
  )
})

test_that("plot_cluster_heatmap succeeds with valid input", {
  skip_if_not_installed("pheatmap")
  source(file.path(viz_dir, "heatmaps.R"))

  df <- data.frame(
    cluster = rep(1:3, each = 20),
    CD4  = rnorm(60, mean = 5),
    CD8  = rnorm(60, mean = 3),
    CD19 = rnorm(60, mean = 4)
  )
  out_path <- file.path(tempdir(), "cluster_hm_test.png")
  expect_no_error(
    plot_cluster_heatmap(df, cluster_col = "cluster",
                         marker_cols = c("CD4", "CD8", "CD19"),
                         output_path = out_path)
  )

  # Cleanup
  base <- sub("\\.[^.]+$", "", out_path)
  unlink(paste0(base, ".png"))
  unlink(paste0(base, ".svg"))
})

# =============================================================================
# Tests for plot_upset()
# =============================================================================
test_that("plot_upset validates missing set_columns", {
  skip_if_not_installed("UpSetR")
  source(file.path(viz_dir, "upset_plots.R"))

  df <- data.frame(A = c(1, 0, 1), B = c(0, 1, 1))
  expect_error(
    plot_upset(df, set_columns = c("A", "B", "C"),
               output_path = file.path(tempdir(), "upset_test.png")),
    "set_columns not found"
  )
})

test_that("plot_upset handles empty data with warning", {
  skip_if_not_installed("UpSetR")
  source(file.path(viz_dir, "upset_plots.R"))

  df <- data.frame(A = integer(0), B = integer(0))
  expect_warning(
    result <- plot_upset(df, set_columns = c("A", "B"),
                         output_path = file.path(tempdir(), "upset_empty.png")),
    "Empty data"
  )
  expect_null(result)
})

test_that("plot_upset succeeds with valid binary matrix", {
  skip_if_not_installed("UpSetR")
  source(file.path(viz_dir, "upset_plots.R"))

  set.seed(42)
  df <- data.frame(
    SetA = sample(0:1, 50, replace = TRUE),
    SetB = sample(0:1, 50, replace = TRUE),
    SetC = sample(0:1, 50, replace = TRUE)
  )
  out_path <- file.path(tempdir(), "upset_valid.png")
  expect_no_error(
    plot_upset(df, set_columns = c("SetA", "SetB", "SetC"),
               output_path = out_path)
  )

  # Check files were created
  base <- sub("\\.[^.]+$", "", out_path)
  expect_true(file.exists(paste0(base, ".png")))

  # Cleanup
  unlink(paste0(base, ".png"))
  unlink(paste0(base, ".svg"))
})

# =============================================================================
# Tests for plot_heatmap()
# =============================================================================
test_that("plot_heatmap handles Inf and NaN values gracefully", {
  skip_if_not_installed("pheatmap")
  source(file.path(viz_dir, "heatmaps.R"))

  mtx <- matrix(rnorm(30), nrow = 5, ncol = 6)
  colnames(mtx) <- paste0("M", 1:6)
  rownames(mtx) <- paste0("C", 1:5)
  mtx[1, 1] <- Inf
  mtx[2, 3] <- NaN

  out_path <- file.path(tempdir(), "heatmap_inf_test.png")
  # Should not error; Inf/NaN are replaced with NA internally
  expect_no_error(
    plot_heatmap(mtx, out_path, scale = "none")
  )

  base <- sub("\\.[^.]+$", "", out_path)
  unlink(paste0(base, ".png"))
  unlink(paste0(base, ".svg"))
})

# =============================================================================
# Tests for save_base_plot_formats()
# =============================================================================
test_that("save_base_plot_formats creates PNG and SVG files", {
  out_base <- file.path(tempdir(), "base_plot_test")
  draw_fn <- function() {
    plot(1:10, 1:10, main = "Test")
  }
  expect_no_error(save_base_plot_formats(out_base, draw_fn, width = 6, height = 4))
  expect_true(file.exists(paste0(out_base, ".png")))
  expect_true(file.exists(paste0(out_base, ".svg")))

  # Cleanup
  unlink(paste0(out_base, ".png"))
  unlink(paste0(out_base, ".svg"))
})

# =============================================================================
# Tests for population_barplot()
# =============================================================================
test_that("population_barplot creates output files", {
  skip_if_not_installed("ggplot2")

  pop_df <- data.frame(
    Population = rep(c("CD4+", "CD8+"), each = 3),
    Count = c(100, 150, 120, 80, 90, 110),
    Sample = rep(c("S1", "S2", "S3"), 2)
  )
  out_file <- file.path(tempdir(), "pop_barplot_test.pdf")
  expect_no_error(population_barplot(pop_df, out_file, "Test Title", "Cell Count"))

  base <- sub("\\.[^.]+$", "", out_file)
  # save_plot_formats creates its own extensions
  expect_true(file.exists(paste0(base, ".pdf")) ||
              file.exists(paste0(base, ".png")) ||
              file.exists(paste0(base, ".svg")))

  # Cleanup
  unlink(paste0(base, ".pdf"))
  unlink(paste0(base, ".png"))
  unlink(paste0(base, ".svg"))
})
