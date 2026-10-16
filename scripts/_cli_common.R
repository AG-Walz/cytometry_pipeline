# =============================================================================
# scripts/_cli_common.R : helpers shared by the CLI wrappers in scripts/
# =============================================================================
# Sourced by the wrappers (sync_fcs_names.R, cytokine_combinations.R,
# umap_from_flowjo.R) from the directory Rscript's --file= argument names.
# Nothing here prints: umap_from_flowjo.R's stdout is parsed by the Dash app.
# =============================================================================

# Path of the running script. Works under Rscript (--file=) and source(); under
# source() the srcfile fallback resolves to this file, which sits in the same
# scripts/ directory, so resolve_src_dir() still lands on repo/src.
find_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1) {
    return(normalizePath(sub("^--file=", "", file_arg), mustWork = FALSE))
  }
  if (!is.null(sys.function(1)) && !is.null(attr(body(sys.function(1)), "srcfile"))) {
    sf <- attr(body(sys.function(1)), "srcfile")
    if (!is.null(sf$filename)) return(normalizePath(sf$filename, mustWork = FALSE))
  }
  NA_character_
}

# The repo 'src' directory: --src_dir if given, else repo/scripts/x.R -> repo/src,
# else src/ under the working directory.
resolve_src_dir <- function(user_src) {
  if (!is.null(user_src)) {
    sd <- normalizePath(user_src, mustWork = FALSE)
    if (!dir.exists(sd)) stop("--src_dir does not exist: ", sd)
    return(sd)
  }
  sp <- find_script_path()
  if (!is.na(sp)) {
    cand <- file.path(dirname(dirname(sp)), "src")
    if (dir.exists(cand)) return(normalizePath(cand))
  }
  cand <- file.path(getwd(), "src")
  if (dir.exists(cand)) return(normalizePath(cand))
  stop("Could not locate the 'src' directory. Pass --src_dir explicitly.")
}

# Split a comma-separated CLI string, trimming and dropping empty items.
# `empty` is what a NULL or all-blank input becomes: NULL for the wrappers
# whose module treats NULL as "not given" (umap, cytokine combinations),
# character(0) for sync_fcs_names(), whose defaults are character().
split_csv <- function(s, empty = NULL) {
  if (is.null(s)) return(empty)
  v <- trimws(strsplit(s, ",")[[1]])
  v <- v[nzchar(v)]
  if (length(v) == 0) empty else v
}

# "png,pdf" -> list(pdf = TRUE, png = TRUE, svg = FALSE), the form
# save_plot_formats() takes.
parse_formats <- function(s) {
  parts <- split_csv(s)
  list(pdf = "pdf" %in% parts, png = "png" %in% parts, svg = "svg" %in% parts)
}
