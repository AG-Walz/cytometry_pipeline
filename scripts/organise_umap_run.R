#!/usr/bin/env Rscript
# =============================================================================
# scripts/organise_umap_run.R: move a flat umap_from_flowjo.R run folder into
# the layout runs have written since 2026-10-15 (run_subdir() in
# src/visualization/umap_from_flowjo.R):
#
#   overview/  markers/  markers/per_sample/<marker>/  metadata/<column>/
#   splits/<column>/  changes/<column>/  data/  model/
#
# Files are moved, not copied, and keep their names; plot_index.json is
# rewritten with the new paths. A file it cannot place stays where it is and
# is listed. Run with --dry_run first to see the plan.
#
# Usage:
#   Rscript scripts/organise_umap_run.R --run DIR [--dry_run]
# =============================================================================

suppressPackageStartupMessages(library(optparse))
opt <- parse_args(OptionParser(option_list = list(
  make_option("--run", type = "character", help = "Run folder to organise"),
  make_option("--dry_run", action = "store_true", default = FALSE,
              help = "Only print what would move")
)))
if (is.null(opt$run) || !dir.exists(opt$run)) {
  stop("--run must be an existing run folder.", call. = FALSE)
}

script_dir <- tryCatch(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), error = function(e) ".")
src_dir <- normalizePath(file.path(script_dir, "..", "src"), mustWork = TRUE)
source(file.path(src_dir, "R_shared", "plot_utils.R"))
source(file.path(src_dir, "visualization", "umap_from_flowjo.R"))

run <- normalizePath(opt$run)
files <- list.files(run, all.files = FALSE, no.. = TRUE)
files <- files[file.info(file.path(run, files))$isdir %in% FALSE]

# Plot kinds and columns from the index, where the name alone is ambiguous
# (a split column with an underscore cannot be told from its level).
index <- NULL
idx_file <- file.path(run, "plot_index.json")
if (file.exists(idx_file)) {
  index <- tryCatch(jsonlite::fromJSON(idx_file), error = function(e) NULL)
  if (is.data.frame(index)) index$stem <- basename(as.character(index$file))
}
split_cols <- sub("^combination_frequencies_(.*)\\.csv$", "\\1",
                  grep("^combination_frequencies_.*\\.csv$", files, value = TRUE))

DATA_FILES <- c("cells.csv.gz", "supercells.csv.gz", "combination_counts.csv",
                "gate_decisions.csv", "separation_metrics.csv",
                "background_pairs.csv", "sample_metadata.csv")
TOP_FILES <- c("run_info.json", "plot_index.json", "README.txt")

# Longest known split column that prefixes `tail` ("<col>_<level>").
column_of <- function(tail) {
  hit <- split_cols[startsWith(tail, paste0(split_cols, "_")) |
                      tail == paste0(split_cols, "_panel")]
  if (length(hit)) hit[which.max(nchar(hit))] else sub("_.*$", "", tail)
}

target_of <- function(f) {
  if (f %in% TOP_FILES) return(NA_character_)
  if (f %in% DATA_FILES) return(run_subdir("data"))
  if (f == "umap_model.uwot") return(run_subdir("model"))
  m <- regmatches(f, regexec("^combination_(frequencies|change)_(.*)\\.csv$", f))[[1]]
  if (length(m)) return(run_subdir("table", column = m[3]))
  if (!grepl("\\.(png|pdf|svg)$", f)) return(NA_character_)
  stem <- sub("\\.(png|pdf|svg)$", "", f)
  if (!is.null(index) && stem %in% index$stem) {
    e <- index[match(stem, index$stem), ]
    if (e$kind %in% c("split", "split_panel", "diff") && !is.na(e$column)) {
      return(run_subdir(e$kind, column = e$column))
    }
    if (e$kind == "combination") return(run_subdir("overview"))
  }
  if (grepl("_(by_combination|by_file|by_workspace|by_subpopulation|by_level[0-9]+)$", stem)) {
    return(run_subdir("overview"))
  }
  hits <- c(split = regexpr("_split_", stem, fixed = TRUE)[1],
            diff = regexpr("_diff_", stem, fixed = TRUE)[1],
            marker = regexpr("_marker_", stem, fixed = TRUE)[1],
            annot = regexpr("_annot_", stem, fixed = TRUE)[1])
  hits <- hits[hits > 0]
  if (!length(hits)) return(NA_character_)
  kind <- names(hits)[which.min(hits)]
  tail <- substring(stem, hits[[kind]] + nchar(kind) + 2L)
  switch(kind,
    split = run_subdir("split", column = column_of(tail)),
    diff = run_subdir("diff", column = column_of(tail)),
    marker = if (grepl("_file_", tail, fixed = TRUE)) {
      run_subdir("marker_file", marker = sub("_file_.*$", "", tail))
    } else run_subdir("marker"),
    annot = run_subdir("annot", column = if (grepl("_level_", tail, fixed = TRUE))
      sub("_level_.*$", "", tail) else sub("_ALL$", "", tail)))
}

target <- vapply(files, target_of, character(1))
moving <- files[!is.na(target) & nzchar(target)]
stay <- setdiff(files[is.na(target)], TOP_FILES)
plan <- table(target[moving])
cat(sprintf("%s: %d file(s) to move into %d folder(s)%s\n", run, length(moving),
            length(plan), if (opt$dry_run) " (dry run)" else ""))
for (d in names(plan)) cat(sprintf("  %-50s %6d\n", d, plan[[d]]))
if (length(stay)) {
  cat("Left where they are (not recognised):", paste(stay, collapse = ", "), "\n")
}
if (opt$dry_run) quit(status = 0)

for (d in unique(target[moving])) dir.create(file.path(run, d), recursive = TRUE,
                                             showWarnings = FALSE)
ok <- file.rename(file.path(run, moving), file.path(run, target[moving], moving))
if (!all(ok)) {
  warning(sum(!ok), " file(s) could not be moved: ",
          paste(utils::head(moving[!ok], 5), collapse = ", "), call. = FALSE)
}

# plot_index.json with the paths the files now have.
if (is.data.frame(index)) {
  stem_ext <- sub("\\.(png|pdf|svg)$", "", moving)
  new_dir <- stats::setNames(target[moving], stem_ext)
  index$file <- ifelse(index$stem %in% names(new_dir),
                       file.path(new_dir[index$stem], index$stem),
                       as.character(index$file))
  index$stem <- NULL
  jsonlite::write_json(index, idx_file, dataframe = "rows", pretty = TRUE,
                       na = "null")
}
cat(sprintf("Moved %d file(s).\n", sum(ok)))
