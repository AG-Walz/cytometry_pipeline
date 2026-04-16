#!/usr/bin/env Rscript
# =============================================================================
# preprocessing/gate_utils.R — Gate Application & Boolean Gating Utilities
# =============================================================================
# Utilities for applying sequences of flowCore gates to flowFrames and
# combining gate results with boolean logic (AND/OR/NOT expressions).
# =============================================================================

if (!exists("shared_dir")) shared_dir <- { .s <- Sys.getenv("CYTOMETRY_SRC", ""); if (nzchar(.s)) file.path(.s, "R_shared") else file.path(dirname(getwd()), "src", "R_shared") }
source(file.path(shared_dir, "io_utils.R"))

#' Apply a sequence of gates to a flowFrame
#'
#' Input from: define_gate() (preprocessing/automated_gating.R),
#'             run_auto_gating() (preprocessing/automated_gating.R)
#' Output to:  quality_control() (preprocessing/quality_control.R),
#'             boolean_gate() (this file),
#'             downstream analysis
#'
#' @param flowframe flowFrame object
#' @param gate_list List of flowCore gate objects (rectangleGate, polygonGate, etc.)
#' @param logical_op "and" (intersection) or "or" (union)
#' @return list(gated_ff, gate_indices, population_stats)
apply_gate_sequence <- function(flowframe, gate_list, logical_op = "and") {
  require(flowCore)

  logical_op <- tolower(logical_op)
  if (!logical_op %in% c("and", "or")) {
    stop("logical_op must be 'and' or 'or', got '", logical_op, "'")
  }
  if (length(gate_list) == 0) {
    stop("gate_list must contain at least one gate.")
  }

  n_cells <- nrow(flowframe)
  cat("Applying", length(gate_list), "gate(s) with '", logical_op,
      "' logic to", n_cells, "events\n")

  # --- Apply each gate and collect logical vectors ---------------------------
  gate_results <- vector("list", length(gate_list))

  for (i in seq_along(gate_list)) {
    gate <- gate_list[[i]]

    filter_result <- tryCatch({
      flowCore::filter(flowframe, gate)
    }, error = function(e) {
      stop("Gate ", i, " failed: ", e$message)
    })

    gate_results[[i]] <- as.logical(filter_result@subSet)

    n_in <- sum(gate_results[[i]])
    cat("  Gate", i, ":", n_in, "/", n_cells,
        "(", round(100 * n_in / n_cells, 1), "%)\n")
  }

  # --- Combine gate results --------------------------------------------------
  if (logical_op == "and") {
    combined <- Reduce(`&`, gate_results)
  } else {
    combined <- Reduce(`|`, gate_results)
  }

  gate_indices <- which(combined)
  n_gated <- length(gate_indices)

  cat("Combined (", logical_op, "):", n_gated, "/", n_cells,
      "(", round(100 * n_gated / n_cells, 1), "%)\n")

  # --- Build gated flowFrame -------------------------------------------------
  gated_ff <- flowframe[gate_indices, ]

  list(
    gated_ff = gated_ff,
    gate_indices = gate_indices,
    population_stats = list(
      n_parent   = n_cells,
      n_gated    = n_gated,
      pct_parent = 100 * n_gated / n_cells,
      n_gates    = length(gate_list),
      logical_op = logical_op
    )
  )
}

#' Boolean gating (combine existing gates with AND/OR/NOT)
#'
#' Evaluates a boolean expression string referencing named gate results.
#' Gate names in the expression are replaced with their logical vectors
#' and the expression is evaluated.
#'
#' Input from: apply_gate_sequence() (this file),
#'             run_auto_gating() (preprocessing/automated_gating.R)
#' Output to:  downstream analysis, visualization
#'
#' @param flowframe flowFrame object
#' @param gate_results Named list of logical vectors from previous gates.
#'   Names should match the gate identifiers used in the expression
#'   (e.g., list("CD3+" = logical_vec, "CD4+" = logical_vec, "CD8+" = logical_vec)).
#' @param expression String expression, e.g., "CD3+ & CD4+ & !CD8+"
#' @return Logical vector of included cells
boolean_gate <- function(flowframe, gate_results, expression) {
  require(flowCore)

  n_cells <- nrow(flowframe)

  # --- Validate inputs -------------------------------------------------------
  if (!is.list(gate_results) || is.null(names(gate_results))) {
    stop("gate_results must be a named list of logical vectors.")
  }

  # Validate all logical vectors have correct length
  for (nm in names(gate_results)) {
    if (!is.logical(gate_results[[nm]])) {
      stop("gate_results[['", nm, "']] must be a logical vector.")
    }
    if (length(gate_results[[nm]]) != n_cells) {
      stop("gate_results[['", nm, "']] has length ", length(gate_results[[nm]]),
           " but flowframe has ", n_cells, " cells.")
    }
  }

  cat("Boolean gating: '", expression, "' on ", n_cells, " events\n", sep = "")

  # --- Parse and evaluate expression -----------------------------------------
  # Gate names may contain characters that are meaningful in regex (e.g. `+`
  # in "CD3+"), so use fixed (literal) replacement. Sort by length so that
  # longer names are replaced before their shorter prefixes (e.g. "CD3+CD4+"
  # before "CD3+").
  gate_names <- names(gate_results)
  gate_names_sorted <- gate_names[order(nchar(gate_names), decreasing = TRUE)]

  # Use a placeholder per gate to avoid overlap between replaced substrings
  # and new gate-name tokens introduced by earlier substitutions.
  placeholders <- paste0("__GATE_", seq_along(gate_names_sorted), "__")
  eval_expr <- expression
  for (i in seq_along(gate_names_sorted)) {
    eval_expr <- gsub(gate_names_sorted[i], placeholders[i], eval_expr, fixed = TRUE)
  }
  for (i in seq_along(gate_names_sorted)) {
    eval_expr <- gsub(placeholders[i],
                      paste0("gate_results[['", gate_names_sorted[i], "']]"),
                      eval_expr, fixed = TRUE)
  }

  # Validate only expected tokens remain (logical operators and gate references)
  check_expr <- eval_expr
  check_expr <- gsub("gate_results\\[\\['.+?'\\]\\]", "", check_expr)
  check_expr <- gsub("[&|!() ]", "", check_expr)
  if (nchar(trimws(check_expr)) > 0) {
    warning("Unrecognized tokens in boolean expression: '", trimws(check_expr),
            "'. Check that all gate names in the expression match names in gate_results.")
  }

  result <- tryCatch({
    eval(parse(text = eval_expr))
  }, error = function(e) {
    stop("Failed to evaluate boolean expression '", expression, "': ", e$message,
         "\nParsed as: ", eval_expr,
         "\nAvailable gates: ", paste(gate_names, collapse = ", "))
  })

  if (!is.logical(result) || length(result) != n_cells) {
    stop("Boolean expression did not produce a logical vector of length ", n_cells)
  }

  n_selected <- sum(result)
  cat("  Result:", n_selected, "/", n_cells,
      "(", round(100 * n_selected / n_cells, 1), "%)\n")

  result
}
