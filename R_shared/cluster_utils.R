# =============================================================================
# R_shared/cluster_utils.R — Cluster Name & Label Utilities
# =============================================================================
# Deduplicated from: 01_split_clusters.R, GLMM.R, Fine_clustering.R, Functions.R
# =============================================================================

#' Clean cluster name by removing common prefixes/suffixes
#' @param cluster_name Character string cluster name
#' @return Cleaned cluster name
clean_cluster_name <- function(cluster_name) {
  if (is.na(cluster_name) || cluster_name == "") return("Unknown")
  cleaned <- gsub("^AI_Cluster_", "", cluster_name)
  cleaned <- gsub(" \\(Clustering annotated\\)$", "", cleaned)
  cleaned
}

#' Extract cell type group from cluster name
#' @param cluster_name Character string cluster name
#' @return Group label string
extract_cluster_group <- function(cluster_name) {
  cleaned_name <- clean_cluster_name(cluster_name)

  group_patterns <- list(
    "BCs"            = c("^BCs(_|$)", "^B_?[Cc]ells?(_|$)"),
    "CD3-"           = c("^CD3-(_|$)", "^CD3_negative(_|$)", "^CD3neg(_|$)"),
    "CD4TCs"         = c("^CD4TCs(_|$)", "^CD4_?[Tt](_|$)", "^CD4_?[Tt][Cc]ells?(_|$)"),
    "CD8TCs"         = c("^CD8TCs(_|$)", "^CD8_?[Tt](_|$)", "^CD8_?[Tt][Cc]ells?(_|$)"),
    "NKCs"           = c("^NKCs(_|$)", "^NK_?[Cc]ells?(_|$)", "^Natural_?[Kk]iller(_|$)"),
    "NKTCells"       = c("^NKTCells(_|$)", "^NKT_?[Cc]ells?(_|$)", "^NK_?[Tt](_|$)"),
    "gdTCs"          = c("^gdTCs(_|$)", "^[gG]amma_?[dD]elta(_|$)", "^gdTCR(_|$)",
                         "^gdTCR-CD56-TCells(_|$)"),
    "TCells"         = c("^TCells(_|$)", "^T_?[Cc]ells?(_|$)"),
    "LymphoidCells"  = c("^LymphoidCells(_|$)", "^Lymphoid_?[Cc]ells?(_|$)")
  )

  for (group_name in names(group_patterns)) {
    for (pattern in group_patterns[[group_name]]) {
      if (grepl(pattern, cleaned_name, perl = TRUE)) {
        if (grepl("_[0-9]+$", cleaned_name)) {
          number_part <- gsub(".*_([0-9]+)$", "\\1", cleaned_name)
          return(paste0(group_name, "_", number_part))
        }
        return(group_name)
      }
    }
  }
  cleaned_name
}

#' Create a short cluster identifier
#' @param cluster_name Full cluster name
#' @param max_length Maximum length for the ID
#' @return Short identifier string
create_cluster_id <- function(cluster_name, max_length = 30) {
  id <- clean_cluster_name(cluster_name)
  id <- gsub("[^A-Za-z0-9_-]", "_", id)
  if (nchar(id) > max_length) id <- substr(id, 1, max_length)
  id
}

#' Flatten a list of label vectors to single labels using gating tree depth
#' @param label_lists List where each element is a character vector of labels
#' @param gating_tree Named numeric vector of label depths
#' @return Character vector of single labels
flatten_label_list <- function(label_lists, gating_tree) {
  sapply(label_lists, function(labels) {
    if (length(labels) == 0) return("root")
    if (length(labels) == 1) return(labels[[1]])
    depths <- gating_tree[labels]
    depths[is.na(depths)] <- -Inf
    labels[which.max(depths)]
  }, USE.NAMES = FALSE)
}

#' Extract labels from a binary (0/1) data.frame
#' @param df Data frame with binary columns
#' @return List of character vectors (selected column names per row)
labels_from_binary_df <- function(df) {
  apply(df, 1, function(row) names(row)[which(row == 1)])
}

#' Extract labels from a binary data.table (first column = index)
#' @param dt data.table with first column as index and remaining as binary
#' @return List of character vectors
labels_from_binary_dt <- function(dt) {
  mat <- as.matrix(dt[, -1])
  col_names <- names(dt)[-1]
  apply(mat, 1, function(row) {
    selected <- col_names[which(row == 1)]
    if (length(selected) == 0) selected <- character(0)
    selected
  }, simplify = FALSE)
}

#' Convert binary cluster columns to a single cluster assignment column
#' @param cells_data Data frame with binary cluster columns
#' @param cluster_cols Character vector of binary cluster column names
#' @return Data frame with added 'cluster_assignment' column
convert_binary_to_single_cluster <- function(cells_data, cluster_cols) {
  message("Converting ", length(cluster_cols),
          " binary cluster columns to single cluster assignment...")
  cells_data$cluster_assignment <- NA_character_

  for (i in seq_len(nrow(cells_data))) {
    memberships <- character(0)
    for (col in cluster_cols) {
      if (!is.na(cells_data[[col]][i]) && cells_data[[col]][i] == 1) {
        memberships <- c(memberships, extract_cluster_group(col))
      }
    }
    cells_data$cluster_assignment[i] <- if (length(memberships) == 0) {
      "Unassigned"
    } else {
      memberships[1]
    }
  }
  cells_data
}

#' Translate label strings using a dictionary
#' @param string_list Character vector of labels
#' @param dictionary Named character vector (names = original, values = standardized)
#' @return Character vector with translated labels
translate_labels <- function(string_list, dictionary) {
  idx <- match(string_list, names(dictionary))
  ifelse(!is.na(idx), dictionary[idx], string_list)
}
