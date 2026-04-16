#!/usr/bin/env Rscript
# =============================================================================
# clustering/cluster_assignment.R — Unified Cluster Assignment Interface
# =============================================================================
# New module providing a single interface for getting cluster labels
# from Leiden, FlowSOM, FlowJo, or CSV-based clustering
# =============================================================================

if (!exists("shared_dir")) shared_dir <- { .s <- Sys.getenv("CYTOMETRY_SRC", ""); if (nzchar(.s)) file.path(.s, "R_shared") else file.path(dirname(getwd()), "src", "R_shared") }
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "cluster_utils.R"))

#' Get cluster assignments from various sources
#'
#' @param source Type of clustering source: "leiden", "flowjo", "csv", "flowsom"
#' @param ... Arguments passed to the specific loader
#' @return Data frame with at minimum cell_id and cluster columns
get_cluster_assignments <- function(source = c("leiden", "flowjo", "csv", "flowsom"),
                                     ...) {
  source <- match.arg(source)

  switch(source,
    "leiden" = get_clusters_from_leiden(...),
    "flowjo" = get_clusters_from_flowjo(...),
    "csv"    = get_clusters_from_csv(...),
    "flowsom" = get_clusters_from_flowsom(...)
  )
}

#' Load Leiden cluster assignments from CSV
#' @param csv_path Path to leiden_clusters.csv
get_clusters_from_leiden <- function(csv_path) {
  require(data.table)
  dt <- fread(csv_path)
  if (!"cluster" %in% names(dt)) stop("No 'cluster' column found")
  data.frame(
    cell_id = seq_len(nrow(dt)),
    cluster = dt$cluster,
    cluster_name = paste0("Leiden_", dt$cluster)
  )
}

#' Load FlowJo cluster assignments
#' @param wsp_path Path to .wsp file
#' @param fcs_dir FCS directory
get_clusters_from_flowjo <- function(wsp_path, fcs_dir = dirname(wsp_path)) {
  source(file.path(dirname(sys.frame(1)$ofile %||% "."), "flowjo_import.R"))
  result <- import_flowjo_gates(wsp_path, fcs_dir)
  # Convert binary columns to single assignment
  pop_cols <- setdiff(names(result$assignments), c("cell_id", "sample"))
  df <- convert_binary_to_single_cluster(result$assignments, pop_cols)
  # Rename cluster_assignment to cluster for consistency with other loaders
  if ("cluster_assignment" %in% names(df) && !"cluster" %in% names(df)) {
    names(df)[names(df) == "cluster_assignment"] <- "cluster"
  }
  df
}

#' Load cluster assignments from a CSV file
#' @param csv_path Path to CSV with cluster column
#' @param cluster_col Name of the cluster column
get_clusters_from_csv <- function(csv_path, cluster_col = "cluster") {
  require(data.table)
  dt <- fread(csv_path)
  if (!cluster_col %in% names(dt)) {
    # Try to find cluster-like columns
    candidates <- grep("cluster|population|gate", names(dt),
                       ignore.case = TRUE, value = TRUE)
    if (length(candidates) > 0) {
      cluster_col <- candidates[1]
      cat("Using column:", cluster_col, "\n")
    } else {
      stop("No cluster column found in: ", csv_path)
    }
  }
  data.frame(
    cell_id = seq_len(nrow(dt)),
    cluster = dt[[cluster_col]],
    cluster_name = clean_cluster_name(as.character(dt[[cluster_col]]))
  )
}

#' Load FlowSOM cluster assignments
#' @param fsom_object FlowSOM result object
get_clusters_from_flowsom <- function(fsom_object) {
  clusters <- FlowSOM::GetMetaclusters(fsom_object)
  data.frame(
    cell_id = seq_len(length(clusters)),
    cluster = as.integer(clusters),
    cluster_name = paste0("FlowSOM_", clusters)
  )
}
