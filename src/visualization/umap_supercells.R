# =============================================================================
# visualization/umap_supercells.R: SuperCellCyto supercells as the UMAP fit set.
#
# Fitting UMAP on every cell of a multi-million-cell cohort is slow (the
# neighbour graph and the initialisation grow faster than the cell count), so
# the markers layout can fit on a smaller set and project every cell onto it
# (run_umap(fit_matrix = ...)). A random subset is one such set; supercells
# are the other. SuperCellCyto (Putri et al. 2023, Genome Biology) groups each
# sample's cells into supercells of about `gamma` similar cells and summarises
# each by its mean, so the fit set is gamma times smaller than the cohort while
# every phenotype keeps a share in proportion to its cells. Supercells are
# formed within one sample at a time and never mix samples.
#
# SuperCellCyto is optional: only --fit_subset supercell needs it, and that
# says how to install it when it is missing.
# =============================================================================

# Stops with how to install SuperCellCyto when it is missing; called before
# the workspaces load, so a missing package costs seconds, not the load.
require_supercellcyto <- function() {
  if (!requireNamespace("SuperCellCyto", quietly = TRUE)) {
    stop("--fit_subset supercell needs the R package SuperCellCyto. Install ",
         "it with remotes::install_github(\"phipsonlab/SuperCellCyto\") (or ",
         "bioconductor-supercellcyto from bioconda on R 4.5), or use ",
         "--fit_subset random.", call. = FALSE)
  }
  invisible(TRUE)
}

# Mean marker values of the SuperCellCyto supercells of `mat`.
#
# @param mat Numeric matrix, cells x markers (the values the embedding uses).
# @param samples Sample of each row of `mat`; supercells never span samples.
# @param gamma Graining level: about this many cells per supercell (> 1).
# @param k_knn Neighbours of SuperCellCyto's within-sample kNN graph.
# @param n_workers Samples processed in parallel (forked; 1 on Windows).
# @return list(matrix = supercells x markers (same columns as `mat`),
#   supercells = one row per fit row: sample, supercell (number within its
#   sample; NA for a passed-through cell) and cells (its size),
#   cell_supercell = the supercell number of every row of `mat` (NA when
#   passed through), n_cells, n_supercells, passed_through = cells of samples
#   too small to be grained, which enter the fit set as themselves).
supercell_fit_matrix <- function(mat, samples, gamma = 20, k_knn = 5L,
                                 n_workers = 1L) {
  require_supercellcyto()
  if (!is.finite(gamma) || gamma <= 1) {
    stop("supercell gamma must be above 1 (cells per supercell).", call. = FALSE)
  }
  k_knn <- as.integer(k_knn)
  if (is.na(k_knn) || k_knn < 1L) {
    stop("supercell k_knn must be a whole number of at least 1.", call. = FALSE)
  }
  stopifnot(length(samples) == nrow(mat))
  markers <- colnames(mat)

  # A sample needs enough cells for a kNN graph and at least two supercells;
  # smaller ones (rare, e.g. a nearly empty tube) join the fit set cell by cell.
  n_per <- table(samples)
  small <- names(n_per)[n_per < max(2 * gamma, k_knn + 2)]
  is_small <- samples %in% small
  passed <- mat[is_small, , drop = FALSE]

  out <- passed[0, , drop = FALSE]
  out_info <- data.frame(sample = character(0), supercell = integer(0),
                         cells = integer(0), stringsAsFactors = FALSE)
  cell_supercell <- rep(NA_integer_, nrow(mat))
  if (any(!is_small)) {
    # Internal column names: a marker called like the sample or cell id column,
    # or one data.table cannot take as a name, cannot collide or break.
    internal <- sprintf("m%d", seq_along(markers))
    df <- as.data.frame(mat[!is_small, , drop = FALSE])
    names(df) <- internal
    df$sc_sample <- as.character(samples[!is_small])
    df$sc_cell <- seq_len(nrow(df))
    dt <- data.table::as.data.table(df)
    parallel <- n_workers > 1L && .Platform$OS.type == "unix"
    # One task per sample, handed to whichever worker is free, largest sample
    # first (load_balancing sorts them). A sample's time grows faster than its
    # cell count (about n^1.8 on MTB: 80k cells 71 s), and BiocParallel's
    # default of one contiguous chunk per worker put the 11 largest samples of
    # the size-sorted list on one worker, which then ran for an hour after the
    # rest were idle. Neither changes the supercells, only who computes them.
    bp <- if (parallel) {
      BiocParallel::MulticoreParam(workers = as.integer(n_workers),
                                   tasks = length(unique(df$sc_sample)))
    } else {
      BiocParallel::SerialParam()
    }
    res <- SuperCellCyto::runSuperCellCyto(
      dt, markers = internal, sample_colname = "sc_sample",
      cell_id_colname = "sc_cell", n_pc = min(10L, length(internal)),
      gam = gamma, k_knn = k_knn, BPPARAM = bp, load_balancing = parallel)
    sc <- as.data.frame(res$supercell_expression_matrix)
    # In sample order, then supercell number: the rows come back in the order
    # the workers finished, and uwot's result depends on the row order, so
    # without this the embedding would depend on the worker count.
    sc_num <- suppressWarnings(as.integer(sub("^SuperCell_([0-9]+)_.*$", "\\1",
                                              sc$SuperCellId)))
    ord <- order(match(sc$sc_sample, unique(df$sc_sample)), sc_num)
    sc <- sc[ord, , drop = FALSE]
    sc_num <- sc_num[ord]
    out <- as.matrix(sc[, internal, drop = FALSE])
    colnames(out) <- markers
    # Which supercell every cell went into, by number within its sample (the
    # sample is already a column of every per-cell table).
    map <- as.data.frame(res$supercell_cell_map)
    map_num <- suppressWarnings(as.integer(sub("^SuperCell_([0-9]+)_.*$", "\\1",
                                               map$SuperCellID)))
    grained <- which(!is_small)
    cell_supercell[grained[as.integer(map$CellId)]] <- map_num
    size <- table(paste(map$Sample, map_num))
    out_info <- data.frame(
      sample = sc$sc_sample, supercell = sc_num,
      cells = as.integer(size[paste(sc$sc_sample, sc_num)]),
      stringsAsFactors = FALSE)
  }
  fit <- rbind(out, passed)
  rownames(fit) <- NULL
  out_info <- rbind(out_info, data.frame(
    sample = as.character(samples[is_small]), supercell = rep(NA_integer_, nrow(passed)),
    cells = rep(1L, nrow(passed)), stringsAsFactors = FALSE))
  list(matrix = fit, supercells = out_info, cell_supercell = cell_supercell,
       n_cells = nrow(mat), n_supercells = nrow(out),
       passed_through = nrow(passed))
}
