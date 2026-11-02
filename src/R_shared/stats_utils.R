# =============================================================================
# R_shared/stats_utils.R — Statistical Analysis Utility Functions
# =============================================================================
# Deduplicated from: 01_split_clusters.R, GLMM.R, Fine_clustering.R, Stats.R
# =============================================================================

#' Fisher's exact test for cluster enrichment between two conditions
#' @param cluster_col Column name with binary (0/1) cluster membership
#' @param cells_data Data frame with cluster_col and 'condition' columns
#' @param cond_levels Character vector of length 2 with condition levels
#' @return tibble with test results
calculate_enrichment <- function(cluster_col, cells_data, cond_levels) {
  inC <- as.logical(cells_data[[cluster_col]])
  condA <- cells_data$condition == cond_levels[1]
  condB <- cells_data$condition == cond_levels[2]

  a <- sum(inC & condA, na.rm = TRUE)
  b <- sum(!inC & condA, na.rm = TRUE)
  cc <- sum(inC & condB, na.rm = TRUE)
  d <- sum(!inC & condB, na.rm = TRUE)

  if (any(c(a, b, cc, d) < 0) || (a + b == 0) || (cc + d == 0)) {
    return(tibble::tibble(
      in_A = a, not_in_A = b, in_B = cc, not_in_B = d,
      odds_ratio = NA, p_value = NA, ci_lower = NA, ci_upper = NA,
      method = "fisher_failed"
    ))
  }

  m <- matrix(c(a, b, cc, d), nrow = 2, byrow = TRUE,
              dimnames = list(c("inCluster", "notInCluster"), cond_levels))
  tryCatch({
    ft <- fisher.test(m)
    tibble::tibble(
      in_A = a, not_in_A = b, in_B = cc, not_in_B = d,
      odds_ratio = unname(ft$estimate),
      p_value = ft$p.value,
      ci_lower = ft$conf.int[1],
      ci_upper = ft$conf.int[2],
      method = "fisher"
    )
  }, error = function(e) {
    tibble::tibble(
      in_A = a, not_in_A = b, in_B = cc, not_in_B = d,
      odds_ratio = NA, p_value = NA, ci_lower = NA, ci_upper = NA,
      method = "fisher_error"
    )
  })
}

#' Mixed-effects logistic regression for cluster enrichment
#' @param cluster_col Column name with binary cluster membership
#' @param cells_data Data frame with cluster_col, condition, file_name columns
#' @param condition_levels Character vector of length 2
#' @return tibble with model results
do_mixed_effects <- function(cluster_col, cells_data, condition_levels) {
  model_data <- cells_data %>%
    dplyr::select(file_name, condition, cluster = dplyr::all_of(cluster_col)) %>%
    dplyr::mutate(
      condition = factor(condition, levels = condition_levels),
      file_name = factor(file_name),
      cluster = as.numeric(cluster)
    ) %>%
    dplyr::filter(!is.na(cluster), !is.na(condition))

  if (length(unique(model_data$condition)) < 2) {
    return(tibble::tibble(
      estimate = NA, std_error = NA, p_value = NA,
      odds_ratio = NA, ci_lower = NA, ci_upper = NA,
      method = "insufficient_conditions"
    ))
  }

  cluster_sum <- sum(model_data$cluster, na.rm = TRUE)
  if (cluster_sum == 0 || cluster_sum == nrow(model_data)) {
    return(tibble::tibble(
      estimate = NA, std_error = NA, p_value = NA,
      odds_ratio = NA, ci_lower = NA, ci_upper = NA,
      method = "no_cluster_variation"
    ))
  }

  tryCatch({
    model <- lme4::glmer(
      cluster ~ condition + (1 | file_name),
      data = model_data, family = binomial,
      control = lme4::glmerControl(optimizer = "bobyqa",
                                   optCtrl = list(maxfun = 20000))
    )
    fixed_effects <- lme4::fixef(model)
    se <- sqrt(diag(vcov(model)))
    cond_idx <- grep("condition", names(fixed_effects))
    if (length(cond_idx) == 0) {
      return(tibble::tibble(
        estimate = NA, std_error = NA, p_value = NA,
        odds_ratio = NA, ci_lower = NA, ci_upper = NA,
        method = "no_condition_effect"
      ))
    }
    est <- fixed_effects[cond_idx]
    se_val <- se[cond_idx]
    z_val <- est / se_val
    p_val <- 2 * (1 - pnorm(abs(z_val)))

    tibble::tibble(
      estimate = est, std_error = se_val, p_value = p_val,
      odds_ratio = exp(est),
      ci_lower = exp(est - 1.96 * se_val),
      ci_upper = exp(est + 1.96 * se_val),
      method = "mixed_effects"
    )
  }, error = function(e) {
    tibble::tibble(
      estimate = NA, std_error = NA, p_value = NA,
      odds_ratio = NA, ci_lower = NA, ci_upper = NA,
      method = paste0("mixed_effects_error: ", e$message)
    )
  })
}

#' Check overdispersion in a Poisson model
#' @param model A fitted glm or glmer model
#' @return List with ratio, p_value, and is_overdispersed
check_overdispersion <- function(model) {
  tryCatch({
    if (requireNamespace("performance", quietly = TRUE)) {
      result <- performance::check_overdispersion(model)
      list(
        ratio = result$dispersion_ratio,
        p_value = result$p_value,
        is_overdispersed = result$p_value < 0.05
      )
    } else {
      # Manual check
      resid_dev <- deviance(model)
      df_resid <- df.residual(model)
      ratio <- resid_dev / df_resid
      p_val <- pchisq(resid_dev, df_resid, lower.tail = FALSE)
      list(ratio = ratio, p_value = p_val, is_overdispersed = ratio > 1.5)
    }
  }, error = function(e) {
    list(ratio = NA, p_value = NA, is_overdispersed = NA)
  })
}

#' Fit a negative binomial GLMM for differential abundance
#' @param count_data Data frame with count, condition, sample_id columns
#' @param formula Model formula
#' @param max_iter Maximum optimizer iterations
#' @return Fitted model object or NULL on failure
fit_nb_glmm <- function(count_data, formula = count ~ condition + (1 | sample_id),
                         max_iter = 20000) {
  tryCatch({
    lme4::glmer.nb(
      formula, data = count_data,
      control = lme4::glmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = max_iter)
      )
    )
  }, error = function(e) {
    warning("NB GLMM failed: ", e$message, ". Falling back to Poisson GLMM. ",
            "Poisson assumes mean=variance; p-values may be anti-conservative.")
    tryCatch({
      message("Falling back to Poisson GLMM...")
      lme4::glmer(
        formula, data = count_data, family = poisson,
        control = lme4::glmerControl(
          optimizer = "bobyqa",
          optCtrl = list(maxfun = max_iter)
        )
      )
    }, error = function(e2) {
      message("Poisson GLMM also failed: ", e2$message)
      NULL
    })
  })
}

#' Run PCA on an expression matrix
#' @param expr_matrix Numeric matrix (rows = cells, columns = markers)
#' @param n_components Number of PCs to compute
#' @param scale Scale variables before PCA
#' @return List with coordinates, variance_explained, loadings
run_pca <- function(expr_matrix, n_components = 2, scale = TRUE) {
  pca_result <- prcomp(expr_matrix, center = TRUE, scale. = scale,
                       rank. = n_components)
  var_explained <- (pca_result$sdev^2 / sum(pca_result$sdev^2))[1:n_components]

  list(
    coordinates = pca_result$x[, 1:n_components, drop = FALSE],
    variance_explained = var_explained,
    loadings = pca_result$rotation[, 1:n_components, drop = FALSE],
    pca_object = pca_result
  )
}

#' Run UMAP on an expression matrix
#' @param expr_matrix Numeric matrix
#' @param n_neighbors Number of neighbors
#' @param min_dist Minimum distance parameter
#' @param n_components Number of output dimensions
#' @param metric Distance metric
#' @param seed Random seed
#' @param y Optional labels for supervised UMAP (a factor for categories).
#'   Either one per row of expr_matrix or one per FITTED row; NULL = unsupervised.
#' @param target_weight uwot's weight of `y` against the data (0 = data only,
#'   1 = labels only). Used only when `y` is given.
#' @param fit_idx Optional sorted, distinct row indices to fit on; the other
#'   rows are projected with umap_transform. Takes precedence over `subsample`.
#' @param transform_threads Threads for the projection's neighbour search
#'   (default: n_threads). The search gives the same neighbours at any thread
#'   count and the projection's optimiser stays serial, so this changes only
#'   the speed, never the coordinates.
#' @param fit_matrix Optional matrix with the same columns as expr_matrix (e.g.
#'   supercell means) to fit on instead of rows of expr_matrix; every row of
#'   expr_matrix is then projected. Unsupervised only; excludes fit_idx.
#' @param model Optional model from an earlier fit with the same data and
#'   settings (uwot::load_uwot): the fit is skipped and the cells are only
#'   projected, as a resumed run does.
#' @param on_fit Optional function(model), called once the model is fitted and
#'   before any cell is projected (e.g. to save it: on a cohort the projection
#'   takes most of an hour, and a failure there must not cost the fit).
#' @return Matrix of UMAP coordinates. When a model was fitted on a subset
#'   (subsample, fit_idx or fit_matrix) it rides along as attr(, "model"), so a
#'   caller can save it (uwot::save_uwot) and project further cells later.
run_umap <- function(expr_matrix, n_neighbors = 15, min_dist = 0.1,
                     n_components = 2, metric = "euclidean", seed = 42,
                     n_threads = 1L, progress = FALSE, verbose = FALSE,
                     init = "auto", subsample = 0L, chunk = 250000L,
                     y = NULL, target_weight = 0.5, fit_idx = NULL,
                     transform_threads = NULL, fit_matrix = NULL,
                     on_fit = NULL, model = NULL) {
  set.seed(seed)
  n_threads <- max(1L, as.integer(n_threads))
  n_obs <- nrow(expr_matrix)

  # A caller-chosen fit set (e.g. one balanced across populations) replaces the
  # random subsample. Checked here because a duplicated or unsorted index would
  # not fail in uwot: it would silently give some cells two coordinates and
  # leave others to the projection.
  if (!is.null(fit_idx)) {
    fit_idx <- as.integer(fit_idx)
    if (!length(fit_idx) || anyNA(fit_idx) || any(fit_idx < 1L) ||
        any(fit_idx > n_obs) || anyDuplicated(fit_idx) ||
        is.unsorted(fit_idx)) {
      stop("run_umap: fit_idx must be sorted, distinct row indices between 1 ",
           "and ", n_obs, ".", call. = FALSE)
    }
  }
  # A separate fit matrix (supercell means) stands in for the cells: the model
  # is fitted on its rows and every cell is projected. Its columns must be the
  # cells' columns in the same order, or the projection would mix markers up.
  if (!is.null(fit_matrix)) {
    if (!is.null(fit_idx) || !is.null(y)) {
      stop("run_umap: fit_matrix cannot be combined with fit_idx or y.",
           call. = FALSE)
    }
    if (!identical(colnames(fit_matrix), colnames(expr_matrix)) ||
        nrow(fit_matrix) < 2L || !all(is.finite(fit_matrix))) {
      stop("run_umap: fit_matrix must be finite, have at least 2 rows and ",
           "the same columns as expr_matrix, in the same order.", call. = FALSE)
    }
  }
  # The rows the model is fitted on, which is what the init rule below is
  # about: a spectral init on a 100k fit set is cheap however many cells are
  # projected afterwards. Equal to n_obs whenever fit_idx is not given.
  n_fit <- if (!is.null(fit_matrix)) nrow(fit_matrix)
           else if (is.null(fit_idx)) n_obs else length(fit_idx)

  # ---------------------------------------------------------------------------
  # WHY THIS FUNCTION IS MORE THAN A CALL TO uwot::umap
  #
  # uwot has two thread pools and three stages, and only some of them scale:
  #
  #   neighbour search   n_threads      parallel
  #   smooth-kNN + graph —              largely serial
  #   INITIALISATION     —              SERIAL, and on millions of points this
  #                                     is the whole run. init = "spectral" (the
  #                                     default) is an eigendecomposition of a
  #                                     graph with one node per cell.
  #   optimisation       n_sgd_threads  parallel
  #
  # Measured on this repo's data: a run over 7.9M cells sat at ~180% CPU with 29
  # threads alive for hours — 28 workers idle while one did the initialisation.
  # More cores cannot fix that; only a different algorithm can. Hence `init`
  # and `subsample` below.
  # ---------------------------------------------------------------------------

  # "auto": spectral is the best initialisation and is affordable up to a few
  # hundred thousand points. Past that it stops paying for itself — and then
  # stops working at all.
  #
  # MEASURED, 500,000 cells x 7 markers, 12 threads:
  #     init = "spectral"   2023 s   30 min of it inside irlba, which then
  #                                  printed "Spectral initialization failed to
  #                                  converge, using random initialization
  #                                  instead" and fell back to random anyway
  #     init = "pca"         166 s
  #     init = "random"      156 s
  #
  # So above the limit spectral is not a better layout bought with time; it is
  # half an hour spent arriving at the random layout. The optimisation itself
  # was 47 s of that 166 s — the neighbour search and the initialisation are
  # the run, which is why more cores did not help.
  #
  # 300,000 is set below the size at which convergence was observed to fail,
  # not at it. The switch is announced, because it changes the embedding.
  SPECTRAL_LIMIT <- 300000L
  # Projection onto a given model fits nothing, so there is no start to pick.
  if (!is.null(model) && is.null(fit_matrix) && is.null(fit_idx)) {
    n_fit <- min(n_fit, SPECTRAL_LIMIT)
  }
  if (identical(init, "auto")) {
    init <- if (n_fit > SPECTRAL_LIMIT) "pca" else "spectral"
    if (init == "pca") {
      message("uwot: ", n_fit, " cells is past the ", SPECTRAL_LIMIT,
              " point where spectral initialisation stops paying for itself ",
              "(it is single-threaded); using init = \"pca\". Pass ",
              "init = \"spectral\" to force it.")
    }
  }

  cb <- NULL
  if (isTRUE(progress) && exists("emit_progress", mode = "function")) {
    cb <- function(epoch, n_epochs, coords) {
      emit_progress("embed", epoch, n_epochs, "optimising")
      NULL
    }
  }

  # THREADS AND REPRODUCIBILITY. n_sgd_threads defaults to asynchronous Hogwild
  # SGD, where the interleaving of updates decides the result — so the same seed
  # gives a different embedding. batch = TRUE updates at the end of each epoch
  # instead, which makes multi-threaded SGD reproducible again. Serial stays the
  # default so an existing analysis does not change under someone's feet.
  base_args <- list(n_neighbors = n_neighbors,
                    min_dist = min_dist,
                    n_components = n_components,
                    metric = metric,
                    init = init,
                    n_threads = n_threads,
                    verbose = isTRUE(verbose))
  if (n_threads > 1L) {
    base_args$n_sgd_threads <- n_threads
    base_args$batch <- TRUE
  }
  if (!is.null(cb)) base_args$epoch_callback <- cb
  # A PCA start is SCALED. Unscaled, it is the data's own principal component
  # scores: on FlowJo values (0-256) the points start ~100 units apart, the
  # gradients are tiny at that distance, and the optimisation ends with a map
  # of specks scattered over a range of 200 (2026-10-14, MTB: 405,521
  # supercells, which crossed SPECTRAL_LIMIT into "pca"). uwot's own advice is
  # to scale it; "range" gives the 0-10 box the spectral start is put in.
  if (identical(init, "pca")) base_args$init_sdev <- "range"

  subsample <- suppressWarnings(as.integer(subsample))
  if (is.na(subsample)) subsample <- 0L

  # Supervision only enters the FIT: umap_transform places new cells by their
  # data alone, which is what lets a supervised model project cells it never
  # saw. `y` may be given per row of the matrix or per fitted row.
  fit_args <- function(rows) {
    if (is.null(y)) return(base_args)
    yy <- if (length(y) == length(rows)) y
          else if (length(y) == n_obs) y[rows]
          else stop("run_umap: y has ", length(y), " entries; expected ",
                    n_obs, " (one per row) or ", length(rows),
                    " (one per fitted row).", call. = FALSE)
    c(base_args, list(y = yy, target_weight = target_weight))
  }

  # --- the whole matrix, when it is small enough to be worth it ---------------
  # (checked for a degenerate spread like every other fit, below)
  # A model given with no fit set of its own is a map from another run: every
  # row is projected onto it, nothing is fitted.
  project_only <- !is.null(model) && is.null(fit_matrix) && is.null(fit_idx)
  fit_all <- if (!is.null(fit_matrix) || project_only) FALSE
             else if (is.null(fit_idx)) subsample <= 0L || subsample >= n_obs
             else n_fit == n_obs
  if (fit_all) {
    coords <- do.call(uwot::umap, c(list(X = expr_matrix),
                                    fit_args(seq_len(n_obs))))
    .warn_spread(coords)
    return(coords)
  }

  # --- fit on a subsample, project the rest ----------------------------------
  #
  # Every stage above is superlinear in the number of points; the projection of
  # an already-fitted model is not. So the expensive part is paid once on a
  # sample big enough to describe the manifold, and the remaining millions are
  # transformed onto it — which IS threaded and is roughly linear.
  #
  # This is a different estimator, not an approximation of the full run: the
  # fitted cells get an exact embedding and the rest get an out-of-sample
  # projection onto it. It is said out loud for that reason.
  # Said only when a fit happens; a resumed run says it reuses the model.
  say <- function(...) if (is.null(model)) message(...)
  if (project_only) {
    idx <- integer(0)
  } else if (!is.null(fit_matrix)) {
    say("uwot: fitting on ", n_fit, " supplied rows (e.g. supercells), ",
        "then projecting all ", n_obs, " cells onto them.")
    idx <- integer(0)
  } else if (is.null(fit_idx)) {
    say("uwot: fitting on a ", subsample, "-cell subsample of ", n_obs,
        ", then projecting the remaining ", n_obs - subsample, ".")
    idx <- sort(sample.int(n_obs, subsample))
  } else {
    say("uwot: fitting on the ", n_fit, " chosen cells of ", n_obs,
        ", then projecting the remaining ", n_obs - n_fit, ".")
    idx <- fit_idx
  }
  if (is.null(model)) {
    model <- do.call(uwot::umap, c(list(
      X = if (is.null(fit_matrix)) expr_matrix[idx, , drop = FALSE] else fit_matrix,
      ret_model = TRUE), fit_args(idx)))
    if (is.function(on_fit)) on_fit(model)
  } else if (!project_only) {
    if (nrow(model$embedding) != (if (length(idx)) length(idx) else n_fit)) {
      stop("run_umap: the given model was fitted on ", nrow(model$embedding),
           " rows, not on this fit set of ", n_fit, ".", call. = FALSE)
    }
    message("uwot: reusing the saved model; projecting only.")
  } else {
    message("uwot: projecting all ", n_obs, " cells onto the given model.")
  }

  .warn_spread(model$embedding)

  coords <- matrix(NA_real_, nrow = n_obs, ncol = n_components)
  if (length(idx)) coords[idx, ] <- model$embedding

  rest <- setdiff(seq_len(n_obs), idx)
  t_threads <- max(1L, as.integer(if (is.null(transform_threads)) n_threads
                                  else transform_threads))
  # Chunked so the transform reports progress and so a 7-million-row projection
  # does not have to be resident all at once.
  chunks <- split(rest, ceiling(seq_along(rest) / max(1L, as.integer(chunk))))
  for (i in seq_along(chunks)) {
    rows <- chunks[[i]]
    coords[rows, ] <- uwot::umap_transform(
      expr_matrix[rows, , drop = FALSE], model,
      n_threads = t_threads,
      # Threaded optimisation only for a model fitted in batch mode: batch
      # updates do not depend on the thread count, Hogwild ones do.
      n_sgd_threads = if (isTRUE(model$batch)) t_threads else 0,
      # Seeded per chunk, so a projection resumed from a saved model gives the
      # coordinates the uninterrupted run would have.
      seed = seed + i,
      verbose = FALSE)
    if (isTRUE(progress) && exists("emit_progress", mode = "function")) {
      emit_progress("embed", i, length(chunks), "projecting")
    }
  }

  # Every plot downstream indexes coordinates by row against the expression
  # matrix, so a hole here would silently mis-colour cells rather than fail.
  if (anyNA(coords)) {
    stop("umap: ", sum(!stats::complete.cases(coords)), " cell(s) were left ",
         "without coordinates by the subsample projection.", call. = FALSE)
  }
  attr(coords, "model") <- model
  coords
}

# A UMAP layout with spread 1 spans a few tens of units. Far more means the
# optimisation never pulled the points together (a badly scaled start), and
# every plot drawn from it would be specks; said loudly rather than plotted
# silently.
.warn_spread <- function(embedding) {
  width <- max(apply(embedding, 2, function(v) diff(range(v, finite = TRUE))))
  if (is.finite(width) && width > 100) {
    warning("umap: the layout spans ", round(width), " units, where a normal ",
            "one spans a few tens: the optimisation did not converge, and the ",
            "plots will show scattered specks. Try init = \"spectral\" on a ",
            "smaller fit set.", call. = FALSE)
  }
  invisible(width)
}

#' Calculate summary statistics for expression data grouped by cluster
#' @param expr_data Data frame with expression columns and a cluster column
#' @param cluster_col Name of the cluster column
#' @param marker_cols Character vector of marker column names
#' @return tibble with per-cluster per-marker statistics
calculate_statistics <- function(expr_data, cluster_col, marker_cols) {
  expr_data %>%
    tidyr::pivot_longer(cols = dplyr::all_of(marker_cols),
                        names_to = "marker", values_to = "expression") %>%
    dplyr::group_by(!!rlang::sym(cluster_col), marker) %>%
    dplyr::summarise(
      mean = mean(expression, na.rm = TRUE),
      median = median(expression, na.rm = TRUE),
      sd = sd(expression, na.rm = TRUE),
      q25 = quantile(expression, 0.25, na.rm = TRUE),
      q75 = quantile(expression, 0.75, na.rm = TRUE),
      n = dplyr::n(),
      .groups = "drop"
    )
}

#' Get top features by variance or mean
#' @param expr_matrix Numeric matrix
#' @param n_top Number of top features to return
#' @param method "variance" or "mean"
#' @return Character vector of top feature names
get_top_features <- function(expr_matrix, n_top = 20, method = "variance") {
  scores <- switch(method,
    "variance" = apply(expr_matrix, 2, var, na.rm = TRUE),
    "mean"     = apply(expr_matrix, 2, mean, na.rm = TRUE),
    apply(expr_matrix, 2, var, na.rm = TRUE)
  )
  names(sort(scores, decreasing = TRUE))[1:min(n_top, length(scores))]
}
