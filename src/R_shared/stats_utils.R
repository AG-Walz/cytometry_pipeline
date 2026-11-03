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
#' @param sgd_threads Threads for the optimiser; more than 1 switches uwot to
#'   its batch optimiser (reproducible, but different coordinates).
#' @param nn_method Neighbour index: NULL (uwot decides), "annoy" or "hnsw".
# -----------------------------------------------------------------------------
# uwot's stage announcements, turned into progress a GUI can draw.
#
# Between "Read N rows" and the first optimisation epoch uwot does the work that
# dominates a large run — the Annoy index, the kNN query, the smooth-kNN
# calibration, the initialisation — and none of it is countable. A watcher
# therefore saw an empty bar for the longest part of the run, which is
# indistinguishable from a hang.
#
# uwot does ANNOUNCE each stage: tsmessage() is message() with a timestamp, so
# the announcements are conditions and a calling handler can translate them. The
# handler deliberately does NOT muffle them, so the log still reads as it did.
#
# Matching on another package's wording is a parser for a format nobody agreed
# to, and it will drift. It is confined to this one function for exactly that
# reason, it matches only on stems that have been stable across uwot releases,
# and an unrecognised message costs nothing — the phase simply does not advance,
# which is the behaviour we had before.
# -----------------------------------------------------------------------------
emit_uwot_stage <- function(text) {
  if (!exists("emit_progress", mode = "function")) return(invisible(NULL))
  # 3 = index, query, calibration. The graph is what these three produce.
  if (grepl("Building Annoy index", text, fixed = TRUE)) {
    emit_progress("embed_nn", 0L, 3L, "building the neighbour index")
  } else if (grepl("Searching Annoy index", text, fixed = TRUE)) {
    emit_progress("embed_nn", 1L, 3L, "searching the neighbour index")
  } else if (grepl("Commencing smooth kNN", text, fixed = TRUE)) {
    emit_progress("embed_nn", 2L, 3L, "calibrating distances")
  } else if (grepl("Initializing from", text, fixed = TRUE)) {
    # The first word of the initialisation is also the last word of the
    # neighbour graph — uwot does not say "finished", so its successor says it.
    emit_progress("embed_nn", 3L, 3L, "neighbour graph complete")
    emit_progress("embed_init", 0L, 1L,
                  trimws(sub(".*Initializing from ", "", text)))
  } else if (grepl("Commencing optimization", text, fixed = TRUE)) {
    emit_progress("embed_init", 1L, 1L, "layout initialised")
  }
  invisible(NULL)
}

#' Run `expr`, translating uwot's stage messages into progress lines.
with_uwot_progress <- function(expr) {
  withCallingHandlers(expr, message = function(m) {
    emit_uwot_stage(conditionMessage(m))
    # No invokeRestart("muffleMessage"): the message is still uwot's to print.
  })
}

run_umap <- function(expr_matrix, n_neighbors = 15, min_dist = 0.1,
                     n_components = 2, metric = "euclidean", seed = 42,
                     n_threads = 1L, progress = FALSE, verbose = FALSE,
                     init = "auto", subsample = 0L, chunk = 250000L,
                     y = NULL, target_weight = 0.5, fit_idx = NULL,
                     transform_threads = NULL, fit_matrix = NULL,
                     on_fit = NULL, model = NULL,
                     sgd_threads = 1L, nn_method = NULL) {
  set.seed(seed)
  n_threads   <- max(1L, as.integer(n_threads))
  sgd_threads <- max(1L, as.integer(sgd_threads))
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

  # The start is NEVER switched behind the user's back. "auto" is the spectral
  # start, uwot's own default and the best one. It used to become PCA above
  # 300,000 fitted cells, silently in all but a log line, and that changed the
  # embedding without anyone choosing it (J445, 2026-11-02: 2.86M cells, an
  # automatic PCA start, a map of specks).
  #
  # What the threshold measured is still true, so it is said instead:
  # MEASURED, 500,000 cells x 7 markers, 12 threads:
  #     init = "spectral"   2023 s   30 min of it inside irlba, which then
  #                                  printed "Spectral initialization failed to
  #                                  converge, using random initialization
  #                                  instead" and fell back to random anyway
  #     init = "pca"         166 s
  # and 300,000 cells with a spectral start took 98-188 s in total
  # (n_neighbors 15-80, 8 threads, 2026-11-03). Past SPECTRAL_LIMIT the run
  # warns, at once, and names the two ways out; which one is the user's call.
  SPECTRAL_LIMIT <- 300000L
  if (identical(init, "auto")) init <- "spectral"
  # Projection onto a given model fits nothing, so there is no start to warn
  # about.
  project_only_start <- !is.null(model) && is.null(fit_matrix) && is.null(fit_idx)
  # A random subsample is a fit set too (it is applied further down).
  sub_n <- suppressWarnings(as.integer(subsample))
  n_start <- if (is.null(fit_matrix) && is.null(fit_idx) && !is.na(sub_n) &&
                 sub_n > 0L) min(sub_n, n_fit) else n_fit
  if (identical(init, "spectral") && !project_only_start &&
      n_start > SPECTRAL_LIMIT) {
    umap_init_warning(n_start, SPECTRAL_LIMIT)
  }

  cb <- NULL
  if (isTRUE(progress) && exists("emit_progress", mode = "function")) {
    cb <- function(epoch, n_epochs, coords) {
      emit_progress("embed", epoch, n_epochs, "optimising")
      NULL
    }
  }

  # THREADS AND REPRODUCIBILITY — TWO POOLS, AND ONLY ONE OF THEM IS A RISK.
  #
  # uwot takes `n_threads` for the neighbour search (Annoy index build, kNN
  # query, smooth-kNN calibration) and `n_sgd_threads` for the optimiser, and
  # `n_sgd_threads` defaults to 0 — serial — INDEPENDENTLY of n_threads. Only
  # the optimiser is order-sensitive: asynchronous Hogwild SGD lets the
  # interleaving of updates decide the result, so the same seed gives a
  # different embedding, and batch = TRUE (update at the end of each epoch)
  # is what makes threaded SGD reproducible again — at a different answer.
  #
  # The neighbour search is deterministic at ANY thread count. This function
  # used to tie the two together, so the only way to ask for a threaded search
  # was to also switch optimisers; the default therefore pinned BOTH to one
  # core, and a 9.8M-cell run spent 22 minutes building an Annoy index on a
  # single thread of a 28-core machine to protect coordinates the search could
  # not have changed. Verified on 20,000 x 6: n_threads 1 vs 8 with the
  # optimiser serial gives max|difference| = 0, while batch = TRUE gives 7.09.
  #
  # So: thread the search freely, and keep the optimiser serial unless someone
  # asks for it by name.
  base_args <- list(n_neighbors = n_neighbors,
                    min_dist = min_dist,
                    n_components = n_components,
                    metric = metric,
                    init = init,
                    n_threads = n_threads,
                    n_sgd_threads = 0L,
                    verbose = isTRUE(verbose))
  if (sgd_threads > 1L) {
    base_args$n_sgd_threads <- sgd_threads
    base_args$batch <- TRUE
    # The batch optimiser applies each epoch's updates at once, so it needs more
    # epochs than Hogwild to converge; uwot's documentation says so, and its
    # default for a large input (200) is set for the serial optimiser. 500 is
    # uwot's own count for small inputs. A 2.86M-cell batch run at 200 epochs
    # (J445, 2026-11-02) ended with clusters that had not finished pulling
    # together.
    base_args$n_epochs <- 500L
  }

  # THE ONE PART OF THE NEIGHBOUR SEARCH THREADS CANNOT REACH.
  #
  # uwot's Annoy path is annoy_build() then annoy_search(), and only the SECOND
  # takes n_threads — building the index is single-threaded however many cores
  # are free. On 9.8M cells that was 22 minutes of a 97-minute neighbour phase,
  # and no thread count will move it.
  #
  # hnsw_build() DOES take n_threads, so nn_method = "hnsw" is the only way to
  # parallelise the whole search. It is a different approximate index, so it
  # gives a different neighbour graph and therefore different coordinates —
  # which is why it is opt-in and never inferred.
  #
  # NULL (the default) is left alone rather than replaced with "annoy": uwot
  # picks an EXACT search for small inputs, and forcing the approximate one
  # would quietly change the answer for every small run.
  if (!is.null(nn_method) && nzchar(nn_method)) {
    if (identical(nn_method, "hnsw") &&
        !requireNamespace("RcppHNSW", quietly = TRUE)) {
      stop("nn_method = \"hnsw\" needs the RcppHNSW package.", call. = FALSE)
    }
    base_args$nn_method <- nn_method
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
    coords <- with_uwot_progress(
      do.call(uwot::umap, c(list(X = expr_matrix), fit_args(seq_len(n_obs)))))
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
    model <- with_uwot_progress(do.call(uwot::umap, c(list(
      X = if (is.null(fit_matrix)) expr_matrix[idx, , drop = FALSE] else fit_matrix,
      ret_model = TRUE), fit_args(idx))))
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

#' Is this embedding fit to look at? Checked once it exists, recorded in
#' run_info.json and shown by the app above the plots, so a failed map is said
#' where it is seen rather than in a log nobody reads (J445, 2026-11-02: a map
#' 550 units wide went unnoticed until the plots were opened).
#'
#' @param coords Two-column matrix of coordinates.
#' @param max_span Central 98% of cells wider than this (units) = the
#'   optimisation did not converge; a UMAP with spread 1 spans a few tens.
#' @param max_stretch Full range over the 1st-99th percentile range above this
#'   = a few far-off cells stretch the axes and squeeze the map.
#' @return list(suspect, span, core_span, stretch, outliers, reasons).
embedding_check <- function(coords, max_span = 100, max_stretch = 4) {
  coords <- as.matrix(coords)[, 1:2, drop = FALSE]
  ok <- stats::complete.cases(coords)
  coords <- coords[ok, , drop = FALSE]
  if (!nrow(coords)) {
    return(list(suspect = TRUE, span = NA, core_span = NA, stretch = NA,
                outliers = 0L, reasons = I("no coordinates")))
  }
  span <- max(apply(coords, 2, function(v) diff(range(v))))
  lo <- apply(coords, 2, stats::quantile, 0.01, names = FALSE)
  hi <- apply(coords, 2, stats::quantile, 0.99, names = FALSE)
  core <- max(hi - lo)
  stretch <- if (core > 0) span / core else NA_real_
  # Cells further than one core width outside the 1-99% box.
  pad <- core
  far <- sum(coords[, 1] < lo[1] - pad | coords[, 1] > hi[1] + pad |
             coords[, 2] < lo[2] - pad | coords[, 2] > hi[2] + pad)
  reasons <- character(0)
  # Judged on the central 98% of cells, so a handful of far-off cells is
  # reported as what it is (below), not as a failed optimisation.
  if (core > max_span) {
    reasons <- c(reasons, sprintf(paste0(
      "The map is %d units wide, where a normal UMAP spans a few tens: the ",
      "optimisation did not converge, and the plots will look like scattered ",
      "specks."), round(core)))
  } else if (is.finite(stretch) && stretch > max_stretch && far > 0) {
    reasons <- c(reasons, sprintf(paste0(
      "%d cell(s) lie far outside the map and stretch the axes %.1f-fold, so ",
      "the map itself is drawn small. Clipping extreme values ",
      "(--clip_quantile 0.001; in the app 'Clip extreme values') usually ",
      "brings them in."), far, stretch))
  }
  list(suspect = length(reasons) > 0, span = round(span, 1),
       core_span = round(core, 1), stretch = round(stretch, 2),
       outliers = far, reasons = I(reasons))
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

# A spectral start on a large fit set: said as soon as it is known, as a
# message (Rscript holds warning() back to the end of the run, an hour too
# late) and as a warning (so a caller can catch it).
umap_init_warning <- function(n_fit, limit = 300000L) {
  text <- paste0(
    "fitting UMAP on ", format(n_fit, big.mark = ","), " cells with a ",
    "spectral start. Above ", format(limit, big.mark = ","), " cells that start ",
    "is single-threaded and slow (500,000 cells: 34 min), and it may fail to ",
    "converge, in which case uwot falls back to a random start. Either fit on ",
    "fewer cells and project the rest (--umap_subsample, e.g. ",
    format(limit, big.mark = ""), "; in the app 'Random cells to fit on'), or ",
    "choose a PCA start yourself (--umap_init pca).")
  message("WARNING: ", text)
  warning(text, call. = FALSE)
  invisible(text)
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
