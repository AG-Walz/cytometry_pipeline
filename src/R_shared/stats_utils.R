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
#' @return Matrix of UMAP coordinates
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
                     sgd_threads = 1L) {
  set.seed(seed)
  n_threads   <- max(1L, as.integer(n_threads))
  sgd_threads <- max(1L, as.integer(sgd_threads))
  n_obs <- nrow(expr_matrix)

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
  if (identical(init, "auto")) {
    init <- if (n_obs > SPECTRAL_LIMIT) "pca" else "spectral"
    if (init == "pca") {
      message("uwot: ", n_obs, " cells is past the ", SPECTRAL_LIMIT,
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
  }
  if (!is.null(cb)) base_args$epoch_callback <- cb

  subsample <- suppressWarnings(as.integer(subsample))
  if (is.na(subsample)) subsample <- 0L

  # --- the whole matrix, when it is small enough to be worth it ---------------
  if (subsample <= 0L || subsample >= n_obs) {
    return(with_uwot_progress(
      do.call(uwot::umap, c(list(X = expr_matrix), base_args))))
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
  message("uwot: fitting on a ", subsample, "-cell subsample of ", n_obs,
          ", then projecting the remaining ", n_obs - subsample, ".")
  idx <- sort(sample.int(n_obs, subsample))
  model <- with_uwot_progress(
    do.call(uwot::umap, c(list(X = expr_matrix[idx, , drop = FALSE],
                               ret_model = TRUE), base_args)))

  coords <- matrix(NA_real_, nrow = n_obs, ncol = n_components)
  coords[idx, ] <- model$embedding

  rest <- setdiff(seq_len(n_obs), idx)
  # Chunked so the transform reports progress and so a 7-million-row projection
  # does not have to be resident all at once.
  chunks <- split(rest, ceiling(seq_along(rest) / max(1L, as.integer(chunk))))
  for (i in seq_along(chunks)) {
    rows <- chunks[[i]]
    coords[rows, ] <- uwot::umap_transform(
      expr_matrix[rows, , drop = FALSE], model, n_threads = n_threads,
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
  coords
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
