#!/usr/bin/env Rscript
# =============================================================================
# R_shared/paths.R — CWD-independent resolution of the src/ tree
# =============================================================================
# Every src module used to locate its shared helpers with
#     file.path(dirname(getwd()), "src", "R_shared")
# which only resolves when the working directory happens to be an immediate
# sub-directory of the repo. Under Nextflow the working directory is
# work/xx/<hash>/, so that path does not exist and sourcing fails.
#
# This file replaces that fragile idiom. Each module loads it once via the
# standard bootstrap header (see any converted module or docs) and then uses:
#     cyto_src_root()    -> absolute path to <repo>/src
#     cyto_repo_root()   -> absolute path to <repo>
#     cyto_shared_dir()  -> <repo>/src/R_shared
#     source_shared(...) -> source one or more files from R_shared
#     source_src(...)    -> source files given relative to src/ ("sub/file.R")
#
# Resolution order (first that works wins), computed once and cached:
#   1. $CYTOMETRY_SRC (set by the CLI wrappers and the Nextflow env scope)
#   2. this file's own on-disk location (paths.R lives in <src>/R_shared)
#   3. walk up from the current working directory looking for an R_shared/ dir
#   4. legacy dirname(getwd())/src fallback
# =============================================================================

# ---- capture this file's location at *source* time (most reliable signal) ----
.cyto_paths_self <- local({
  of <- NULL
  # a) when source()d, the source frame carries $ofile
  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    f <- frames[[i]]$ofile
    if (!is.null(f) && nzchar(f)) { of <- f; break }
  }
  # b) when run via `Rscript paths.R`, recover from --file=
  if (is.null(of)) {
    a <- commandArgs(FALSE)
    m <- grep("^--file=", a)
    if (length(m)) of <- sub("^--file=", "", a[m[[1]]])
  }
  if (!is.null(of)) normalizePath(of, mustWork = FALSE) else NA_character_
})

# ---- resolve and cache the src/ root ----
.cyto_src_root_cached <- NULL

cyto_src_root <- function() {
  if (!is.null(.cyto_src_root_cached)) return(.cyto_src_root_cached)

  has_shared <- function(d) !is.na(d) && nzchar(d) && dir.exists(file.path(d, "R_shared"))

  root <- NA_character_

  # 1. explicit env var
  env <- Sys.getenv("CYTOMETRY_SRC", "")
  if (nzchar(env) && has_shared(env)) root <- env

  # 2. this file lives in <src>/R_shared -> src root is its grandparent dir
  if (is.na(root) && !is.na(.cyto_paths_self)) {
    cand <- dirname(dirname(.cyto_paths_self))
    if (has_shared(cand)) root <- cand
  }

  # 3. walk up from the current working directory for a dir containing R_shared
  if (is.na(root)) {
    d <- normalizePath(getwd(), mustWork = FALSE)
    while (!has_shared(d) && d != dirname(d)) d <- dirname(d)
    if (has_shared(d)) root <- d
  }

  # 4. legacy fallback (kept so nothing is worse than before)
  if (is.na(root)) root <- file.path(dirname(getwd()), "src")

  root <- normalizePath(root, mustWork = FALSE)
  .cyto_src_root_cached <<- root
  root
}

cyto_repo_root  <- function() dirname(cyto_src_root())
cyto_shared_dir <- function() file.path(cyto_src_root(), "R_shared")

#' Source one or more helper files from src/R_shared
source_shared <- function(...) {
  for (f in c(...)) source(file.path(cyto_shared_dir(), f))
  invisible(TRUE)
}

#' Source one or more files given relative to src/ (e.g. "harmonization/x.R")
source_src <- function(...) {
  for (p in c(...)) source(file.path(cyto_src_root(), p))
  invisible(TRUE)
}

# marker so the bootstrap header can detect a prior load cheaply
.cyto_paths_loaded <- TRUE
