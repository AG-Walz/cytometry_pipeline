# scripts/drift/common.R — shared bootstrap for the drift CLIs
# Sourced by every NN_*.R in this folder. Locates the repo, sources the drift
# modules, parses `--key value` arguments without optparse (so steps 01/02 run
# on a plain R install), and reads the cohort config given with --config.

drift_repo_root <- function() {
  env <- Sys.getenv("CYTOMETRY_REPO", "")
  if (nzchar(env) && dir.exists(file.path(env, "src"))) return(normalizePath(env))
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
  cands <- c(if (length(f)) normalizePath(file.path(dirname(f), c("..", "../.."))), getwd())
  for (c in cands) if (dir.exists(file.path(c, "src", "harmonization", "drift"))) return(c)
  stop("cannot locate the repository root; set CYTOMETRY_REPO", call. = FALSE)
}

REPO <- drift_repo_root()
source(file.path(REPO, "src", "R_shared", "fcs_keywords.R"))
for (m in list.files(file.path(REPO, "src", "harmonization", "drift"), "\\.R$", full.names = TRUE)) source(m)

#' Parse `--name value` and `--flag` arguments. `spec` = list(name = default);
#' a logical default marks a flag. Unknown arguments are an error.
parse_args <- function(spec, args = commandArgs(trailingOnly = TRUE), usage = "") {
  out <- spec
  i <- 1L
  while (i <= length(args)) {
    a <- args[i]
    if (a %in% c("-h", "--help")) { cat(usage, "\n"); quit(status = 0) }
    if (!startsWith(a, "--")) stop("unexpected argument: ", a, "\n", usage, call. = FALSE)
    key <- substring(a, 3)
    if (!key %in% names(spec)) stop("unknown option --", key, "\n", usage, call. = FALSE)
    if (is.logical(spec[[key]])) { out[[key]] <- TRUE; i <- i + 1L; next }
    if (i == length(args)) stop("--", key, " needs a value", call. = FALSE)
    val <- args[i + 1L]
    out[[key]] <- if (is.numeric(spec[[key]])) as.numeric(val) else val
    i <- i + 2L
  }
  for (k in names(spec)) if (identical(spec[[k]], "") && !nzchar(out[[k]]))
    stop("--", k, " is required\n", usage, call. = FALSE)
  out
}

read_drift_config <- function(path) {
  if (!requireNamespace("yaml", quietly = TRUE)) stop("package 'yaml' is required", call. = FALSE)
  if (!file.exists(path)) stop("config file not found: ", path, call. = FALSE)
  cfg <- yaml::read_yaml(path)
  for (k in c("tubes", "instrument_labels", "markers", "reference_instrument"))
    if (is.null(cfg[[k]])) stop("config lacks the `", k, "` block; see config/drift_mtb.yaml", call. = FALSE)
  # configs written before the design block get the MTB layout without protected covariates
  if (is.null(cfg$design)) cfg$design <- default_design()
  else cfg$design <- modifyList(default_design(), cfg$design)
  cfg
}

#' Every column named in cfg$design must exist in tubes.csv (`tubes`).
check_design <- function(cfg, tubes) {
  miss <- setdiff(design_columns(cfg), names(tubes))
  if (length(miss)) stop("design columns missing from tubes.csv: ", paste(miss, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

#' File names excluded by step 3, or none when it has not run.
excluded_tubes <- function(out) {
  f <- file.path(out, "excluded_tubes.csv")
  if (file.exists(f)) read.csv(f)$filename else character(0)
}

#' Every FlowJo workspace under `wsp_dir`.
list_workspaces <- function(wsp_dir) {
  wsps <- list.files(wsp_dir, "\\.wsp$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  if (!length(wsps)) stop("no .wsp under ", wsp_dir, call. = FALSE)
  wsps
}

ensure_dir <- function(d) { dir.create(d, recursive = TRUE, showWarnings = FALSE); d }

write_table <- function(df, path) {
  ensure_dir(dirname(path))
  write.csv(df, path, row.names = FALSE)
  cat("wrote", path, "(", nrow(df), "rows )\n")
}

need_packages <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("missing R packages: ", paste(miss, collapse = ", "),
                         ". Run scripts/drift/00_setup_env.sh", call. = FALSE)
}
