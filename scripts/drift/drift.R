#!/usr/bin/env Rscript
# scripts/drift/drift.R: one entry point for the drift analysis (docs/plans/drift_usability_plan.md, section 4).
#   Rscript scripts/drift/drift.R <init|doctor|check|run|status> [options]      (--help on every sub-command)
# The step scripts in this folder are unchanged; this wrapper decides which of them need to run
# (an output is missing or older than an input), runs them with the same Rscript that runs this
# file, logs to <out>/run_all.log and writes <out>/.done/<step>.json after every step.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
HERE <- file.path(REPO, "scripts", "drift")
RSCRIPT <- file.path(R.home("bin"), "Rscript")

top_usage <- "drift.R <sub-command> [options]

  init    --fcs_dir DIR --out DIR [--config CFG] [--cohort NAME]      pre-fill sample_sheet.csv and a cohort config
  doctor  [--bf_env NAME]                                             check R packages, CmdStan and the BayesFlow env
  check   --out DIR --config CFG [--n_rep 10] [--wsp_dir DIR --fcs_dir DIR]   separability of week and acquisition day
  run     --out DIR --config CFG [--sample_sheet CSV] [--fcs_dir DIR] [--wsp_dir DIR] ...   run the steps that are stale
  status  --out DIR                                                   what exists, what is stale, when it last ran

Run `drift.R <sub-command> --help` for the options of one sub-command."

# ---------------------------------------------------------------- step table
# needs/makes are file names under --out (globs allowed in needs); the config is an input of every
# step. optional = TRUE: a missing script is reported as skipped instead of stopping the run.
STEPS <- list(
  list(id = "1",  name = "metadata",        script = "01_extract_metadata.R", needs = character(0),
       makes = c("tubes.csv", "bridge_report.csv"), dirs = "fcs_dir"),
  list(id = "1b", name = "check_design",    script = "check_design.R", needs = "tubes.csv",
       makes = c("separability.csv", "separability_summary.csv"), optional = TRUE),
  list(id = "2",  name = "panel",           script = "02_panel_reconcile.R", needs = "tubes.csv",
       makes = c("panel_map.csv", "panel_signatures.csv", "panel_gaps.csv")),
  list(id = "3",  name = "tube QC",         script = "03_tube_qc.R", needs = c("tubes.csv", "panel_map.csv"),
       makes = c("qc.csv", "excluded_tubes.csv")),
  list(id = "4",  name = "neg summaries",   script = "04_neg_summaries.R",
       needs = c("tubes.csv", "panel_map.csv", "panel_signatures.csv", "excluded_tubes.csv", "excluded_manual.csv"),
       makes = c("neg_summaries.csv", "neg_emd.csv", "gate_bounds.csv"), dirs = c("wsp_dir", "fcs_dir")),
  list(id = "5",  name = "frequencies",     script = "05_gated_frequencies.R",
       needs = c("tubes.csv", "panel_signatures.csv", "excluded_tubes.csv", "excluded_manual.csv"),
       makes = c("gated_frequencies.csv", "all_gate_counts.csv"), dirs = c("wsp_dir", "fcs_dir")),
  list(id = "6a", name = "brms",            script = "06_yardstick_brms.R",
       needs = c("tubes.csv", "neg_summaries.csv", "neg_emd.csv", "gated_frequencies.csv"),
       makes = c("effects_brms_A.csv", "effects_brms_BG.csv", "effects_brms_B.csv", "covariates_brms.csv", "fits_brms_summary.csv")),
  list(id = "6b", name = "stan",            script = "06_yardstick_stan.R",
       needs = c("tubes.csv", "neg_summaries.csv", "neg_emd.csv", "gated_frequencies.csv"),
       makes = c("effects_stan_A.csv", "effects_stan_BG.csv", "effects_stan_B.csv", "covariates_stan.csv")),
  list(id = "6c", name = "lme4",            script = "06_lme4.R",
       needs = c("tubes.csv", "neg_summaries.csv", "gated_frequencies.csv"),
       makes = c("lme4_variance_components.csv", "assumption_check.csv")),
  list(id = "6d", name = "biology",         script = "06_biology.R", needs = c("tubes.csv", "gated_frequencies.csv"),
       makes = c("week_slope_sensitivity.csv", "week_slope_lodo.csv"), optional = TRUE),
  list(id = "7",  name = "bayesflow v1",    script = "07_bayesflow_v1.py", python = TRUE,
       needs = c("tubes.csv", "neg_summaries.csv", "neg_emd.csv"),
       makes = c("effects_bayesflow_v1.csv", "bf_v1_recovery.csv", "bf_v1_sbc.csv", "bf_v1_timing.json")),
  list(id = "8",  name = "bayesflow v2",    script = "08_bayesflow_v2.py", python = TRUE,
       needs = c("tubes.csv", "neg_summaries.csv", "neg_cells/*.csv.gz"),
       makes = c("effects_bayesflow_v2.csv", "bf_v2_recovery.csv", "bf_v2_sbc.csv", "bf_v2_ppc.csv", "bf_v2_timing.json")),
  list(id = "9",  name = "decision",        script = "09_decision.R",
       needs = c("tubes.csv", "gated_frequencies.csv", "neg_summaries.csv", "effects_*.csv"),
       makes = c("decisions.csv", "controls_summary.csv")),
  list(id = "10", name = "report",          script = "10_report.Rmd",
       needs = c("tubes.csv", "decisions.csv", "effects_*.csv", "lme4_variance_components.csv"),
       makes = "report.html")
)
STEP_IDS <- vapply(STEPS, `[[`, "", "id")
BF_STEPS <- c("7", "8")

step_by_id <- function(id) STEPS[[match(id, STEP_IDS)]]

#' Resolve a user-typed step id: "6" means "6a", "6a" stays, anything else must be a known id.
resolve_id <- function(x) {
  x <- trimws(x)
  if (x %in% STEP_IDS) return(x)
  if (paste0(x, "a") %in% STEP_IDS) return(paste0(x, "a"))
  stop("unknown step '", x, "'; known steps: ", paste(STEP_IDS, collapse = ", "), call. = FALSE)
}

input_files <- function(step, out, extra = character(0)) {
  f <- unlist(lapply(step$needs, function(p) if (grepl("[*?]", p)) Sys.glob(file.path(out, p)) else file.path(out, p)))
  c(f, extra)
}

#' "missing", "stale" or "ok" from the outputs and the existing inputs of a step.
step_state <- function(step, out, extra_inputs = character(0)) {
  outs <- file.path(out, step$makes)
  if (!all(file.exists(outs))) return("missing")
  ins <- input_files(step, out, extra_inputs); ins <- ins[file.exists(ins)]
  if (length(ins) && min(file.mtime(outs)) < max(file.mtime(ins))) return("stale")
  "ok"
}

# ---------------------------------------------------------------- .done markers (no jsonlite needed)
json_str <- function(x) paste0('"', gsub('"', '\\\\"', gsub("\\\\", "\\\\\\\\", x)), '"')
write_done <- function(out, step, started, cmd, config) {
  d <- ensure_dir(file.path(out, ".done"))
  commit <- tryCatch(suppressWarnings(system2("git", c("-C", shQuote(REPO), "rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)),
                     error = function(e) character(0))
  commit <- if (length(commit) == 1 && nzchar(commit)) commit else NA
  fields <- c(step = json_str(step$id), name = json_str(step$name),
              time = json_str(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
              duration_s = sprintf("%.1f", as.numeric(difftime(Sys.time(), started, units = "secs"))),
              git_commit = if (is.na(commit)) "null" else json_str(commit),
              config_md5 = json_str(unname(tools::md5sum(config))),
              command = json_str(cmd))
  writeLines(c("{", paste0('  "', names(fields), '": ', fields, c(rep(",", length(fields) - 1), "")), "}"),
             file.path(d, paste0(step$id, ".json")))
}
read_done_time <- function(out, id) {
  f <- file.path(out, ".done", paste0(id, ".json"))
  if (!file.exists(f)) return("")
  m <- regmatches(readLines(f), regexpr('"time": "[^"]*"', readLines(f)))
  if (length(m)) sub('"time": "([^"]*)"', "\\1", m[1]) else ""
}

# ---------------------------------------------------------------- running commands
LOG <- NULL
say <- function(...) { msg <- paste0(...); cat(msg, "\n", sep = ""); if (!is.null(LOG)) cat(msg, "\n", sep = "", file = LOG, append = TRUE) }

#' Run a shell command, echoing its output to the terminal and appending it to the log (like tee).
run_cmd <- function(cmd) {
  full <- if (is.null(LOG)) cmd else paste(cmd, "2>&1 | tee -a", shQuote(LOG))
  system2("bash", c("-o", "pipefail", "-c", shQuote(full)))
}

#' conda, mamba or micromamba, whichever is on the PATH (or the MAMBA_EXE a micromamba shell exports).
conda_bin <- function() {
  for (b in c("conda", "micromamba", "mamba")) if (nzchar(Sys.which(b))) return(b)
  if (file.exists(Sys.getenv("MAMBA_EXE"))) return(Sys.getenv("MAMBA_EXE"))
  "conda"
}

cohort_name <- function(config) {
  l <- grep("^cohort:", readLines(config, warn = FALSE), value = TRUE)
  if (!length(l)) return("cohort")
  trimws(gsub('^cohort:\\s*"?([^"#]*?)"?\\s*(#.*)?$', "\\1", l[1]))
}

#' The shell command that runs one step.
step_command <- function(step, a, out, config) {
  s <- file.path(HERE, step$script)
  q <- shQuote
  common <- paste("--out", q(out), "--config", q(config))
  if (isTRUE(step$python)) {
    py <- paste(conda_bin(), "run --no-capture-output -n", q(a$bf_env), "python")
    extra <- if (step$id == "7") "--epochs 30 --num_sims 30000 --batch_size 128"
             else "--epochs 15 --batches_per_epoch 100 --batch_size 8 --n_cells 800 --chunk 4 --n_sbc 100 --posterior_samples 1000"
    return(paste(py, q(s), common, extra))
  }
  if (step$id == "10") {
    expr <- sprintf("rmarkdown::render(%s, params = list(out = %s, cohort = %s, config = %s), output_file = 'report.html', output_dir = %s)",
                    q(s), q(out), q(cohort_name(config)), q(config), q(out))
    return(paste(q(RSCRIPT), "-e", q(expr)))
  }
  args <- switch(step$id,
    "1"  = paste("--fcs_dir", q(a$fcs_dir), if (a$sample_sheet != "none") paste("--sample_sheet", q(a$sample_sheet))),
    "1b" = if (a$wsp_dir != "none" && a$fcs_dir != "none") paste("--wsp_dir", q(a$wsp_dir), "--fcs_dir", q(a$fcs_dir)) else "",
    "3"  = paste("--workers", a$workers),
    "4"  = paste("--wsp_dir", q(a$wsp_dir), "--fcs_dir", q(a$fcs_dir), "--all_tubes"),
    "5"  = paste("--wsp_dir", q(a$wsp_dir), "--fcs_dir", q(a$fcs_dir)),
    "6a" = paste("--cores", a$cores), "6b" = paste("--cores", a$cores),
    "")
  paste(q(RSCRIPT), q(s), common, args)
}

# ---------------------------------------------------------------- sub-commands
forward <- function(script, usage, args) {
  if (any(args %in% c("-h", "--help"))) { cat(usage, "\n"); quit(status = 0) }
  s <- file.path(HERE, script)
  if (!file.exists(s)) stop(script, " is not in ", HERE, " yet", call. = FALSE)
  quit(status = system2(RSCRIPT, c(shQuote(s), shQuote(args))))
}

cmd_init <- function(args) forward("00_init_sheet.R",
  "drift.R init --fcs_dir DIR --out DIR [--config CFG] [--cohort NAME]\n  writes <out>/sample_sheet.csv and <out>/drift_<cohort>.yaml (forwards to 00_init_sheet.R)", args)

cmd_check <- function(args) forward("check_design.R",
  "drift.R check --out DIR --config CFG [--n_rep 10] [--wsp_dir DIR --fcs_dir DIR]\n  writes separability.csv and separability_summary.csv and prints a verdict line (forwards to check_design.R)", args)

cmd_doctor <- function(args) {
  usage <- "drift.R doctor [--bf_env NAME]\n  checks the R packages, CmdStan and (with --bf_env) the BayesFlow conda env; hints point to 00_setup_env.sh"
  a <- parse_args(list(bf_env = "none"), args, usage)
  setup <- "bash scripts/drift/00_setup_env.sh (or conda env create -f environment_drift.yml)"
  rows <- list()
  add <- function(what, ok, detail = "", hint = "") rows[[length(rows) + 1L]] <<- data.frame(
    component = what, status = if (is.na(ok)) "skipped" else if (ok) "OK" else "missing", detail = detail,
    hint = if (isTRUE(ok) || is.na(ok)) "" else hint, stringsAsFactors = FALSE)
  add("Rscript", TRUE, paste(RSCRIPT, R.version$major, R.version$minor, sep = " "))
  needed <- c(yaml = "steps 2 to 10", lme4 = "steps 1b, 6c, 6d", brms = "step 6a", cmdstanr = "step 6b",
              CytoML = "steps 4, 5", flowWorkspace = "steps 4, 5", flowCore = "steps 3 to 5",
              rmarkdown = "step 10", ggplot2 = "step 10", knitr = "step 10")
  for (p in names(needed)) {
    ok <- requireNamespace(p, quietly = TRUE)
    add(paste("R package", p), ok, if (ok) as.character(packageVersion(p)) else needed[[p]], setup)
  }
  if (requireNamespace("cmdstanr", quietly = TRUE)) {
    v <- tryCatch(cmdstanr::cmdstan_version(), error = function(e) NULL)
    add("CmdStan", !is.null(v), if (is.null(v)) "step 6b" else v,
        "Rscript -e 'cmdstanr::install_cmdstan(cores = 8)' or set CMDSTAN to the conda cmdstan")
  } else add("CmdStan", FALSE, "cmdstanr not installed", setup)
  if (a$bf_env != "none") {
    cb <- conda_bin()
    if (nzchar(Sys.which(cb))) {
      v <- suppressWarnings(system2(cb, c("run", "-n", shQuote(a$bf_env), "python", "-c",
            shQuote("import bayesflow as bf, keras; print(bf.__version__, keras.backend.backend())")), stdout = TRUE, stderr = FALSE))
      ok <- is.null(attr(v, "status")) && length(v) >= 1
      add(paste0("BayesFlow env ", a$bf_env), ok, if (ok) paste("bayesflow", tail(v, 1)) else "steps 7, 8 (opt-in)", setup)
    } else add("conda", FALSE, "needed for `conda run -n <bf_env> python`", "install conda or micromamba")
  } else add("BayesFlow env", NA, "not checked (pass --bf_env NAME); steps 7 and 8 are opt-in")
  tab <- do.call(rbind, rows)
  cat(sprintf("%-26s %-8s %s\n", "component", "status", "detail"))
  for (i in seq_len(nrow(tab))) cat(sprintf("%-26s %-8s %s\n", tab$component[i], tab$status[i], tab$detail[i]))
  n_miss <- sum(tab$status == "missing")
  if (n_miss) {
    cat(sprintf("\n%d component(s) missing. Install hints:\n", n_miss))
    for (h in unique(tab$hint[tab$status == "missing"])) cat("  ", h, "\n")
  } else cat("\nall components present\n")
  quit(status = as.integer(n_miss > 0))
}

run_usage <- "drift.R run --out DIR --config CFG [--sample_sheet CSV] [--fcs_dir DIR] [--wsp_dir DIR]
             [--cores 4] [--workers 8] [--from N] [--only 5,9,10] [--with_bayesflow] [--bf_env bayesflow_env] [--force]
  Runs every step whose outputs are missing or older than its inputs (or the config). --from N and --only
  restrict the set of steps considered, --force reruns the considered steps regardless. Steps 7 and 8
  (BayesFlow, `conda run -n <bf_env> python`) run only with --with_bayesflow. Steps: 1 1b 2 3 4 5 6a 6b 6c 6d 7 8 9 10.
  Output goes to the terminal and to <out>/run_all.log; <out>/.done/<step>.json records each finished step."

cmd_run <- function(args) {
  a <- parse_args(list(out = "", config = "", sample_sheet = "none", fcs_dir = "none", wsp_dir = "none",
                       cores = 4, workers = 8, from = "none", only = "none",
                       with_bayesflow = FALSE, bf_env = "bayesflow_env", force = FALSE), args, run_usage)
  if (!file.exists(a$config)) stop("config not found: ", a$config, call. = FALSE)
  out <- normalizePath(ensure_dir(a$out)); config <- normalizePath(a$config)
  for (k in c("sample_sheet", "fcs_dir", "wsp_dir")) if (a[[k]] != "none") a[[k]] <- normalizePath(a[[k]], mustWork = TRUE)
  LOG <<- file.path(out, "run_all.log")

  selected <- STEP_IDS
  if (a$from != "none") selected <- selected[seq(match(resolve_id(a$from), STEP_IDS), length(selected))]
  if (a$only != "none") selected <- intersect(selected, vapply(strsplit(a$only, ",")[[1]], resolve_id, ""))
  if (!a$with_bayesflow) selected <- setdiff(selected, BF_STEPS)

  say(sprintf("\n===== %s drift.R run: out=%s config=%s =====", format(Sys.time(), "%F %T"), out, config))
  say("steps considered: ", paste(selected, collapse = " "), if (a$force) " (forced)" else "")
  extra <- c(config, if (a$sample_sheet != "none") a$sample_sheet)
  n_run <- 0L
  for (id in selected) {
    st <- step_by_id(id)
    script <- file.path(HERE, st$script)
    if (!file.exists(script)) {
      if (isTRUE(st$optional)) { say(sprintf("step %-3s %-14s skipped: %s not present yet", id, st$name, st$script)); next }
      stop("step ", id, ": script not found: ", script, call. = FALSE)
    }
    state <- step_state(st, out, extra)
    if (state == "ok" && !a$force) { say(sprintf("step %-3s %-14s up to date", id, st$name)); next }
    for (d in st$dirs) if (a[[d]] == "none") stop("step ", id, " needs --", d, call. = FALSE)
    cmd <- step_command(st, a, out, config)
    say(sprintf("\n===== %s step %s %s (%s) =====", format(Sys.time(), "%F %T"), id, st$name, if (a$force) "forced" else state))
    if (!is.null(LOG)) cat(cmd, "\n", sep = "", file = LOG, append = TRUE)   # the command line goes to the log only
    started <- Sys.time()
    status <- run_cmd(cmd)
    if (status != 0) { say(sprintf("step %s %s FAILED (exit %d); see %s", id, st$name, status, LOG)); quit(status = 1) }
    missing_out <- st$makes[!file.exists(file.path(out, st$makes))]
    if (length(missing_out)) say("warning: step ", id, " finished without writing ", paste(missing_out, collapse = ", "))
    write_done(out, st, started, cmd, config)
    n_run <- n_run + 1L
    say(sprintf("step %s %s done in %.0f s", id, st$name, as.numeric(difftime(Sys.time(), started, units = "secs"))))
  }
  say(sprintf("\n===== %s done: %d step(s) run, %d up to date or skipped =====", format(Sys.time(), "%F %T"), n_run, length(selected) - n_run))
}

cmd_status <- function(args) {
  a <- parse_args(list(out = "", config = "none"), args,
    "drift.R status --out DIR [--config CFG]\n  one line per step: outputs present, stale or not, and the time of the last .done marker")
  out <- a$out
  if (!dir.exists(out)) stop("no such directory: ", out, call. = FALSE)
  extra <- if (a$config != "none") normalizePath(a$config) else character(0)
  cat(sprintf("%-4s %-14s %-8s %-9s %-24s %s\n", "step", "name", "outputs", "state", "last done", "outputs"))
  for (st in STEPS) {
    present <- file.exists(file.path(out, st$makes))
    state <- step_state(st, out, extra)
    script_ok <- file.exists(file.path(HERE, st$script))
    note <- if (!script_ok) " (script not present)" else if (st$id %in% BF_STEPS) " (opt-in)" else ""
    cat(sprintf("%-4s %-14s %-8s %-9s %-24s %s%s\n", st$id, st$name, sprintf("%d/%d", sum(present), length(present)),
                state, read_done_time(out, st$id), paste(st$makes, collapse = ","), note))
  }
}

# ---------------------------------------------------------------- dispatch
argv <- commandArgs(trailingOnly = TRUE)
if (!length(argv) || argv[1] %in% c("-h", "--help")) { cat(top_usage, "\n"); quit(status = 0) }
sub <- argv[1]; rest <- argv[-1]
switch(sub,
  init = cmd_init(rest), doctor = cmd_doctor(rest), check = cmd_check(rest),
  run = cmd_run(rest), status = cmd_status(rest),
  stop("unknown sub-command '", sub, "'\n", top_usage, call. = FALSE))
