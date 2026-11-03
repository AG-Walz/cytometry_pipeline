# =============================================================================
# test_data/test_run_record.R — the settings record every UMAP run writes
# =============================================================================
# The record exists so a run can be repeated and debugged later. What it must
# guarantee, and what these tests pin:
#   * every option is in the command, defaults included, so a later change
#     of a default does not change what an old command repeats;
#   * the command parses back to the same options (negative numbers, empty
#     strings, flags, values with spaces and quotes);
#   * it is written before the run, kept per run, and records the outcome.
#
# Run from the project root:
#   Rscript test_data/test_run_record.R
# =============================================================================

cat("=== Run settings record tests ===\n\n")

.ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(.ofile)) .ofile <- "test_data/test_run_record.R"
base_dir <- normalizePath(file.path(dirname(.ofile), ".."))
if (!dir.exists(file.path(base_dir, "src"))) base_dir <- getwd()

suppressPackageStartupMessages(library(optparse))
source(file.path(base_dir, "src", "R_shared", "run_record.R"))

# --- tiny test harness -------------------------------------------------------
.pass <- 0L
.fail <- 0L
ok <- function(label, cond) {
  if (isTRUE(cond)) {
    .pass <<- .pass + 1L
    cat("  PASS  ", label, "\n", sep = "")
  } else {
    .fail <<- .fail + 1L
    cat("  FAIL  ", label, "\n", sep = "")
  }
}

option_list <- list(
  make_option("--subsample", type = "integer", default = -1L),
  make_option("--name", type = "character", default = "markers"),
  make_option("--empty", type = "character", default = ""),
  make_option("--map", type = "character", default = NULL),
  make_option("--size", type = "double", default = 0.4),
  make_option("--no_resume", action = "store_true", default = FALSE),
  make_option("--output_dir", type = "character", default = "out"))

given <- c("--map=Asuwertung J445 m,TAA=/A/B C;It's=/D", "--no_resume",
           "--size=0.1")
opt <- parse_args(OptionParser(option_list = option_list), args = given)
opts <- run_record_options(opt, option_list, given)
cmd <- run_record_command("script.R", opts)

ok("every option is in the command, defaults included",
   all(vapply(c("--subsample=", "--name=", "--empty=", "--map=", "--size=",
                "--no_resume", "--output_dir="),
              function(f) grepl(f, cmd, fixed = TRUE), logical(1))))
ok("options typed by the user are marked as given",
   setequal(vapply(Filter(function(o) o$given, opts), function(o) o$dest, ""),
            c("map", "no_resume", "size")))

# Round trip through the shell: the command's arguments, as bash splits them,
# parse back to exactly the same options.
args_line <- sub("^'Rscript' \\\\\n  'script.R' \\\\\n  ", "", cmd)
script <- tempfile(fileext = ".sh")
writeLines(c("for a in \\", paste0("  ", args_line, "; do printf '%s\\n' \"$a\"; done")),
           script)
argv_back <- system2("bash", script, stdout = TRUE)
back <- parse_args(OptionParser(option_list = option_list), args = argv_back)
ok("the command parses back to the same options",
   identical(back[names(opt)], opt[names(opt)]))
ok("...including a negative default (-1), not read as a flag",
   identical(back$subsample, -1L))
ok("...an empty string, and a value with spaces, ';' and a quote",
   identical(back$empty, "") &&
     identical(back$map, "Asuwertung J445 m,TAA=/A/B C;It's=/D"))

# A flag left FALSE is left out, and a NULL option is left out.
opt2 <- parse_args(OptionParser(option_list = option_list), args = character(0))
cmd2 <- run_record_command("script.R", run_record_options(opt2, option_list))
ok("a FALSE flag and a NULL option are left out",
   !grepl("--no_resume", cmd2, fixed = TRUE) && !grepl("--map", cmd2, fixed = TRUE))

# --- the record on disk ------------------------------------------------------
out <- tempfile("run_")
inputs_dir <- tempfile("fcs_")
dir.create(inputs_dir)
writeLines("x", file.path(inputs_dir, "a.fcs"))
wsp <- tempfile(fileext = ".wsp")
writeLines("<Workspace/>", wsp)

rec <- suppressMessages(write_run_record(
  out, "script.R", opt, option_list,
  inputs = list(wsp = wsp, fcs_dir = inputs_dir, missing = "/no/such/path"),
  repo = base_dir, packages = c("jsonlite", "not_a_package")))
ok("the record is written before the run, with status 'running'",
   file.exists(file.path(out, "_run", "settings.json")) &&
     identical(jsonlite::fromJSON(file.path(out, "_run", "settings.json"))$status,
               "running"))
ok("command.sh sits beside it and in the dated folder",
   file.exists(file.path(out, "_run", "command.sh")) &&
     file.exists(file.path(rec$dir, "command.sh")))
s <- jsonlite::fromJSON(file.path(rec$dir, "settings.json"))
ok("inputs are fingerprinted: the workspace hashed, the folder listed",
   nchar(s$inputs$wsp$md5) == 32L && identical(s$inputs$fcs_dir$files, 1L) &&
     file.exists(file.path(rec$dir, "fcs_dir_files.tsv")) &&
     is.null(s$inputs$missing))
ok("package versions are recorded, a missing package as null",
   is.character(s$packages$jsonlite) && is.null(s$packages$not_a_package))
ok("the code version is recorded", nchar(s$git$commit) == 40L)

finish_run_record(rec, "failed", simpleError("population not found"))
s <- jsonlite::fromJSON(file.path(out, "_run", "settings.json"))
ok("the outcome and the error are recorded at the end",
   identical(s$status, "failed") && identical(s$error, "population not found") &&
     is.numeric(s$minutes))

rec2 <- suppressMessages(write_run_record(out, "script.R", opt, option_list,
                                          repo = base_dir))
ok("a second run into the same folder keeps the first run's record",
   length(list.dirs(file.path(out, "_run", "settings"), recursive = FALSE)) == 2L)
ok("...and the latest copy is the second run's",
   identical(jsonlite::fromJSON(file.path(out, "_run", "settings.json"))$status,
             "running"))

ok("a record that cannot be written does not stop the run",
   is.null(suppressMessages(write_run_record("/proc/no/such/dir", "script.R", opt,
                                             option_list))))

cat("\n=== ", .pass, " passed, ", .fail, " failed ===\n", sep = "")
if (.fail > 0) quit(status = 1)
