#!/usr/bin/env Rscript
# =============================================================================
# test_data/test_harmonize_marker_names.R : the $PnS harmoniser
# =============================================================================
# Covers src/metadata_qc/harmonize_marker_names.R: census and summary, plan
# precedence (alias > canonical > spelling, transpositions only suggested,
# scatter and Time untouched, CD4 and CD40 kept apart), dry run, in-place
# apply with byte-identical events, idempotence, revert and copies.
#
# Needs neither flowCore nor optparse nor yaml: fixtures come from
# test_data/fcs_fixtures.R. Run from the project root:
#   Rscript test_data/test_harmonize_marker_names.R
# =============================================================================

cat("=== marker name harmoniser tests ===\n\n")

.ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(.ofile)) .ofile <- "test_data/test_harmonize_marker_names.R"
base_dir <- normalizePath(file.path(dirname(.ofile), ".."))
if (!dir.exists(file.path(base_dir, "src"))) base_dir <- getwd()

shared_dir <- file.path(base_dir, "src", "R_shared")
source(file.path(shared_dir, "fcs_keywords.R"))
source(file.path(base_dir, "src", "harmonization", "panel_homogenization.R"))
source(file.path(base_dir, "src", "harmonization", "drift", "drift_panel.R"))
source(file.path(base_dir, "src", "metadata_qc", "harmonize_marker_names.R"))
source(file.path(base_dir, "test_data", "fcs_fixtures.R"))

# --- tiny test harness -------------------------------------------------------
.pass <- 0L; .fail <- 0L
ok <- function(label, cond) {
  if (isTRUE(cond)) { .pass <<- .pass + 1L; cat("  PASS  ", label, "\n", sep = "") }
  else { .fail <<- .fail + 1L; cat("  FAIL  ", label, "\n", sep = "") }
}
eq <- function(label, actual, expected) {
  same <- isTRUE(all.equal(actual, expected))
  if (!same) {
    cat("  FAIL  ", label, "\n        expected: ", paste(format(expected), collapse = " | "),
        "\n        actual:   ", paste(format(actual), collapse = " | "), "\n", sep = "")
    .fail <<- .fail + 1L
  } else { .pass <<- .pass + 1L; cat("  PASS  ", label, "\n", sep = "") }
}
throws <- function(label, expr, pattern = NULL) {
  msg <- tryCatch({ force(expr); NULL }, error = function(e) conditionMessage(e))
  ok(label, !is.null(msg) && (is.null(pattern) || grepl(pattern, msg)))
}
quiet <- function(expr) {
  out <- textConnection(NULL, "w", local = TRUE)
  sink(out); on.exit({ sink(); close(out) }, add = TRUE)
  force(expr)
}

root <- file.path(tempdir(), paste0("markers_", Sys.getpid()))
unlink(root, recursive = TRUE); dir.create(root, recursive = TRUE)

# Six parameters: scatter, three fluorochromes, one extra, Time.
CHANNELS <- c("FSC-A", "PE-A", "Pacific Blue-A", "APC-Cy7-A", "FITC-A", "Time")
make_file <- function(dir, name, markers = list()) {
  kw <- list()
  for (i in seq_along(CHANNELS)) kw[[paste0("$P", i, "N")]] <- CHANNELS[i]
  for (i in names(markers)) kw[[paste0("$P", i, "S")]] <- markers[[i]]
  p <- file.path(dir, name)
  if (!dir.exists(dirname(p))) dir.create(dirname(p), recursive = TRUE)
  fcs_fixture(p, keywords = kw, par = length(CHANNELS), tot = 7L, slack = 0L)
  p
}
# A: INFy (transposed), TNFa, CD4.  B: IFN, TNF, CD4.  C: no $PnS at all.
# D: IFNg, TNFa, CD4, and CD40 on FITC-A, which must never merge with CD4.
make_cohort <- function(dir) {
  unlink(dir, recursive = TRUE)
  make_file(dir, "A.fcs", list("2" = "INFy", "3" = "TNFa", "4" = "CD4"))
  make_file(dir, "B.fcs", list("2" = "IFN", "3" = "TNF", "4" = "CD4"))
  make_file(dir, "C.fcs", list())
  make_file(dir, "sub/D.fcs", list("2" = "IFNg", "3" = "TNFa", "4" = "CD4", "5" = "CD40"))
  dir
}
files_of <- function(d) sort(list.files(d, "\\.fcs$", recursive = TRUE, full.names = TRUE))
data_bytes <- function(f) {
  i <- fcs_inspect(f); fx_bytes(f, i$seg$datastart, i$seg$dataend - i$seg$datastart + 1)
}
pns <- function(f, i) fcs_keyword(f, paste0("$P", i, "S"))

CANON <- list(FSC_A = "^FSC-A$", LiveDead = c("live", "dead", "^ALD$"),
              CD4 = c("^CD4$", "^APC-Cy7-A$"), IFNg = c("IFN", "INF", "^PE-A$"),
              TNFa = c("TNF", "^Pacific Blue-A$"))

# --- census -------------------------------------------------------------------
cat("\n-- census ----------------------------------------------------------\n")
d <- make_cohort(file.path(root, "census"))
files <- files_of(d)
md5 <- vapply(files, fx_md5, character(1))
cen <- marker_census(files, d)
eq("one row per file x parameter", nrow(cen), 4L * length(CHANNELS))
eq("filenames are relative to fcs_dir", sort(unique(cen$filename)),
   c("A.fcs", "B.fcs", "C.fcs", "sub/D.fcs"))
eq("channels come from $PnN", cen$channel[cen$filename == "A.fcs"], CHANNELS)
ok("absent $PnS is NA", all(is.na(cen$marker[cen$filename == "C.fcs"])))
eq("present $PnS is read", cen$marker[cen$filename == "A.fcs" & cen$index == 2L], "INFy")
ok("census writes nothing", identical(vapply(files, fx_md5, character(1)), md5))
cs <- census_summary(cen)
pe <- cs[cs$channel == "PE-A", ]
eq("summary lists every spelling of a channel", sort(pe$marker, na.last = TRUE),
   c("IFN", "IFNg", "INFy", NA))
eq("...with file counts", pe$n_files, rep(1L, 4))
eq("summary counts files, not rows", cs$n_files[cs$channel == "APC-Cy7-A" & cs$marker %in% "CD4"], 3L)
rep_txt <- format_marker_report(cen)
# C has none at all (6); scatter, Time and, except in D, FITC-A carry none (8).
ok("the report counts blank parameters", grepl("without \\$PnS: 14", rep_txt))

# --- plan: spelling rule alone ------------------------------------------------
cat("\n-- plan: spelling only ---------------------------------------------\n")
pl <- plan_marker_names(cen)
p <- pl$plan
row <- function(f, i) p[p$filename == f & p$index == i, ]
eq("plan columns", names(p), c("filename", "index", "channel", "marker_old", "marker_new", "source"))
eq("IFN -> IFNg is a spelling fix", row("B.fcs", 2)$marker_new, "IFNg")
eq("...flagged as such", row("B.fcs", 2)$source, "spelling")
eq("TNF -> TNFa likewise", row("B.fcs", 3)$marker_new, "TNFa")
eq("INFy is a transposition and is NOT applied", row("A.fcs", 2)$marker_new, "INFy")
eq("...its row is unchanged", row("A.fcs", 2)$source, "unchanged")
ok("...but it is suggested", nrow(pl$suggestions) == 1L &&
     pl$suggestions$confidence == "transposition" && grepl("INFy", pl$suggestions$names))
ok("no spelling group is reported as a suggestion", !any(pl$suggestions$confidence == "spelling"))
eq("CD4 is left alone", unique(p$marker_new[p$channel == "APC-Cy7-A" & !is.na(p$marker_old)]), "CD4")
eq("CD40 does not merge with CD4", row("sub/D.fcs", 5)$marker_new, "CD40")
ok("scatter and Time are untouched",
   all(p$source[p$channel %in% c("FSC-A", "Time")] == "unchanged"))
ok("the skip rule covers prefixed scatter and Time, not fluorochromes",
   all(grepl(MARKER_SKIP_CHANNEL_REGEX, c("FSC-A", "SSC-H", "SP SSC-A", "Time"), ignore.case = TRUE)) &&
     !any(grepl(MARKER_SKIP_CHANNEL_REGEX, c("PE-A", "BV480-A", "Aqua LiveDead-A"), ignore.case = TRUE)))
ok("blank $PnS stays blank without a canonical map", all(is.na(p$marker_new[p$filename == "C.fcs"])))
sp <- summarise_marker_plan(p)
eq("the grouped summary has one row per rename", nrow(sp), 2L)
eq("...with file counts", sp$n_files, c(1L, 1L))

# --- plan: canonical map ------------------------------------------------------
cat("\n-- plan: canonical map ---------------------------------------------\n")
pl <- plan_marker_names(cen, canonical = CANON)
p <- pl$plan
eq("blank $PnS is filled from the channel", row("C.fcs", 2)$marker_new, "IFNg")
eq("...as canonical", row("C.fcs", 2)$source, "canonical")
eq("INFy IS renamed once a canonical pattern says so", row("A.fcs", 2)$marker_new, "IFNg")
eq("IFN too", row("B.fcs", 2)$marker_new, "IFNg")
eq("TNF -> TNFa via canonical, not spelling", row("B.fcs", 3)$source, "canonical")
eq("a name already canonical is unchanged", row("sub/D.fcs", 2)$source, "unchanged")
ok("FSC-A is not renamed to FSC_A even though the map has a pattern for it",
   is.na(row("A.fcs", 1)$marker_new) && row("A.fcs", 1)$source == "unchanged")
ok("Time keeps no $PnS", is.na(row("A.fcs", 6)$marker_new))
eq("CD40 on an unmapped channel is left alone", row("sub/D.fcs", 5)$marker_new, "CD40")
eq("nothing is left to suggest", nrow(pl$suggestions), 0L)

# --- plan: aliases win over everything ----------------------------------------
cat("\n-- plan: alias precedence ------------------------------------------\n")
pl <- plan_marker_names(cen, canonical = CANON, aliases = c(INFy = "IFN-gamma"))
p <- pl$plan
eq("an explicit alias beats the canonical map", row("A.fcs", 2)$marker_new, "IFN-gamma")
eq("...and is labelled alias", row("A.fcs", 2)$source, "alias")
eq("other rows still follow the canonical map", row("B.fcs", 2)$marker_new, "IFNg")
pl <- plan_marker_names(cen, aliases = parse_marker_aliases("INFy=IFNg;TNF=TNFa"))
p <- pl$plan
eq("the CLI alias form is accepted", row("A.fcs", 2)$marker_new, "IFNg")
eq("an alias to the current spelling is unchanged",
   plan_marker_names(cen, aliases = c(CD4 = "CD4"))$plan$source[
     cen$filename == "A.fcs" & cen$index == 4L], "unchanged")

# --- dry run --------------------------------------------------------------------
cat("\n-- dry run ---------------------------------------------------------\n")
res <- quiet(apply_marker_plan(pl$plan, d))
ok("dry run visits every file", nrow(res) == 4L)
eq("dry run projects a strategy", sort(unique(res$status)), c("ok", "unchanged"))
ok("dry run changes no bytes", identical(vapply(files, fx_md5, character(1)), md5))
ok("dry run writes no ledger", !file.exists(file.path(d, MARKER_LEDGER_DEFAULT)))
throws("apply without in_place or out_dir is refused",
       quiet(apply_marker_plan(pl$plan, d, dry_run = FALSE)), "refusing to guess")
throws("in_place and out_dir together are refused",
       quiet(apply_marker_plan(pl$plan, d, dry_run = FALSE, in_place = TRUE,
                               out_dir = file.path(root, "x"))), "mutually exclusive")
ok("...and still nothing was written", identical(vapply(files, fx_md5, character(1)), md5))

# --- apply in place -----------------------------------------------------------
cat("\n-- apply in place --------------------------------------------------\n")
d <- make_cohort(file.path(root, "apply"))
files <- files_of(d)
before <- lapply(files, data_bytes)
cen <- marker_census(files, d)
pl <- plan_marker_names(cen, canonical = CANON)
res <- quiet(apply_marker_plan(pl$plan, d, in_place = TRUE, dry_run = FALSE))
eq("every file with changes is written", res$status[res$filename != "sub/D.fcs"], rep("ok", 3))
eq("a file already canonical is reported unchanged", res$status[res$filename == "sub/D.fcs"], "unchanged")
fA <- file.path(d, "A.fcs"); fB <- file.path(d, "B.fcs"); fC <- file.path(d, "C.fcs")
eq("A reads back the new $P2S", pns(fA, 2), "IFNg")
eq("B reads back both renames", c(pns(fB, 2), pns(fB, 3)), c("IFNg", "TNFa"))
eq("C gained $PnS from the channel map", c(pns(fC, 2), pns(fC, 3), pns(fC, 4)),
   c("IFNg", "TNFa", "CD4"))
ok("C's scatter and Time still have none", is.na(pns(fC, 1)) && is.na(pns(fC, 6)))
eq("$PnN is never touched", fcs_keyword(fC, "$P2N"), "PE-A")
ok("event bytes are identical in every file",
   all(mapply(function(f, b) identical(data_bytes(f), b), files, before)))
ok("every file still validates", all(vapply(files, function(f) fcs_inspect(f)$ok, logical(1))))
ok("no temp or sidecar files survive",
   !length(list.files(d, pattern = "fcstmp|fcsold|fcsprefix", recursive = TRUE)))
ledger <- file.path(d, MARKER_LEDGER_DEFAULT)
ok("the ledger sits next to the data", file.exists(ledger))
led <- utils::read.csv(ledger, stringsAsFactors = FALSE, colClasses = "character")
eq("the ledger records the replaced value", led$old_value[led$keyword == "$P2S" &
                                                           basename(led$path) == "A.fcs"], "INFy")
eq("...and marks a created keyword as absent",
   unique(led$old_value[basename(led$path) == "C.fcs" & nzchar(led$keyword)]), MARKER_LEDGER_ABSENT)
eq("one run id", length(unique(led$run_id)), 1L)

# idempotence
cen2 <- marker_census(files, d)
pl2 <- plan_marker_names(cen2, canonical = CANON)
eq("a fresh plan after the run changes nothing", unique(pl2$plan$source), "unchanged")
md5 <- vapply(files, fx_md5, character(1))
res2 <- quiet(apply_marker_plan(pl2$plan, d, in_place = TRUE, dry_run = FALSE))
eq("a second apply is a no-op", unique(res2$status), "unchanged")
ok("...byte for byte", identical(vapply(files, fx_md5, character(1)), md5))
eq("...and adds no ledger rows", nrow(utils::read.csv(ledger, colClasses = "character")), nrow(led))

# --- revert ---------------------------------------------------------------------
cat("\n-- revert ----------------------------------------------------------\n")
quiet(revert_marker_plan(ledger, verbose = FALSE))
eq("revert restores A's transposed spelling", pns(fA, 2), "INFy")
eq("revert restores B", c(pns(fB, 2), pns(fB, 3)), c("IFN", "TNF"))
ok("revert removes the $PnS C never had", all(is.na(c(pns(fC, 2), pns(fC, 3), pns(fC, 4)))))
ok("event bytes survive the revert",
   all(mapply(function(f, b) identical(data_bytes(f), b), files, before)))
ok("files are still valid", all(vapply(files, function(f) fcs_inspect(f)$ok, logical(1))))
throws("a second revert has nothing to undo", quiet(revert_marker_plan(ledger, verbose = FALSE)),
       "already been reverted")

# --- copies ---------------------------------------------------------------------
cat("\n-- copies ----------------------------------------------------------\n")
d <- make_cohort(file.path(root, "copysrc"))
files <- files_of(d)
md5 <- vapply(files, fx_md5, character(1))
outd <- file.path(root, "copyout")
cen <- marker_census(files, d)
pl <- plan_marker_names(cen, canonical = CANON)
res <- quiet(apply_marker_plan(pl$plan, d, out_dir = outd, dry_run = FALSE))
ok("originals are untouched", identical(vapply(files, fx_md5, character(1)), md5))
eq("copies mirror the layout", sort(list.files(outd, "\\.fcs$", recursive = TRUE)),
   c("A.fcs", "B.fcs", "C.fcs", "sub/D.fcs"))
eq("copies carry the new names", pns(file.path(outd, "B.fcs"), 2), "IFNg")
ok("an unchanged file is copied verbatim",
   identical(fx_md5(file.path(outd, "sub", "D.fcs")), fx_md5(file.path(d, "sub", "D.fcs"))))
ok("the ledger lands with the copies", file.exists(file.path(outd, MARKER_LEDGER_DEFAULT)))
throws("out_dir inside fcs_dir is refused",
       quiet(apply_marker_plan(pl$plan, d, out_dir = file.path(d, "out"), dry_run = FALSE)),
       "must not be")

# --- symlinks -------------------------------------------------------------------
cat("\n-- symlinks --------------------------------------------------------\n")
d <- make_cohort(file.path(root, "links"))
target <- file.path(root, "links", "A.fcs")
link <- file.path(d, "L.fcs")
if (isTRUE(file.symlink(target, link))) {
  files <- files_of(d)
  cen <- marker_census(files, d)
  pl <- plan_marker_names(cen, canonical = CANON)
  res <- quiet(apply_marker_plan(pl$plan, d, in_place = TRUE, dry_run = FALSE))
  eq("a symlink is refused", res$status[res$filename == "L.fcs"], "refused")
  eq("...while the real files are written", res$status[res$filename == "B.fcs"], "ok")
} else {
  cat("  SKIP  symlinks not supported here\n")
}

# --- entry point ----------------------------------------------------------------
cat("\n-- entry point -----------------------------------------------------\n")
d <- make_cohort(file.path(root, "entry"))
files <- files_of(d)
md5 <- vapply(files, fx_md5, character(1))
out <- quiet(harmonize_marker_names(d, probe_only = TRUE))
ok("probe returns the census and no plan", !is.null(out$census) && is.null(out$plan))
out <- quiet(harmonize_marker_names(d, canonical = CANON))
ok("the default is a dry run", identical(vapply(files, fx_md5, character(1)), md5))
ok("...that reports the plan grouped by channel", grepl("Planned \\$PnS changes", out$report))
throws("apply needs a destination", quiet(harmonize_marker_names(d, canonical = CANON, apply = TRUE)),
       "refusing to guess")
out <- quiet(harmonize_marker_names(d, canonical = CANON, apply = TRUE, in_place = TRUE))
eq("apply in place writes", pns(file.path(d, "B.fcs"), 2), "IFNg")
quiet(harmonize_marker_names(d, revert = TRUE))
eq("revert through the entry point", pns(file.path(d, "B.fcs"), 2), "IFN")

unlink(root, recursive = TRUE)
cat("\n=== ", .pass, " passed, ", .fail, " failed ===\n", sep = "")
if (.fail > 0) quit(status = 1)
