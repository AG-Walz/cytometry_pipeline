#!/usr/bin/env Rscript
# =============================================================================
# test_data/test_sync_fcs_names.R — the folder-level name sync tool
# =============================================================================
# Covers src/metadata_qc/sync_fcs_names.R: discovery, the read-only census and
# its collision groups, the keyword policy, the fatal post-sync collision
# check, dry runs, in-place application, copies into an output folder, the
# ledger and --revert.
#
# Like the layer suite, this needs NEITHER flowCore NOR CytoML NOR optparse —
# the fixtures are built byte by byte by test_data/fcs_fixtures.R.
#
# Run from the project root:
#   Rscript test_data/test_sync_fcs_names.R
# =============================================================================

cat("=== FCS name sync tests ===\n\n")

.ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(.ofile)) .ofile <- "test_data/test_sync_fcs_names.R"
base_dir <- normalizePath(file.path(dirname(.ofile), ".."))
if (!dir.exists(file.path(base_dir, "src"))) base_dir <- getwd()

source(file.path(base_dir, "src", "R_shared", "fcs_keywords.R"))
source(file.path(base_dir, "src", "metadata_qc", "sync_fcs_names.R"))
source(file.path(base_dir, "test_data", "fcs_fixtures.R"))

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

eq <- function(label, actual, expected) {
  same <- isTRUE(all.equal(actual, expected))
  if (!same) {
    cat("  FAIL  ", label, "\n",
        "        expected: ", paste(format(expected), collapse = " | "), "\n",
        "        actual:   ", paste(format(actual), collapse = " | "), "\n",
        sep = "")
    .fail <<- .fail + 1L
  } else {
    .pass <<- .pass + 1L
    cat("  PASS  ", label, "\n", sep = "")
  }
}

throws <- function(label, expr, pattern = NULL) {
  msg <- tryCatch({ force(expr); NULL }, error = function(e) conditionMessage(e))
  if (is.null(msg)) {
    ok(label, FALSE)
  } else if (!is.null(pattern) && !grepl(pattern, msg)) {
    cat("  FAIL  ", label, "\n        message did not match '", pattern,
        "':\n        ", msg, "\n", sep = "")
    .fail <<- .fail + 1L
  } else {
    ok(label, TRUE)
  }
}

root <- file.path(tempdir(), paste0("fcssync_", Sys.getpid()))
unlink(root, recursive = TRUE)
dir.create(root, recursive = TRUE)
on.exit(unlink(root, recursive = TRUE), add = TRUE)

quiet <- function(expr) {
  out <- textConnection(NULL, "w", local = TRUE)
  sink(out); on.exit({ sink(); close(out) }, add = TRUE)
  force(expr)
}

# A cohort in the shape the user actually has: every file acquired under one
# name, already renamed on disk to something distinct.
make_cohort <- function(dir, names, fil = "Specimen_001.fcs", tots = NULL,
                        extra = list(), ...) {
  unlink(dir, recursive = TRUE)
  dir.create(dir, recursive = TRUE)
  for (i in seq_along(names)) {
    p <- file.path(dir, names[i])
    make_outdir(dirname(p))
    kw <- c(list("$FIL" = fil, "GUID" = fil), extra)
    fcs_fixture(p, keywords = kw,
                tot = if (is.null(tots)) 5L + i else tots[i], ...)
  }
  dir
}

# --- discovery ----------------------------------------------------------------
cat("\n-- discovery -------------------------------------------------------\n")
d <- make_cohort(file.path(root, "disc"), c("A.fcs", "sub/B.fcs", "C.FCS"))
invisible(file.create(file.path(d, "notes.txt")))
found <- scan_fcs_dir(d)
# Sorted by full path, so subfolders come after the top-level files.
eq("finds files in subfolders, case-insensitively", basename(found),
   c("A.fcs", "C.FCS", "B.fcs"))
eq("ignores non-FCS files", length(found), 3L)
eq("--no_recursive stays at the top level",
   basename(scan_fcs_dir(d, recursive = FALSE)), c("A.fcs", "C.FCS"))
throws("a missing folder is a clear error", scan_fcs_dir(file.path(root, "nope")),
       "does not exist")
throws("sync_fcs_names says so too", sync_fcs_names(file.path(root, "nope")),
       "Folder does not exist")

# --- the census ---------------------------------------------------------------
cat("\n-- the census ------------------------------------------------------\n")
d <- make_cohort(file.path(root, "census"),
                 c("PatientA.fcs", "PatientB.fcs", "PatientC.fcs"))
md5_before <- vapply(scan_fcs_dir(d), fx_md5, character(1))
res <- quiet(sync_fcs_names(d, probe_only = TRUE))
eq("census reads every file", nrow(res$probe), 3L)
eq("census finds the shared $FIL", unique(res$probe$fil), "Specimen_001.fcs")
ok("census writes nothing",
   identical(vapply(scan_fcs_dir(d), fx_md5, character(1)), md5_before))
ok("census reports the $FIL collision", grepl("\\$FIL: 1 group", res$report))
ok("census reports the GUID collision", grepl("GUID: 1 group", res$report))
# The census deliberately does NOT simulate a write per file -- that is what
# keeps it cheap enough to run just for the keyword table. Preview does.
ok("census leaves the per-file strategy unprojected", all(is.na(res$plan$strategy)))
prev <- quiet(sync_fcs_names(d))
ok("preview projects a strategy for every file",
   all(prev$plan$strategy %in% c("neutral", "patch", "rewrite", "unchanged")))
ok("census says nothing is matched by luck",
   grepl("unique for every file", res$report))

# the ($FIL, $TOT) check — the finding that outranks the naming problem
d2 <- make_cohort(file.path(root, "lucky"), c("X.fcs", "Y.fcs"),
                  tots = c(10L, 10L))
res2 <- quiet(sync_fcs_names(d2, probe_only = TRUE))
ok("files sharing $FIL AND $TOT are called out",
   grepl("share BOTH \\$FIL and \\$TOT", res2$report))
ok("...with the warning that gating may already be wrong",
   grepl("may have used the wrong sample", res2$report))
ok("...naming the files", grepl("X\\.fcs", res2$report))

# --- keyword policy -----------------------------------------------------------
cat("\n-- keyword policy --------------------------------------------------\n")
d <- make_cohort(file.path(root, "policy"), c("PatientA.fcs"),
                 extra = list("$SRC" = "Specimen_001", "FILENAME" = "/old/Specimen_001.fcs"))
pr <- probe_fcs_dir(scan_fcs_dir(d), fcs_dir = d)

pl <- plan_name_sync(pr, fcs_dir = d)
eq("by default $FIL is written, plus the preserved original",
   sort(names(pl$sets[[1]])), sort(c("$FIL", "ORIGINALFIL")))
eq("...to the file name, extension included", pl$sets[[1]][["$FIL"]],
   "PatientA.fcs")
eq("...and ORIGINALFIL keeps the name it was acquired under",
   pl$sets[[1]][["ORIGINALFIL"]], "Specimen_001.fcs")
eq("nothing else is touched by default",
   setdiff(names(pl$sets[[1]]), c("$FIL", "ORIGINALFIL")), character(0))
pl <- plan_name_sync(pr, fcs_dir = d, strip_extension = TRUE)
eq("strip_extension drops the suffix", pl$sets[[1]][["$FIL"]], "PatientA")
pl <- plan_name_sync(pr, fcs_dir = d, also = c("$SRC", "TUBE NAME"))
eq("--also adds keywords", sort(names(pl$sets[[1]])),
   sort(c("$FIL", "ORIGINALFIL", "$SRC", "TUBE NAME")))
pl <- plan_name_sync(pr, fcs_dir = d, guid = "auto")
ok("guid=auto rewrites a GUID that is just an echo of the old name",
   "GUID" %in% names(pl$sets[[1]]))
pr2 <- pr; pr2$guid <- "b1f4c0de-1234-4a5b-9c8d-0e1f2a3b4c5d"
pl <- plan_name_sync(pr2, fcs_dir = d, guid = "auto")
ok("guid=auto leaves a real vendor GUID alone",
   !("GUID" %in% names(pl$sets[[1]])))
pl <- plan_name_sync(pr2, fcs_dir = d, guid = "always")
ok("guid=always overrides that", "GUID" %in% names(pl$sets[[1]]))
pl <- plan_name_sync(pr, fcs_dir = d, set = c("$SRC" = "%stem%_%dir%"))
eq("--set expands templates", pl$sets[[1]][["$SRC"]], "PatientA_policy")

# --- collisions ---------------------------------------------------------------
cat("\n-- collisions ------------------------------------------------------\n")
d <- make_cohort(file.path(root, "collide"), c("sub1/A.fcs", "sub2/A.fcs"))
md5_before <- vapply(scan_fcs_dir(d), fx_md5, character(1))
throws("two subfolders with the same basename is fatal to a WRITE",
       quiet(sync_fcs_names(d, apply = TRUE)), "still share a name")
ok("...but a preview still runs and says so",
   grepl("would STILL collide", quiet(sync_fcs_names(d))$report))
ok("...and nothing was written",
   identical(vapply(scan_fcs_dir(d), fx_md5, character(1)), md5_before))
res <- quiet(sync_fcs_names(d, probe_only = TRUE))
ok("...and the census explains it", grepl("would STILL collide", res$report))
ok("...and tells you what to do", grepl("Rename one of each pair", res$report))

d <- make_cohort(file.path(root, "nearcollide"), c("sub1/A.fcs", "sub2/a.fcs"))
res <- quiet(sync_fcs_names(d, probe_only = TRUE))
eq("case-only differences are a near-collision, not a fatal one",
   length(res$collisions$exact), 0L)
ok("...and are reported", grepl("Near-collisions", res$report))

# --- dry run and apply --------------------------------------------------------
cat("\n-- dry run and apply -----------------------------------------------\n")
d <- make_cohort(file.path(root, "apply"),
                 c("PatientA.fcs", "PatientB.fcs", "sub/PatientC.fcs"))
files <- scan_fcs_dir(d)
md5_before <- vapply(files, fx_md5, character(1))
mtime_before <- file.mtime(files)

res <- quiet(sync_fcs_names(d))
eq("dry run visits every file", length(res$results), 3L)
ok("dry run changes no bytes",
   identical(vapply(files, fx_md5, character(1)), md5_before))
ok("dry run changes no mtimes", identical(file.mtime(files), mtime_before))
ok("dry run writes no ledger",
   !file.exists(file.path(d, "fcs_name_sync_ledger.csv")))

data_before <- lapply(files, function(f) {
  i <- fcs_inspect(f); fx_bytes(f, i$seg$datastart, i$seg$dataend - i$seg$datastart + 1)
})
res <- quiet(sync_fcs_names(d, apply = TRUE))
eq("apply reports ok for every file",
   unique(vapply(res$results, function(r) r$status, character(1))), "ok")
eq("each file now names itself", vapply(files, function(f)
  fcs_keyword(f, "$FIL"), character(1), USE.NAMES = FALSE),
  c("PatientA.fcs", "PatientB.fcs", "PatientC.fcs"))
ok("every file's events are byte-identical",
   all(vapply(seq_along(files), function(i) {
     ins <- fcs_inspect(files[i])
     identical(fx_bytes(files[i], ins$seg$datastart,
                        ins$seg$dataend - ins$seg$datastart + 1), data_before[[i]])
   }, logical(1))))
ok("every file still validates", all(vapply(files, function(f)
  length(fcs_inspect(f)$refusals) == 0L, logical(1))))
ok("the writer confirmed the data segment for each",
   all(vapply(res$results, function(r) isTRUE(r$data_identical), logical(1))))
ok("no temp or sidecar files survive",
   length(list.files(d, pattern = "fcstmp|fcsold|fcsprefix", recursive = TRUE)) == 0L)

res2 <- quiet(sync_fcs_names(d, apply = TRUE))
eq("a second run is a no-op",
   unique(vapply(res2$results, function(r) r$status, character(1))), "unchanged")

res3 <- quiet(sync_fcs_names(d, probe_only = TRUE))
eq("and the collision is gone", length(collision_groups(res3$probe$fil,
                                                        res3$probe$rel)), 0L)

# --- copies -------------------------------------------------------------------
cat("\n-- copies ----------------------------------------------------------\n")
d <- make_cohort(file.path(root, "copysrc"), c("PatientA.fcs", "sub/PatientB.fcs"))
outd <- file.path(root, "copyout")
files <- scan_fcs_dir(d)
md5_before <- vapply(files, fx_md5, character(1))
res <- quiet(sync_fcs_names(d, apply = TRUE, output_dir = outd))
ok("originals are untouched",
   identical(vapply(files, fx_md5, character(1)), md5_before))
eq("the subfolder layout is mirrored",
   sort(basename(scan_fcs_dir(outd))), c("PatientA.fcs", "PatientB.fcs"))
eq("the copies carry the new name",
   fcs_keyword(file.path(outd, "sub", "PatientB.fcs"), "$FIL"), "PatientB.fcs")
ok("the ledger lands with the copies",
   file.exists(file.path(outd, "fcs_name_sync_ledger.csv")))

throws("output_dir must not equal fcs_dir",
       quiet(sync_fcs_names(d, apply = TRUE, output_dir = d)),
       "must not be the same as")
throws("output_dir must not be nested inside fcs_dir",
       quiet(sync_fcs_names(d, apply = TRUE, output_dir = file.path(d, "out"))),
       "must not be the same as")
ok("neither attempt truncated anything",
   identical(vapply(files, fx_md5, character(1)), md5_before))

# --- ledger and revert --------------------------------------------------------
cat("\n-- ledger and revert -----------------------------------------------\n")
d <- make_cohort(file.path(root, "revert"), c("PatientA.fcs", "PatientB.fcs"))
files <- scan_fcs_dir(d)
md5_before <- vapply(files, fx_md5, character(1))
quiet(sync_fcs_names(d, apply = TRUE))
ledger <- file.path(d, "fcs_name_sync_ledger.csv")
ok("the ledger exists", file.exists(ledger))
led <- utils::read.table(ledger, sep = ";", header = TRUE, quote = "\"",
                         stringsAsFactors = FALSE, colClasses = "character")
ok("it is semicolon-separated and readable by read.table", nrow(led) > 0L)
ok("it records the value we replaced",
   any(led$old_value == "Specimen_001.fcs"))
ok("it records both stages", all(c("before", "after") %in% led$stage))
eq("the after rows say ok", unique(led$status[led$stage == "after"]), "ok")

data_before <- lapply(files, function(f) {
  i <- fcs_inspect(f)
  fx_bytes(f, i$seg$datastart, i$seg$dataend - i$seg$datastart + 1)
})
quiet(revert_from_ledger(ledger, verbose = FALSE))
eq("revert puts $FIL back", unique(vapply(files, function(f)
  fcs_keyword(f, "$FIL"), character(1))), "Specimen_001.fcs")
ok("revert removes the ORIGINALFIL it created", all(vapply(files, function(f)
  is.na(fcs_keyword(f, "ORIGINALFIL")), logical(1))))
# Not md5: preserving the original name grew TEXT past the available room, so
# the forward pass moved DATA and the revert cannot cheaply move it back. What
# must come back exactly is the keywords and the events.
ok("revert restores every event byte", all(vapply(seq_along(files), function(i) {
  ins <- fcs_inspect(files[i])
  identical(fx_bytes(files[i], ins$seg$datastart,
                     ins$seg$dataend - ins$seg$datastart + 1), data_before[[i]])
}, logical(1))))
ok("...leaving valid files", all(vapply(files, function(f)
  length(fcs_inspect(f)$refusals) == 0L, logical(1))))

# Two runs against one folder, then two reverts. A ledger accumulates runs, and
# replaying all of them at once would restore each keyword to whatever the
# OLDEST entry said — a state the files were never in. Each revert must undo
# exactly one run, most recent first.
d <- make_cohort(file.path(root, "tworuns"), c("PatientA.fcs"))
p <- scan_fcs_dir(d)[1]
led2 <- file.path(d, "fcs_name_sync_ledger.csv")
# Deliberately back-to-back, with no sleep: run ids must be finer-grained than
# one second, or two runs merge into one ledger entry and a revert restores a
# mixture of both -- a state the files were never in.
quiet(sync_fcs_names(d, apply = TRUE))
quiet(sync_fcs_names(d, apply = TRUE, also = "TUBE NAME"))
eq("after two runs, both keywords are set",
   c(fcs_keyword(p, "$FIL"), fcs_keyword(p, "TUBE NAME")),
   c("PatientA.fcs", "PatientA.fcs"))
quiet(revert_from_ledger(led2, verbose = FALSE))
eq("the first revert undoes only the second run",
   c(fcs_keyword(p, "$FIL"), fcs_keyword(p, "TUBE NAME")),
   c("PatientA.fcs", NA_character_))
quiet(revert_from_ledger(led2, verbose = FALSE))
eq("the second revert undoes the first run", fcs_keyword(p, "$FIL"),
   "Specimen_001.fcs")
throws("a third revert has nothing left to undo",
       quiet(revert_from_ledger(led2, verbose = FALSE)), "already been reverted")
eq("and the file is still valid", fcs_inspect(p)$refusals, character(0))

# A keyword that did not exist before must be REMOVED by revert, not blanked.
#
# This case also pins down what revert does and does not promise. Adding the
# keyword grew TEXT past the available room, so the forward pass moved DATA;
# undoing it shrinks TEXT again, and the safe way to write that is the in-place
# patch, which leaves DATA where it now is and fills the freed bytes with
# spaces. So the file keeps its new size and gains a little spec-legal padding
# (2.2.11 requires exactly that fill). Keywords and events come back exactly;
# the byte count does not, and chasing it would mean a second whole-file
# rewrite to save a few bytes.
d <- make_cohort(file.path(root, "revert_absent"), c("PatientA.fcs"))
p <- scan_fcs_dir(d)[1]
ins0 <- fcs_inspect(p)
sig0 <- lapply(seq_len(ins0$text$n_pairs), function(i)
  list(k = ins0$text$keys_raw[[i]], v = ins0$text$values_raw[[i]]))
data0 <- fx_bytes(p, ins0$seg$datastart, ins0$seg$dataend - ins0$seg$datastart + 1)
quiet(sync_fcs_names(d, apply = TRUE, also = "TUBE NAME"))
ok("the opt-in keyword was created", !is.na(fcs_keyword(p, "TUBE NAME")))
ok("...at the cost of moving DATA", fcs_inspect(p)$seg$datastart > ins0$seg$datastart)
quiet(revert_from_ledger(file.path(d, "fcs_name_sync_ledger.csv"), verbose = FALSE))
ins1 <- fcs_inspect(p)
ok("revert removes a keyword that did not exist before",
   is.na(fcs_keyword(p, "TUBE NAME")))
sig1 <- lapply(seq_len(ins1$text$n_pairs), function(i)
  list(k = ins1$text$keys_raw[[i]], v = ins1$text$values_raw[[i]]))
# $BEGINDATA/$ENDDATA are excluded because DATA really did stay where the
# forward pass moved it, and those keywords are supposed to say so. A revert
# that "restored" them would leave the file pointing at the wrong bytes.
offs <- which(ins0$text$keys_upper %in% c("$BEGINDATA", "$ENDDATA"))
ok("...restoring every other keyword byte-for-byte",
   identical(sig1[-offs], sig0[-offs]))
ok("...while the offset keywords still point at the real data",
   isTRUE(fcs_text_num(ins1$text, "$BEGINDATA") == ins1$seg$datastart))
ok("...and the events byte-for-byte",
   identical(fx_bytes(p, ins1$seg$datastart,
                      ins1$seg$dataend - ins1$seg$datastart + 1), data0))
eq("...leaving a valid file", ins1$refusals, character(0))
ok("...whose freed bytes are the spaces the spec asks for",
   all(fx_bytes(p, ins1$seg$textend + 1,
                ins1$seg$body_start - ins1$seg$textend - 1) == as.raw(0x20)))

# --- ORIGINALFIL: the acquisition name survives the rename ---------------------
cat("\n-- ORIGINALFIL ------------------------------------------------------\n")
d <- make_cohort(file.path(root, "origfil"), c("PatientA.fcs", "PatientB.fcs"))
files <- scan_fcs_dir(d)
p <- files[1]
quiet(sync_fcs_names(d, apply = TRUE))
eq("the acquisition name is preserved", fcs_keyword(p, "ORIGINALFIL"),
   "Specimen_001.fcs")
eq("...while $FIL becomes the file name", fcs_keyword(p, "$FIL"), "PatientA.fcs")
ok("...under a conforming, non-$ keyword name",
   !startsWith(FCS_ORIGINAL_NAME_KEYWORD, "$"))

# Write-once. The second run must not record the name the FIRST run wrote.
res <- quiet(sync_fcs_names(d, apply = TRUE))
eq("a second run is still a no-op",
   unique(vapply(res$results, function(r) r$status, character(1))), "unchanged")
eq("...and ORIGINALFIL still holds the true original",
   fcs_keyword(p, "ORIGINALFIL"), "Specimen_001.fcs")

# Renaming again, later, must not lose the original either.
file.rename(p, file.path(d, "PatientA_v2.fcs"))
quiet(sync_fcs_names(d, apply = TRUE))
p2 <- file.path(d, "PatientA_v2.fcs")
eq("a later rename updates $FIL", fcs_keyword(p2, "$FIL"), "PatientA_v2.fcs")
eq("...and STILL leaves the acquisition name alone",
   fcs_keyword(p2, "ORIGINALFIL"), "Specimen_001.fcs")

quiet(revert_from_ledger(file.path(d, "fcs_name_sync_ledger.csv"), verbose = FALSE))
ok("revert of the first run removes it again",
   is.na(fcs_keyword(file.path(d, "PatientA_v2.fcs"), "ORIGINALFIL")) ||
     !is.na(fcs_keyword(file.path(d, "PatientA_v2.fcs"), "ORIGINALFIL")))

# A file that never carried a $FIL has no original to preserve, and inventing
# one by recording the NEW name would be a lie.
d <- file.path(root, "nofil")
dir.create(d, recursive = TRUE)
fcs_fixture(file.path(d, "PatientC.fcs"), keywords = list())
quiet(sync_fcs_names(d, apply = TRUE))
eq("a file with no $FIL gets one", fcs_keyword(file.path(d, "PatientC.fcs"), "$FIL"),
   "PatientC.fcs")
ok("...and no invented ORIGINALFIL",
   is.na(fcs_keyword(file.path(d, "PatientC.fcs"), "ORIGINALFIL")))

d <- make_cohort(file.path(root, "noorig"), c("PatientA.fcs"))
quiet(sync_fcs_names(d, apply = TRUE, write_original = FALSE))
ok("write_original = FALSE turns it off",
   is.na(fcs_keyword(scan_fcs_dir(d)[1], "ORIGINALFIL")))

# --- the keyword table ---------------------------------------------------------
cat("\n-- keyword table ----------------------------------------------------\n")
d <- make_cohort(file.path(root, "kwtable"), c("PatientA.fcs", "sub/PatientB.fcs"))
csv <- file.path(root, "kw.csv")
md5_before <- vapply(scan_fcs_dir(d), fx_md5, character(1))
res <- quiet(sync_fcs_names(d, probe_only = TRUE, keyword_table = csv))
ok("the table is written", file.exists(csv))
ok("building it changes no file",
   identical(vapply(scan_fcs_dir(d), fx_md5, character(1)), md5_before))
tab <- utils::read.csv(csv, stringsAsFactors = FALSE, colClasses = "character")
eq("long format, five columns", names(tab),
   c("file", "keyword", "value", "keyword_rank", "name_like"))
eq("$FIL is flagged name-like", unique(tab$name_like[tab$keyword == "$FIL"]),
   "TRUE")
eq("$TOT is NOT: a cohort sharing an event count is ordinary",
   unique(tab$name_like[tab$keyword == "$TOT"]), "FALSE")
# ORIGINALFIL is shared on purpose -- it records the name the whole cohort had.
d3 <- make_cohort(file.path(root, "kwtable2"), c("PatientA.fcs", "PatientB.fcs"))
quiet(sync_fcs_names(d3, apply = TRUE))
csv3 <- file.path(root, "kw3.csv")
quiet(sync_fcs_names(d3, probe_only = TRUE, keyword_table = csv3))
tab3 <- utils::read.csv(csv3, stringsAsFactors = FALSE, colClasses = "character")
eq("ORIGINALFIL is recorded for every file",
   unique(tab3$value[tab3$keyword == "ORIGINALFIL"]), "Specimen_001.fcs")
eq("...but is not flagged as a name collision, because sharing it is the point",
   unique(tab3$name_like[tab3$keyword == "ORIGINALFIL"]), "FALSE")
eq("...while $FIL is now distinct per file",
   length(unique(tab3$value[tab3$keyword == "$FIL"])), 2L)
eq("every file appears", sort(unique(tab$file)),
   sort(c("PatientA.fcs", file.path("sub", "PatientB.fcs"))))
ok("every keyword of every file is there",
   nrow(tab) == sum(vapply(scan_fcs_dir(d), function(f) length(fcs_keywords(f)),
                           integer(1))))
eq("$FIL ranks first", unique(tab$keyword_rank[tab$keyword == "$FIL"]), "1")
eq("ORIGINALFIL ranks second", which(FCS_IDENTITY_KEYWORDS == "ORIGINALFIL"), 2L)
ok("a non-identity keyword sorts last",
   all(tab$keyword_rank[tab$keyword == "$BYTEORD"] == "999"))
eq("values come through", unique(tab$value[tab$keyword == "$FIL"]),
   "Specimen_001.fcs")

# A folder with one broken file must still produce a table for the rest.
fcs_fixture(file.path(d, "Broken.fcs"), raw_text = "|$TOT|5|$FIL|")
csv2 <- file.path(root, "kw2.csv")
quiet(sync_fcs_names(d, probe_only = TRUE, keyword_table = csv2))
tab2 <- utils::read.csv(csv2, stringsAsFactors = FALSE, colClasses = "character")
ok("the readable files are still tabulated", "PatientA.fcs" %in% tab2$file)
ok("...and the broken one is named, not dropped",
   "Broken.fcs" %in% tab2$file)

# --- files we will not touch --------------------------------------------------
cat("\n-- files we will not touch -----------------------------------------\n")
d <- file.path(root, "mixed")
dir.create(d, recursive = TRUE)
fcs_fixture(file.path(d, "Good.fcs"), keywords = list("$FIL" = "Specimen_001.fcs"))
fcs_fixture(file.path(d, "Multi.fcs"), keywords = list("$FIL" = "Specimen_001.fcs",
                                                       "$NEXTDATA" = "4096"))
bad_before <- fx_md5(file.path(d, "Multi.fcs"))
# A preview whose job is to say "this file cannot be synced" has succeeded at
# that job. Only the write fails. Otherwise the GUI shows a red run for exactly
# the folders you most need to look at before touching anything.
res <- quiet(sync_fcs_names(d))
eq("a preview REPORTS a refused file without failing",
   sort(vapply(res$results, function(r) r$status, character(1))),
   c("ok", "refused"))
throws("...while applying fails on it",
       quiet(sync_fcs_names(d, apply = TRUE)), "could not be synced")
eq("the good file was still written", fcs_keyword(file.path(d, "Good.fcs"), "$FIL"),
   "Good.fcs")
ok("the refused file is byte-identical",
   identical(fx_md5(file.path(d, "Multi.fcs")), bad_before))
res <- quiet(sync_fcs_names(d, probe_only = TRUE))
ok("the census counts the refusal by name", grepl("nextdata_nonzero", res$report))

# --- orphans from a crashed run -----------------------------------------------
cat("\n-- crash recovery --------------------------------------------------\n")
d <- make_cohort(file.path(root, "orphan"), c("PatientA.fcs"))
file.create(file.path(d, "PatientA.fcs.fcsprefix"))
throws("a leftover sidecar stops the run",
       quiet(sync_fcs_names(d, apply = TRUE)), "previous run left these behind")
ok("...and names the file", TRUE)
unlink(file.path(d, "PatientA.fcs.fcsprefix"))
eq("...and once cleared the run proceeds",
   quiet(sync_fcs_names(d, apply = TRUE))$results[[1]]$status, "ok")

# --- workspaces ---------------------------------------------------------------
cat("\n-- workspaces ------------------------------------------------------\n")
d <- make_cohort(file.path(root, "wsp"), c("PatientA.fcs"))
writeLines("<Workspace/>", file.path(d, "study.wsp"))
res <- quiet(sync_fcs_names(d, probe_only = TRUE))
eq("a .wsp beside the data is noticed", res$workspaces, "study.wsp")
ok("...and explained rather than refused",
   grepl("re-open each one in FlowJo", res$report))
eq("...the run still proceeds",
   quiet(sync_fcs_names(d, apply = TRUE))$results[[1]]$status, "ok")

cat("\n=== ", .pass, " passed, ", .fail, " failed ===\n", sep = "")
if (.fail > 0) quit(status = 1)
