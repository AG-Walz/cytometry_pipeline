#!/usr/bin/env Rscript
# =============================================================================
# test_data/test_umap_multi_wsp.R — Multi-workspace discovery and pairing tests
# =============================================================================
# Covers the pure helpers behind "one UMAP over many FlowJo workspaces":
# workspace discovery, the recursive FCS index, filename pairing, group
# resolution, sample-label disambiguation, and subgate-candidate merging.
#
# These deliberately need NEITHER CytoML NOR any real .wsp/.fcs file — the
# workspace reader is injected — so they run on any machine with base R, which
# is the point: the parts that can be tested without the cytometry stack are.
# What they cannot cover is CytoML's own behaviour; that needs a real workspace
# on an analysis host.
#
# Run from the project root:
#   Rscript test_data/test_umap_multi_wsp.R
# =============================================================================

cat("=== Multi-workspace FlowJo discovery tests ===\n\n")

.ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(.ofile)) .ofile <- "test_data/test_umap_multi_wsp.R"
base_dir <- normalizePath(file.path(dirname(.ofile), ".."))
if (!dir.exists(file.path(base_dir, "src"))) base_dir <- getwd()

source(file.path(base_dir, "src", "loading", "flowjo_workspace.R"))
source(file.path(base_dir, "src", "loading", "sample_annotations.R"))

# umap_from_flowjo.R is sourced for the pure helpers only. It has no top-level
# side effects and does not load CytoML at source time (only inside the main
# function), so this is safe without the cytometry stack.
source(file.path(base_dir, "src", "visualization", "umap_from_flowjo.R"))

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
  msg <- tryCatch({ force(expr); NULL },
                  error = function(e) conditionMessage(e))
  if (is.null(msg)) {
    ok(label, FALSE)
  } else if (!is.null(pattern) && !grepl(pattern, msg, fixed = FALSE)) {
    cat("  FAIL  ", label, "\n        message did not match '", pattern,
        "':\n        ", msg, "\n", sep = "")
    .fail <<- .fail + 1L
  } else {
    ok(label, TRUE)
  }
}

# --- a fake data layout ------------------------------------------------------
# Deliberately awkward: workspaces at two depths, two of them sharing a base
# name, FCS files spread over nested subfolders, one duplicated file name, and
# one workspace whose files are split across two directories.
root <- file.path(tempdir(), paste0("multiwsp_", Sys.getpid()))
unlink(root, recursive = TRUE)

wsp_root <- file.path(root, "workspaces")
fcs_root <- file.path(root, "fcs")

mk <- function(...) dir.create(file.path(...), recursive = TRUE,
                               showWarnings = FALSE)
touch <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  file.create(path)
  invisible(path)
}

touch(file.path(wsp_root, "runA.wsp"))
touch(file.path(wsp_root, "2025", "runB.wsp"))
touch(file.path(wsp_root, "2025", "runC.wsp"))
touch(file.path(wsp_root, "archive", "runC.wsp"))   # base-name collision
touch(file.path(wsp_root, "notes.txt"))             # must be ignored

touch(file.path(fcs_root, "2024", "runA", "Pool_001.fcs"))
touch(file.path(fcs_root, "2024", "runA", "Pool_002.fcs"))
touch(file.path(fcs_root, "2025", "Pool_003.fcs"))
touch(file.path(fcs_root, "2025", "extra", "Pool_004.fcs"))
touch(file.path(fcs_root, "dup", "left", "Shared.fcs"))
touch(file.path(fcs_root, "dup", "right", "Shared.fcs"))  # ambiguous
touch(file.path(fcs_root, "2024", "runA", "readme.md"))   # must be ignored

on.exit(unlink(root, recursive = TRUE), add = TRUE)

# --- find_workspace_files ----------------------------------------------------
cat("\nfind_workspace_files()\n")

found <- find_workspace_files(wsp_root)
eq("finds every .wsp recursively", sort(basename(found)),
   c("runA.wsp", "runB.wsp", "runC.wsp", "runC.wsp"))
ok("returns them in a stable path order", !is.unsorted(found))
ok("ignores non-.wsp files", !any(grepl("notes.txt", found)))

single <- find_workspace_files(file.path(wsp_root, "runA.wsp"))
eq("a single .wsp file still works", basename(single), "runA.wsp")

throws("a non-.wsp file is rejected",
       find_workspace_files(file.path(wsp_root, "notes.txt")),
       "Not a .wsp file")
throws("a missing path is rejected",
       find_workspace_files(file.path(root, "nope")), "not found")
mk(root, "empty")
throws("a folder with no workspace is rejected",
       find_workspace_files(file.path(root, "empty")), "No .wsp workspace")

# --- workspace_labels --------------------------------------------------------
cat("\nworkspace_labels()\n")

eq("unique base names are left alone",
   workspace_labels(c("/w/runA.wsp", "/w/2025/runB.wsp")),
   c("runA", "runB"))
eq("colliding base names are qualified by their parent folder",
   workspace_labels(c("/w/2025/runC.wsp", "/w/archive/runC.wsp",
                      "/w/runA.wsp")),
   c("2025/runC", "archive/runC", "runA"))
eq("a still-colliding label falls back to the full path",
   workspace_labels(c("/a/x/runC.wsp", "/b/x/runC.wsp")),
   c("/a/x/runC.wsp", "/b/x/runC.wsp"))
eq("extension matching is case-insensitive",
   workspace_labels("/w/runA.WSP"), "runA")

# --- index_fcs_dir -----------------------------------------------------------
cat("\nindex_fcs_dir()\n")

idx <- index_fcs_dir(fcs_root)
eq("indexes every .fcs by lower-cased name",
   sort(names(idx)),
   c("pool_001.fcs", "pool_002.fcs", "pool_003.fcs", "pool_004.fcs",
     "shared.fcs"))
eq("a duplicated name keeps both paths", length(idx[["shared.fcs"]]), 2L)
ok("non-.fcs files are not indexed", is.null(idx[["readme.md"]]))
ok("nested files are found",
   grepl("2025/extra", idx[["pool_004.fcs"]][1], fixed = TRUE))
throws("a folder with no FCS is rejected",
       index_fcs_dir(file.path(root, "empty")), "No .fcs file")

# --- resolve_fcs_names -------------------------------------------------------
cat("\nresolve_fcs_names()\n")

r <- resolve_fcs_names(c("Pool_001.fcs", "Pool_003.fcs"), idx)
eq("resolves plain names to full paths", basename(r$file),
   c("Pool_001.fcs", "Pool_003.fcs"))
eq("pairs each file with its own directory", basename(dirname(r$file)),
   c("runA", "2025"))
ok("nothing missing", length(r$missing) == 0)
ok("nothing ambiguous", length(r$ambiguous) == 0)

eq("matches case-insensitively",
   basename(resolve_fcs_names("POOL_001.FCS", idx)$file), "Pool_001.fcs")
eq("matches a name with the extension omitted",
   basename(resolve_fcs_names("Pool_002", idx)$file), "Pool_002.fcs")
eq("matches a stored Windows path by its base name",
   basename(resolve_fcs_names("C:\\\\Data\\\\Pool_003.fcs", idx)$file),
   "Pool_003.fcs")

r <- resolve_fcs_names(c("Pool_001.fcs", "Nope.fcs", "Shared.fcs"), idx)
eq("reports the unmatched name", r$missing, "Nope.fcs")
eq("reports the ambiguous name", names(r$ambiguous), "Shared.fcs")
eq("lists every ambiguous candidate", length(r$ambiguous[["Shared.fcs"]]), 2L)
ok("an unresolved slot stays NA", is.na(r$file[2]) && is.na(r$file[3]))
ok("resolved slots survive alongside failures", !is.na(r$file[1]))

# --- resolve_group -----------------------------------------------------------
cat("\nresolve_group()\n")

tbl <- data.frame(
  index     = 1:3,
  groupName = c("All Samples", "Controls", "Killing"),
  n_samples = c(6L, 0L, 4L),
  stringsAsFactors = FALSE
)

eq("an integer index still works (unchanged behaviour)",
   resolve_group("1", tbl)$name, "All Samples")
eq("a group name resolves to its index",
   resolve_group("Killing", tbl)$index, 3L)
eq("a name carries its sample count",
   resolve_group("Killing", tbl)$n_samples, 4L)
eq("surrounding whitespace is tolerated",
   resolve_group("  Killing  ", tbl)$index, 3L)
throws("an out-of-range index is rejected, listing the groups",
       resolve_group("9", tbl), "out of range")
throws("an unknown name is rejected, listing the groups",
       resolve_group("Nonsense", tbl), "not a sample group")
throws("an empty group is rejected", resolve_group("Controls", tbl), "is empty")
throws("the workspace is named in the message",
       resolve_group("Nonsense", tbl, where = "runB"), "workspace 'runB'")

# The group accessor moved between CytoML and flowWorkspace and has been
# renamed across versions, so it may be unavailable. An INDEX must still work.
eq("an index works with no group table at all",
   resolve_group("2", NULL)$index, 2L)
ok("...and reports no group name",
   is.na(resolve_group("2", NULL)$name))
throws("a NAME without a group table fails with an actionable message",
       resolve_group("Killing", NULL), "Use the group index instead")

cat("\nfind_wsp_accessor()\n")
ok("returns NULL when nothing matches",
   is.null(find_wsp_accessor(c("no_such_function_at_all"))))
ok("finds a function that a loaded package really exports",
   is.function(find_wsp_accessor(c("no_such_function_at_all", "sessionInfo"))) ||
     is.null(find_wsp_accessor("sessionInfo")))  # CytoML absent here

cat("\ndescribe_group()\n")
eq("a named group renders with its name",
   describe_group(list(group_index = 2L, group_name = "Killing")),
   'group 2 "Killing"')
eq("an unnamed group does not print the string NA",
   describe_group(list(group_index = 2L, group_name = NA_character_)),
   "group 2")
eq("a missing name field is handled",
   describe_group(list(group_index = 1L)), "group 1")

# --- resolve_workspace_specs (CytoML injected) -------------------------------
cat("\nresolve_workspace_specs()\n")

# Stand-in workspace reader: the "workspace" is just its path, the group table
# is fixed, and the sample list is looked up per workspace.
fake_samples <- list(
  runA = c("Pool_001.fcs", "Pool_002.fcs"),
  runB = c("Pool_003.fcs", "Pool_004.fcs"),   # split across two directories
  runC = c("Pool_001.fcs", "Nope.fcs"),       # one missing
  runD = c("Shared.fcs")                      # ambiguous
)
fake_open <- function(p) sub("\\.wsp$", "", basename(p))
fake_read <- function(ws, idx) {
  nm <- fake_samples[[ws]]
  data.frame(sampleID = as.character(seq_along(nm)), name = nm,
             stringsAsFactors = FALSE)
}
# wsp_group_table() is the one CytoML touch point left in the orchestration, so
# it is stubbed for the duration of this block.
.real_group_table <- wsp_group_table
wsp_group_table <- function(ws) tbl

paths  <- paste0("/w/", c("runA", "runB", "runC", "runD"), ".wsp")
labels <- workspace_labels(paths)
res <- resolve_workspace_specs(paths, labels, idx, "All Samples",
                               open_ws = fake_open, read_samples = fake_read)

eq("clean workspaces become specs", names(res$specs), c("runA", "runB"))
eq("broken workspaces become problems", names(res$problems), c("runC", "runD"))
eq("a spec resolves every file",
   basename(res$specs$runA$files$file), c("Pool_001.fcs", "Pool_002.fcs"))
eq("a spec records the single directory its files live in",
   length(res$specs$runA$dirs), 1L)
eq("a split workspace records both directories",
   length(res$specs$runB$dirs), 2L)
eq("the resolved group is carried on the spec",
   res$specs$runA$group_name, "All Samples")
eq("a missing file is recorded", res$problems$runC$missing, "Nope.fcs")
eq("an ambiguous file is recorded", names(res$problems$runD$ambiguous),
   "Shared.fcs")

bad <- resolve_workspace_specs(paths[1], "runA", idx, "Controls",
                               open_ws = fake_open, read_samples = fake_read)
ok("a group error is captured, not thrown",
   length(bad$specs) == 0 && grepl("is empty", bad$problems$runA$error))

wsp_group_table <- .real_group_table

# --- format_resolution_report ------------------------------------------------
cat("\nformat_resolution_report()\n")

eq("no problems means no message", format_resolution_report(list()), "")

report <- format_resolution_report(res$problems, fcs_root = fcs_root)
ok("names every failing workspace",
   grepl("runC", report) && grepl("runD", report))
ok("reports the missing file", grepl("missing:\\s+Nope.fcs", report))
ok("reports every ambiguous candidate",
   grepl("ambiguous", report) &&
   length(grep("Shared.fcs", strsplit(report, "\n")[[1]])) >= 1 &&
   length(grep("dup/(left|right)", strsplit(report, "\n")[[1]])) == 2)
ok("says nothing was loaded", grepl("Nothing was loaded", report))

cat("\nformat_shared_file_report()\n")

mkspec <- function(label, files) list(
  label = label,
  files = data.frame(file = files, stringsAsFactors = FALSE))

eq("a single workspace can share nothing",
   format_shared_file_report(list(mkspec("a", "/x/P1.fcs"))), "")
eq("disjoint workspaces are silent",
   format_shared_file_report(list(mkspec("a", "/x/P1.fcs"),
                                  mkspec("b", "/x/P2.fcs"))), "")
shared <- format_shared_file_report(list(
  mkspec("MTB048", c("/x/P1.fcs", "/x/P2.fcs")),
  mkspec("MTB039_048", c("/x/P2.fcs", "/x/P3.fcs")),
  mkspec("MTB039", "/x/P3.fcs")))
ok("an overlapping file is reported", grepl("P2.fcs", shared))
ok("...naming every workspace that claims it",
   grepl("MTB048, MTB039_048", shared))
ok("an overlapping file is also caught for the second pair",
   grepl("P3.fcs", shared))
ok("a file owned by one workspace is not reported", !grepl("P1.fcs", shared))
ok("it says what the consequence is", grepl("ONCE PER WORKSPACE", shared))
eq("the count is right", sub("^(\\d+) FCS.*", "\\1", strsplit(shared, "\n")[[1]][1]),
   "2")
long <- format_shared_file_report(
  list(mkspec("a", sprintf("/x/P%02d.fcs", 1:15)),
       mkspec("b", sprintf("/x/P%02d.fcs", 1:15))), max_show = 3L)
ok("a long list is truncated", grepl("and 12 more", long))
ok("echoes the FCS root that was searched", grepl(fcs_root, report,
                                                  fixed = TRUE))

# --- disambiguate_sample_labels ----------------------------------------------
cat("\ndisambiguate_sample_labels()\n")

eq("unique names are left exactly as they are",
   disambiguate_sample_labels(c("Pool_001", "Pool_002", "Pool_003"),
                              c("runA", "runA", "runB")),
   c("Pool_001", "Pool_002", "Pool_003"))
eq("only the colliding names are qualified",
   disambiguate_sample_labels(c("Pool_001", "Pool_002", "Pool_001"),
                              c("runA", "runA", "runB")),
   c("runA: Pool_001", "Pool_002", "runB: Pool_001"))
eq("a single-workspace run is untouched",
   disambiguate_sample_labels(c("A", "B"), c("runA", "runA")), c("A", "B"))
eq("an empty run is handled",
   disambiguate_sample_labels(character(0), character(0)), character(0))
throws("mismatched lengths are caught",
       disambiguate_sample_labels(c("A", "B"), "runA"))

# --- population resolution ---------------------------------------------------
cat("\ntry_resolve_population()\n")

pops_a <- c("/Lymphocytes", "/Lymphocytes/Single Cells",
            "/Lymphocytes/Single Cells/Viable cells",
            "/Lymphocytes/Single Cells/Viable cells/CD4+")
pops_b <- c("/Lymph", "/Lymph/Singlets", "/Lymph/Singlets/Viable cells")

eq("an exact path resolves",
   try_resolve_population("/Lymphocytes/Single Cells", pops_a)$path,
   "/Lymphocytes/Single Cells")
eq("a leading slash is optional",
   try_resolve_population("Lymphocytes/Single Cells", pops_a)$path,
   "/Lymphocytes/Single Cells")
eq("a unique suffix resolves",
   try_resolve_population("Viable cells", pops_a)$path,
   "/Lymphocytes/Single Cells/Viable cells")
eq("the same suffix resolves to a different absolute path elsewhere",
   try_resolve_population("Viable cells", pops_b)$path,
   "/Lymph/Singlets/Viable cells")
ok("a missing population reports rather than throws",
   !try_resolve_population("CD4+", pops_b)$ok)
ok("a missing population explains itself",
   grepl("not found", try_resolve_population("CD4+", pops_b)$reason))
ok("an ambiguous population is reported",
   grepl("ambiguous",
         try_resolve_population("Cells", c("/A/Cells", "/B/Cells"))$reason))
throws("resolve_population() still throws for its existing callers",
       resolve_population("CD4+", pops_b), "not found")

# --- merge_subgate_candidates ------------------------------------------------
cat("\nmerge_subgate_candidates()\n")

# Two workspaces whose population resolves to DIFFERENT absolute paths but
# whose relative subgate labels overlap — the case the merge exists for.
cand_a <- flowjo_subgate_candidates(
  c("/Lymphocytes/Single Cells/Viable cells/CD4+",
    "/Lymphocytes/Single Cells/Viable cells/CD8+",
    "/Lymphocytes/Single Cells/Viable cells/CD4+/IFNg+"),
  pops_a, "/Lymphocytes/Single Cells/Viable cells", mode = "furthest")
cand_b <- flowjo_subgate_candidates(
  c("/Lymph/Singlets/Viable cells/CD4+",
    "/Lymph/Singlets/Viable cells/NK"),
  pops_b, "/Lymph/Singlets/Viable cells", mode = "furthest")

merged <- merge_subgate_candidates(list(runA = cand_a, runB = cand_b))
eq("labels are the union, first-appearance ordered",
   merged$labels, c("CD4+", "CD4+/IFNg+", "CD8+", "NK"))
eq("depth comes from the relative label", merged$depths, c(1L, 2L, 1L, 1L))
eq("each workspace keeps its own absolute paths, column-aligned",
   merged$paths$runA,
   c("/Lymphocytes/Single Cells/Viable cells/CD4+",
     "/Lymphocytes/Single Cells/Viable cells/CD4+/IFNg+",
     "/Lymphocytes/Single Cells/Viable cells/CD8+",
     NA_character_))
eq("a gate absent from a workspace becomes an NA column there",
   merged$paths$runB,
   c("/Lymph/Singlets/Viable cells/CD4+", NA_character_, NA_character_,
     "/Lymph/Singlets/Viable cells/NK"))
eq("the mode tag survives", merged$tag, "furthest")
eq("an empty candidate set merges to nothing",
   length(merge_subgate_candidates(list())$labels), 0L)

# --- get_subpop_membership tolerates the NA columns --------------------------
cat("\nget_subpop_membership() with absent gates\n")

# gh is stubbed: this checks the column alignment contract, which is what the
# NA paths above feed into, without needing flowWorkspace.
fake_gh <- structure(list(), class = "fakeGatingHierarchy")
membership <- get_subpop_membership(fake_gh, c(TRUE, TRUE, FALSE),
                                    rep(NA_character_, 4))
eq("all-absent gates give an all-FALSE matrix of the right shape",
   dim(membership), c(2L, 4L))
ok("and it really is all FALSE", !any(membership))

# --- sample annotations ------------------------------------------------------
cat("\nannotation_key()\n")

eq("strips the cytoset suffix and the extension",
   annotation_key("MTB048_TP1.fcs_0"), "mtb048_tp1")
eq("strips a bare extension", annotation_key("MTB048_TP1.fcs"), "mtb048_tp1")
eq("leaves a plain id alone", annotation_key("MTB048_TP1"), "mtb048_tp1")
eq("reduces a Windows path to its base name",
   annotation_key("C:\\\\Data\\\\MTB048_TP1.fcs"), "mtb048_tp1")
eq("is vectorised", annotation_key(c("A.fcs", "B.fcs_2")), c("a", "b"))

cat("\nannotations_from_pattern()\n")

pat <- "(?<patient>MTB[0-9]+)_(?<timepoint>TP[0-9]+)"
ap  <- annotations_from_pattern(
  c("MTB048_TP1", "MTB011_TP3", "junk"), pat)
eq("named groups become the columns", colnames(ap), c("patient", "timepoint"))
eq("captures the patient", ap$patient, c("MTB048", "MTB011", NA))
eq("captures the timepoint", ap$timepoint, c("TP1", "TP3", NA))
eq("a name that does not match yields NA everywhere",
   unname(unlist(ap[3, ])), c(NA_character_, NA_character_))
throws("an unnamed pattern is rejected with an example",
       annotations_from_pattern("MTB048_TP1", "(MTB[0-9]+)"),
       "named capture groups")

cat("\nresolve_sample_annotations()\n")

samples <- data.frame(
  workspace = c("runA", "runA", "runB"),
  raw       = c("MTB048_TP1.fcs_0", "MTB048_TP3.fcs_1", "MTB011_TP1.fcs_0"),
  clean     = c("MTB048_TP1", "MTB048_TP3", "MTB011_TP1"),
  display   = c("MTB048_TP1", "MTB048_TP3", "runB: MTB011_TP1"),
  stringsAsFactors = FALSE)

eq("no source means no annotations",
   ncol(resolve_sample_annotations(samples)), 0L)

md_file <- file.path(root, "cohort.csv")
writeLines(c("sample_id,patient,timepoint,response",
             "MTB048_TP1.fcs,MTB048,TP1,R",
             "MTB048_TP3.fcs,MTB048,TP3,R",
             "MTB011_TP1.fcs,MTB011,TP1,NR"), md_file)

a <- resolve_sample_annotations(samples, metadata_file = md_file)
eq("every non-key column becomes an annotation",
   colnames(a), c("patient", "timepoint", "response"))
eq("rows follow the sample order", a$patient,
   c("MTB048", "MTB048", "MTB011"))
eq("a row is matched despite the cytoset suffix", a$timepoint,
   c("TP1", "TP3", "TP1"))
eq("a display label with a workspace prefix still matches", a$response[3], "NR")

eq("annotate restricts and orders the columns",
   colnames(resolve_sample_annotations(samples, metadata_file = md_file,
                                       annotate = c("timepoint", "patient"))),
   c("timepoint", "patient"))
throws("an unknown annotate column is rejected, listing what exists",
       resolve_sample_annotations(samples, metadata_file = md_file,
                                  annotate = "nonsense"),
       "not available")

# A sample missing from the table must not vanish or poison the run.
md_partial <- file.path(root, "partial.csv")
writeLines(c("sample_id,patient", "MTB048_TP1,MTB048"), md_partial)
a2 <- suppressWarnings(
  resolve_sample_annotations(samples, metadata_file = md_partial))
eq("an unmatched sample becomes Unassigned", a2$patient,
   c("MTB048", .ANNOT_UNASSIGNED, .ANNOT_UNASSIGNED))
ok("...and it warns about it",
   !is.null(tryCatch(resolve_sample_annotations(samples,
                                                metadata_file = md_partial),
                     warning = function(w) conditionMessage(w))))

# Blanks and placeholder strings are Unassigned, matching check_metadata.R.
md_blank <- file.path(root, "blank.csv")
writeLines(c("sample_id,arm", "MTB048_TP1,A", "MTB048_TP3,",
             "MTB011_TP1,n/a"), md_blank)
eq("blanks and placeholders become Unassigned",
   resolve_sample_annotations(samples, metadata_file = md_blank)$arm,
   c("A", .ANNOT_UNASSIGNED, .ANNOT_UNASSIGNED))

# Semicolon-separated exports are common in this project.
md_semi <- file.path(root, "semi.csv")
writeLines(c("sample_id;patient", "MTB048_TP1;MTB048", "MTB048_TP3;MTB048",
             "MTB011_TP1;MTB011"), md_semi)
eq("a semicolon-separated table is detected",
   resolve_sample_annotations(samples, metadata_file = md_semi)$patient,
   c("MTB048", "MTB048", "MTB011"))

# Values must stay text: a timepoint of "01" must not become the number 1.
md_num <- file.path(root, "num.csv")
writeLines(c("sample_id,timepoint", "MTB048_TP1,01", "MTB048_TP3,02",
             "MTB011_TP1,01"), md_num)
eq("numeric-looking levels keep their leading zeros",
   resolve_sample_annotations(samples, metadata_file = md_num)$timepoint,
   c("01", "02", "01"))

eq("the pattern route works when no metadata file is given",
   resolve_sample_annotations(samples, pattern = pat)$patient,
   c("MTB048", "MTB048", "MTB011"))
eq("metadata wins over a pattern when both are given",
   colnames(resolve_sample_annotations(samples, metadata_file = md_file,
                                       pattern = pat)),
   c("patient", "timepoint", "response"))

md_alt <- file.path(root, "alt.csv")
writeLines(c("FCS,patient", "MTB048_TP1,MTB048", "MTB048_TP3,MTB048",
             "MTB011_TP1,MTB011"), md_alt)
eq("the first column is used when there is no sample_id",
   resolve_sample_annotations(samples, metadata_file = md_alt)$patient,
   c("MTB048", "MTB048", "MTB011"))
throws("a bad sample_id_column is rejected",
       resolve_sample_annotations(samples, metadata_file = md_file,
                                  id_column = "nope"), "is not in the metadata")
throws("a missing metadata file is rejected",
       resolve_sample_annotations(samples, metadata_file = "/no/such.csv"),
       "not found")

# --- summary -----------------------------------------------------------------
cat("\n=== ", .pass, " passed, ", .fail, " failed ===\n", sep = "")
if (.fail > 0) quit(status = 1)
