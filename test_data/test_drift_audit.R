#!/usr/bin/env Rscript
# test_data/test_drift_audit.R — unit tests for src/harmonization/drift/drift_audit.R (acquisition audit). Base R.
#   Rscript test_data/test_drift_audit.R
cat("=== drift acquisition audit tests ===\n\n")
.ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(.ofile)) .ofile <- "test_data/test_drift_audit.R"
base_dir <- normalizePath(file.path(dirname(.ofile), ".."))
if (!dir.exists(file.path(base_dir, "src"))) base_dir <- getwd()
source(file.path(base_dir, "src", "harmonization", "drift", "drift_audit.R"))
n_pass <- 0L; n_fail <- 0L
check <- function(label, ok) { if (isTRUE(ok)) { n_pass <<- n_pass + 1L; cat("  ok  ", label, "\n") } else { n_fail <<- n_fail + 1L; cat("  FAIL", label, "\n") } }

kw <- data.frame(filename = rep(c("a.fcs", "b.fcs"), each = 6),
                 key = rep(c("$P1N", "$P1V", "$P2N", "$P2V", "SPILL", "CST PERFORMANCE EXPIRED"), 2),
                 value = c("PE-A", "520", "FITC-A", "500", "m1", "x", "PE-A", "552", "FITC-A", "500", "m2", NA), stringsAsFactors = FALSE)
kw <- kw[!is.na(kw$value), ]
pm <- data.frame(filename = rep(c("a.fcs", "b.fcs"), each = 2), channel = rep(c("PE-A", "FITC-A"), 2), canonical = rep(c("IFNg", "CD107a"), 2))
acq <- acquisition_settings(kw, pm)
check("voltages per canonical marker, voltage set, compensation id, CS&T expired flag",
      identical(acq$V_IFNg, c(520, 552)) && identical(acq$V_CD107a, c(500, 500)) && length(unique(acq$voltage_set)) == 2 &&
        length(unique(acq$compensation_id)) == 2 && identical(acq$cst_expired, c(TRUE, FALSE)))
gc <- data.frame(filename = c("a.fcs", "a.fcs", "b.fcs"), Population = c("/L/Single Cells", "/L/Single Cells/Live cells", "/L/Single Cells/Live cells"),
                 Parent = c("/L", "/L/Single Cells", "/L/Single Cells"), Count = c(900, 720, 300), ParentCount = c(1000, 900, 1000))
q <- sample_quality(gc, "Live cells")
check("viability = live / parent of the live gate, live-cell count", isTRUE(all.equal(q$viability, c(.8, .3))) && identical(q$live_cells, c(720, 300)))
d <- data.frame(day = rep(c("d1", "d2", "d3"), each = 4), op = rep(c("A", "A", "B"), each = 4), pat = rep(1:4, 3))
M <- determination_matrix(d, c(day = "day", operator = "op", patient = "pat"))
check("day determines operator fully; operator does not determine day; patient independent of day",
      isTRUE(all.equal(M["day", "operator"], 1)) && M["operator", "day"] < 1 && abs(M["patient", "day"]) < 1e-8)
check("eta squared", isTRUE(all.equal(eta_squared(c(1, 1, 3, 3), c("a", "a", "b", "b")), 1)) && abs(eta_squared(c(1, 3, 1, 3), c("a", "a", "b", "b"))) < 1e-8)
x <- tempfile(fileext = ".wsp")
writeLines(c('<Workspace><SampleList><Sample><DataSet uri="file:/x/MTB%20a.fcs" /><Transformations>',
             '<transforms:biex transforms:length="256" transforms:neg="0" transforms:width="-100" transforms:pos="4.4" >',
             '<data-type:parameter data-type:name="Comp-PE-A" /></transforms:biex>',
             '<transforms:linear transforms:minRange="0" transforms:maxRange="262144" >',
             '<data-type:parameter data-type:name="FSC-A" /></transforms:linear></Transformations></Sample></SampleList></Workspace>'), x)
ft <- flowjo_transformations(x)
check("FlowJo transformations read per sample and parameter from the workspace XML",
      nrow(ft) == 2 && ft$filename[1] == "MTB a.fcs" && ft$parameter[ft$transform == "biex"] == "PE-A" && grepl("pos=4.4", ft$settings[ft$transform == "biex"]))
cat(sprintf("\n%d passed, %d failed\n", n_pass, n_fail))
if (n_fail > 0) quit(status = 1)
