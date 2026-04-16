# Ensure src scripts can locate R_shared regardless of how the test is invoked.
# The src scripts check the CYTOMETRY_SRC env var first when resolving shared_dir.
local({
  d <- getwd()
  while (nchar(d) > 1) {
    if (file.exists(file.path(d, "R_shared", "io_utils.R"))) {
      Sys.setenv(CYTOMETRY_SRC = normalizePath(d, mustWork = FALSE))
      return(invisible())
    }
    d <- dirname(d)
  }
  warning("helper-setup.R could not locate an R_shared sibling; src scripts may fail")
})
