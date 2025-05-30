library("irr")
library(data.table)

# Required Libraries: flowWorkspace, dplyr, openxlsx
# !For RStudio only: Set working directory to the same path as the current file
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

#================================ Parameters ===================================
# Output directory
output_dir <- paste0(getwd(), "/venn_diagrams")

# Create a directory to save the Venn diagrams
dir.create(output_dir, showWarnings = FALSE)

# Load the label lists into data #data <- label_lists
data <- data.frame(fread(paste0(output_dir,
                                "/label_listsICK-043-02_Unmixed.fcs.csv")))

data_subsets <- list(data[, c(1:3, 4)],
                     data[, c(1:3, 5)], data, data[, 1:2], data[, 1:3],
                     data[, 1:4], data[, 1:5], data[, 2:3], data[, 2:4],
                     data[, 2:5], data[, 3:4], data[, 3:5], data[, 4:5])
subset_labels <- c("FlowJo vs Broad", "FlowJo vs Narrow", "Full data", "LJ",
                   "LS", "LB", "LN", "JS", "JB", "JN", "SB", "SN", "BN")

#================================= Main code ===================================

summary(data)

log_file <- file.path(output_dir, "fleiss_kappa_log.txt")
log_con <- file(log_file, open = "wt")

fleiss_kappa <- list()
for (i in seq_along(data_subsets)){
  step_start <- Sys.time()
  cat("====================================\n", file = log_con, append = TRUE)
  cat(paste0("Analysis ", i, ": ", subset_labels[[i]], "\n"), file = log_con,
      append = TRUE)
  print(i)
  print(paste0("Analysing ", subset_labels[[i]]))
  fleiss_kappa[[i]] <- irr::kappam.fleiss(data_subsets[[i]])
  print(fleiss_kappa[[i]])
  capture.output(fleiss_kappa[[i]], file = log_con, append = TRUE)
  step_end <- Sys.time()
  step_elapsed <- step_end - step_start
  cat(sprintf("Elapsed time for this step: %s\n", step_elapsed), file = log_con,
      append = TRUE)
}
close(log_con)
