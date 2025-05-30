# Required Libraries
library(ggplot2)
library(venn)
library(VennDiagram)
library(grid)
library(eulerr)
library(venneuler)
library(rJava)

# Example Input: Two lists of the same size
# Generate two lists of 35,000 entries each
set.seed(123)  # Set seed for reproducibility

source("Functions.R")

#================== Parameters ===================

# Taget FCS filename
target_fcs <- "ICK-043-02_Unmixed.fcs"

# Custom name for the workspace (keep identifier as short as possible)
custom_name <- c("Lea", "Jonas", "Samu", "Broad", "Narrow")

# User provided dictionary for translating labels
# This dictionary should map the original labels to the desired labels
label_dict <- c(
  "root" = "Removed", #"Erys"

  "/no erys" = "Removed", #"Doublets"
  "/No Erythrocytes" = "Removed", #"Doublets"
  "/No_erys" = "Removed", #"Doublets"

  "/no erys/Single Cells" = "Removed", #"Myelos"
  "/No Erythrocytes/Single Cells" = "Removed", #"Myelos"
  "/No_erys/singlets"  = "Removed", #"Myelos"

  "/no erys/Single Cells/lymphocytes" = "Removed", #"ncLymphos"
  "/No Erythrocytes/Single Cells/Lymphocytes"  = "Removed", #"ncLymphos"
  "/No_erys/singlets/lymphocytes"  = "Removed", #"ncLymphos"

  "/no erys/Single Cells/lymphocytes/CD3+ TCRyd-" = "ncT cells",
  "/no erys/Single Cells/lymphocytes/CD3+ TCRyd-/CD3+ CD56-"  = "ncT cells",
  "/No Erythrocytes/Single Cells/Lymphocytes/CD3+gdTCR-" = "ncT cells",
  "/No Erythrocytes/Single Cells/Lymphocytes/CD3+gdTCR-/CD3+CD56-" =
    "ncT cells",
  "/No Erythrocytes/Single Cells/Lymphocytes/CD3+ TCRgd-" = "ncT cells",
  "/No_erys/singlets/lymphocytes/CD3+/gd-/TCs (CD3+CD56-)" = "ncT cells",
  "/No_erys/singlets/lymphocytes/CD3+/gd-" = "ncT cells",
  "/No_erys/singlets/lymphocytes/CD3+" = "ncT cells",
  "/No Erythrocytes/Single Cells/Lymphocytes/CD3+ TCRgd-/CD3+ CD56-" = "ncT cells",

  "/no erys/Single Cells/lymphocytes/CD3+ TCRyd-/CD3+ CD56-/CD4 T cells" =
    "CD4T cells",
  "/No Erythrocytes/Single Cells/Lymphocytes/CD3+gdTCR-/CD3+CD56-/CD4+" =
    "CD4T cells",
  "/No_erys/singlets/lymphocytes/CD3+/gd-/TCs (CD3+CD56-)/CD4TCs" = "CD4T cells",
  "/No Erythrocytes/Single Cells/Lymphocytes/CD3+ TCRgd-/CD3+ CD56-/CD4+" =
    "CD4T cells",

  "/no erys/Single Cells/lymphocytes/CD3+ TCRyd-/CD3+ CD56-/CD8 Tcells" =
    "CD8T cells",
  "/No Erythrocytes/Single Cells/Lymphocytes/CD3+gdTCR-/CD3+CD56-/CD8+" =
    "CD8T cells",
  "/No_erys/singlets/lymphocytes/CD3+/gd-/TCs (CD3+CD56-)/CD8TCs" = "CD8T cells",
  "/No Erythrocytes/Single Cells/Lymphocytes/CD3+ TCRgd-/CD3+ CD56-/CD8+" =
    "CD8T cells",

  "/no erys/Single Cells/lymphocytes/CD3+ TCRyd-/NKTcells" = "NKT cells",
  "/No Erythrocytes/Single Cells/Lymphocytes/CD3+gdTCR-/CD3+CD56+" = "NKT cells",
  "/No_erys/singlets/lymphocytes/CD3+/gd-/NKTCs" = "NKT cells",
  "/No Erythrocytes/Single Cells/Lymphocytes/CD3+ TCRgd-/CD3+ CD56+ " = "NKT cells",

  "/no erys/Single Cells/lymphocytes/yd Tcells" = "gdT cells",
  "/No Erythrocytes/Single Cells/Lymphocytes/CD3+ TCRgd+" = "gdT cells",
  "/No_erys/singlets/lymphocytes/CD3+/gdTCs" = "gdT cells",
  "/No_erys/singlets/lymphocytes/gdTCs" = "gdT cells",

  "/no erys/Single Cells/lymphocytes/CD3-" = "ncCD3-",
  "/No Erythrocytes/Single Cells/Lymphocytes/CD3-" = "ncCD3-",
  "/No_erys/singlets/lymphocytes/CD3-" = "ncCD3-",
  "/No_erys/singlets/lymphocytes/CD3-/CD19-CD56-" = "ncCD3-",
  "/No_erys/singlets/lymphocytes/CD3-/CD19-CD56-/CD3-CD56-CD19-CD11b-" =
    "ncCD3-",
  "/No_erys/singlets/lymphocytes/CD3-/CD19-CD56-/DCs" = "ncCD3-",

  "/no erys/Single Cells/lymphocytes/CD3-/B cells" = "B cells",
  "/No Erythrocytes/Single Cells/Lymphocytes/CD3-/CD19+ CD56-" = "B cells",
  "/No_erys/singlets/lymphocytes/CD3-/BCs" = "B cells",

  "/no erys/Single Cells/lymphocytes/CD3-/NK cells" = "NK cells",
  "/No Erythrocytes/Single Cells/Lymphocytes/CD3-/CD19- CD56+" = "NK cells",
  "/No_erys/singlets/lymphocytes/CD3-/NKCs" = "NK cells"
)
# Output directory
output_dir <- paste0(getwd(), "/venn_diagrams")

# Create a directory to save the Venn diagrams
dir.create(output_dir, showWarnings = FALSE)

# ========================= Load flowjo population lists ===========================

# Check whether the target FCS is present in all of the FlowJo workspaces
missing_fcs <- sapply(pop_strings_list, function(x) {
  !any(grepl(target_fcs, names(x), fixed = TRUE))
})
# If any workspaces are missing target FCS, stop the script and print a message
if (any(missing_fcs)) {
  stop(
    sprintf(
      "The following workspaces are missing the target FCS '%s': %s",
      target_fcs,
      paste(names(pop_strings_list)[missing_fcs], collapse = ", ")
    )
  )
} else {
  message(sprintf("Target FCS '%s' found in all workspaces.", target_fcs))
}
# Set of label lists
fj_label_lists <- sapply(pop_strings_list, function(x) {
  # Find the name in x that contains target_fcs as a substring
  idx <- grep(target_fcs, names(x), fixed = TRUE)
  if (length(idx) == 1) {
    x[[idx]]
  } else if (length(idx) > 1) {
    warning(sprintf("Multiple matches for '%s' in names(x), using the first.", 
                    target_fcs))
    x[[idx[1]]]
  } else {
    stop(sprintf("No match found for '%s' in names(x)", target_fcs))
  }
})
# Translate the labels in the lists
for (i in seq_along(fj_label_lists)) {
  fj_label_lists[[i]] <- translate_labels(fj_label_lists[[i]], label_dict)
}

# Save fj_label_lists into CSV files using fwriter for efficiency
fj_label_lists_df <- as.data.frame(fj_label_lists, stringsAsFactors = FALSE)
fwrite(fj_label_lists_df, file = file.path(output_dir, paste0("fj_label_lists",
                                                              target_fcs,
                                                              ".csv")))

# ========================= Load csv string lists ==============================

# Normalize the target FCS base name (remove .fcs and replace hyphens/underscores)
target_base <- gsub("\\.fcs$", "", target_fcs)
target_base_norm <- gsub("[-_]", "", target_base)

# Normalize the names in csv_label_lists
csv_names_norm <- gsub("[-_]", "", names(idx_labels))

# Find the index of the matching label
match_idx <- which(grepl(target_base_norm, csv_names_norm, fixed = TRUE))
if (length(match_idx) == 0) {
  stop(sprintf("No match found for FCS file for CSV '%s'", target_base))
} else if (length(match_idx) > 1) {
  message(sprintf("Found %d CSVs matching FCS '%s'. Using all matches.", length(match_idx), target_fcs))
  cy_label_lists <- lapply(match_idx, function(idx) {
    tmp <- idx_labels[[idx]]
    tmp[tmp == ""] <- "Removed"
    tmp
  })
  # Combine with fj_label_lists
  label_lists <- cbind(fj_label_lists[, c(1, 2, 3)], do.call(cbind, cy_label_lists))
} else {
  cy_label_lists <- idx_labels[[match_idx]]
  cy_label_lists[cy_label_lists == ""] <- "Removed"
  label_lists <- cbind(fj_label_lists[, c(1, 2, 3)], cy_label_lists)
}

# Save fj_label_lists into CSV files using fwriter for efficiency
cy_label_lists_df <- as.data.frame(cy_label_lists, stringsAsFactors = FALSE)
fwrite(cy_label_lists_df, file = file.path(output_dir, paste0("cy_label_lists",
                                                              target_fcs,
                                                              ".csv")))

# Save fj_label_lists into CSV files using fwriter for efficiency
label_lists_df <- as.data.frame(label_lists, stringsAsFactors = FALSE)
fwrite(label_lists_df, file = file.path(output_dir, paste0("label_lists",
                                                           target_fcs,
                                                           ".csv")))

# ============================= Main Code ======================================

# Ensure all lists are of the same size
list_lengths <- sapply(label_lists, length)
if (length(unique(list_lengths)) != 1) {
  stop("All lists must have the same size.")
}
unique_combinations <- unique(label_lists)
print(unique_combinations)

# Get all unique strings across all lists (concatenate, don't unlist)
unique_populations <- unique(as.vector(unlist(label_lists)))

# Generate a Venn diagram for each unique string
for (population in unique_populations) {
  print(paste("Plotting", population))
  if (is.na(population) || population == "") {
    warning("Skipping empty or NA population for Venn diagram.")
    next
  }
  plot_and_save_venn(lists = label_lists, string = population, 
                     output_dir = output_dir, scale_size = TRUE,
                     custom_name = custom_name)
}

# Print completion message
print(paste("Venn diagrams saved to:", output_dir))
