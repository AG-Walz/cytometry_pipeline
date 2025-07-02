library(ggplot2)
library(data.table)
# !For RStudio only: Set working directory to the same path as the current file
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Output directory
output_dir <- paste0(getwd(), "/venn_diagrams")

# Taget FCS filename
target_fcs <- "ICK-043-02_Unmixed.fcs"

# Custom name for the workspace (keep identifier as short as possible)
custom_name <- c("Lea", "Jonas", "Samu", "Broad", "Narrow")

data <- label_lists_df
data <- data.frame(fread(paste0(output_dir,
                                "/label_lists", target_fcs,".csv")))

bin_size <- 1000
n_bins <- ceiling(nrow(data) / bin_size)

# Collect all unique labels across all columns
all_labels <- unique(unlist(lapply(data, unique)))

# Define a color palette for all labels (grey for Removed, rest from Set2)
label_palette <- setNames(
  c(RColorBrewer::brewer.pal(length(all_labels) - 1, "Spectral"),"#d3d3d3"),
  c(setdiff(all_labels, "Removed"), "Removed")
)

for (wsp in seq_along(colnames(data))) {
  col <-colnames(data)[[wsp]]
  label_counts <- data.frame()
  removed_bins <- c()
  for (i in seq_len(n_bins)) {
    start_idx <- (i - 1) * bin_size + 1
    end_idx <- min(i * bin_size, nrow(data))
    bin_data <- data[start_idx:end_idx, col]
    bin_table <- as.data.frame(table(bin_data), stringsAsFactors = FALSE)
    bin_table$bin <- i
    label_counts <- rbind(label_counts, bin_table)
    if (all(bin_data == "Removed")) {
      removed_bins <- c(removed_bins, i)
    }
  }

  # Plot histogram
  p <- ggplot(label_counts, aes(x = factor(bin), y = Freq, fill = bin_data)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(title = paste("Label distribution for", custom_name[wsp], "in file", 
                       target_fcs),
         x = paste0("Bins ordered by time as in FCS-file (", bin_size, 
                    " cells each)"), y = "Cell Count",
         fill = "Label") +
    theme_minimal()+
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    scale_fill_manual(values = label_palette, na.value = "black")
  # Annotate bins where all are "Removed"
  if (length(removed_bins) > 0) {
    p <- p + geom_text(
      data = subset(label_counts, 
                    bin %in% removed_bins & bin_data == "Removed"),
      aes(label = "R"), vjust = -0.5, color = "darkgrey", size = 3
    )
  }
  print(p)
  # Optionally save the plot
  ggsave(file.path(output_dir, paste0("hist_", custom_name[wsp], "_", target_fcs, ".png")),
         plot = p, width = 8, height = 4)
}

all_removed_segments <- list()

for (col in colnames(data)) {
  removed_bins <- c()
  for (i in seq_len(n_bins)) {
    start_idx <- (i - 1) * bin_size + 1
    end_idx <- min(i * bin_size, nrow(data))
    bin_data <- data[start_idx:end_idx, col]
    if (all(bin_data == "Removed")) {
      removed_bins <- c(removed_bins, i)
    }
  }
  # Convert bin numbers to row ranges and merge consecutive bins
  if (length(removed_bins) > 0) {
    ranges <- list()
    current_start <- removed_bins[1]
    current_end <- removed_bins[1]
    for (j in 2:length(removed_bins)) {
      if (removed_bins[j] == current_end + 1) {
        current_end <- removed_bins[j]
      } else {
        ranges[[length(ranges) + 1]] <- c(current_start, current_end)
        current_start <- removed_bins[j]
        current_end <- removed_bins[j]
      }
    }
    ranges[[length(ranges) + 1]] <- c(current_start, current_end)
    # Print row ranges
    cat(sprintf("Column '%s':\n", col))
    for (r in ranges) {
      row_start <- (r[1] - 1) * bin_size + 1
      row_end <- min(r[2] * bin_size, nrow(data))
      cat(sprintf("  All 'Removed' from row %d to %d\n", row_start, row_end))
    }
    all_removed_segments[[col]] <- ranges
  }
}

