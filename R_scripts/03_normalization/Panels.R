# Create an empty list to store marker names and channel names from each flow frame
marker_channel_list <- list()

# Iterate over each flow frame to extract marker names and channel names
for (i in seq_along(flow_frames)) {
  marker_names <- markernames(flow_frames[[i]])
  marker_channel_list[[i]] <- data.frame(FlowFrame_ID = i,
                                         Marker_Names = marker_names)
}

# Combine marker names and channel names from all flow frames into a single table
marker_channel_table <- do.call(rbind, marker_channel_list)

# Export the table as a CSV file
write.csv(marker_channel_table, "marker_channel_table.csv", row.names = FALSE)

table <- read.csv("/home/ubuntu/CellCnn/Data/Liver/all/marker_final.csv", header = TRUE, row.names = 1)
colnames(table)<- sample_names

list_of_vectors <- as.list(table)

# Create a hash table to store vectors and their corresponding indices
hash_table <- list()

# Iterate through each vector in the list
for (i in 1:length(list_of_vectors)) {
  # Calculate hash value for the vector
  hash_value <- paste(list_of_vectors[[i]], collapse = ",")
  
  # Check if the hash value already exists in the hash table
  if (hash_value %in% names(hash_table)) {
    # If the hash value exists, append the index to the corresponding list
    hash_table[[hash_value]] <- c(hash_table[[hash_value]], i)
  } else {
    # If the hash value does not exist, create a new list with the index
    hash_table[[hash_value]] <- list(i)
  }
}

# Filter out unique hash values
unique_hashes <- hash_table[names(hash_table) != ""]

# Print the indices of vectors that are equal to each other
for (hash_value in unique_hashes) {
  print(paste0("Equal vectors:", hash_value, "\n"))
}

df_clean <- as.data.frame(lapply(ss_table, function(x) gsub("[ -]", "", x)))

# Function to compare pairwise columns
compare_columns <- function(df) {
  n <- ncol(df)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      result <- df[, i] == df[, j]
      if (!all(result)) {
        print(paste("Columns", names(df)[i], "and", names(df)[j], "differ:"))
        print(df[!result, c(names(df)[i], names(df)[j])])
      }
    }
  }
}

compare_columns(df_clean)
