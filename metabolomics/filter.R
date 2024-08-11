library(KEGGREST)
library(dplyr)
library(yaml)

# Load parameters from config file
config <- yaml::read_yaml("/users/2875659p/sharedscratch/rna_seq/metabolomics/meta/config_2.yaml")

# Read the input files
if (config$mode == "positive") {
  input_csv <- config$input_csv_positive
  peaks_csv <- config$pos_peaks_csv
} else if (config$mode == "negative") {
  input_csv <- config$input_csv_negative
  peaks_csv <- config$neg_peaks_csv
} else {
  stop("Unknown ionization mode specified.")
}

data <- read.csv(file.path(config$results, input_csv), header = TRUE)
peak_table <- read.csv(file.path(config$results, peaks_csv))

# Define a function to get the compound name for a given KEGG ID
get_compound_name <- function(kegg_id) {
  compound_info <- tryCatch({
    keggGet(kegg_id)
  }, error = function(e) {
    return(NULL)  # Return NULL if there is an error
  })

  if (!is.null(compound_info) && length(compound_info) > 0) {
    return(compound_info[[1]]$NAME[1])  # Return the first name
  } else {
    return(NA)  # Return NA if no information is found or there's an error
  }
}

# Extract KEGG IDs from the specified column and clean them
kegg_ids <- trimws(as.character(data[[config$column_index]]))

# Apply the function to the list of KEGG IDs
compound_names <- sapply(kegg_ids, get_compound_name)

# Add the compound names to the original data frame
data$Compound_Name <- compound_names

# Remove rows with NA in Compound_Name
data <- data %>%
  filter(!is.na(Compound_Name))

# Process peak_table
peak_table <- peak_table[-1,]
split_column <- strsplit(as.character(peak_table$Sample), "__")
peak_table$Query.Mass <- sapply(split_column, function(x) x[1])
peak_table$Retention.Time <- sapply(split_column, function(x) x[2])
peak_table$Query.Mass <- as.numeric(peak_table$Query.Mass)
peak_table$Retention.Time <- as.numeric(peak_table$Retention.Time)

# Process compound_table
data$Query.Mass <- as.numeric(data$Query.Mass)
data$Retention.Time <- as.numeric(data$Retention.Time)

# Merge the tables
merged_table <- data %>%
  left_join(peak_table, by = c("Retention.Time", "Query.Mass"))

# Select specific columns
selected_columns <- merged_table[, c(7, 9:16)]

# Write the filtered table to a CSV file
write.csv(selected_columns, file.path(config$results, config$output_filtered_csv), row.names = FALSE)

