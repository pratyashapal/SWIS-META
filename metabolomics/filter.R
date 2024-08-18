# Loading necessary libraries
library(KEGGREST)
library(dplyr)
library(yaml)

# Loading parameters from config file
config <- yaml::read_yaml("/SWIS-META/metabolomics/config_2.yaml")

# Determining the output directory based on ionization mode
output_dir <- file.path(config$results, config$mode)

# Setting the working directory to the output directory
setwd(output_dir)

# Reading the input files
if (config$mode == "positive") {
  input_csv <- file.path(config$input_csv_positive)
  peaks_csv <- file.path(config$pos_peaks_csv)
} else if (config$mode == "negative") {
  input_csv <- file.path(config$input_csv_negative)
  peaks_csv <- file.path(config$neg_peaks_csv)
} else {
  stop("Unknown ionization mode specified.")
}

data <- read.csv(input_csv, header = TRUE)
peak_table <- read.csv(peaks_csv)

# Defining a function to get the compound name for a given KEGG ID
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

# Extracting KEGG IDs from the specified column and clean them
kegg_ids <- trimws(as.character(data[[config$column_index]]))

# Applying the function to the list of KEGG IDs
compound_names <- sapply(kegg_ids, get_compound_name)

# Adding the compound names to the original data frame
data$Compound_Name <- compound_names

# Removing rows with NA in Compound_Name
data <- data %>%
  filter(!is.na(Compound_Name))

# Processing peak_table
peak_table <- peak_table[-1,]
split_column <- strsplit(as.character(peak_table$Sample), "__")
peak_table$Query.Mass <- sapply(split_column, function(x) x[1])
peak_table$Retention.Time <- sapply(split_column, function(x) x[2])
peak_table$Query.Mass <- as.numeric(peak_table$Query.Mass)
peak_table$Retention.Time <- as.numeric(peak_table$Retention.Time)

# Processing compound_table
data$Query.Mass <- as.numeric(data$Query.Mass)
data$Retention.Time <- as.numeric(data$Retention.Time)

# Merging the tables
merged_table <- data %>%
  left_join(peak_table, by = c("Retention.Time", "Query.Mass"))

# Selecting specific columns
selected_columns <- merged_table[, c(7, 9:16)]

# Write the filtered table to a CSV file in the output directory
output_file <- file.path(output_dir, paste0(config$output_prefix, ".csv"))
write.csv(selected_columns, output_file, row.names = FALSE)

