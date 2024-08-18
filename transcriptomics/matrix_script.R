# Loading necessary libraries
library(yaml)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Loading configuration from the YAML file
config <- yaml::read_yaml("/SWIS-META/transcriptomics/config.yaml")
count_dir <- config$COUNT_DIR
out_final <- config$OUT_FINAL

read_featureCounts <- function(filename) {
  # Reading the first few lines to determine the number of lines to skip
  lines <- readLines(filename, warn = FALSE)
  skip_lines <- sum(grepl("^#", lines))  # Count comment lines to skip
  
  # Reading the file, skipping comment lines
  dat <- read.delim(filename, header = TRUE, sep = "\t", stringsAsFactors = FALSE, skip = skip_lines)
  
  # Verifying that 'Geneid' is among the column names
  if (!"Geneid" %in% colnames(dat)) {
    stop("Geneid column not found in ", filename)
  }
  
  # Extracting 'Geneid' column and the counts column (assumed to be the last column)
  # Assuming the counts are in the last column. Adjust if the file structure is different.
  gene_id_col <- which(colnames(dat) == "Geneid")
  counts_col <- ncol(dat)
  dat_subset <- dat[, c(gene_id_col, counts_col)]
  
  # Renaming the columns: Geneid and the sample name based on the file name
  colnames(dat_subset) <- c("Geneid", gsub("_counts\\.txt$", "", basename(filename)))
  return(dat_subset)
}

count_files <- list.files(path = count_dir, pattern = "_counts.txt$", full.names = TRUE)
count_list <- lapply(count_files, read_featureCounts)

# Merging all count data into a single data frame by gene ID
gene_counts <- Reduce(function(x, y) merge(x, y, by = "Geneid", all = TRUE), count_list)

# Mapping ENSEMBL gene IDs to gene symbols
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = gene_counts$Geneid,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

# Replacing ENSEMBL IDs with gene symbols in the data frame
gene_counts$Geneid <- gene_symbols[match(gene_counts$Geneid, names(gene_symbols))]

# Aggregating counts by gene symbol, summing counts for duplicated symbols
aggregate_cols <- setdiff(names(gene_counts), "Geneid")
gene_counts_aggregated <- aggregate(. ~ Geneid, data = gene_counts, FUN = sum)

# Setting Geneid column as row names
rownames(gene_counts_aggregated) <- gene_counts_aggregated$Geneid
gene_counts_aggregated$Geneid <- NULL  # Remove the now redundant 'Geneid' column

write.csv(gene_counts_aggregated, file.path(out_final, "gene_symbol_counts_matrix.csv"), row.names = TRUE)

