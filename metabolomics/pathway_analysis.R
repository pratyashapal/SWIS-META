# Load necessary packages
library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(patchwork)
library(pathview)
library(dplyr)
library(KEGGREST)
library(tibble)
library(DOSE)
library(ggplot2)
library(DESeq2)
library(edgeR)
library(limma)
library(fgsea)
library(pheatmap)
library(tidyr)
library(yaml)

# Read config file
config <- yaml::read_yaml("/SWIS-META/metabolomics/config_main.yaml")
results_dir <- config$results

# Define file paths
filtered_combined_path <- file.path(results_dir, "filtered_combined.csv")
filtered_combined_log_path <- file.path(results_dir, "filtered_combined_log.csv")

# Load data
file <- read.csv(filtered_combined_path)
file_log <- read.csv(filtered_combined_log_path, header = TRUE, stringsAsFactors = FALSE)

compound_data <- file_log[, c("X", "log2.FC.")]
colnames(compound_data) <- c("CompoundName", "LogFoldChange")

# Function to get KEGG pathway mappings for compounds
get_kegg_pathways_for_compounds <- function(compounds) {
  pathway_list <- list()
  for (compound in compounds) {
    tryCatch({
      kegg_compound <- keggFind("compound", compound)
      if (length(kegg_compound) > 0) {
        compound_id <- names(kegg_compound)[1]
        pathways <- keggLink("pathway", compound_id)
        pathway_ids <- gsub("path:", "", pathways)
        pathway_list[[compound]] <- pathway_ids
      }
    }, error = function(e) {})
  }
  return(pathway_list)
}

# Get pathway mappings for compounds
compound_pathways <- get_kegg_pathways_for_compounds(compound_data$CompoundName)
compound_pathway_df <- do.call(rbind, lapply(names(compound_pathways), function(compound) {
  data.frame(CompoundName = compound, PathwayID = compound_pathways[[compound]], stringsAsFactors = FALSE)
}))

# Merge compound data with pathway mappings
compound_data_with_pathways <- merge(compound_data, compound_pathway_df, by.x = 'CompoundName', by.y = 'CompoundName')
compound_data_with_pathways <- na.omit(compound_data_with_pathways)

# Function to get pathway names from KEGG
get_pathway_names <- function(pathway_ids) {
  pathway_names <- sapply(pathway_ids, function(id) {
    tryCatch({
      info <- keggGet(paste0("path:", id))
      if (length(info) > 0) {
        return(info[[1]]$NAME)
      } else {
        return(NA)
      }
    }, error = function(e) {
      return(NA)
    })
  })
  return(pathway_names)
}

# Retrieve common names for the pathways
pathway_names <- get_pathway_names(unique(compound_data_with_pathways$PathwayID))

# Create a data frame for PathwayID and PathwayName
pathway_names_df <- data.frame(
  PathwayID = names(pathway_names),
  PathwayName = pathway_names,
  stringsAsFactors = FALSE
)

# Merge with the original compound data
compound_data_with_names <- merge(compound_data_with_pathways, pathway_names_df, by = "PathwayID", all.x = TRUE)

# Calculate the average log fold change for each pathway
pathway_avg_logFC <- compound_data_with_names %>%
  group_by(PathwayName) %>%
  summarise(avg_logFC = mean(LogFoldChange, na.rm = TRUE),
            .groups = 'drop')

# Separate pathways into upregulated and downregulated based on average log fold change
pathway_avg_logFC <- pathway_avg_logFC %>%
  mutate(Direction = ifelse(avg_logFC > 0, "Upregulated", "Downregulated"))

# Specify pathways to exclude
excluded_pathways <- c("Autophagy - animal", "Autophagy - other", "Autophagy - yeast", "Ubiquinone and other terpenoid-quinone biosynthesis", "Tuberculosis", "Systemic lupus erythematosus", "Pathogenic Escherichia coli infection")

# Filter out specified pathways
pathway_avg_logFC <- pathway_avg_logFC %>%
  filter(!PathwayName %in% excluded_pathways)

# Filter top 10 upregulated and top 10 downregulated pathways
top_upregulated <- pathway_avg_logFC %>%
  filter(Direction == "Upregulated") %>%
  arrange(desc(avg_logFC)) %>%
  head(10)

top_downregulated <- pathway_avg_logFC %>%
  filter(Direction == "Downregulated") %>%
  arrange(avg_logFC) %>%
  head(10)

# Combine top pathways
top_pathways <- bind_rows(top_upregulated, top_downregulated)

# Plot the data
ggplot(top_pathways, aes(x = reorder(PathwayName, avg_logFC), y = avg_logFC, fill = Direction)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Top 10 Upregulated and Downregulated Pathways",
       x = "Pathway",
       y = "Average Log Fold Change") +
  scale_fill_manual(values = c("Upregulated" = "steelblue", "Downregulated" = "darkorange")) +
  theme_minimal() +
  ggsave(filename = file.path(results_dir, "plots/top_pathways.pdf"), device = "pdf")

# Count the number of occurrences of each PathwayName
pathway_counts <- compound_data_with_names %>%
  group_by(PathwayName) %>%
  summarise(Count = n(), .groups = 'drop')

# Order pathways by compound count (descending) and select top 20
top_20_pathways <- pathway_counts[order(pathway_counts$Count, decreasing = TRUE), ][1:20, ]

# Reverse the order to have the highest counts at the top
top_20_pathways <- top_20_pathways[rev(seq_len(nrow(top_20_pathways))), ]

# Set up plot margins to accommodate long names and shift y-axis label
par(mar = c(7, 18, 4, 3) + 0.1)  # Increase left margin to fit long names and shift y-axis label

# Create the horizontal bar plot
barplot(
  top_20_pathways$Count,
  names.arg = top_20_pathways$PathwayName,
  horiz = TRUE,         # Bars are horizontal
  col = "lightblue",
  main = "Number of Compounds per KEGG Pathway (Top 20)",
  xlab = "Number of Compounds",
  ylab = "",            # Remove default y-axis label
  cex.names = 0.7,      # Adjust text size for readability
  las = 1,              # Rotate axis labels (1 for horizontal text)
  cex.axis = 0.9        # Adjust axis text size for better fit
)
mtext("Pathways", side = 2, line = 0.1, cex = 1)
ggsave(filename = file.path(results_dir, "plots/top_20_pathways.pdf"), device = "pdf")

# Save the compound-pathway mapping to CSV
write.csv(compound_data_with_names, file.path(results_dir, "compound_pathway_with_names.csv"), row.names = FALSE)

