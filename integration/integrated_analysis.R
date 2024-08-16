# Load necessary packages
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(patchwork)
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
library(igraph)
library(ggraph)
library(yaml)

# Read config file
config <- yaml::read_yaml("config.yaml")
metabolomics_results_dir <- config$metabolomics_results_dir
transcriptomics_results_dir <- config$transcriptomics_results_dir
results_dir <- config$results_dir

# Create the results directory if it doesn't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# Define file paths
compound_pathway_file <- file.path(metabolomics_results_dir)
gene_pathway_file <- file.path(transcriptomics_results_dir)

# Read compound and gene data
compound_data_with_names <- read.csv(compound_pathway_file, header = TRUE, stringsAsFactors = FALSE)
gene_data <- read.csv(gene_pathway_file, header = TRUE, stringsAsFactors = FALSE)

# Extract relevant columns

gene_data <- gene_data[, c(1, 3)]
colnames(gene_data) <- c("GeneSymbol", "LogFoldChange")

# Map gene symbols to Entrez IDs using org.Hs.eg.db
gene_data$EntrezID <- mapIds(org.Hs.eg.db, 
                             keys = gene_data$GeneSymbol, 
                             column = "ENTREZID", 
                             keytype = "SYMBOL", 
                             multiVals = "first")

# Remove rows with missing Entrez IDs
gene_data <- na.omit(gene_data)

# Perform KEGG pathway enrichment analysis for genes
kegg_gene_enrichment <- enrichKEGG(gene = gene_data$EntrezID, 
                                   organism = 'hsa', 
                                   keyType = 'kegg', 
                                   pAdjustMethod = 'BH', 
                                   qvalueCutoff = 0.05)

# Extract pathways from KEGG enrichment results
pathway_data <- as.data.frame(kegg_gene_enrichment@result)

# Filter to include only pathways related to "Metabolism" and "Environmental Information Processing"
filtered_pathways <- pathway_data %>%
  filter(grepl("Metabolism", category, ignore.case = TRUE) | 
           grepl("Environmental Information Processing", category, ignore.case = TRUE))

# Split geneID into separate rows
gene_pathway_df <- filtered_pathways %>%
  dplyr::select(Description, geneID) %>%
  separate_rows(geneID, sep = "/") %>%
  rename(PathwayName = Description, EntrezID = geneID)

# Merge gene data with pathways
gene_data_with_pathways <- merge(gene_data, gene_pathway_df, by = 'EntrezID')
gene_data_with_pathways <- na.omit(gene_data_with_pathways)


# Combine gene and compound data by pathway
combined_data <- inner_join(gene_data_with_pathways, compound_data_with_names, by = "PathwayName")

# Create summary table of integrated pathways
integrated_pathways <- combined_data %>% 
  group_by(PathwayName) %>% 
  summarise(genes = paste(GeneSymbol, collapse = ", "),
            compounds = paste(CompoundName, collapse = ", "),
            avg_gene_lfc = mean(LogFoldChange.x),
            avg_compound_lfc = mean(LogFoldChange.y))

# Save the integrated pathways to a CSV file
write.csv(integrated_pathways, file.path(results_dir, "integrated_pathways.csv"), row.names = FALSE)

# Plot average log fold changes for genes and compounds in all pathways
p1 <- ggplot(integrated_pathways, aes(x = reorder(PathwayName, avg_gene_lfc))) +
  geom_bar(aes(y = avg_gene_lfc, fill = "Gene LFC"), stat = "identity", position = "dodge", width = 0.4) +
  geom_bar(aes(y = avg_compound_lfc, fill = "Compound LFC"), stat = "identity", position = "dodge", width = 0.4) +
  labs(title = "Average Log Fold Change of Genes and Compounds in Pathways",
       x = "Pathway",
       y = "Average Log Fold Change") +
  coord_flip() +
  theme_minimal() +
  scale_fill_manual(values = c("Gene LFC" = "steelblue", "Compound LFC" = "darkorange")) +
  theme(legend.title = element_blank())

ggsave(file.path(results_dir, "avg_log_fold_changes.pdf"), plot = p1, device = "pdf", width = 10, height = 8)

# Bar plot of average log fold changes for genes
p2 <- ggplot(integrated_pathways, aes(x = reorder(PathwayName, avg_gene_lfc), y = avg_gene_lfc)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Average Gene Log Fold Changes per Pathway", x = "Pathway", y = "Average Gene Log Fold Change") +
  theme_minimal()

ggsave(file.path(results_dir, "avg_gene_log_fold_changes.pdf"), plot = p2, device = "pdf", width = 10, height = 8)

# Bar plot of average log fold changes for compounds
p3 <- ggplot(integrated_pathways, aes(x = reorder(PathwayName, avg_compound_lfc), y = avg_compound_lfc)) +
  geom_bar(stat = "identity", fill = "darkorange") +
  coord_flip() +
  labs(title = "Average Compound Log Fold Changes per Pathway", x = "Pathway", y = "Average Compound Log Fold Change") +
  theme_minimal()

ggsave(file.path(results_dir, "avg_compound_log_fold_changes.pdf"), plot = p3, device = "pdf", width = 10, height = 8)

# Prepare data for heatmap
heatmap_data <- integrated_pathways %>%
  select(PathwayName, avg_gene_lfc, avg_compound_lfc) %>%
  column_to_rownames('PathwayName') %>%
  as.matrix()

# Plot heatmap
p4 <- pheatmap(heatmap_data, 
               cluster_rows = TRUE, 
               cluster_cols = FALSE, 
               display_numbers = FALSE,
               main = "Heatmap of Average Log Fold Changes",
               fontsize_row = 10,
               color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

# Save heatmap as PDF
ggsave(file.path(results_dir, "heatmap_log_fold_changes.pdf"), plot = p4, device = "pdf", width = 10, height = 8)

# Create a network graph of pathways, genes, and compounds
edges <- combined_data %>%
  select(PathwayName, GeneSymbol, CompoundName)

nodes <- data.frame(
  name = unique(c(edges$PathwayName, edges$GeneSymbol, edges$CompoundName)),
  type = c(rep("Pathway", length(unique(edges$PathwayName))),
           rep("Gene", length(unique(edges$GeneSymbol))),
           rep("Compound", length(unique(edges$CompoundName))))
)

network <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)

p5 <- ggraph(network, layout = "fr") +
  geom_edge_link(aes(edge_alpha = 0.2)) +
  geom_node_point(aes(color = type), size = 5) +
  geom_node_text(aes(label = name), repel = TRUE) +
  scale_color_manual(values = c("Pathway" = "blue", "Gene" = "green", "Compound" = "red")) +
  theme_void() +
  labs(title = "Network of Pathways, Genes, and Compounds")

ggsave(file.path(results_dir, "network_graph.pdf"), plot = p5, device = "pdf", width = 20, height = 15)


# Separate analysis for upregulated and downregulated genes and compounds
upregulated_genes <- gene_data_with_pathways %>% filter(LogFoldChange > 1)
downregulated_genes <- gene_data_with_pathways %>% filter(LogFoldChange < -1)
upregulated_compounds <- compound_data_with_names %>% filter(LogFoldChange > 1)
downregulated_compounds <- compound_data_with_names %>% filter(LogFoldChange < -1)

# Combined analysis for upregulated genes and compounds
combined_upregulated <- inner_join(upregulated_genes, upregulated_compounds, by = "PathwayName")
combined_upregulated_pathways <- combined_upregulated %>% 
  group_by(PathwayName) %>% 
  summarise(genes = paste(GeneSymbol, collapse = ", "),
            compounds = paste(CompoundName, collapse = ", "),
            avg_gene_lfc = mean(LogFoldChange.x),
            avg_compound_lfc = mean(LogFoldChange.y))

# Combined analysis for downregulated genes and compounds
combined_downregulated <- inner_join(downregulated_genes, downregulated_compounds, by = "PathwayName")
combined_downregulated_pathways <- combined_downregulated %>% 
  group_by(PathwayName) %>% 
  summarise(genes = paste(GeneSymbol, collapse = ", "),
            compounds = paste(CompoundName, collapse = ", "),
            avg_gene_lfc = mean(LogFoldChange.x),
            avg_compound_lfc = mean(LogFoldChange.y))

# Combine the data for gene and compound log fold changes into a single data frame
combined_upregulated_pathways_long <- combined_upregulated_pathways %>%
  pivot_longer(cols = c(avg_gene_lfc, avg_compound_lfc), 
               names_to = "Type", 
               values_to = "Log_Fold_Change")

# Create the faceted bar plot for upregulated pathways
p6 <- ggplot(combined_upregulated_pathways_long, aes(x = reorder(PathwayName, Log_Fold_Change), y = Log_Fold_Change, fill = Type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  coord_flip() +
  labs(title = "Average Upregulated Log Fold Changes per Pathway",
       x = "Pathway", 
       y = "Average Log Fold Change",
       fill = "Type") +
  scale_fill_manual(values = c("avg_gene_lfc" = "steelblue", "avg_compound_lfc" = "darkorange"),
                    labels = c("Gene Log Fold Change", "Compound Log Fold Change")) +
  theme_minimal() +
  theme(
    text = element_text(size = 12),  # Adjust text size
    axis.text.x = element_text(angle = 0, vjust = 0.5),  # Align x-axis text
    axis.text.y = element_text(size = 10),  # Adjust y-axis text size
    axis.title = element_text(size = 14, face = "bold"),  # Bold axis titles
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # Center and bold plot title
    legend.title = element_text(size = 12),  # Adjust legend title size
    legend.text = element_text(size = 10),  # Adjust legend text size
    panel.grid.major = element_line(color = "gray90"),  # Light grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )

# Save the plot
ggsave(file.path(results_dir, "upregulated_pathways_log_fold_changes.pdf"), plot = p6, device = "pdf", width = 10, height = 8)

# Combine the data for gene and compound log fold changes into a single data frame
combined_downregulated_pathways_long <- combined_downregulated_pathways %>%
  pivot_longer(cols = c(avg_gene_lfc, avg_compound_lfc), 
               names_to = "Type", 
               values_to = "Log_Fold_Change")

# Create the faceted bar plot for downregulated pathways
p7 <- ggplot(combined_downregulated_pathways_long, aes(x = reorder(PathwayName, Log_Fold_Change), y = Log_Fold_Change, fill = Type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  coord_flip() +
  labs(title = "Average Downregulated Log Fold Changes per Pathway",
       x = "Pathway", 
       y = "Average Log Fold Change",
       fill = "Type") +
  scale_fill_manual(values = c("avg_gene_lfc" = "steelblue", "avg_compound_lfc" = "darkorange"),
                    labels = c("Gene Log Fold Change", "Compound Log Fold Change")) +
  theme_minimal() +
  theme(
    text = element_text(size = 12),  # Adjust text size
    axis.text.x = element_text(angle = 0, vjust = 0.5),  # Align x-axis text
    axis.text.y = element_text(size = 10),  # Adjust y-axis text size
    axis.title = element_text(size = 14, face = "bold"),  # Bold axis titles
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # Center and bold plot title
    legend.title = element_text(size = 12),  # Adjust legend title size
    legend.text = element_text(size = 10),  # Adjust legend text size
    panel.grid.major = element_line(color = "gray90"),  # Light grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )

# Save the plot
ggsave(file.path(results_dir, "downregulated_pathways_log_fold_changes.pdf"), plot = p7, device = "pdf", width = 10, height = 8)

