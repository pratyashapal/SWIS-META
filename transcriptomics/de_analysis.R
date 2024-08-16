# Load required libraries
library(edgeR)
library(limma)
library(fgsea)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DESeq2)
library(dplyr)
library(enrichplot)
library(patchwork)
library(cowplot)
library(ggrepel)
library(DOSE)
library(tibble)
library(yaml)  # For reading YAML files

# Load configuration file
config <- yaml.load_file("config_test2.yaml")

# Load the count matrix and metadata from config paths
count_matrix <- read.csv(config$count_matrix_path, row.names = 1, header = TRUE)
meta_data <- read.csv(config$metadata_path, header = TRUE)
meta_data = meta_data[-c(9,10),]
rownames(meta_data) = meta_data$sample

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = meta_data, design = ~ condition)
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)

# Access normalized counts, if needed
normalized_counts <- counts(dds, normalized = TRUE)
log_norm_counts <- log2(normalized_counts + 1)
pca_res <- prcomp(t(log_norm_counts))

# Create a data frame for plotting
pca_data <- data.frame(Sample = rownames(pca_res$x),
                       PC1 = pca_res$x[,1],
                       PC2 = pca_res$x[,2],
                       Condition = meta_data$condition[rownames(pca_res$x)])
pca_data$Condition <- meta_data[match(pca_data$Sample, rownames(meta_data)), "condition"]

# Create the PCA plot with bold sample labels
pca_log_norm_counts <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition, label = Sample)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3) + 
  theme_minimal() +
  ggtitle("PCA of Gene Expression") +
  xlab(paste("PC1 -", round(summary(pca_res)$importance[2,1] * 100, 2), "% Variance")) +
  ylab(paste("PC2 -", round(summary(pca_res)$importance[2,2] * 100, 2), "% Variance")) +
  scale_color_manual(values = c("bm" = "royalblue3", "cns" = "firebrick2")) +
  theme(
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  )

ggsave(filename = file.path(config$output_dir, "pca_plot.pdf"), plot = pca_log_norm_counts, device = "pdf")

# Significant genes
res = results(dds)
sig_res <- subset(res, padj < 0.05 & abs(log2FoldChange) >1)
sig_genes_df <- as.data.frame(sig_res)
write.csv(sig_genes_df, file.path(config$output_dir, "sig_genes.csv"))

# Sample-to-sample correlation matrix
sample_cor_matrix <- cor(t(normalized_counts))
color_palette <- colorRampPalette(c("seagreen", "khaki", "coral1"))(100)

sample_corr_matrix <- pheatmap(
  sample_cor_matrix,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  scale = "none",
  color = color_palette,
  main = "Sample-to-Sample Correlation",
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_legend = TRUE,
  fontsize_row = 10,
  fontsize_col = 10,
  border_color = NA,
  legend_breaks = seq(-1, 1, by = 0.5),
  legend_labels = c("-1", "-0.5", "0", "0.5", "1")
)

ggsave(filename = file.path(config$output_dir, "correlation_matrix.pdf"), plot = sample_corr_matrix$gtable, device = "pdf")

# Heatmap of top 25 significant genes
top25_genes <- head(rownames(sig_genes_df[order(sig_genes_df$padj),]), 25)
top25_counts <- counts(dds, normalized = TRUE)[top25_genes,]
log_top25_counts <- log2(top25_counts + 1)

heatmap_top25 <- pheatmap(log_top25_counts,
                          clustering_distance_rows = "euclidean",
                          clustering_distance_cols = "euclidean",
                          scale = "row",
                          annotation_col = meta_data,
                          show_rownames = TRUE,
                          show_colnames = TRUE,
                          color = colorRampPalette(c("seagreen", "khaki", "coral3"))(100),
                          main = "Heatmap of Top 25 Differentially Expressed Genes")

ggsave(filename = file.path(config$output_dir, "heatmap_top25_genes.pdf"), plot = heatmap_top25$gtable, device = "pdf")

# Volcano plot
top_upregulated <- sig_genes_df %>%
  filter(log2FoldChange > 1) %>%
  arrange(padj) %>%
  head(5)

top_downregulated <- sig_genes_df %>%
  filter(log2FoldChange < -1) %>%
  arrange(padj) %>%
  head(5)

top_genes <- rbind(top_upregulated, top_downregulated)
top_genes$label <- rownames(top_genes)
top_genes$color <- ifelse(top_genes$log2FoldChange > 1, "Upregulated", "Downregulated")

sig_genes_df$color <- ifelse(sig_genes_df$log2FoldChange > 1, "Upregulated",
                             ifelse(sig_genes_df$log2FoldChange < -1, "Downregulated", "Not Significant"))

volcano_plot <- ggplot(sig_genes_df, aes(x = log2FoldChange, y = -log10(padj), color = color)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "firebrick2", "Downregulated" = "royalblue3", "Not Significant" = "grey")) +
  geom_point(data = top_genes, aes(x = log2FoldChange, y = -log10(padj), color = color), size = 3) +
  geom_text_repel(data = top_genes, aes(label = label), size = 3.5, fontface = "bold", nudge_x = 0.5, nudge_y = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray") +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted p-value") +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

ggsave(filename = file.path(config$output_dir, "volcano_plot.pdf"), plot = volcano_plot, device = "pdf")

# KEGG analysis
perform_kegg_analysis <- function(gene_list, title_suffix, n_top_pathways = 10, n_top_metabolism = 20) {
  gene_entrez <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  gene_entrez <- gene_entrez[!is.na(gene_entrez$ENTREZID),]
  kegg_result <- enrichKEGG(gene = gene_entrez$ENTREZID, organism = 'hsa', keyType = 'kegg', pAdjustMethod = "BH", qvalueCutoff = 0.05)
  kegg_results_df <- as.data.frame(kegg_result@result)
  top_kegg_res <- kegg_results_df %>%
    arrange(p.adjust) %>%
    slice_head(n = n_top_pathways)
  metabolism_kegg_res <- kegg_results_df %>%
    filter(grepl("Metabolism", Description)) %>%
    arrange(p.adjust) %>%
    slice_head(n = n_top_metabolism)

  # Save KEGG results as CSV
  write.csv(kegg_results_df, file.path(config$output_dir, paste0("kegg_results_", title_suffix, ".csv")))
  write.csv(top_kegg_res, file.path(config$output_dir, paste0("top_kegg_pathways_", title_suffix, ".csv")))
  write.csv(metabolism_kegg_res, file.path(config$output_dir, paste0("top_metabolism_kegg_", title_suffix, ".csv")))

  # Plot KEGG
  if (nrow(top_kegg_res) > 0) {
    kegg_plot <- ggplot(top_kegg_res, aes(x = reorder(Description, -log10(p.adjust)), y = -log10(p.adjust))) +
      geom_bar(stat = "identity", fill = "steelblue") +
      coord_flip() +
      labs(title = paste("Top KEGG Pathways", title_suffix), x = "Pathway", y = "-Log10 Adjusted p-value") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)
      )
    ggsave(filename = file.path(config$output_dir, paste0("kegg_pathways_", title_suffix, ".pdf")), plot = kegg_plot, device = "pdf")
  }
  
  if (nrow(metabolism_kegg_res) > 0) {
    metabolism_kegg_plot <- ggplot(metabolism_kegg_res, aes(x = reorder(Description, -log10(p.adjust)), y = -log10(p.adjust))) +
      geom_bar(stat = "identity", fill = "darkorange") +
      coord_flip() +
      labs(title = paste("Top KEGG Metabolism Pathways", title_suffix), x = "Pathway", y = "-Log10 Adjusted p-value") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)
      )
    ggsave(filename = file.path(config$output_dir, paste0("kegg_metabolism_", title_suffix, ".pdf")), plot = metabolism_kegg_plot, device = "pdf")
  }
}

# Perform KEGG analysis
perform_kegg_analysis(rownames(sig_genes_df), "differentially_expressed_genes")

