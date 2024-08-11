#!/usr/bin/env Rscript

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
library(cowplot)
library(ggrepel)
library(DOSE)
library(tibble)
library(yaml)

# Read the configuration file
config <- yaml::read_yaml("config_main.yaml")

# Extract file paths from config
count_matrix_file <- config$count_matrix_file
results_dir <- config$results_dir

# Define output file names
output_files <- list(
  volcano_plot = file.path(results_dir, "volcano_plot.pdf"),
  pca_plot = file.path(results_dir, "pca_plot.pdf"),
  sample_correlation_heatmap = file.path(results_dir, "sample_correlation_heatmap.pdf"),
  significant_genes_heatmap = file.path(results_dir, "significant_genes_heatmap.pdf"),
  top25_heatmap = file.path(results_dir, "top25_heatmap.pdf"),
  kegg_enrichment_csv = file.path(results_dir, "kegg_enrichment.csv"),
  kegg_enrichment_gene_names_csv = file.path(results_dir, "kegg_pathway_enrichment_gene_names.csv"),
  sig_genes_csv = file.path(results_dir, "sig_genes_df.csv"),
  gsea_upregulated_top10 = file.path(results_dir, "gsea_upregulated_top10.pdf"),
  gsea_downregulated_top10 = file.path(results_dir, "gsea_downregulated_top10.pdf"),
  kegg_upregulated_top10 = file.path(results_dir, "kegg_upregulated_top10.pdf"),
  kegg_downregulated_top10 = file.path(results_dir, "kegg_downregulated_top10.pdf"),
  metabolism_upregulated_top10 = file.path(results_dir, "metabolism_upregulated_top10.pdf"),
  metabolism_downregulated_top10 = file.path(results_dir, "metabolism_downregulated_top10.pdf")
)

# Load the count matrix
count_matrix <- read.csv(count_matrix_file, row.names = 1, header = TRUE)

# Define sample information
col_data <- data.frame(
  sample = colnames(count_matrix),
  condition = factor(c(rep("bm", 4), rep("cns", 4)))
)
rownames(col_data) <- col_data$sample

# Create DESeq2 object and perform differential analysis
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = col_data, design = ~ condition)
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
res <- results(dds)

# Extract significant results
sig_res <- subset(res, padj < 0.05)
sig_genes_df <- as.data.frame(sig_res)

# Save significant genes to CSV
write.csv(sig_genes_df, output_files$sig_genes_csv)

# Identify top 5 upregulated and downregulated genes
top_upregulated <- sig_genes_df %>%
  filter(log2FoldChange > 1) %>%
  arrange(padj) %>%
  head(5)

top_downregulated <- sig_genes_df %>%
  filter(log2FoldChange < -1) %>%
  arrange(padj) %>%
  head(5)

# Combine top genes for labeling
top_genes <- rbind(top_upregulated, top_downregulated)
top_genes$label <- rownames(top_genes)

# Create volcano plot
volcano_plot <- ggplot(sig_genes_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_point(data = top_genes, aes(x = log2FoldChange, y = -log10(padj)), color = "red", size = 3) +
  geom_text_repel(data = top_genes, aes(label = label)) +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted p-value") +
  theme_minimal()
ggsave(output_files$volcano_plot, plot = volcano_plot, device = "pdf", width = 10, height = 10, dpi = 300)

# Access normalized counts and perform PCA
normalized_counts <- counts(dds, normalized = TRUE)
log_norm_counts <- log2(normalized_counts + 1)

# PCA analysis
pca_res <- prcomp(t(log_norm_counts))
pca_data <- data.frame(
  Sample = rownames(pca_res$x),
  PC1 = pca_res$x[,1],
  PC2 = pca_res$x[,2],
  Condition = col_data$condition[rownames(pca_res$x)]
)

# PCA Plot
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition, label = Sample)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3) +
  theme_minimal() +
  ggtitle("PCA of Gene Expression") +
  xlab(paste("PC1 -", round(summary(pca_res)$importance[2,1] * 100, 2), "% Variance")) +
  ylab(paste("PC2 -", round(summary(pca_res)$importance[2,2] * 100, 2), "% Variance")) +
  scale_color_manual(values = c("bm" = "blue", "cns" = "red")) +
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))
ggsave(output_files$pca_plot, plot = pca_plot, device = "pdf", width = 10, height = 10, dpi = 300)

# Sample-to-Sample Correlation Heatmap
sample_cor_matrix <- cor(log_norm_counts, method = "pearson")
pheatmap(sample_cor_matrix,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         scale = "none",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Sample-to-Sample Correlation",
         show_rownames = TRUE,
         show_colnames = TRUE,
         annotation_legend = TRUE,
         filename = output_files$sample_correlation_heatmap)

# Heatmap of significant genes
sig_gene_names <- rownames(sig_genes_df)
norm_counts <- counts(dds, normalized = TRUE)[sig_gene_names,]
log_norm_counts <- log2(norm_counts + 1)

pheatmap(log_norm_counts,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         scale = "row",
         annotation_col = col_data,
         show_rownames = FALSE,
         show_colnames = FALSE,
         color = colorRampPalette(c("navyblue", "white", "firebrick3"))(100),
         main = "Heatmap of Significant Genes",
         filename = output_files$significant_genes_heatmap)

# Subset top DEGs for heatmap
top_degs <- sig_genes_df[order(sig_genes_df$padj), ]
top_degs_subset <- head(top_degs, 25)
top_degs_genes <- rownames(top_degs_subset)

log_norm_counts_subset <- log_norm_counts[top_degs_genes, ]

pheatmap(log_norm_counts_subset,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         scale = "row",
         annotation_col = col_data,
         show_rownames = TRUE,
         show_colnames = FALSE,
         color = colorRampPalette(c("navyblue", "white", "firebrick3"))(100),
         main = "Heatmap of Top 25 DEGs",
         filename = output_files$top25_heatmap)

# KEGG enrichment analysis
degs_final <- rownames(sig_genes_df)
degs_final_entrez <- bitr(degs_final, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
degs_final_entrez <- degs_final_entrez[!is.na(degs_final_entrez$ENTREZID),]

kegg <- enrichKEGG(gene = degs_final_entrez$ENTREZID, organism = 'hsa', keyType = 'kegg')
write.csv(as.data.frame(kegg), file = output_files$kegg_enrichment_csv)
write.csv(as.data.frame(kegg) %>% mutate(GeneNames = mapIds(org.Hs.eg.db, as.character(.$geneID), 'SYMBOL', 'ENTREZID')),
          file = output_files$kegg_enrichment_gene_names_csv)

# Plot top 10 KEGG pathways (Upregulated and Downregulated)
top10_kegg_upregulated <- head(kegg@result[order(kegg@result$p.adjust),], 10)
top10_kegg_downregulated <- head(kegg@result[order(kegg@result$p.adjust, decreasing = TRUE)], 10)

dotplot_kegg_upregulated <- dotplot(as.data.frame(top10_kegg_upregulated), showCategory = 10) +
  ggtitle("Top 10 KEGG Upregulated Pathways") +
  theme_minimal()
ggsave(output_files$kegg_upregulated_top10, plot = dotplot_kegg_upregulated, device = "pdf", width = 10, height = 10, dpi = 300)

dotplot_kegg_downregulated <- dotplot(as.data.frame(top10_kegg_downregulated), showCategory = 10) +
  ggtitle("Top 10 KEGG Downregulated Pathways") +
  theme_minimal()
ggsave(output_files$kegg_downregulated_top10, plot = dotplot_kegg_downregulated, device = "pdf", width = 10, height = 10, dpi = 300)

# GSEA analysis
gsea_upregulated <- gseKEGG(geneList = sort(sig_genes_df$log2FoldChange, decreasing = TRUE), organism = 'hsa')
gsea_downregulated <- gseKEGG(geneList = sort(sig_genes_df$log2FoldChange, decreasing = FALSE), organism = 'hsa')

# Filter for top 10 pathways
top10_gsea_upregulated <- head(gsea_upregulated@result[order(gsea_upregulated@result$p.adjust),], 10)
top10_gsea_downregulated <- head(gsea_downregulated@result[order(gsea_downregulated@result$p.adjust),], 10)

# Plot GSEA results
dotplot_gsea_upregulated <- dotplot(as.data.frame(top10_gsea_upregulated), showCategory = 10) +
  ggtitle("Top 10 GSEA Upregulated Pathways") +
  theme_minimal()
ggsave(output_files$gsea_upregulated_top10, plot = dotplot_gsea_upregulated, device = "pdf", width = 10, height = 10, dpi = 300)

dotplot_gsea_downregulated <- dotplot(as.data.frame(top10_gsea_downregulated), showCategory = 10) +
  ggtitle("Top 10 GSEA Downregulated Pathways") +
  theme_minimal()
ggsave(output_files$gsea_downregulated_top10, plot = dotplot_gsea_downregulated, device = "pdf", width = 10, height = 10, dpi = 300)

# Metabolism pathways enrichment analysis
metabolism_terms <- c("Metabolism", "Metabolic pathways")
metabolism_pathways <- enrichKEGG(gene = degs_final_entrez$ENTREZID, organism = 'hsa', keyType = 'kegg')
metabolism_pathways_filtered <- metabolism_pathways[grep(paste(metabolism_terms, collapse = "|"), metabolism_pathways$Description),]

# Filter for top 10 pathways
top10_metabolism_upregulated <- head(metabolism_pathways_filtered@result[order(metabolism_pathways_filtered@result$p.adjust),], 10)
top10_metabolism_downregulated <- head(metabolism_pathways_filtered@result[order(metabolism_pathways_filtered@result$p.adjust, decreasing = TRUE)], 10)

# Plot metabolism pathways
dotplot_metabolism_upregulated <- dotplot(as.data.frame(top10_metabolism_upregulated), showCategory = 10) +
  ggtitle("Top 10 Upregulated Metabolism Pathways") +
  theme_minimal()
ggsave(output_files$metabolism_upregulated_top10, plot = dotplot_metabolism_upregulated, device = "pdf", width = 10, height = 10, dpi = 300)

dotplot_metabolism_downregulated <- dotplot(as.data.frame(top10_metabolism_downregulated), showCategory = 10) +
  ggtitle("Top 10 Downregulated Metabolism Pathways") +
  theme_minimal()
ggsave(output_files$metabolism_downregulated_top10, plot = dotplot_metabolism_downregulated, device = "pdf", width = 10, height = 10, dpi = 300)

# End of script
message("Analysis complete. Results are saved in the results directory.")

