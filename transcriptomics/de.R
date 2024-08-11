# scripts/analyze.R
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

# Load configuration
config <- list(
  paths = list(
    count_matrix = "/users/2875659p/sharedscratch/rna_seq/transcriptomics/snake_out_final/gene_symbol_counts_matrix (1).csv",
    output_directory = "output/",
    plots_directory = "plots/"
  )
)

# Define paths from config
count_matrix_path <- config$paths$count_matrix
output_directory <- config$paths$output_directory
plots_directory <- config$paths$plots_directory

# Ensure directories exist
dir.create(output_directory, showWarnings = FALSE)
dir.create(plots_directory, showWarnings = FALSE)

# Load the count matrix
count_matrix <- read.csv(count_matrix_path, row.names = 1, header = TRUE)

# Define sample information
col_data <- data.frame(
  sample = colnames(count_matrix),
  condition = c(rep("bm", 4), rep("cns", 4))
)
rownames(col_data) <- col_data$sample
col_data$condition <- factor(col_data$condition)

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = col_data, design = ~ condition)
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)

# Get results
res <- results(dds)
sig_res <- subset(res, padj < 0.05)
sig_genes_df <- as.data.frame(sig_res)

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
  geom_text(data = top_genes, aes(label = label), vjust = 1.5, color = "black") +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted p-value") +
  theme_minimal()
ggsave(filename = file.path(plots_directory, "volcano_plot.pdf"), plot = volcano_plot, device = "pdf",  width = 10, height = 6, dpi = 300)

# Access normalized counts
normalized_counts <- counts(dds, normalized = TRUE)
log_norm_counts <- log2(normalized_counts + 1)

# PCA analysis
pca_res <- prcomp(t(log_norm_counts))
pca_data <- data.frame(Sample = rownames(pca_res$x),
                       PC1 = pca_res$x[,1],
                       PC2 = pca_res$x[,2],
                       Condition = col_data$condition[rownames(pca_res$x)])

pca_log_norm_counts <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle("PCA of Gene Expression") +
  xlab(paste("PC1 -", round(summary(pca_res)$importance[2,1] * 100, 2), "% Variance")) +
  ylab(paste("PC2 -", round(summary(pca_res)$importance[2,2] * 100, 2), "% Variance")) +
  scale_color_manual(values = c("bm" = "blue", "cns" = "red")) +
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))
ggsave(filename = file.path(plots_directory, "pca_log_norm_counts.pdf"), plot = pca_log_norm_counts, device = "pdf",  width = 10, height = 6, dpi = 300)

# PCA plot with labels
pca_log_norm_counts_with_labels <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition, label = Sample)) +
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
ggsave(filename = file.path(plots_directory, "pca_log_norm_counts_with_labels.pdf"), plot = pca_log_norm_counts_with_labels, device = "pdf",  width = 12, height = 8, dpi = 300)

# Calculate correlation matrix and plot
sample_cor_matrix <- cor(log_norm_counts, method = "pearson")
sample_corr_matrix <- pheatmap(sample_cor_matrix,
                               clustering_distance_rows = "euclidean",
                               clustering_distance_cols = "euclidean",
                               scale = "none",
                               color = colorRampPalette(c("blue", "white", "red"))(100),
                               main = "Sample-to-Sample Correlation",
                               show_rownames = TRUE,
                               show_colnames = TRUE,
                               annotation_legend = TRUE)

ggsave(filename = file.path(plots_directory, "sample_corr_matrix.pdf"), plot = sample_corr_matrix, device = "pdf", width = 10, height = 10, dpi = 300)

# VST PCA analysis
vsd <- vst(dds, blind = FALSE)
pca_vst_res <- prcomp(t(assay(vsd)))
pca_vst_data <- data.frame(Sample = rownames(pca_vst_res$x),
                            PC1 = pca_vst_res$x[,1],
                            PC2 = pca_vst_res$x[,2],
                            Condition = col_data$condition[rownames(pca_vst_res$x)])

pca_var_stab <- ggplot(pca_vst_data, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle("PCA of VST Normalized Counts") +
  xlab(paste("PC1 -", round(summary(pca_vst_res)$importance[2,1] * 100, 2), "% Variance")) +
  ylab(paste("PC2 -", round(summary(pca_vst_res)$importance[2,2] * 100, 2), "% Variance")) +
  scale_color_manual(values = c("bm" = "blue", "cns" = "red")) +
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))
ggsave(filename = file.path(plots_directory, "pca_var_stab.pdf"), plot = pca_var_stab, device = "pdf",  width = 10, height = 6, dpi = 300)

# Heatmap of significant genes
heatmap <- pheatmap(assay(vsd)[rownames(sig_genes_df), ],
                    clustering_distance_rows = "euclidean",
                    clustering_distance_cols = "euclidean",
                    scale = "row",
                    color = colorRampPalette(c("blue", "white", "red"))(100),
                    main = "Heatmap of Significant Genes",
                    show_rownames = TRUE,
                    show_colnames = TRUE,
                    annotation_legend = TRUE)

ggsave(filename = file.path(plots_directory, "heatmap_sig_genes.pdf"), plot = heatmap, device = "pdf",  width = 12, height = 10, dpi = 300)

# Heatmap of top 25 DEGs
top25_genes <- head(sig_genes_df[order(sig_genes_df$padj), ], 25)
heatmap_top25 <- pheatmap(assay(vsd)[rownames(top25_genes), ],
                          clustering_distance_rows = "euclidean",
                          clustering_distance_cols = "euclidean",
                          scale = "row",
                          color = colorRampPalette(c("blue", "white", "red"))(100),
                          main = "Heatmap of Top 25 DEGs",
                          show_rownames = TRUE,
                          show_colnames = TRUE,
                          annotation_legend = TRUE)

ggsave(filename = file.path(plots_directory, "heatmap_top25_degs.pdf"), plot = heatmap_top25, device = "pdf",  width = 12, height = 10, dpi = 300)

# GO BP Enrichment
go_bp <- enrichGO(gene = rownames(sig_genes_df),
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.05)

go_bp_dotplot <- dotplot(go_bp, showCategory=30) + ggtitle("GO Biological Process Enrichment")
ggsave(filename = file.path(plots_directory, "go_bp_dotplot.pdf"), plot = go_bp_dotplot, device = "pdf",  width = 10, height = 6, dpi = 300)

# GO MF Enrichment
go_mf <- enrichGO(gene = rownames(sig_genes_df),
                  OrgDb = org.Hs.eg.db,
                  ont = "MF",
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.05)

go_mf_dotplot <- dotplot(go_mf, showCategory=30) + ggtitle("GO Molecular Function Enrichment")
ggsave(filename = file.path(plots_directory, "go_mf_dotplot.pdf"), plot = go_mf_dotplot, device = "pdf",  width = 10, height = 6, dpi = 300)

# GO CC Enrichment
go_cc <- enrichGO(gene = rownames(sig_genes_df),
                  OrgDb = org.Hs.eg.db,
                  ont = "CC",
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.05)

go_cc_dotplot <- dotplot(go_cc, showCategory=30) + ggtitle("GO Cellular Component Enrichment")
ggsave(filename = file.path(plots_directory, "go_cc_dotplot.pdf"), plot = go_cc_dotplot, device = "pdf", width = 10, height = 6, dpi = 300)

# KEGG Enrichment
kegg <- enrichKEGG(gene = rownames(sig_genes_df),
                    organism = 'hsa',
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.05)

kegg_dotplot <- dotplot(kegg, showCategory=30) + ggtitle("KEGG Pathway Enrichment")
ggsave(filename = file.path(plots_directory, "kegg_dotplot.pdf"), plot = kegg_dotplot, device = "pdf", width = 10, height = 6, dpi = 300)


