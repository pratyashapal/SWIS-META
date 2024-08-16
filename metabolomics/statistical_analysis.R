library(ggplot2)
library(iheatmapr)
library(randomForest)
library(pls)
library(MetaboAnalystR)
library(yaml)

# Configuration
config <- yaml::read_yaml("/SWIS-META/metabolomics/config_main.yaml")

output_dir <- file.path(config$results, config$mode)
setwd(output_dir)


results_dir <- config$results

# Define input and output paths based on ionization mode
if (config$mode == "positive") {
  ionization_dir <- file.path(results_dir, "positive")
  filter_table <- config$filter_table_positive
} else if (config$mode == "negative") {
  ionization_dir <- file.path(results_dir, "negative")
  filter_table <- config$filter_table_negative
} else {
  stop("Unknown ionization mode specified.")
}

# Read the filter table
filter_table <- read.csv(filter_table)

# Analysis and plotting
mSet <- InitDataObjects("pktable", "stat", FALSE)
mSet <- Read.TextData(mSet, filter_table, "rowu", "disc")
mSet <- SanityCheckData(mSet)
mSet <- ReplaceMin(mSet)
mSet <- SanityCheckData(mSet)
mSet <- FilterVariable(mSet, "F", 10, "iqr", 25, "mean", 0)
mSet <- PreparePrenormData(mSet)
mSet <- Normalization(mSet, "NULL", "LogNorm", "NULL", ratio=FALSE, ratioNum=20)
mSet <- PlotNormSummary(mSet, file.path(ionization_dir, "norm_0_.png"), "png", 72, width=NA)
mSet <- PlotSampleNormSummary(mSet, file.path(ionization_dir, "snorm_0_.png"), "png", 72, width=NA)
mSet <- FC.Anal(mSet, 2.0, 0, FALSE)
mSet <- PlotFC(mSet, file.path(ionization_dir, "fc_0_.png"), "png", 72, width=NA)
mSet <- Ttests.Anal(mSet, F, 0.05, FALSE, TRUE, "fdr", FALSE)
mSet <- PlotTT(mSet, file.path(ionization_dir, "tt_0_.png"), "png", 72, width=NA)
mSet <- Volcano.Anal(mSet, FALSE, 2.0, 0, F, 0.1, TRUE, "raw")
mSet <- PlotVolcano(mSet, file.path(ionization_dir, "volcano_0_.png"), 1, 0, "png", 72, width=NA, -1)
mSet <- Volcano.Anal(mSet, FALSE, 2.0, 0, F, 0.1, TRUE, "raw")
mSet <- PlotVolcano(mSet, file.path(ionization_dir, "volcano_0_.png"), 1, 0, "png", 72, width=NA, -1)
mSet <- PlotCorrHeatMap(mSet, file.path(ionization_dir, "corr_1_.png"), "png", 72, width=NA, "row", "pearson", "bwm", F, F, 6, 10, 0.0)
mSet <- PlotCorrHeatMap(mSet, file.path(ionization_dir, "corr_2_.png"), "png", 72, width=NA, "col", "pearson", "bwm", F, F, 6, 10, 0.0)
mSet <- PCA.Anal(mSet)
mSet <- PlotPCAPairSummary(mSet, file.path(ionization_dir, "pca_pair_0_.png"), "png", 72, width=NA, 5)
mSet <- PlotPCAScree(mSet, file.path(ionization_dir, "pca_scree_0_.png"), "png", 72, width=NA, 5)
mSet <- PlotPCA2DScore(mSet, file.path(ionization_dir, "pca_score2d_0_.png"), "png", 72, width=NA, 1, 2, 0.95, 0, 0, "na")
mSet <- PlotPCALoading(mSet, file.path(ionization_dir, "pca_loading_0_.png"), "png", 72, width=NA, 1, 2)
mSet <- PlotPCABiplot(mSet, file.path(ionization_dir, "pca_biplot_0_.png"), "png", 72, width=NA, 1, 2)
mSet <- PlotPCA3DLoading(mSet, file.path(ionization_dir, "pca_loading3d_0_.json"), "json", 1, 2, 3)
mSet <- PLSR.Anal(mSet, reg=TRUE)
mSet <- PlotPLSPairSummary(mSet, file.path(ionization_dir, "pls_pair_0_.png"), "png", 72, width=NA, 5)
mSet <- PlotPLS2DScore(mSet, file.path(ionization_dir, "pls_score2d_0_.png"), "png", 72, width=NA, 1, 2, 0.95, 0, 0, "na")
mSet <- PlotPLS3DScoreImg(mSet, file.path(ionization_dir, "pls_score3d_0_.png"), "png", 72, width=NA, 1, 2, 3, 40)
mSet <- PlotPLSLoading(mSet, file.path(ionization_dir, "pls_loading_0_.png"), "png", 72, width=NA, 1, 2)
mSet <- PlotPLS3DLoading(mSet, file.path(ionization_dir, "pls_loading3d_0_.json"), "json", 1, 2, 3)
mSet <- PlotPLS.Imp(mSet, file.path(ionization_dir, "pls_imp_0_.png"), "png", 72, width=NA, "vip", "Comp. 1", 15, FALSE)
mSet <- RF.Anal(mSet, 500, 7, 1)
mSet <- PlotRF.Classify(mSet, file.path(ionization_dir, "rf_cls_0_.png"), "png", 72, width=NA)
mSet <- PlotRF.VIP(mSet, file.path(ionization_dir, "rf_imp_0_.png"), "png", 72, width=NA)
mSet <- PlotRF.Outlier(mSet, file.path(ionization_dir, "rf_outlier_0_.png"), "png", 72, width=NA)
mSet <- PlotHeatMap(mSet, file.path(ionization_dir, "heatmap_1_.png"), "png", 72, width=NA, "norm", "row", "euclidean", "ward.D", "bwm", 8, 8, 10.0, 0.02, 10, 10, T, T, NULL, T, F, T, T, T, T)
mSet <- PlotSubHeatMap(mSet, file.path(ionization_dir, "heatmap_2_.png"), "png", 72, width=NA, "norm", "row", "euclidean", "ward.D", "bwm", 8, 8, 10.0, 0.02, 10, 10, "tanova", 25, F, T, T, F)

