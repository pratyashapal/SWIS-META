# Load required libraries
library(MetaboAnalystR)
library(OptiLCMS)
library(yaml)

# Define the path to the config file
config_file <- "config/config_pre.yaml"

# Read the configuration file
config <- yaml::read_yaml(config_file)

# Extract input files
input_dirs <- config$input_dirs
bm_pos <- list.files(input_dirs$bm_pos, full.names = TRUE)
bm_neg <- list.files(input_dirs$bm_neg, full.names = TRUE)
brain_pos <- list.files(input_dirs$brain_pos, full.names = TRUE)
brain_neg <- list.files(input_dirs$brain_neg, full.names = TRUE)
qc_pos <- list.files(input_dirs$qc_pos, full.names = TRUE)
qc_neg <- list.files(input_dirs$qc_neg, full.names = TRUE)

# Combine all files
all_files <- c(bm_pos, bm_neg, brain_pos, brain_neg, qc_pos, qc_neg)

# Extract parameters
roi_params <- config$roi_params
rt_idx <- roi_params$rt_idx
rm_conts <- roi_params$rm_conts

annotation_params <- config$annotation_params
polarity <- annotation_params$polarity
mz_abs_add <- annotation_params$mz_abs_add

optimization_params <- config$optimization_params
method <- optimization_params$method
ncore <- optimization_params$ncore

# ROI Extraction
mSet <- PerformROIExtraction(datapath = all_files, 
                             rt.idx = rt_idx, 
                             rmConts = rm_conts)

# Save intermediate result if necessary
saveRDS(mSet, "results/roi_extraction.rds")

# Parameter Optimization
best_params <- PerformParamsOptimization(mSet, 
                                         param = NULL, 
                                         method = method, 
                                         ncore = ncore)

# Save intermediate result if necessary
saveRDS(best_params, "results/best_params.rds")

# Import Raw Data
mSet <- ImportRawMSData(path = all_files, 
                        plotSettings = SetPlotParam(Plot = TRUE))

# Save intermediate result if necessary
saveRDS(mSet, "results/raw_data_import.rds")

# Peak Profiling
mSet <- PerformPeakProfiling(mSet, Params = best_params, plotSettings = SetPlotParam(Plot=TRUE))

# Save intermediate result if necessary
saveRDS(mSet, "results/peak_profiling.rds")

# Peak Annotation
annParams <- SetAnnotationParam(polarity = polarity, mz_abs_add = mz_abs_add)
mSet <- PerformPeakAnnotation(mSet, annParams)

# Export annotated peak lists
Export.Annotation(mSet, file="results/annotated_peaklist_positive_raw.csv")

# Adjust polarity for negative annotation if needed
annParams_neg <- SetAnnotationParam(polarity = "negative", mz_abs_add = mz_abs_add)
mSet <- PerformPeakAnnotation(mSet, annParams_neg)
Export.Annotation(mSet, file="results/annotated_peaklist_negative_raw.csv")

# Format and Export Peak List for MetaboAnalyst
mSet <- FormatPeakList(mSet, annParams, filtIso = FALSE, filtAdducts = FALSE, missPercent = 1)
Export.PeakTable(mSet, file="results/metaboanalyst_input_positive_raw.csv")

mSet <- FormatPeakList(mSet, annParams, filtIso = FALSE, filtAdducts = FALSE, missPercent = 1)
Export.PeakTable(mSet, file="results/metaboanalyst_input_negative_raw.csv")

