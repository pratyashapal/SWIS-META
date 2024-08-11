# Load required libraries
library(MetaboAnalystR)
library(yaml)

# Load parameters from config file
config <- yaml::read_yaml("/users/2875659p/sharedscratch/rna_seq/metabolomics/meta/config.yaml")

# Create objects for storing processed data
mSet <- InitDataObjects("mass_table", "mummichog", FALSE)

# Set parameters from config file
mSet <- SetPeakFormat(mSet, "colu")
mSet <- UpdateInstrumentParameters(mSet, config$ppm, config$mode, config$include_internal_standards, config$retention_time_window)
mSet <- SetRTincluded(mSet, config$retention_time_unit)

# Read Peak table
mSet <- Read.TextData(mSet, config$input_csv, "mpt", "disc")
mSet <- SanityCheckMummichogData(mSet)

# Replace minimum value and data filtration
mSet <- ReplaceMin(mSet)
mSet <- SanityCheckMummichogData(mSet)
mSet <- FilterVariable(mSet, config$filter_method, config$filter_threshold, config$filter_by_fold_change, config$top_percentage, config$filter_by_stat)

# Perform data normalization
mSet <- PreparePrenormData(mSet)
mSet <- Normalization(mSet, "NULL", config$normalization_method, "NULL", ratio=config$ratio, ratioNum=config$ratio_num)
mSet <- PlotNormSummary(mSet, "norm_0_", config$output_format, config$output_dpi, width=NA)
mSet <- PlotSampleNormSummary(mSet, "snorm_0_", config$output_format, config$output_dpi, width=NA)

# Perform functional analysis with mummichog algorithm
mSet <- SetPeakEnrichMethod(mSet, "mum", "v2")
mSet <- PreparePeakTable4PSEA(mSet)
mSet <- SetMummichogPval(mSet, config$pval_cutoff)
mSet <- PerformPSEA(mSet, config$organism, config$database_version, config$num_permutations, config$random_seed)

# Plot and save results
output_prefix <- config$output_prefix
output_format <- config$output_format
output_dpi <- config$output_dpi

mSet <- PlotPeaks2Paths(mSet, output_prefix, output_format, output_dpi, width=NA)

