# Loading required libraries
library(MetaboAnalystR)
library(yaml)

# Loading parameters from config file
config <- yaml::read_yaml("/SWIS-META/metabolomics/config_1.yaml")

# Determining the output directory based on ionization mode
output_dir <- file.path(config$results, config$mode)

# Creating the output directory if it doesn't exist
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# Setting the working directory to the output directory
setwd(output_dir)

# Creating objects for storing processed data
mSet <- InitDataObjects("mass_table", "mummichog", FALSE)

# Setting parameters from config file
mSet <- SetPeakFormat(mSet, "colu")
mSet <- UpdateInstrumentParameters(mSet, config$ppm, config$mode, config$include_internal_standards, config$retention_time_window)
mSet <- SetRTincluded(mSet, config$retention_time_unit)

# Reading Peak table
mSet <- Read.TextData(mSet, config$input_csv, "mpt", "disc")
mSet <- SanityCheckMummichogData(mSet)

# Replacing minimum value and data filtration
mSet <- ReplaceMin(mSet)
mSet <- SanityCheckMummichogData(mSet)
mSet <- FilterVariable(mSet, config$filter_method, config$filter_threshold, config$filter_by_fold_change, config$top_percentage, config$filter_by_stat)

# Performing data normalization
mSet <- PreparePrenormData(mSet)
mSet <- Normalization(mSet, "NULL", config$normalization_method, "NULL", ratio=config$ratio, ratioNum=config$ratio_num)
mSet <- PlotNormSummary(mSet, paste0("norm_0_", config$output_format), config$output_format, config$output_dpi, width=NA)
mSet <- PlotSampleNormSummary(mSet, paste0("snorm_0_", config$output_format), config$output_format, config$output_dpi, width=NA)

# Performing functional analysis with mummichog algorithm
mSet <- SetPeakEnrichMethod(mSet, "mum", "v2")
mSet <- PreparePeakTable4PSEA(mSet)
mSet <- SetMummichogPval(mSet, config$pval_cutoff)
mSet <- PerformPSEA(mSet, config$organism, config$database_version, config$num_permutations, config$random_seed)

# Plotting and saving results
mSet <- PlotPeaks2Paths(mSet, paste0(config$output_prefix, ".", config$output_format), config$output_format, config$output_dpi, width=NA)
