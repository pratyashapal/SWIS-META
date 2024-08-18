# Loading necessary libraries
library(dplyr)
library(yaml)

# Loading configuration
config <- yaml::read_yaml("/SWIS-META/metabolomics/config_main.yaml")

output_dir <- file.path(config$results, config$mode)
setwd(output_dir)

results_dir <- config$results

# Defining file paths for t-test results
pos_ttest_file <- config$pos_ttest
neg_ttest_file <- config$neg_ttest

# Processing t-test results
pos <- read.csv(pos_ttest_file)
neg <- read.csv(neg_ttest_file)

pos <- pos[which(pos$p.value < 0.05),]
neg <- neg[which(neg$p.value < 0.05),]

comp <- pos[,1]
comp_n <- neg[,1]

strip <- function(com) {
  sub("[.;].*", "", com)
}

comp_p_stip <- sapply(comp, strip)
comp_n_stip <- sapply(comp_n, strip)

pos$X <- comp_p_stip
neg$X <- comp_n_stip

combined <- bind_rows(pos, neg)

filtered_combined <- combined %>%
  group_by(X) %>%
  filter(p.value == max(p.value)) %>%
  ungroup()

write.csv(filtered_combined, file.path(results_dir, "filtered_combined.csv"), row.names = FALSE)

# Defining file paths for volcano results
pos_volcano_file <- config$pos_volcano
neg_volcano_file <- config$neg_volcano

# Processing volcano results
pos_vol <- read.csv(pos_volcano_file)
neg_vol <- read.csv(neg_volcano_file)

pos_vol <- pos_vol[which(pos_vol$raw.pval < 0.05),]
neg_vol <- neg_vol[which(neg_vol$raw.pval < 0.05),]

comp_vol <- pos_vol[,1]
comp_n_vol <- neg_vol[,1]

vol_p_stip <- sapply(comp_vol, strip)
vol_n_stip <- sapply(comp_n_vol, strip)

pos_vol$X <- vol_p_stip
neg_vol$X <- vol_n_stip

combined_vol <- bind_rows(pos_vol, neg_vol)

filtered_combined_vol <- combined_vol %>%
  group_by(X) %>%
  filter(raw.pval == max(raw.pval)) %>%
  slice(1) %>%
  ungroup()

write.csv(filtered_combined_vol, file.path(results_dir, "filtered_combined_log.csv"), row.names = FALSE)

