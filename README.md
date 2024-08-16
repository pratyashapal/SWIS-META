# SWIS-META
SWIS-META is a comprehensive analysis pipeline designed for integrating transcriptomic and metabolomic data using Snakemake. This document provides detailed instructions for installation, setup, and execution of the SWIS-META pipeline.

Installation and Setup
To get started with SWIS-META, follow these steps:

Clone the Repository:
git clone https://github.com/pratyashapal/SWIS-META.git

Navigate to the Repository Directory:
cd SWIS-META

Create and Activate the Conda Environment:
conda env create -f environment.yml
conda activate swis-meta

This setup creates a dedicated Conda environment with all necessary dependencies for running SWIS-META. Ensure that you have the latest R packages and MetaboAnalyst installed separately for compatibility.

Implementation
The SWIS-META workflow is managed by Snakemake, which structures the process into discrete, traceable steps. The workflow is divided into three main components:

Transcriptomics
Metabolomics
Integrated Analysis
Each component has its own set of Snakefiles and configuration files, which can be edited to customize the analysis.

Transcriptomics Pipeline
The transcriptomics component processes raw RNA-seq data through several steps:

Quality Control:

Snakefile_1: Performs FastQC and MultiQC.
Snakefile_2: Trims reads and regenerates QC reports.
Xenograft Check:

Snakefile_3: Filters xenograft reads.
Alignment and Counting:

Snakefile_4: Aligns reads and generates gene count matrix.
Differential Expression Analysis:

Snakefile_5: Performs DE analysis, including PCA, heatmaps, volcano plots, and pathway enrichment.
Steps Overview

Quality Control: FastQC and MultiQC reports are generated to assess read quality. Optional trimming is applied if specified.
Xenograft Check: Xengsort tool is used to filter human reads from xenografts.
Alignment and Counting: Reads are aligned using HISAT2, and counts are quantified with FeatureCounts.
Differential Expression: DE analysis is conducted using DESeq2, providing various visualizations and enrichment results.


Metabolomics Pipeline
The metabolomics component automates several analytical steps:

Data Preprocessing:

Snakefile_pre: Processes raw LC-MS data and generates a feature table.
Functional Analysis:

Snakefile_1: Performs functional analysis using MetaboAnalystR.
Filtering:

Snakefile_2: Filters compound data.
Statistical Analysis:

Snakefile_3: Conducts statistical analyses, including fold change and volcano plots.
Pathway Analysis:

Snakefile_5: Analyzes pathways and generates relevant visualizations.

Steps Overview

Data Preprocessing: Extracts and aligns peaks from LC-MS data.
Functional Analysis: Uses Mummichog for biological inference of mass peaks.
Filtering: Maps KEGG IDs to compound names and filters data.
Statistical Analysis: Performs various statistical tests and visualizations.
Pathway Analysis: Identifies significant pathways and generates plots.


Integrated Pipeline
The integrated pipeline combines transcriptomic and metabolomic data for comprehensive pathway analysis:

Data Preparation:

Gene symbols are mapped to Entrez IDs, and metabolomics compound data is integrated.
Pathway Enrichment Analysis:

KEGG pathway enrichment is performed using clusterProfiler.
Data Integration and Visualization:

Merges gene and compound data, visualizing pathway activity through bar plots and heatmaps.
Steps Overview
Data Preparation: Maps gene symbols and integrates compound data.
Pathway Enrichment: Identifies significant pathways linking genes and compounds.
Data Integration: Merges and visualizes average log fold changes across pathways.
