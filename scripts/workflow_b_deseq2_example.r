# R SCRIPT: DESeq2 Analysis Example (Workflow B)

# --- 1. Packages and Setup ---
# Bioconductor package for Differential Gene Expression (DGE) analysis
library(DESeq2) 
# Standard package for data visualization
library(ggplot2) 

# --- 2. Core Functions ---

# Load raw count data and sample metadata
# countData <- read.table("featurecounts_results.txt", ...)
# colData <- data.frame(condition = c("C34", "C38", ...))

# 1. Create the DESeq2 object
# This combines the counts, metadata, and the experimental design (~ condition)
dds <- DESeqDataSetFromMatrix(
    countData = countData,
    colData = colData,
    design = ~ condition
)

# 2. Run the main DESeq2 pipeline
dds <- DESeq(dds)

# 3. Extract the results table (38C vs 34C)
res <- results(dds, contrast=c("condition", "C38", "C34"), alpha=0.05)

# --- 3. Visualization Example ---

# 4. Prepare data for PCA plotting 
# Note: Data must be transformed (e.g., using rlog or normalized log counts) 
#       to ensure variance is not dependent on the mean.
# rld <- rlog(dds, blind=FALSE) 

# plotPCA(rld, intgroup="condition")