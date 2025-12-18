# R SCRIPT: Sleuth Analysis Example (Workflow C)

# --- 1. Packages and Setup ---
# Primary package for Differential Transcript Expression (DTE) analysis
library(sleuth)

# --- 2. Core Functions ---

# Define a data frame mapping sample IDs to Kallisto output files
# s2c <- data.frame(sample=..., condition=..., path=...) 

# 1. Prepare the Sleuth object
# Reads Kallisto abundance files (.h5) and prepares the data for modeling.
so <- sleuth_prep(
    s2c, 
    target_mapping = NULL, 
    read_bootstrap_summary = TRUE
)

# 2. Fit the 'full' model (includes condition effect)
so <- sleuth_fit(so, ~ condition, 'full')

# 3. Fit the 'reduced' model (baseline without condition)
so <- sleuth_fit(so, ~ 1, 'reduced')

# 4. Perform the Likelihood Ratio Test (LRT)
# This compares the two models to identify transcripts where 'condition' is significant.
so <- sleuth_lrt(so, 'reduced', 'full')

# 5. Get the final results table
sleuth_table <- sleuth_results(so, 'reduced:full', test_type = 'lrt')