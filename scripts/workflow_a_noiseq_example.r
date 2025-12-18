# R SCRIPT: NOISeq Analysis Example (Workflow A)

# --- 1. Packages ---
# Bioconductor package for non-parametric DGE analysis
library(NOISeq) 
# library(edgeR) # Often needed for TMM normalization

# --- 2. Core Functions ---

# Data preparation steps (loading counts and metadata)
# countData <- read.table("htseq_results.txt", ...) 
# colData <- data.frame(condition = c("C34", "C38", ...))

# 1. Create the NOISeq object
mydata <- readData(
    data = countData, 
    factors = colData, 
    type = "count"
)

# 2. Filter low-count features (e.g., keep genes with CPM > 1 in at least 1 sample)
filt_data <- filtered.data(
    mydata, 
    factor = 1, 
    norm = FALSE, 
    method = 1, 
    cpm = 1 
)

# 3. Normalize the data (TMM - Trimmed Mean of M-values)
# TMM was specified in your report as the normalization method for this pipeline.
mydata_norm <- tmm(filt_data, long = 1000, lc = 0, k = 0)

# 4. Run the NOISeq differential expression test
# Uses the non-parametric distribution ("n")
mynoiseq <- noiseq(
    mydata_norm, 
    k = 0.5, 
    norm = "n", 
    factor = "condition", 
    replicates = "technical" 
)

# 5. Extract results with custom cutoffs
# q = Probability cutoff (PROB > 0.7, as per your documentation)
# M = Log2(Fold Change) cutoff (M value threshold, set to 1.0 for a 2-fold change)
noiseq_results <- degenes(
    mynoiseq, 
    q = 0.7,  
    M = 1.0   
) 
# write.csv(noiseq_results, file="noiseq_DGE_results.csv")