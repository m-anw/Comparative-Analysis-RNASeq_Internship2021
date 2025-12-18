# R SCRIPT: visualization_volcano_plot.R
# DESCRIPTION: Generates a Volcano Plot from DESeq2 output, highlighting significant genes.

# --- 1. Data Preparation ---

# Load data. Assume the first column contains the Gene IDs (used as row names).
res <- read.csv("DESeq2_results.csv", row.names = 1) 

# Remove any rows where the adjusted p-value (padj) is missing (NA), 
# as these genes cannot be plotted on the Y-axis.
res <- res[!is.na(res$padj), ]

# Calculate -log10(padj) for the Y-axis of the volcano plot.
res$neg_log10_padj <- -log10(res$padj)

# --- 2. Define Groups and Plot ---

# Define cutoffs used for coloring.
P_CUTOFF <- 0.05
FC_CUTOFF <- 0.0 # Log2(1) = 0. We are using the sign to define direction.

# Identify significant genes based on the cutoffs.
# Up-regulated genes (in 38C vs 34C): padj < 0.05 AND negative log2FC.
up_idx   <- res$padj < P_CUTOFF & res$log2FoldChange < FC_CUTOFF 
# Down-regulated genes: padj < 0.05 AND positive log2FC.
down_idx <- res$padj < P_CUTOFF & res$log2FoldChange > FC_CUTOFF 

# 1. Start plot: Plot all genes as gray background points.
plot(res$log2FoldChange, res$neg_log10_padj, 
     pch = 20, 
     col = "gray", 
     xlab = "Log2(Fold Change)", 
     ylab = "-Log10(Adjusted P-Value)",
     main = "Volcano Plot (C38 vs C34)") # Added main title for clarity

# 2. Overlay significant up-regulated genes (Red).
points(res$log2FoldChange[up_idx], res$neg_log10_padj[up_idx], 
       pch = 20, col = "red")

# 3. Overlay significant down-regulated genes (Light Blue).
points(res$log2FoldChange[down_idx], res$neg_log10_padj[down_idx], 
       pch = 20, col = "lightskyblue")

# 4. Add legend showing total counts.
legend("topright",
       legend = c(paste("up =", sum(up_idx)), 
                  paste("down =", sum(down_idx))),
       col = c("red", "lightskyblue"),
       pch = 20, 
       bty = "n") # 'n' means no box around the legend.