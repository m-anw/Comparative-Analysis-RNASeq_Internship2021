# MINIMAL BASH WORKFLOW (HISAT2, FeatureCounts, MultiQC)

# --- 1. Alignment (HISAT2) ---
# Aligns paired-end reads to the genome index.
# Example: HISAT2 alignment for a single sample (SRR12759173).

# hisat2 -x <index_path> -1 <R1_fastq> -2 <R2_fastq> -p 8 -S <output.sam>

echo "Running HISAT2 Alignment..."
hisat2 -x /path/to/genome_index \
       -1 /path/to/SRR12759173_1.fastq.gz \
       -2 /path/to/SRR12759173_2.fastq.gz \
       -p 8 \
       -S SRR12759173.sam 

# SAM to BAM conversion and sorting
# samtools view -bS SRR12759173.sam | samtools sort -o SRR12759173.sorted.bam

# --- 2. Quantification (FeatureCounts for Workflow B) ---
# Counts reads assigned to genomic features (genes/exons) using BAM files.

# featureCounts -p -T 8 -a <gtf_file> -o <counts.txt> <*.sorted.bam>

echo "Running FeatureCounts Quantification (Workflow B)..."
featureCounts -p \
              -T 8 \
              -a /path/to/annotation.gtf \
              -o featurecounts_results.txt \
              /path/to/all_samples/*.sorted.bam

# --- 3. Quality Control (MultiQC) ---
# Aggregates all log files (HISAT2, FeatureCounts, etc.) into a single report.

echo "Running MultiQC Aggregation..."
multiqc /path/to/logs_and_summaries -o multiqc_report_dir