#!/bin/bash
# --- ATAC-seq Pipeline for REAL DATA (hg38) ---
# This script executes all steps from FASTQ to Peak Calling using the installed tools.

# Fail if any command fails
set -e
# Print commands as they are executed
set -x

# --- CONFIGURATION ---
SAMPLE_LIST=("ERR15487529") # CHANGED to your real sample accession!
# The name of the Bowtie2 index files built earlier
GENOME_INDEX_PATH="/mnt/d/NGS/reference/hg38" 
# Effective genome size for MACS2 on hg38
GENOME_SIZE="2.7e9"
# Number of threads to use for parallel processing
THREADS=4

# --- DIRECTORY SETUP ---
RAW_DATA_DIR="./raw_data"
QC_DIR="./qc_reports"
TRIM_DIR="./trimmed_data"
ALIGN_DIR="./aligned_data"
FILTER_DIR="./filtered_data"
PEAKS_DIR="./peaks"
BIGWIG_DIR="./bigwig_tracks"

# Create all necessary output directories
mkdir -p $QC_DIR $TRIM_DIR $ALIGN_DIR $FILTER_DIR $PEAKS_DIR $BIGWIG_DIR multiqc_report

# --- STEP 0: CHECK DATA ---
echo "--- STEP 0: VERIFYING RAW DATA FILES ---"
if [ ! -f "$RAW_DATA_DIR/${SAMPLE_LIST[0]}_1.fastq.gz" ]; then
    echo "ERROR: Raw data file $RAW_DATA_DIR/${SAMPLE_LIST[0]}_1.fastq.gz not found."
    echo "Please download the real FASTQ files and name them: ${SAMPLE_LIST[0]}_1.fastq.gz and ${SAMPLE_LIST[0]}_2.fastq.gz"
    exit 1
fi

# --- MAIN PIPELINE LOOP ---
for SAMPLE in "${SAMPLE_LIST[@]}"; do
    echo "Processing sample: ${SAMPLE}..."

    R1="$RAW_DATA_DIR/${SAMPLE}_1.fastq.gz"
    R2="$RAW_DATA_DIR/${SAMPLE}_2.fastq.gz"

    # --- STEP 1: INITIAL QC (FASTQC) ---
    echo "--- STEP 1: INITIAL QC ---"
    fastqc -o $QC_DIR $R1 $R2

    # --- STEP 2: TRIMMING (Trim Galore!) ---
    echo "--- STEP 2: TRIMMING ADAPTERS ---"
    trim_galore --paired --nextera \
        --output_dir $TRIM_DIR \
        $R1 $R2

    TRIM_R1="$TRIM_DIR/${SAMPLE}_1_val_1.fq.gz"
    TRIM_R2="$TRIM_DIR/${SAMPLE}_2_val_2.fq.gz"

    # --- STEP 3: ALIGNMENT (Bowtie2) ---
    echo "--- STEP 3: ALIGNMENT TO HG38 ---"
    # Aligns and pipes SAM output directly to BAM conversion/sort
    bowtie2 -x $GENOME_INDEX_PATH -1 $TRIM_R1 -2 $TRIM_R2 \
        --very-sensitive-local \
        -p $THREADS | samtools view -@ $THREADS -bS - > $ALIGN_DIR/${SAMPLE}.bam

    # --- STEP 4: FILTERING AND SORTING ---
    echo "--- STEP 4: FILTERING AND SORTING ---"
    
    # Sort BAM file by read name (for filtering mitochondrial reads later)
    samtools sort -@ $THREADS -n -o $ALIGN_DIR/${SAMPLE}_sorted_name.bam $ALIGN_DIR/${SAMPLE}.bam

    # Filter: -F 1804 (Removes unmapped, secondary, PCR/Optical duplicates, failed QC)
    # Then remove chrM reads, and finally coordinate sort
    samtools view -@ $THREADS -b -F 1804 $ALIGN_DIR/${SAMPLE}_sorted_name.bam | \
    samtools view -@ $THREADS -h - | grep -v 'chrM' | \
    samtools view -@ $THREADS -b - | \
    samtools sort -@ $THREADS -o $FILTER_DIR/${SAMPLE}_clean.bam -
    
    # Index the final clean BAM file
    samtools index -@ $THREADS $FILTER_DIR/${SAMPLE}_clean.bam
    CLEAN_BAM="$FILTER_DIR/${SAMPLE}_clean.bam"

    # --- STEP 5: PEAK CALLING (MACS2) ---
    echo "--- STEP 5: PEAK CALLING ---"
    macs2 callpeak \
        -t $CLEAN_BAM \
        -f BAMPE \
        -g $GENOME_SIZE \
        -n ${SAMPLE}_ATAC \
        --outdir $PEAKS_DIR \
        -q 0.05 \
        --nomodel \
        --keep-dup all

    # --- STEP 6: BIGWIG GENERATION (deepTools) ---
    echo "--- STEP 6: GENERATING BIGWIG TRACKS (RPKM normalized) ---"
    # Normalizing by RPKM (Restored normalization for real data)
    bamCoverage -b $CLEAN_BAM \
                -o $BIGWIG_DIR/${SAMPLE}_coverage.bw \
                --normalizeUsing RPKM \
                --binSize 10 \
                -p $THREADS \
                --ignoreForNormalization chrX chrY

done

# --- STEP 7: FINAL QC REPORT (MultiQC) ---
echo "--- STEP 7: COMPILING FINAL MULTIQC REPORT ---"
multiqc . -o multiqc_report

echo "--- PIPELINE COMPLETE! ---"
echo "Outputs are in: peaks/, bigwig_tracks/, and multiqc_report/"
