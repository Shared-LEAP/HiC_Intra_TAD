#!/bin/bash
#SBATCH --job-name=HiC_step1_fastqc
#SBATCH --output=fastqc.log
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100G
#SBATCH --nodelist=n001

# ============================================================
# STEP 1: QUALITY CONTROL WITH FASTQC
# ============================================================
# Purpose: Assess quality of raw FASTQ files before processing
# Input: Raw .fastq.gz files
# Output: FastQC HTML reports for quality assessment
# ============================================================

echo "=========================================="
echo "STEP 1: FastQC Quality Control"
echo "Started at: $(date)"
echo "=========================================="

# ============================================================
# DIRECTORY SETUP
# ============================================================
WORK_DIR="/YOUR DIRECTORY HERE/4DN_pipeline"
FASTQ_DIR="/YOUR DIRECTORY HERE/Fastq.gz"
PROCESS_DIR="${WORK_DIR}/HiC_processing"
QC_DIR="${PROCESS_DIR}/qc/fastqc_raw"

# Tool paths
CONDA_BASE="/YOUR DIRECTORY HERE/anaconda3"
FASTQC="${CONDA_BASE}/envs/fastqc/bin/fastqc"

# Create directories
mkdir -p ${QC_DIR}
mkdir -p ${PROCESS_DIR}/logs

# Change to working directory
cd ${PROCESS_DIR}

echo ""
echo "Working directory: ${PROCESS_DIR}"
echo "FASTQ source: ${FASTQ_DIR}"
echo "QC output: ${QC_DIR}"
echo "Threads: ${SLURM_CPUS_PER_TASK}"

# ============================================================
# VERIFY FASTQ FILES
# ============================================================
echo ""
echo "=== Verifying FASTQ Files ==="

if [[ ! -d ${FASTQ_DIR} ]]; then
    echo "ERROR: FASTQ directory not found: ${FASTQ_DIR}"
    exit 1
fi

# Count available FASTQ files
FASTQ_COUNT=$(ls -1 ${FASTQ_DIR}/*.fastq.gz 2>/dev/null | wc -l)

if [[ ${FASTQ_COUNT} -eq 0 ]]; then
    echo "ERROR: No .fastq.gz files found in ${FASTQ_DIR}"
    exit 1
fi

echo "Found ${FASTQ_COUNT} FASTQ files to process"
echo ""
echo "Files to analyze:"
ls -lh ${FASTQ_DIR}/*.fastq.gz | awk '{print "  "$9" ("$5")"}'

# ============================================================
# RUN FASTQC
# ============================================================
echo ""
echo "=== Running FastQC ==="
echo "This may take 15-30 minutes per file depending on size..."
echo ""

# Process each FASTQ file
SUCCESS_COUNT=0
FAIL_COUNT=0

for FASTQ in ${FASTQ_DIR}/*.fastq.gz; do
    BASENAME=$(basename ${FASTQ} .fastq.gz)
    
    echo "-----------------------------------"
    echo "Processing: ${BASENAME}"
    echo "File: ${FASTQ}"
    
    # Check if file is readable and not empty
    if [[ ! -s ${FASTQ} ]]; then
        echo "  ✗ SKIPPED: File is empty or not readable"
        ((FAIL_COUNT++))
        continue
    fi
    
    # Run FastQC
    echo "  Running FastQC..."
    ${FASTQC} \
        --outdir ${QC_DIR} \
        --threads ${SLURM_CPUS_PER_TASK} \
        --quiet \
        ${FASTQ}
    
    # Check if successful
    if [[ -f "${QC_DIR}/${BASENAME}_fastqc.html" ]]; then
        echo "  ✓ SUCCESS: Report generated"
        ((SUCCESS_COUNT++))
    else
        echo "  ✗ FAILED: No report generated"
        ((FAIL_COUNT++))
    fi
done

# ============================================================
# SUMMARY
# ============================================================
echo ""
echo "=========================================="
echo "FastQC SUMMARY"
echo "=========================================="
echo "Total files processed: ${FASTQ_COUNT}"
echo "  Successful: ${SUCCESS_COUNT}"
echo "  Failed: ${FAIL_COUNT}"
echo ""

if [[ ${SUCCESS_COUNT} -gt 0 ]]; then
    echo "FastQC reports available in:"
    echo "  ${QC_DIR}"
    echo ""
    echo "Generated reports:"
    ls -1 ${QC_DIR}/*.html | awk '{print "  "$0}'
    echo ""
    echo "=== NEXT STEPS ==="
    echo "1. Review FastQC HTML reports in: ${QC_DIR}/"
    echo "2. Look for:"
    echo "   - Per base sequence quality (should be >20, ideally >30)"
    echo "   - Adapter content (if high, trimming is needed)"
    echo "   - Sequence duplication levels (expected to be high for Hi-C)"
    echo "   - Overrepresented sequences"
    echo ""
    echo "3. After reviewing reports, run: sbatch step2_trim.sh"
else
    echo "ERROR: No FastQC reports were generated!"
    exit 1
fi

echo ""
echo "Completed at: $(date)"
echo "=========================================="
