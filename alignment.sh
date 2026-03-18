#!/bin/bash
#SBATCH --job-name=HiC_step2_align
#SBATCH --output=HiC_processing/logs/step2_alignment.log
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100G
#SBATCH --nodelist=n001

# ============================================================
# STEP 3: BWA ALIGNMENT
# ============================================================
# Purpose: Align trimmed Hi-C reads to reference genome
# Input: Trimmed FASTQ files from trimmed/ directory
# Output: BAM files in aligned/ directory
# ============================================================

echo "=========================================="
echo "STEP 3: BWA Alignment"
echo "Started at: $(date)"
echo "=========================================="

# ============================================================
# DIRECTORY SETUP
# ============================================================
WORK_DIR="/YOUR DIRECTORY HERE/"
PROCESS_DIR="${WORK_DIR}/HiC_processing"
FASTQ_DIR="/YOUR DIRECTORY HERE/Fastq.gz"
ALIGNED_DIR="${PROCESS_DIR}/aligned"
REFERENCE="/YOUR DIRECTORY HERE/references/hg38.fa"

# Tool paths
CONDA_BASE="/YOUR DIRECTORY HERE/anaconda3"
BWA="${CONDA_BASE}/envs/bwa/bin/bwa"
SAMTOOLS="${CONDA_BASE}/envs/samtools/bin/samtools"

# Create directories
mkdir -p ${ALIGNED_DIR}
mkdir -p ${PROCESS_DIR}/logs/alignment

cd ${PROCESS_DIR}

echo ""
echo "Working directory: ${PROCESS_DIR}"
echo "Input: ${FASTQ_DIR}"
echo "Output: ${ALIGNED_DIR}"
echo "Reference: ${REFERENCE}"
echo "Threads: ${SLURM_CPUS_PER_TASK}"

# ============================================================
# VERIFY REFERENCE GENOME
# ============================================================
echo ""
echo "=== Verifying Reference Genome ==="

if [[ ! -f ${REFERENCE} ]]; then
    echo "ERROR: Reference genome not found: ${REFERENCE}"
    echo "Please ensure hg38.fa is in ${PROCESS_DIR}/references/"
    exit 1
fi

# Check for BWA index files
INDEX_FILES=("${REFERENCE}.amb" "${REFERENCE}.ann" "${REFERENCE}.bwt" "${REFERENCE}.pac" "${REFERENCE}.sa")
INDEX_MISSING=0

for IDX_FILE in "${INDEX_FILES[@]}"; do
    if [[ ! -f ${IDX_FILE} ]]; then
        echo "  ✗ Missing index file: ${IDX_FILE}"
        INDEX_MISSING=1
    fi
done

if [[ ${INDEX_MISSING} -eq 1 ]]; then
    echo ""
    echo "ERROR: BWA index files are missing!"
    echo "Please create index with: bwa index ${REFERENCE}"
    exit 1
fi

echo "  ✓ Reference genome and index files verified"

# ============================================================
# DEFINE SAMPLE PAIRS
# ============================================================
echo ""
echo "=== Setting Up Sample Pairs ==="

# Replicate 1 pairs
REP1_PAIRS=(
    "rep1_pair1:4DNFI49T8A9Y:4DNFIJYA9SGN"
    "rep1_pair2:4DNFI4TDJG6B:4DNFISHZEA1D"
    "rep1_pair3:4DNFIC7PJLUM:4DNFIV7BFANE"
)

# Replicate 2 pairs
REP2_PAIRS=(
    "rep2_pair1:4DNFIC1P6NXF:4DNFIX2DEJKT"
    "rep2_pair2:4DNFINHCOYVE:4DNFIC8NUJ51"
)

ALL_PAIRS=("${REP1_PAIRS[@]}" "${REP2_PAIRS[@]}")

echo "Total pairs to align: ${#ALL_PAIRS[@]}"

# ============================================================
# ALIGN EACH PAIR
# ============================================================
echo ""
echo "=== Starting Alignment ==="
echo "Using BWA-MEM with -SP5M flags (optimized for Hi-C)"
echo ""

SUCCESS_COUNT=0
FAIL_COUNT=0

for PAIR in "${ALL_PAIRS[@]}"; do
    IFS=':' read -r SAMPLE R1_BASE R2_BASE <<< "${PAIR}"
    
    # Input files (raw FASTQ - no trimming per 4DN pipeline)
    R1_INPUT="${FASTQ_DIR}/${R1_BASE}.fastq.gz"
    R2_INPUT="${FASTQ_DIR}/${R2_BASE}.fastq.gz"
    
    # Output files
    BAM_OUTPUT="${ALIGNED_DIR}/${SAMPLE}.bam"
    LOG_FILE="logs/alignment/${SAMPLE}_bwa.log"
    
    echo "-----------------------------------"
    echo "Sample: ${SAMPLE}"
    echo "  R1: ${R1_BASE}"
    echo "  R2: ${R2_BASE}"
    
    # Verify input files exist
    if [[ ! -f ${R1_INPUT} ]]; then
        echo "  ✗ ERROR: R1 not found: ${R1_INPUT}"
        ((FAIL_COUNT++))
        continue
    fi
    
    if [[ ! -f ${R2_INPUT} ]]; then
        echo "  ✗ ERROR: R2 not found: ${R2_INPUT}"
        ((FAIL_COUNT++))
        continue
    fi
    
    # Show input file sizes
    R1_SIZE=$(ls -lh ${R1_INPUT} | awk '{print $5}')
    R2_SIZE=$(ls -lh ${R2_INPUT} | awk '{print $5}')
    echo "  Input sizes: R1=${R1_SIZE}, R2=${R2_SIZE}"
    
    # Run BWA-MEM
    echo "  Running BWA-MEM alignment..."
    echo "  Started at: $(date)"
    
    # BWA flags for Hi-C:
    # -S: Skip pairing (Hi-C reads from ligation junctions)
    # -P: Mark secondary alignments
    # -5: For split alignment (chimeric reads at ligation junction)
    # -M: Mark shorter split hits as secondary
    ${BWA} mem \
        -SP5M \
        -t ${SLURM_CPUS_PER_TASK} \
        ${REFERENCE} \
        ${R1_INPUT} \
        ${R2_INPUT} \
        2> ${LOG_FILE} \
        | ${SAMTOOLS} view -bS - > ${BAM_OUTPUT}
    
    # Check if successful
    if [[ -s ${BAM_OUTPUT} ]]; then
        BAM_SIZE=$(ls -lh ${BAM_OUTPUT} | awk '{print $5}')
        echo "  ✓ Alignment successful"
        echo "  Output BAM size: ${BAM_SIZE}"
        
        # Get alignment statistics
        echo "  Generating alignment statistics..."
        ${SAMTOOLS} flagstat ${BAM_OUTPUT} > logs/alignment/${SAMPLE}_flagstat.txt
        
        # Show key stats
        echo "  Alignment stats:"
        grep "mapped (" logs/alignment/${SAMPLE}_flagstat.txt | head -1 | awk '{print "    "$0}'
        grep "properly paired" logs/alignment/${SAMPLE}_flagstat.txt | awk '{print "    "$0}'
        
        ((SUCCESS_COUNT++))
    else
        echo "  ✗ Alignment failed - check log: ${LOG_FILE}"
        ((FAIL_COUNT++))
    fi
    
    echo "  Finished at: $(date)"
    echo ""
done

# ============================================================
# SUMMARY
# ============================================================
echo ""
echo "=========================================="
echo "ALIGNMENT SUMMARY"
echo "=========================================="
echo "Total pairs processed: ${#ALL_PAIRS[@]}"
echo "  Successful: ${SUCCESS_COUNT}"
echo "  Failed: ${FAIL_COUNT}"
echo ""

if [[ ${SUCCESS_COUNT} -gt 0 ]]; then
    echo "Aligned BAM files available in:"
    echo "  ${ALIGNED_DIR}/"
    echo ""
    echo "Files created:"
    ls -lh ${ALIGNED_DIR}/*.bam | awk '{print "  "$9" ("$5")"}'
    echo ""
    echo "Alignment statistics available in:"
    echo "  logs/alignment/*_flagstat.txt"
    echo ""
    echo "=== QUALITY CHECK ==="
    echo "Review mapping rates in flagstat files:"
    for FLAGSTAT in logs/alignment/*_flagstat.txt; do
        if [[ -f ${FLAGSTAT} ]]; then
            SAMPLE=$(basename ${FLAGSTAT} _flagstat.txt)
            MAPPED=$(grep "mapped (" ${FLAGSTAT} | head -1 | awk '{print $5}')
            echo "  ${SAMPLE}: ${MAPPED} mapped"
        fi
    done
    echo ""
    echo "Expected: >80% mapping rate for good Hi-C data"
    echo ""
    echo "=== NEXT STEPS ==="
    echo "Run: sbatch step4_parse_filter.sh"
else
    echo "ERROR: No files were successfully aligned!"
    exit 1
fi

echo ""
echo "Completed at: $(date)"
echo "=========================================="
