#!/bin/bash
#SBATCH --job-name=HiC_step4_parse
#SBATCH --output=parse.log
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100G
#SBATCH --nodelist=n001

# ============================================================
# STEP 4: PAIRTOOLS PARSING AND FILTERING
# ============================================================
# Purpose: Convert BAM to pairs, filter, deduplicate
# Input: BAM files from aligned/ directory
# Output: Filtered, deduplicated pairs files
# ============================================================

echo "=========================================="
echo "STEP 4: Pairtools Parse & Filter"
echo "Started at: $(date)"
echo "=========================================="

# ============================================================
# DIRECTORY SETUP
# ============================================================
WORK_DIR="/home/jnavarrete/3D_Chromatin/4DN_pipeline"
PROCESS_DIR="${WORK_DIR}/HiC_processing"
ALIGNED_DIR="${PROCESS_DIR}/aligned"
PAIRS_DIR="${PROCESS_DIR}/pairs"
CHROMSIZES="/home/jnavarrete/3D_Chromatin/4DN_pipeline/HiC_processing/references/hg38.chrom.sizes"


# Tool paths
CONDA_BASE="/home/jnavarrete/anaconda3"

# CRITICAL: Activate conda environment so pairtools can find samtools/pbgzip
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate pairtools

# Now tools are in PATH
PAIRTOOLS="pairtools"
PAIRIX="pairix"


# Create directories
mkdir -p ${PAIRS_DIR}
mkdir -p ${PROCESS_DIR}/logs/pairtools
mkdir -p ${PROCESS_DIR}/temp

cd ${PROCESS_DIR}

echo ""
echo "Working directory: ${PROCESS_DIR}"
echo "Input: ${ALIGNED_DIR}"
echo "Output: ${PAIRS_DIR}"
echo "Chrom sizes: ${CHROMSIZES}"
echo "Threads: ${SLURM_CPUS_PER_TASK}"

# ============================================================
# VERIFY CHROMOSOME SIZES FILE
# ============================================================
echo ""
echo "=== Verifying Chromosome Sizes ==="

if [[ ! -f ${CHROMSIZES} ]]; then
    echo "ERROR: Chromosome sizes file not found: ${CHROMSIZES}"
    echo ""
    echo "Create it with:"
    echo "  samtools faidx ${PROCESS_DIR}/references/hg38.fa"
    echo "  cut -f1,2 ${PROCESS_DIR}/references/hg38.fa.fai > ${CHROMSIZES}"
    exit 1
fi

echo "  ✓ Chromosome sizes file verified"
echo "  Chromosomes:"
head -5 ${CHROMSIZES} | awk '{print "    "$1": "$2" bp"}'
echo "    ..."

# ============================================================
# DEFINE SAMPLE PAIRS
# ============================================================
echo ""
echo "=== Setting Up Sample Pairs ==="

# Individual samples to process
SAMPLES=(
    "rep1_pair1"
    "rep1_pair2"
    "rep1_pair3"
    "rep2_pair1"
    "rep2_pair2"
)

echo "Total samples to process: ${#SAMPLES[@]}"

# ============================================================
# PROCESS EACH SAMPLE
# ============================================================
echo ""
echo "=== Processing Individual Samples ==="
echo ""

SUCCESS_COUNT=0
FAIL_COUNT=0

for SAMPLE in "${SAMPLES[@]}"; do
    BAM_INPUT="${ALIGNED_DIR}/${SAMPLE}.bam"
    
    echo "=========================================="
    echo "Sample: ${SAMPLE}"
    echo "=========================================="
    
    # Verify BAM exists
    if [[ ! -f ${BAM_INPUT} ]]; then
        echo "  ✗ ERROR: BAM file not found: ${BAM_INPUT}"
        ((FAIL_COUNT++))
        continue
    fi
    
    BAM_SIZE=$(ls -lh ${BAM_INPUT} | awk '{print $5}')
    echo "  Input BAM: ${BAM_INPUT} (${BAM_SIZE})"
    
    # Create temp directory for this sample
    TEMP_DIR="temp/${SAMPLE}"
    mkdir -p ${TEMP_DIR}
    
    # ----------------------------------------------------------
    # STEP 4A: PARSE BAM TO PAIRSAM
    # ----------------------------------------------------------
    echo ""
    echo "  [1/4] Parsing BAM to pairsam format..."
    
    PAIRSAM="${TEMP_DIR}/${SAMPLE}.pairsam.gz"
    
    ${PAIRTOOLS} parse \
        --chroms-path ${CHROMSIZES} \
        --min-mapq 30 \
        --walks-policy 5unique \
        --max-inter-align-gap 30 \
        --output ${PAIRSAM} \
        ${BAM_INPUT} \
        2> logs/pairtools/${SAMPLE}_parse.log
    
    if [[ ! -s ${PAIRSAM} ]]; then
        echo "    ✗ Parse failed - check log"
        ((FAIL_COUNT++))
        continue
    fi
    
    # Count contacts
    PARSE_COUNT=$(zcat ${PAIRSAM} | grep -v "^#" | wc -l)
    echo "    ✓ Parsed: ${PARSE_COUNT} contacts"
    
    # ----------------------------------------------------------
    # STEP 4B: SORT PAIRSAM
    # ----------------------------------------------------------
    echo ""
    echo "  [2/4] Sorting pairsam..."
    
    SORTED="${TEMP_DIR}/${SAMPLE}.sorted.pairsam.gz"
    
    ${PAIRTOOLS} sort \
        --nproc ${SLURM_CPUS_PER_TASK} \
        --tmpdir ${TEMP_DIR} \
        --output ${SORTED} \
        ${PAIRSAM} \
        2> logs/pairtools/${SAMPLE}_sort.log
    
    if [[ ! -s ${SORTED} ]]; then
        echo "    ✗ Sort failed - check log"
        ((FAIL_COUNT++))
        continue
    fi
    
    echo "    ✓ Sorted successfully"
    
    # ----------------------------------------------------------
    # STEP 4C: DEDUPLICATE
    # ----------------------------------------------------------
    echo ""
    echo "  [3/4] Deduplicating..."
    
    DEDUP="${TEMP_DIR}/${SAMPLE}.dedup.pairsam.gz"
    STATS="${PAIRS_DIR}/${SAMPLE}_dedup_stats.txt"
    
    ${PAIRTOOLS} dedup \
        --max-mismatch 1 \
        --mark-dups \
        --output ${DEDUP} \
        --output-stats ${STATS} \
        ${SORTED} \
        2> logs/pairtools/${SAMPLE}_dedup.log
    
    if [[ ! -s ${DEDUP} ]]; then
        echo "    ✗ Dedup failed - check log"
        ((FAIL_COUNT++))
        continue
    fi
    
    # Show dedup statistics
    if [[ -f ${STATS} ]]; then
        DUP_RATE=$(grep "duplicate" ${STATS} | head -1 | awk '{print $NF}')
        echo "    ✓ Deduplicated: ${DUP_RATE} duplicate rate"
    fi
    
    # ----------------------------------------------------------
    # STEP 4D: SELECT VALID PAIRS AND SPLIT
    # ----------------------------------------------------------
    echo ""
    echo "  [4/4] Selecting valid Hi-C pairs..."
    
    VALID="${TEMP_DIR}/${SAMPLE}.valid.pairsam.gz"
    FINAL_PAIRS="${PAIRS_DIR}/${SAMPLE}.pairs.gz"
    
    # Select only valid unique pairs (UU, UR, RU)
    ${PAIRTOOLS} select \
        '(pair_type == "UU") or (pair_type == "UR") or (pair_type == "RU")' \
        --output ${VALID} \
        ${DEDUP} \
        2> logs/pairtools/${SAMPLE}_select.log
    
    if [[ ! -s ${VALID} ]]; then
        echo "    ✗ Select failed - check log"
        ((FAIL_COUNT++))
        continue
    fi
    
    # Split pairsam to pairs (remove SAM)
    ${PAIRTOOLS} split \
        --output-pairs ${FINAL_PAIRS} \
        ${VALID} \
        2> logs/pairtools/${SAMPLE}_split.log
    
    if [[ ! -s ${FINAL_PAIRS} ]]; then
        echo "    ✗ Split failed - check log"
        ((FAIL_COUNT++))
        continue
    fi
    
    # ----------------------------------------------------------
    # STEP 4E: INDEX PAIRS FILE
    # ----------------------------------------------------------
    echo ""
    echo "  [5/4] Indexing pairs file..."
    
    ${PAIRIX} -f -p pairs ${FINAL_PAIRS} \
        2> logs/pairtools/${SAMPLE}_pairix.log
    
    if [[ -f "${FINAL_PAIRS}.px2" ]]; then
        echo "    ✓ Indexed successfully"
    else
        echo "    ⚠ Warning: Indexing may have failed"
    fi
    
    # ----------------------------------------------------------
    # SUMMARY FOR THIS SAMPLE
    # ----------------------------------------------------------
    VALID_COUNT=$(zcat ${FINAL_PAIRS} | grep -v "^#" | wc -l)
    PAIRS_SIZE=$(ls -lh ${FINAL_PAIRS} | awk '{print $5}')
    
    echo ""
    echo "  =========================================="
    echo "  SAMPLE SUMMARY: ${SAMPLE}"
    echo "  =========================================="
    echo "  Final valid pairs: ${VALID_COUNT}"
    echo "  Output file: ${PAIRS_SIZE}"
    echo "  Pairs file: ${FINAL_PAIRS}"
    echo "  Index file: ${FINAL_PAIRS}.px2"
    
    # Show pair type distribution
    echo ""
    echo "  Pair type distribution:"
    zcat ${FINAL_PAIRS} | grep -v "^#" | cut -f8 | sort | uniq -c | \
        awk '{printf "    %s: %s\n", $2, $1}'
    
    # Cleanup temp files
    rm -rf ${TEMP_DIR}
    
    ((SUCCESS_COUNT++))
    echo ""
done

# ============================================================
# MERGE REPLICATES
# ============================================================
echo ""
echo "=========================================="
echo "MERGING BIOLOGICAL REPLICATES"
echo "=========================================="

# Merge Replicate 1
REP1_FILES=()
for SAMPLE in rep1_pair1 rep1_pair2 rep1_pair3; do
    PAIRS_FILE="${PAIRS_DIR}/${SAMPLE}.pairs.gz"
    if [[ -f ${PAIRS_FILE} ]]; then
        REP1_FILES+=("${PAIRS_FILE}")
    fi
done

if [[ ${#REP1_FILES[@]} -gt 0 ]]; then
    echo ""
    echo "Merging Replicate 1 (${#REP1_FILES[@]} files)..."
    
    ${PAIRTOOLS} merge \
        --nproc ${SLURM_CPUS_PER_TASK} \
        --output ${PAIRS_DIR}/rep1_merged.pairs.gz \
        ${REP1_FILES[@]} \
        2> logs/pairtools/merge_rep1.log
    
    if [[ -f "${PAIRS_DIR}/rep1_merged.pairs.gz" ]]; then
        ${PAIRIX} -f -p pairs ${PAIRS_DIR}/rep1_merged.pairs.gz
        
        REP1_COUNT=$(zcat ${PAIRS_DIR}/rep1_merged.pairs.gz | grep -v "^#" | wc -l)
        REP1_SIZE=$(ls -lh ${PAIRS_DIR}/rep1_merged.pairs.gz | awk '{print $5}')
        echo "  ✓ Replicate 1 merged: ${REP1_COUNT} pairs (${REP1_SIZE})"
    else
        echo "  ✗ Replicate 1 merge failed"
    fi
fi

# Merge Replicate 2
REP2_FILES=()
for SAMPLE in rep2_pair1 rep2_pair2; do
    PAIRS_FILE="${PAIRS_DIR}/${SAMPLE}.pairs.gz"
    if [[ -f ${PAIRS_FILE} ]]; then
        REP2_FILES+=("${PAIRS_FILE}")
    fi
done

if [[ ${#REP2_FILES[@]} -gt 0 ]]; then
    echo ""
    echo "Merging Replicate 2 (${#REP2_FILES[@]} files)..."
    
    ${PAIRTOOLS} merge \
        --nproc ${SLURM_CPUS_PER_TASK} \
        --output ${PAIRS_DIR}/rep2_merged.pairs.gz \
        ${REP2_FILES[@]} \
        2> logs/pairtools/merge_rep2.log
    
    if [[ -f "${PAIRS_DIR}/rep2_merged.pairs.gz" ]]; then
        ${PAIRIX} -f -p pairs ${PAIRS_DIR}/rep2_merged.pairs.gz
        
        REP2_COUNT=$(zcat ${PAIRS_DIR}/rep2_merged.pairs.gz | grep -v "^#" | wc -l)
        REP2_SIZE=$(ls -lh ${PAIRS_DIR}/rep2_merged.pairs.gz | awk '{print $5}')
        echo "  ✓ Replicate 2 merged: ${REP2_COUNT} pairs (${REP2_SIZE})"
    else
        echo "  ✗ Replicate 2 merge failed"
    fi
fi

# ============================================================
# FINAL SUMMARY
# ============================================================
echo ""
echo "=========================================="
echo "FINAL SUMMARY"
echo "=========================================="
echo "Individual samples processed: ${SUCCESS_COUNT}/${#SAMPLES[@]}"
echo ""

if [[ -f "${PAIRS_DIR}/rep1_merged.pairs.gz" ]]; then
    echo "Replicate 1 merged pairs:"
    echo "  File: ${PAIRS_DIR}/rep1_merged.pairs.gz"
    echo "  Size: $(ls -lh ${PAIRS_DIR}/rep1_merged.pairs.gz | awk '{print $5}')"
    echo "  Pairs: $(zcat ${PAIRS_DIR}/rep1_merged.pairs.gz | grep -v "^#" | wc -l)"
fi

if [[ -f "${PAIRS_DIR}/rep2_merged.pairs.gz" ]]; then
    echo ""
    echo "Replicate 2 merged pairs:"
    echo "  File: ${PAIRS_DIR}/rep2_merged.pairs.gz"
    echo "  Size: $(ls -lh ${PAIRS_DIR}/rep2_merged.pairs.gz | awk '{print $5}')"
    echo "  Pairs: $(zcat ${PAIRS_DIR}/rep2_merged.pairs.gz | grep -v "^#" | wc -l)"
fi

echo ""
echo "Deduplication statistics available in:"
ls -1 ${PAIRS_DIR}/*_dedup_stats.txt 2>/dev/null | awk '{print "  "$0}'

echo ""
echo "=== NEXT STEPS ==="
echo "Run: sbatch step5_matrices.sh"

echo ""
echo "Completed at: $(date)"
echo "=========================================="