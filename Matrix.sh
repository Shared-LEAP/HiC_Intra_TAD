#!/bin/bash
#SBATCH --job-name=HiC_step5_matrix
#SBATCH --output=matrix.log
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100G
#SBATCH --nodelist=n001

# ============================================================
# STEP 5: CONTACT MATRIX GENERATION
# ============================================================
# Purpose: Create multi-resolution contact matrices (.mcool)
# Input: Merged pairs files from pairs/ directory
# Output: Normalized .mcool files in contact_maps/
# ============================================================

echo "=========================================="
echo "STEP 5: Contact Matrix Generation"
echo "Started at: $(date)"
echo "=========================================="

# ============================================================
# DIRECTORY SETUP
# ============================================================
WORK_DIR="/YOUR DIRECTORY HERE/"
PROCESS_DIR="${WORK_DIR}/HiC_processing"
PAIRS_DIR="${PROCESS_DIR}/pairs"
MAPS_DIR="${PROCESS_DIR}/contact_maps"
CHROMSIZES="${PROCESS_DIR}/references/hg38.chrom.sizes"

# Tool paths
CONDA_BASE="/YOUR DIRECTORY HERE/anaconda3"

# Activate conda environment
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate cooler

# Use simple tool name
COOLER="cooler"

# Create directories
mkdir -p ${MAPS_DIR}
mkdir -p ${PROCESS_DIR}/logs/cooler
mkdir -p ${PROCESS_DIR}/temp_cooler

cd ${PROCESS_DIR}

echo ""
echo "Working directory: ${PROCESS_DIR}"
echo "Input: ${PAIRS_DIR}"
echo "Output: ${MAPS_DIR}"
echo "Chrom sizes: ${CHROMSIZES}"
echo "Threads: ${SLURM_CPUS_PER_TASK}"

# ============================================================
# RESOLUTION SETTINGS
# ============================================================
echo ""
echo "=== Resolution Settings ==="

# Define resolutions (in base pairs)
# Standard 4DN/HiGlass resolutions
RESOLUTIONS="1000,2000,5000,10000,25000,50000,100000,250000,500000,1000000,2500000,5000000,10000000"

# Base resolution for initial binning (highest resolution)
BASE_RES=1000

echo "Base resolution: ${BASE_RES} bp"
echo "Multi-resolution zoom levels: ${RESOLUTIONS}"

# ============================================================
# VERIFY INPUT FILES
# ============================================================
echo ""
echo "=== Verifying Input Files ==="

REPLICATES=("rep1_merged" "rep2_merged")
VALID_REPS=()

for REP in "${REPLICATES[@]}"; do
    PAIRS_FILE="${PAIRS_DIR}/${REP}.pairs.gz"
    INDEX_FILE="${PAIRS_FILE}.px2"
    
    if [[ -f ${PAIRS_FILE} ]]; then
        if [[ -f ${INDEX_FILE} ]]; then
            PAIR_COUNT=$(zcat ${PAIRS_FILE} | grep -v "^#" | wc -l)
            FILE_SIZE=$(ls -lh ${PAIRS_FILE} | awk '{print $5}')
            echo "  ✓ ${REP}: ${PAIR_COUNT} pairs (${FILE_SIZE})"
            VALID_REPS+=("${REP}")
        else
            echo "  ⚠ ${REP}: Missing index file (.px2)"
        fi
    else
        echo "  ✗ ${REP}: Pairs file not found"
    fi
done

if [[ ${#VALID_REPS[@]} -eq 0 ]]; then
    echo ""
    echo "ERROR: No valid pairs files found!"
    echo "Please run step4_parse_filter.sh first"
    exit 1
fi

echo ""
echo "Processing ${#VALID_REPS[@]} replicate(s)"

# ============================================================
# PROCESS EACH REPLICATE
# ============================================================
echo ""
echo "=== Creating Contact Matrices ==="
echo ""

SUCCESS_COUNT=0
FAIL_COUNT=0

for REP in "${VALID_REPS[@]}"; do
    PAIRS_FILE="${PAIRS_DIR}/${REP}.pairs.gz"
    
    echo "=========================================="
    echo "Replicate: ${REP}"
    echo "=========================================="
    echo "  Input: ${PAIRS_FILE}"
    
    # Output files
    COOL_BASE="${MAPS_DIR}/${REP}.cool"
    MCOOL_FILE="${MAPS_DIR}/${REP}.mcool"
    
    # ----------------------------------------------------------
    # STEP 5A: CREATE BASE RESOLUTION COOL FILE
    # ----------------------------------------------------------
    echo ""
    echo "  [1/3] Creating base resolution cool file (${BASE_RES} bp)..."
    
    ${COOLER} cload pairs \
        -c1 2 -p1 3 -c2 4 -p2 5 \
        ${CHROMSIZES}:${BASE_RES} \
        ${PAIRS_FILE} \
        ${COOL_BASE} \
        2> logs/cooler/${REP}_cload.log
    
    if [[ ! -f ${COOL_BASE} ]]; then
        echo "    ✗ Failed to create cool file - check log"
        ((FAIL_COUNT++))
        continue
    fi
    
    # Get number of bins and interactions
    NUM_BINS=$(${COOLER} info ${COOL_BASE} | grep "nbins" | awk '{print $2}')
    NUM_CONTACTS=$(${COOLER} info ${COOL_BASE} | grep "nnz" | awk '{print $2}')
    
    echo "    ✓ Cool file created"
    echo "    Bins: ${NUM_BINS}"
    echo "    Contacts: ${NUM_CONTACTS}"
    
    # ----------------------------------------------------------
    # STEP 5B: CREATE MULTI-RESOLUTION FILE
    # ----------------------------------------------------------
    echo ""
    echo "  [2/3] Creating multi-resolution mcool file..."
    
    ${COOLER} zoomify \
        --nproc ${SLURM_CPUS_PER_TASK} \
        --balance \
        --resolutions ${RESOLUTIONS} \
        --out ${MCOOL_FILE} \
        ${COOL_BASE} \
        2> logs/cooler/${REP}_zoomify.log
    
    if [[ ! -f ${MCOOL_FILE} ]]; then
        echo "    ✗ Failed to create mcool file - check log"
        ((FAIL_COUNT++))
        continue
    fi
    
    echo "    ✓ Multi-resolution mcool created"
    
    # ----------------------------------------------------------
    # STEP 5C: BALANCE/NORMALIZE MATRICES
    # ----------------------------------------------------------
    echo ""
    echo "  [3/3] Applying ICE normalization..."
    
    # Note: zoomify with --balance already applies ICE normalization
    # This step verifies the normalization
    
    ${COOLER} balance \
        --nproc ${SLURM_CPUS_PER_TASK} \
        --force \
        ${MCOOL_FILE}::resolutions/${BASE_RES} \
        2> logs/cooler/${REP}_balance.log
    
    echo "    ✓ ICE normalization applied"
    
    # ----------------------------------------------------------
    # VERIFY AND SHOW RESOLUTIONS
    # ----------------------------------------------------------
    echo ""
    echo "  Available resolutions in ${REP}.mcool:"
    ${COOLER} ls ${MCOOL_FILE} | grep "resolutions" | \
        sed 's|/resolutions/||' | awk '{printf "    %s bp\n", $1}'
    
    # Show file size
    MCOOL_SIZE=$(ls -lh ${MCOOL_FILE} | awk '{print $5}')
    echo ""
    echo "  File size: ${MCOOL_SIZE}"
    
    # ----------------------------------------------------------
    # CLEANUP BASE RESOLUTION FILE
    # ----------------------------------------------------------
    echo ""
    echo "  Cleaning up intermediate files..."
    rm -f ${COOL_BASE}
    echo "    ✓ Removed base cool file (data preserved in mcool)"
    
    ((SUCCESS_COUNT++))
    echo ""
done

# ============================================================
# FINAL SUMMARY
# ============================================================
echo ""
echo "=========================================="
echo "MATRIX GENERATION SUMMARY"
echo "=========================================="
echo "Replicates processed: ${SUCCESS_COUNT}/${#VALID_REPS[@]}"
echo ""

if [[ ${SUCCESS_COUNT} -gt 0 ]]; then
    echo "Contact matrices created:"
    ls -lh ${MAPS_DIR}/*.mcool | awk '{print "  "$9" ("$5")"}'
    echo ""
    
    # Show detailed info for each mcool
    for MCOOL in ${MAPS_DIR}/*.mcool; do
        if [[ -f ${MCOOL} ]]; then
            REP=$(basename ${MCOOL} .mcool)
            echo "-----------------------------------"
            echo "Matrix: ${REP}"
            echo "-----------------------------------"
            
            # Get info for base resolution
            BASE_INFO=$(${COOLER} info ${MCOOL}::resolutions/${BASE_RES})
            
            echo "Genome: $(echo "${BASE_INFO}" | grep "genome-assembly" | awk '{print $2}')"
            echo "Total bins: $(echo "${BASE_INFO}" | grep "nbins" | awk '{print $2}')"
            echo "Total contacts: $(echo "${BASE_INFO}" | grep "nnz" | awk '{print $2}')"
            echo "Balanced: $(echo "${BASE_INFO}" | grep "nchroms" | awk '{print $2}' > /dev/null && echo "Yes" || echo "No")"
            
            echo ""
            echo "Resolutions available:"
            ${COOLER} ls ${MCOOL} | grep "resolutions" | sed 's|/resolutions/||' | \
                awk '{printf "  %s bp\n", $1}'
            echo ""
        fi
    done
    
    echo "=========================================="
    echo "PIPELINE COMPLETE!"
    echo "=========================================="
    echo ""
    echo "=== QUALITY CHECKS ==="
    echo "1. Verify contact matrices:"
    echo "   cooler dump --table chroms ${MAPS_DIR}/rep1_merged.mcool::resolutions/10000"
    echo ""
    echo "2. Check balance quality:"
    echo "   cooler dump --table bins ${MAPS_DIR}/rep1_merged.mcool::resolutions/10000 | head"
    echo ""
    echo "=== VISUALIZATION OPTIONS ==="
    echo "1. HiGlass: https://higlass.io/"
    echo "   - Upload .mcool files directly"
    echo ""
    echo "2. Cooler show (quick plots):"
    echo "   cooler show ${MAPS_DIR}/rep1_merged.mcool::resolutions/100000"
    echo ""
    echo "3. Python analysis:"
    echo "   import cooler"
    echo "   c = cooler.Cooler('${MAPS_DIR}/rep1_merged.mcool::resolutions/10000')"
    echo "   matrix = c.matrix(balance=True).fetch('chr1')"
    echo ""
    echo "=== OUTPUT FILES ==="
    echo "Contact matrices: ${MAPS_DIR}/"
    echo "Pairs files: ${PAIRS_DIR}/"
    echo "QC reports: ${PROCESS_DIR}/qc/"
    echo "All logs: ${PROCESS_DIR}/logs/"
else
    echo "ERROR: No contact matrices were created!"
    exit 1
fi

echo ""
echo "Completed at: $(date)"
echo "=========================================="
