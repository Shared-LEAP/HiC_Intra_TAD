# Interactive 4DN Hi-C Pipeline
## Cardiomyocyte Differentiation Hi-C Analysis

**Author:** Josue Navarrete  
**Date:** January 2026  
**Purpose:** Step-by-step processing of Hi-C paired-end data to generate contact matrices

---

## Overview

This pipeline processes Hi-C data through 5 separate, interactive steps:

1. **Quality Control** - FastQC on raw reads
2. **Alignment** - BWA-MEM mapping to reference genome
3. **Parsing & Filtering** - Convert to pairs, deduplicate, filter
4. **Matrix Generation** - Create multi-resolution contact matrices

---

## Dataset Information

**Source:** 4DN Data Portal  
**Experiment:** Hi-C on WTC-11 cardiac differentiation timecourse  
**Sample:** Day 23 ventricular cardiac myocytes  
**Sequencing:** In situ Hi-C, paired-end Illumina

**Files:**
- Replicate 1: 3 pairs (6 FASTQ files)
- Replicate 2: 2 pairs (4 FASTQ files)
- **Total:** 10 FASTQ files

---

## Directory Structure

```
/(Your Home Directory)
├── Fastq/                          # Raw FASTQ files (downloaded)
└── HiC_processing/
    ├── qc/
    │   ├── fastqc_raw/            # QC reports on raw data
    │   └── fastqc_trimmed/        # QC reports on trimmed data
    ├── aligned/                    # BAM files from BWA
    ├── pairs/                      # Pairs files (.pairs.gz)
    ├── contact_maps/               # Final .mcool matrices
    ├── references/
    │   ├── hg38.fa                # Reference genome
    │   ├── hg38.fa.* (BWA index)  # BWA index files
    │   └── hg38.chrom.sizes       # Chromosome sizes
    ├── logs/                       # All log files
    └── temp*/                      # Temporary files (auto-cleaned)
```

---

## Prerequisites

### 1. Reference Genome Setup

**Download hg38 reference:**
```bash
cd /(Your Home Directory)/references

# Option A: From UCSC
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz

# Option B: From GENCODE
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
mv GRCh38.primary_assembly.genome.fa hg38.fa
```

**Create BWA index:**
```bash
bwa index hg38.fa
# This creates: .amb, .ann, .bwt, .pac, .sa files
```

**Create chromosome sizes file:**
```bash
samtools faidx hg38.fa
cut -f1,2 hg38.fa.fai > hg38.chrom.sizes
```

### 2. Conda Environments

Create separate environments for each tool:

```bash
# FastQC
conda create -n fastqc -c bioconda fastqc
conda activate fastqc
which fastqc  # Verify path

# BWA
conda create -n bwa -c bioconda bwa samtools
conda activate bwa
which bwa

# Pairtools
conda create -n pairtools -c bioconda pairtools pairix
conda activate pairtools
which pairtools

# Cooler
conda create -n cooler -c bioconda cooler
conda activate cooler
which cooler
```

**Update tool paths in scripts if needed!**

### 3. Verify FASTQ Files

```bash
ls -lh /(Your Home Directory)/Fastq/
# Should show all 10 .fastq.gz files
```

---

## Running the Pipeline

### STEP 1: Quality Control (FastQC)

**Purpose:** Assess quality of raw sequencing data

**Run:**
```bash
cd /(Your Home Directory)/HiC_processing
sbatch step1_fastqc.sh
```

**Check progress:**
```bash
tail -f logs/step1_fastqc_*.log
squeue -u $USER
```

**Review results:**
```bash
# Open HTML reports in browser
ls qc/fastqc_raw/*.html

# Look for:
# - Per base sequence quality (should be >Q30)
# - Adapter content (if high, trimming needed)
# - Sequence duplication (high is expected for Hi-C)
```

**Expected runtime:** 1-2 hours

---

### STEP 2: Alignment with BWA-MEM

**Purpose:** Map reads to reference genome

**Run:**
```bash
sbatch step3_align.sh
```

**Check progress:**
```bash
tail -f logs/step3_align_*.log

# Monitor individual alignments
ls -lh aligned/*.bam
```

**Review mapping statistics:**
```bash
# Check mapping rates (should be >80%)
cat logs/alignment/*_flagstat.txt

# Example output:
# 50000000 + 0 in total (QC-passed reads + QC-failed reads)
# 45000000 + 0 mapped (90.00% : N/A)
```

**BWA flags used:**
- `-S`: Skip mate rescue (Hi-C specific)
- `-P`: Mark secondary alignments
- `-5`: Split alignment for chimeric reads
- `-M`: Mark shorter split hits as secondary

**Expected runtime:** 4-8 hours

---

### STEP 3: Parse, Filter, and Generate Pairs

**Purpose:** Convert BAM to pairs format, deduplicate, filter valid contacts

**Run:**
```bash
sbatch step4_parse_filter.sh
```

**Check progress:**
```bash
tail -f logs/step4_parse_*.log

# Monitor pairs generation
ls -lh pairs/*.pairs.gz
```

**Review statistics:**
```bash
# Check deduplication stats
cat pairs/*_dedup_stats.txt

# View first few pairs
zcat pairs/rep1_pair1.pairs.gz | head -20

# Check merged replicates
zcat pairs/rep1_merged.pairs.gz | grep -v "^#" | wc -l
```

**Processing steps:**
1. Parse BAM → pairsam (min MAPQ 30)
2. Sort by genomic coordinates
3. Mark duplicates
4. Select valid pairs (UU, UR, RU)
5. Split to pairs format
6. Index with pairix
7. Merge technical replicates

**Expected runtime:** 6-12 hours

---

### STEP 4: Generate Contact Matrices

**Purpose:** Create multi-resolution contact matrices with normalization

**Run:**
```bash
sbatch step5_matrices.sh
```

**Check progress:**
```bash
tail -f logs/step5_matrices_*.log

# Monitor matrix creation
ls -lh contact_maps/*.mcool
```

**Review matrices:**
```bash
# List resolutions
cooler ls contact_maps/rep1_merged.mcool

# Get matrix info
cooler info contact_maps/rep1_merged.mcool::resolutions/10000

# Quick visualization
cooler show contact_maps/rep1_merged.mcool::resolutions/100000

# Export to text (for custom analysis)
cooler dump contact_maps/rep1_merged.mcool::resolutions/10000 > matrix_10kb.txt
```

**Resolutions created:**
- 1kb, 2kb, 5kb, 10kb, 25kb, 50kb, 100kb
- 250kb, 500kb, 1Mb, 2.5Mb, 5Mb, 10Mb

**Normalization:** ICE (Iterative Correction and Eigenvector decomposition)

**Expected runtime:** 2-4 hours

---

## Output Files

### Final Products

**Contact Matrices (mcool format):**
```
contact_maps/rep1_merged.mcool  # Biological replicate 1
contact_maps/rep2_merged.mcool  # Biological replicate 2
```

**Pairs Files:**
```
pairs/rep1_merged.pairs.gz      # All valid pairs, replicate 1
pairs/rep2_merged.pairs.gz      # All valid pairs, replicate 2
pairs/rep1_merged.pairs.gz.px2  # Pairix index
```

### Intermediate Files

**Individual sample pairs:**
```
pairs/rep1_pair1.pairs.gz
pairs/rep1_pair2.pairs.gz
pairs/rep1_pair3.pairs.gz
pairs/rep2_pair1.pairs.gz
pairs/rep2_pair2.pairs.gz
```

**QC Reports:**
```
qc/fastqc_raw/*.html           # QC
```

---

## Quality Control Checkpoints

### After Each Step

**Step 1 (FastQC):**
- [ ] All 10 files have HTML reports
- [ ] Per-base quality mostly >Q30
- [ ] Note any adapter contamination

**Step 2 (Alignment):**
- [ ] Mapping rate >80% for all samples
- [ ] BAM files are reasonable size
- [ ] No errors in BWA logs

**Step 3 (Pairs):**
- [ ] Valid pairs extracted for all samples
- [ ] Duplicate rate 20-40% (typical for Hi-C)
- [ ] Merged replicates have millions of pairs
- [ ] Pairs indexed successfully (.px2 files exist)

**Step 4 (Matrices):**
- [ ] .mcool files created for both replicates
- [ ] All resolutions present
- [ ] ICE normalization applied
- [ ] File sizes reasonable (several GB each)

---

## Visualization

### Option 1: HiGlass (Recommended)

**Upload to HiGlass:**
```bash
# Use HiGlass web interface: https://higlass.io/
# Upload .mcool files directly
```

**Or install locally:**
```bash
pip install higlass-manage
higlass-manage start
higlass-manage ingest contact_maps/rep1_merged.mcool
```

### Option 2: Python Analysis

```python
import cooler
import matplotlib.pyplot as plt
import numpy as np

# Load matrix at 10kb resolution
c = cooler.Cooler('contact_maps/rep1_merged.mcool::resolutions/10000')

# Extract chromosome 1
matrix = c.matrix(balance=True).fetch('chr1')

# Plot
plt.figure(figsize=(10, 10))
plt.imshow(np.log10(matrix + 1), cmap='YlOrRd')
plt.title('Chr1 Hi-C Contact Matrix (10kb)')
plt.colorbar(label='log10(contacts + 1)')
plt.savefig('chr1_heatmap.png', dpi=300)
```

### Option 3: Juicebox

**Convert to .hic format (if needed):**
```bash
# This requires additional tools not in this pipeline
# See 4DN documentation for hic file generation
```

---

## Performance Tips

1. **Use fast local storage** for temp directories
2. **Increase threads** if you have CPUs available
3. **Use node-local storage** on cluster if available
4. **Run steps in parallel** for different samples if resources permit

---

## Next Steps After Pipeline

1. **Replicate comparison** - Compare rep1 vs rep2
2. **TAD calling** - Use tools like TADbit, Arrowhead
3. **Loop calling** - Use HiCCUPS, Mustache
4. **Compartment analysis** - A/B compartments
5. **Differential analysis** - Compare across time points

---

## Citations

**4DN Hi-C Pipeline:**
- Abdennur N, Mirny LA (2020). Cooler: scalable storage for Hi-C data and other genomically labeled arrays. Bioinformatics.

**Tools Used:**
- BWA: Li H, Durbin R (2009). Fast and accurate short read alignment. Bioinformatics.
- Pairtools: Open2C (2021). https://github.com/open2c/pairtools
- Cooler: https://github.com/open2c/cooler

---

## Support

For issues with:
- **Pipeline scripts:** Check logs in `logs/` directory
- **4DN documentation:** https://data.4dnucleome.org/
- **Cooler:** https://cooler.readthedocs.io/
- **Pairtools:** https://pairtools.readthedocs.io/

:P
