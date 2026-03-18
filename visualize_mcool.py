#!/usr/bin/env python3
"""Simple Hi-C visualization for both replicates"""

import cooler
import matplotlib.pyplot as plt
import numpy as np

# Configuration
RESOLUTION = 10000  # 100kb
CHROMOSOMES = ['chr1', 'chr2']

# Load replicates
print("Loading matrices...")
rep1 = cooler.Cooler(f'/YOUR DIRECTORY/contact_maps/rep1_merged.mcool::resolutions/{RESOLUTION}')
rep2 = cooler.Cooler(f'/YOUR DIRECTORY/contact_maps/rep2_merged.mcool::resolutions/{RESOLUTION}')

# Visualize each chromosome
for chrom in CHROMOSOMES:
    print(f"Plotting {chrom}...")
    
    # Fetch data
    m1 = rep1.matrix(balance=True).fetch(chrom)
    m2 = rep2.matrix(balance=True).fetch(chrom)
    
    # Create side-by-side plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 9))
    
    # Rep1
    im1 = ax1.imshow(np.log10(m1 + 1), cmap='YlOrRd')
    ax1.set_title(f'Replicate 1 - {chrom}', fontsize=16)
    plt.colorbar(im1, ax=ax1, label='log10(contacts)')
    
    # Rep2
    im2 = ax2.imshow(np.log10(m2 + 1), cmap='YlOrRd')
    ax2.set_title(f'Replicate 2 - {chrom}', fontsize=16)
    plt.colorbar(im2, ax=ax2, label='log10(contacts)')
    
    # Save
    plt.tight_layout()
    plt.savefig(f'{chrom}_comparison.png', dpi=300)
    print(f"  Saved: {chrom}_comparison.png")
    plt.close()

print("\nDone! Created:")
for chrom in CHROMOSOMES:
    print(f"  - {chrom}_comparison.png")
