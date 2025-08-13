# -*- coding: utf-8 -*-
"""
Created on Fri Jul 11 12:45:30 2025

@author: Isabella Buck

This script permutes a percentage of bases in enhancer sequences.
You can adjust the percentage and output directory as needed.
"""

import numpy as np
import os

amount = 5  # Percentage of bases to permute
output_dir = r"C:\Users\izzye\OneDrive - Johns Hopkins\Documents\BRIGHT\seq\enhancers"  # Directory to save permuted FASTA files

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)


# Maps of enhancer and window sequences (FASTA format)
enhancers = ['IntP2B', 'IntP2A', 'st10', 'st10R', 'pCRM', 'Peak9820',
             'Peak9821', 'AEL4', '5P3', 'Peak9817', 'Peak9823']

enhancer_regions = {}

for enhancer in enhancers:

    # Read enhancer sequence
    enhancer_path = os.path.join(output_dir, enhancer + ".txt") 
    try:
        with open(enhancer_path, "r") as file:
            content = file.read()
            # Remove newlines and split if needed
            enhancer_seq = content.strip()
        
        # Pair enhancer sequence with random and window sequences
        enhancer_regions[enhancer] = enhancer_seq

    except Exception as e:
        print(f"Error reading {enhancer_path}: {str(e)}")
        continue
    
def permute(seq, x):
    """
    Permute x% of the bases in seq, leaving the rest unchanged.
    """
    seq = np.array(list(seq))
    n = len(seq)
    num_bases = int(n * x / 100)
    indices = np.random.choice(n, num_bases, replace=False)
    bases = seq[indices]
    np.random.shuffle(bases)
    seq[indices] = bases
    return ''.join(seq)

# Permute each enhancer and save to FASTA file
for name, dna in enhancer_regions.items():
    new_dna = permute(dna, amount)
    print(f"Original ({name}):", dna)
    print(f"Permuted ({name}, {amount}%):", new_dna)
    output_path = os.path.join(output_dir, f"{name}_{amount}.fa")
    with open(output_path, "w") as f:
        f.write(f">{name}_{amount}\n")
        f.write(new_dna + "\n")