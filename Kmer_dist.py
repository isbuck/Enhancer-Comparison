# -*- coding: utf-8 -*-
"""
Created on Mon Jun  9 15:28:33 2025

@author: Isabella Buck

This script computes k-mer distances (through various Alfpy metrics) between enhancer sequences and window/random sequences,
saving results to bedGraph files.
"""

from alfpy.utils import seqrecords
from alfpy import word_pattern, word_vector, word_distance
from alfpy.utils import distmatrix
import os
from Bio.Seq import Seq
import subprocess

# Modify these values as needed
k = 3                    # k-mer size (e.g., 3 for 3-mers)
distance = 'canberra'    # Alfpy distance metric (e.g., 'canberra', 'euclidean', 'cosine')
pval = 0.05             # P-value threshold
batch_size = 10         # Number of sequences per batch
count_max = 2000        # Maximum number of batches to process (ensure batch_size * count_max > 10,000)
input_dir = r"C:\Users\izzye\OneDrive - Johns Hopkins\Documents\BRIGHT\seq\windows_randoms"  # Directory containing input FASTA files
base_data_dir = r"C:\Users\izzye\OneDrive - Johns Hopkins\Documents\BRIGHT\seq"  # Base directory for output files
enhancer_dir = r"C:\Users\izzye\OneDrive - Johns Hopkins\Documents\BRIGHT\seq\enhancers"  # Base directory for enhancer files


group_name = f"{distance}"

# Set output subfolder path based on p-value
if pval == 0.05:
    subfolder_path = os.path.join(base_data_dir, group_name, f"k={k}")
else:
    subfolder_path = os.path.join(base_data_dir, group_name, f"k={k}", str(pval))

os.makedirs(subfolder_path, exist_ok=True)

# Maps of enhancer and window sequences (FASTA format)
enhancers = ['IntP2B', 'IntP2A', 'st10', 'st10R', 'pCRM', 'Peak9820',
             'Peak9821', 'AEL4', '5P3', 'Peak9817', 'Peak9823']

random_map = {}
window_map = {}

for enhancer in enhancers:
    random = "rand_" + enhancer + ".fa"
    window = "window_" + enhancer + ".fa"

    # Read enhancer sequence
    enhancer_path = os.path.join(enhancer_dir, enhancer + ".txt") 
    try:
        with open(enhancer_path, "r") as file:
            content = file.read()
            # Remove newlines and split if needed
            enhancer_seq = content.strip()
            if "_" in enhancer_seq:
                enhancer_seq = enhancer_seq.split("_")[0]
        
        # Pair enhancer sequence with random and window sequences
        random_map[random] = enhancer_seq
        window_map[window] = enhancer_seq
    except Exception as e:
        print(f"Error reading {enhancer_path}: {str(e)}")
        continue

base_data_dir = r"C:\Users\izzye\OneDrive - Johns Hopkins\Documents\BRIGHT\seq"
group_name = f"{distance}"

# Set output subfolder path based on p-value
if pval == 0.05:
    subfolder_path = os.path.join(base_data_dir, group_name, f"k={k}")
else:
    subfolder_path = os.path.join(base_data_dir, group_name, f"k={k}", str(pval))

os.makedirs(subfolder_path, exist_ok=True)

def check_files(seq_dict):
    """
    Check if all files in seq_dict exist in the input directory.
    Prints status for each file.
    """
    for file in seq_dict.keys():  
        file_path = os.path.join(input_dir, file)
        if os.path.exists(file_path):  
            print(f"File found: {file_path}")
        else:
            print(f"XX File missing: {file_path}")

def save_bedgraph_line(out_file, header, score):
    """
    Write one line to a BedGraph file.
    header: "chr:start-end" or "chr_start_end"
    score: computed distance value (string or float)
    """
    if ":" in header and "-" in header:
        chrom, coords = header.split(":", 1)
        start, end = coords.split("-")
        start = int(start)
        end = int(end)
        out_file.write(f"{chrom}\t{start}\t{end}\t{score}\n")
    elif "_" in header:
        chrom, coords = header.split("_",1)
        start, end = coords.split("_")
        out_file.write(f"{chrom}\t{start}\t{end}\t{score}\n")
    else:
        print(f"Skipping header without valid coordinates: {header}")

def compute_distance(seq_dict, label, input_dir, subfolder_path, batch_size):
    """
    Compute k-mer distance between enhancer and window/random sequences.
    Saves results to bedGraph files.
    Args:
        seq_dict (dict): Mapping of filenames to enhancer sequences.
        label (str): 'window' or 'random'.
        input_dir (str): Directory containing input FASTA files.
        subfolder_path (str): Directory to save output files.
        batch_size (int): Number of sequences per batch.
    """
    check_files(seq_dict)

    for file, enhancer_seq in seq_dict.items():
        print(f"\nComparing enhancer to sequences in: {file}")

        # Read all sequences from FASTA file
        with open(os.path.join(input_dir, file)) as fh:
            recs_full = seqrecords.read_fasta(fh)

        # Extract enhancer name and type from filename
        enhancer_name = file.split("_", 1)[1].replace(".fa", "")
        rand_win = file.split("_",1)[0].replace("rand","Random").replace("window","Window")
        enhancer_seq_obj = Seq(enhancer_seq)
        rev_complement_seq = str(enhancer_seq_obj.reverse_complement())

        total = len(recs_full.id_list)
        print(f"Total records in {file}: {total}")

        count = 0
        start = 0
        end = batch_size

        # Output files for enhancer and reverse complement distances
        enhancer_outfile = os.path.join(subfolder_path, f"BedResults{rand_win}_{enhancer_name}.bedGraph")
        reverse_outfile = os.path.join(subfolder_path, f"BedResults{rand_win}_{enhancer_name}_rev.bedGraph")

        with open(enhancer_outfile, 'w') as enhancer_file, open(reverse_outfile, 'w') as reverse_file:
            print(f"Saving enhancer distances to: {enhancer_outfile}")
            print(f"Saving reverse complement distances to: {reverse_outfile}")
            
            while count < count_max and start < total:
                # Prepare batch of sequences
                recs = recs_full.__class__()
                recs.id_list = recs_full.id_list[start:end].copy()
                recs.seq_list = recs_full.seq_list[start:end].copy()
                recs.add(enhancer_name, seq=enhancer_seq)
                recs.add(f"{enhancer_name}_rev", seq=rev_complement_seq)

                # Compute k-mer patterns and distances
                pattern = word_pattern.create(recs.seq_list, word_size=k)
                counts_obj = word_vector.Counts(recs.length_list, pattern)
                dist_obj = word_distance.Distance(counts_obj, distance)
                matrix = distmatrix.create(recs.id_list, dist_obj)

                last_seq_idx = len(recs) - 1
                second_last_seq_idx = len(recs) - 2
                distance_values = matrix.data

                # Save enhancer distances
                for idx, target_id in enumerate(recs.id_list):
                    if idx == second_last_seq_idx or idx == last_seq_idx:
                        continue
                    d = distance_values[second_last_seq_idx][idx]
                    header = f"{target_id}"
                    save_bedgraph_line(enhancer_file, header, f"{d:.6f}")

                # Save reverse complement distances
                for idx, target_id in enumerate(recs.id_list):
                    if idx == last_seq_idx or idx == second_last_seq_idx:
                        continue
                    d = distance_values[last_seq_idx][idx]
                    header = f"{target_id}"
                    save_bedgraph_line(reverse_file, header, f"{d:.6f}")

                count += 1
                start += batch_size
                end = start + batch_size

if __name__ == "__main__":
    
    # Run distance computation for both window and random sequences
    compute_distance(window_map, 'window', input_dir, subfolder_path, batch_size)
    compute_distance(random_map, 'random', input_dir, subfolder_path, batch_size)

    # Optionally run visualization script
    results_vis = subprocess.run(["python", "Kmer_visualizer.py"], capture_output=True, text=True, check=True)
    print(results_vis.stdout)

    print("finished")
