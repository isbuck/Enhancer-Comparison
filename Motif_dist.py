# -*- coding: utf-8 -*-
"""
Created on Wed Jun 25 13:30:18 2025

@author: Isabella Buck

This script compares motif content between enhancer regions and windows/randoms using Jaccard or Summed Canberra distance.
It outputs bedGraph files with motif similarity scores.
"""

import pandas as pd
import os
from scipy.spatial.distance import canberra


# Modify these values as needed
cutoff_metric = 'q-value'    # Column name for motif significance cutoff
cutoff_qval = 0.05          # Q-value threshold for motif filtering
pval = 0.01                 # P-value threshold
path = 1                    # 1 for Jaccard, 2 for Canberra distance
perm = ""                   # Prefix for permuted files (e.g., "" for normal, "perm_" for permuted)
base_dir = r"C:\Users\izzye\OneDrive - Johns Hopkins\Documents\BRIGHT\motif\full"  # Base directory for motif files

# Set up directories based on path and mutation
mut = "" if perm == "" else "\mutation"
if path == 1:
    data_dir = rf"{base_dir}{mut}\jaccard\qval={cutoff_qval}\pval={pval}\bedgraph"
    file_place = "jaccard"
elif path == 2:
    data_dir = rf"{base_dir}{mut}\canberra\qval={cutoff_qval}\pval={pval}\bedgraph"
    file_place = "canberra"
os.makedirs(data_dir, exist_ok=True)
print(f"Data directory: {data_dir}")

window_files = []
random_files = []
enhancer_files = []

# File lists
enhancers = [
    "st10", "st10R", "IntP2A", "IntP2B", "pCRM", "Peak9820",
    "Peak9821", "Peak9823", "Peak9817", "AEL4", "5P3"
]

for enhancer in enhancers:
    window_files.append(enhancer + "_window_fimo.tsv")
    random_files.append(enhancer + "_random_fimo.tsv")
    enhancer_files.append(perm + enhancer + "_fimo.tsv")

# Directories for enhancer, window, and random motif files
enhancer_dir = os.path.join(base_dir, "enhancers")
window_dir = os.path.join(base_dir, "windows")
random_dir = os.path.join(base_dir, "randoms")

os.makedirs(data_dir, exist_ok=True)

def extract_motifs(fimo_file):
    """
    Extracts unique motif IDs from a FIMO file, filtered by cutoff_metric and cutoff_qval.
    Returns a set of motif IDs.
    """
    try:
        df = pd.read_csv(fimo_file, sep='\t', comment='#')
        if df.empty: 
            print(f"XX{fimo_file} is empty. No motif hits.")
            return set()
        if cutoff_metric in df.columns:
            df = df[df[cutoff_metric] < cutoff_qval]
        print(f"{fimo_file}: {len(df)} motif hits after q-value < {cutoff_qval}")
        return set(df['motif_id'].unique())
    except pd.errors.EmptyDataError:
        print(f"XX {fimo_file} is empty or has no columns.")
        return set()

def jaccard_distance(set1, set2):
    """
    Computes Jaccard distance between two sets of motif IDs.
    """
    intersection = len(set1 & set2)
    union = len(set1 | set2)
    print(f"Jaccard intersection: {intersection}, union: {union}")
    return 1 - (intersection / union) if union else 0.0

def write_bedgraph_line(out_file, region_id, score):
    """
    Writes a single line to a bedGraph file for a region and its motif similarity score.
    """
    try:
        chrom, start, end = region_id.split("_")
        out_file.write(f"{chrom}\t{start}\t{end}\t{score:.6f}\n")
    except Exception as e:
        print(f"XX Invalid region ID format: {region_id} – {e}")

def extract_motif_matrix(fimo_file):
    """
    Extracts a motif score matrix from a FIMO file for Canberra distance calculation.
    Returns a DataFrame with sequence_name as rows and motif_id as columns.
    """
    df = pd.read_csv(fimo_file, sep='\t', comment='#')
    if cutoff_metric in df.columns:
        df = df[df[cutoff_metric] < cutoff_qval]
    print(f"{fimo_file}: {len(df)} motif hits after q-value < {cutoff_qval}")
    if not {'sequence_name', 'motif_id', 'score'}.issubset(df.columns):
        raise ValueError("FIMO file must contain 'sequence_name', 'motif_id', and 'score' columns")
    df['sequence_name'] = df['sequence_name'].astype(str).str.strip()
    df['motif_id'] = df['motif_id'].astype(str).str.strip()
    grouped = df.groupby(['sequence_name', 'motif_id'], as_index=False)['score'].sum()
    matrix = grouped.pivot(index='sequence_name', columns='motif_id', values='score').fillna(0)
    matrix = matrix.sort_index().sort_index(axis=1)
    return matrix

def canberra_distance(enhancer_matrix, window_matrix):
    """
    Computes Summed Canberra distance between enhancer motif profile and each window/random region.
    Returns a dictionary mapping region_id to distance.
    """
    all_motifs = sorted(set(enhancer_matrix.columns) | set(window_matrix.columns))
    enhancer_matrix = enhancer_matrix.reindex(columns=all_motifs, fill_value=0)
    window_matrix = window_matrix.reindex(columns=all_motifs, fill_value=0)
    enhancer_profile = enhancer_matrix.sum(axis=0).values
    distances = {}
    for region_id, row in window_matrix.iterrows():
        window_profile = row.values
        distance = canberra(window_profile, enhancer_profile)
        distances[region_id] = distance
    return distances

if __name__ == "__main__":
    # Main processing loop for each enhancer
    for enhancer_file in enhancer_files:
        enhancer_id = enhancer_file.split("_fimo")[0]
        print("Processing enhancer:", enhancer_id)
        enhancer_path = os.path.join(enhancer_dir, enhancer_file)

        try:
            if path == 2:
                enhancer_matrix = extract_motif_matrix(enhancer_path)
            else:
                enhancer_motifs = extract_motifs(enhancer_path)
        except FileNotFoundError:
            print(f"XX Missing enhancer file: {enhancer_path}")
            continue

        window_out = os.path.join(data_dir, f"{enhancer_id}_windows.bedGraph")
        random_out = os.path.join(data_dir, f"{enhancer_id}_randoms.bedGraph")

        with open(window_out, 'w') as wf, open(random_out, 'w') as rf:
            # Process window files
            for window_file in window_files:
                windows_id = window_file.split("_")[0]
                window_id = perm + windows_id
                if window_id == enhancer_id:
                    path_win = os.path.join(window_dir, window_file)
                    try:
                        if path == 2:
                            window_matrix = extract_motif_matrix(path_win)
                            distances = canberra_distance(enhancer_matrix, window_matrix)
                            for region_id, score in distances.items():
                                write_bedgraph_line(wf, region_id, score)
                        else:
                            df = pd.read_csv(path_win, sep='\t', comment="#")
                            for region_id, group in df.groupby("sequence_name"):
                                region_motifs = set(group["motif_id"])
                                dist = jaccard_distance(enhancer_motifs, region_motifs)
                                write_bedgraph_line(wf, region_id, dist)
                    except Exception as e:
                        print(f"XX Error processing {path_win}: {e}")

            # Process random files
            for random_file in random_files:
                randoms_id = random_file.split("_")[0]
                random_id = perm + randoms_id
                if random_id == enhancer_id:
                    path_rand = os.path.join(random_dir, random_file)
                    try:
                        if path == 2:
                            random_matrix = extract_motif_matrix(path_rand)
                            distances = canberra_distance(enhancer_matrix, random_matrix)
                            for region_id, score in distances.items():
                                write_bedgraph_line(rf, region_id, score)
                        else:
                            df = pd.read_csv(path_rand, sep='\t', comment="#")
                            for region_id, group in df.groupby("sequence_name"):
                                region_motifs = set(group["motif_id"])
                                dist = jaccard_distance(enhancer_motifs, region_motifs)
                                write_bedgraph_line(rf, region_id, dist)
                    except Exception as e:
                        print(f"XX Error processing {path_rand}: {e}")

        print(f"Saved: {window_out}")
        print(f"Saved: {random_out}")

        
