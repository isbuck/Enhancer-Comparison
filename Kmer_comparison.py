# -*- coding: utf-8 -*-
"""
Created on Mon Jun 16 14:01:33 2025

@author: Isabella Buck

This script compares merged BED files to reference regions, generates summary tables and heatmaps,
and visualizes enhancer/overlap regions. 
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os
import numpy as np
import seaborn as sns

from Kmer_dist import k, distance, enhancers

# Direct parameter configuration - edit these values as needed
data_dir = rf"C:\Users\izzye\OneDrive - Johns Hopkins\Documents\BRIGHT\seq\{distance}\k={k}"
output_dir = r"C:\Users\izzye\OneDrive - Johns Hopkins\Documents\BRIGHT\seq\csv"
comp_dir = r"C:\Users\izzye\OneDrive - Johns Hopkins\Documents\BRIGHT\seq\images\overlap_plots"
heat_dir = r"C:\Users\izzye\OneDrive - Johns Hopkins\Documents\BRIGHT\seq\images\heatmaps"

run_name = f"{distance}_k={k}"

# List of merged BED files to process

files = []

# Load bed files
for enhancer in enhancers:
    files.append(enhancer + "_merged_sorted.bed")
    
coords = (13056267, 13086627)

# Reference regions for overlap
overlaps = {
    "sim_A1.0": (13056591, 13057759),
    "sim_st10": (13060273, 13060969),
    "sim_D2.1": (13064360, 13066407),
    "sim_E2.3": (13066277, 13068610),
    "sim_VT040842": (13067711, 13069801),
    "sim_1.6MLE": (13070722, 13072353)
}

def comparison(bed_file):
    """
    Summarizes p-values and ranks for enhancer overlaps, outputs CSV tables and heatmap.
    """
    results = []
    
    # Loop through each BED file
    for file in files:
        name = file.split('_')[0]
        file_path = os.path.join(data_dir, file)
        
        # Read BED file into DataFrame
        df = pd.read_csv(file_path, sep="\t", header=None, names=["chrom", "start", "end", "pvalue"])
        
        # Rank p-values (lower is better)
        df["rank"] = df["pvalue"].rank(method="dense", ascending=True)
        label = file.split("_",2)[0]
        
        # Iterate through each row to check overlaps
        for _, row in df.iterrows():
            pvalue = row["pvalue"]
            rank = row["rank"]
            start = row["start"]
            end = row["end"]
            
            # Print the best (lowest) p-value region for each file
            if int(rank) == 1:
                print(f"{name} {start}, {end}, {pvalue}")
            
            # Check if region overlaps any reference region
            for sim_label, (start_ref, end_ref) in overlaps.items():
                if start_ref <= start <= end_ref or start_ref <= end <= end_ref or start <= start_ref <= end or start <= end_ref <= end:
                    results.append([label, sim_label, pvalue, rank])
                    
    # Create DataFrame of results
    output_df = pd.DataFrame(results, columns=["Enhancer", "Sequence", "p-value", "Rank"])
    
    # Group duplicate entries & store lists of values per Enhancer+Sequence
    grouped_df = output_df.groupby(["Enhancer", "Sequence"]).agg(list).reset_index()
    
    # Enforce order for Enhancer & Sequence
    grouped_df["Enhancer"] = pd.Categorical(grouped_df["Enhancer"], categories=[f.split("_", 2)[0] for f in files], ordered=True)
    grouped_df["Sequence"] = pd.Categorical(grouped_df["Sequence"], categories=list(overlaps.keys()), ordered=True)
    
    # Convert lists into comma-separated strings
    grouped_df["p-value"] = grouped_df["p-value"].apply(lambda x: ', '.join(map(str, x)))
    grouped_df["Rank"] = grouped_df["Rank"].apply(lambda x: ', '.join(map(str, x)))
    
    # Sort the grouped DataFrame
    grouped_df = grouped_df.sort_values(["Enhancer", "Sequence"])
    
    # Pivot tables for p-values and ranks
    pivot_pval = grouped_df.pivot(index="Enhancer", columns="Sequence", values="p-value")
    pivot_rank = grouped_df.pivot(index="Enhancer", columns="Sequence", values="Rank")
    
    # Save tables to CSV
    pivot_pval.to_csv(os.path.join(output_dir, f"{run_name}_pvalues_output.csv"))
    pivot_rank.to_csv(os.path.join(output_dir, f"{run_name}_ranks_output.csv"))
    
    # Prepare numeric p-values for heatmap coloring
    pval_numeric = grouped_df.copy()
    pval_numeric["p-value"] = pval_numeric["p-value"].apply(
        lambda x: min(map(float, str(x).split(','))) if x else np.nan
    )
    pivot_pval_numeric = pval_numeric.pivot(index="Enhancer", columns="Sequence", values="p-value")
    
    # Ensure all rows/columns are represented
    pivot_pval_numeric = pivot_pval_numeric.reindex(
        index=[f.split("_", 2)[0] for f in files],
        columns=["sim_A1.0", "sim_st10", "sim_D2.1", "sim_E2.3", "sim_VT040842", "sim_1.6MLE"]
    )
    
    # Align annotation array to numeric data shape
    pivot_pval = pivot_pval.reindex_like(pivot_pval_numeric)
    pivot_rank = pivot_rank.reindex_like(pivot_pval_numeric)
    annot_array = pivot_pval.copy()
    
    # Build annotation text for heatmap
    for enhancer in annot_array.index:
        for sequence in annot_array.columns:
            pval_text = pivot_pval.loc[enhancer, sequence]
            rank_text = pivot_rank.loc[enhancer, sequence]
            if pd.isna(pval_text) or pd.isna(rank_text):
                annot_array.loc[enhancer, sequence] = ""
            else:
                annot_array.loc[enhancer, sequence] = f"{pval_text}\n({rank_text})"
                
    # Plot heatmap
    plt.figure(figsize=(10, 6))
    sns.heatmap(
        pivot_pval_numeric,  # values for color mapping
        annot=annot_array,   # text annotations
        fmt="",              # raw text, not a number format
        cmap="coolwarm_r",
        cbar_kws={'label': 'p-value'},
        linewidths=0.5,
        linecolor="gray",
        mask=pivot_pval_numeric.isnull()
    )
    plt.title(f"Heatmap with Rank ({run_name})")
    plt.tight_layout()
    combo_path = os.path.join(heat_dir, f"{run_name}_pval_rank_heatmap.png")
    plt.savefig(combo_path, dpi=300)
    plt.show()

def compare(files):
    """
    Visualizes enhancer regions and reference overlaps as colored boxes.
    """
    fig, ax = plt.subplots(figsize=(12, 4))
    y_position = 0
    
    # Plot enhancer regions as blue horizontal boxes
    for file in files:
        file_path = os.path.join(data_dir, file)
        df = pd.read_csv(file_path, sep="\t", header=None, names=["chrom", "start", "end", "pvalue"])
        num = 0
        for _, row in df.iterrows():
            offset = .4 if num % 2 == 0 else .2
            ax.add_patch(patches.Rectangle((row["start"], y_position), row["end"] - row["start"], 0.5, color="blue", alpha=0.6))
            
            # Add p-value below each box
            ax.text((row["start"] + row["end"]) / 2, y_position - offset, f"p = {row['pvalue']:.4f}",
                        ha="center", fontsize=4, color="black")
            num += 1
        y_position += 1  # Move to next row for each BED file
    
    # Plot overlap regions as purple horizontal boxes
    for i, (name, (start, end)) in enumerate(overlaps.items()):
        y_top = y_position + i + 0.5
        ax.add_patch(patches.Rectangle((start, y_top - 0.5), end - start, 0.5, color="purple", alpha=0.6))
       
        # Add vertical dashed lines at region boundaries
        ax.plot([start, start], [y_top,-1], linestyle="dashed", color="navy", alpha=0.9)
        ax.plot([end, end], [y_top, -1], linestyle="dashed", color="navy", alpha=0.9)
        ax.text(start + 100, y_position + i + 0.2, name, fontsize=5, color="indigo", weight="bold")
    
    # Adjust plot settings
    ax.set_xlim(coords)
    ax.set_ylim(-1, y_position + len(overlaps) + 1)
    y_positions = [y + 0.25 for y in range(len(files) + len(overlaps))]
    ax.set_yticks(y_positions)
    ax.set_yticklabels([bed_file.split("_", 1)[0] for bed_file in files] + list(overlaps.keys()))
    ax.set_xlabel("Genomic Coordinates")
    ax.set_title(f"{run_name} Enhancer vs. Overlap Regions")
    plot_path = os.path.join(comp_dir, f"{run_name}_enhancer_overlap_plot.png")
    plt.savefig(plot_path, dpi=300, bbox_inches="tight")
    plt.show()

if __name__ == "__main__":
    # Run comparison and visualization
    compare(files)

    comparison(files)
