# -*- coding: utf-8 -*-
"""
Created on Tue Jun 10 11:36:29 2025

@author: Isabella Buck

This script visualizes motif similarity scores and empirical p-values for enhancer windows.
It generates scatterplots, overlays significant regions, and saves results to disk.
You can adjust directories and thresholds via command-line arguments.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib.patches as mpatches
import argparse
from statsmodels.stats.multitest import multipletests
from scipy.signal import find_peaks

from Motif_dist import data_dir, pval, file_place, cutoff_qval, perm, mut

# Argument parsing for flexible output and thresholds
parser = argparse.ArgumentParser(description="Visualize motif similarity and empirical p-values.")
parser.add_argument("--scatter_dir", type=str, default=rf"C:\\Users\\izzye\\OneDrive - Johns Hopkins\\Documents\\BRIGHT\\motif{mut}\full\images\scatterplots", help="Directory to save scatterplots.")
parser.add_argument("--pval", type=float, default=pval, help="Empirical p-value threshold for significance.")
parser.add_argument("--cutoff_qval", type=float, default=cutoff_qval, help="Q-value threshold for motif filtering.")
args = parser.parse_args()

scatter_dir = args.scatter_dir
pval = args.pval
cutoff_qval = args.cutoff_qval
os.makedirs(scatter_dir, exist_ok=True)

# List of bedGraph files to visualize
bedgraph_files = [
    f"{perm}st10_windows.bedgraph", f"{perm}st10_randoms.bedgraph",
    f"{perm}st10R_windows.bedgraph", f"{perm}st10R_randoms.bedgraph",
    f"{perm}IntP2A_windows.bedgraph", f"{perm}IntP2A_randoms.bedgraph",
    f"{perm}IntP2B_windows.bedgraph", f"{perm}IntP2B_randoms.bedgraph",
    f"{perm}pCRM_windows.bedgraph", f"{perm}pCRM_randoms.bedgraph",
    f"{perm}Peak9820_windows.bedgraph", f"{perm}Peak9820_randoms.bedgraph",
    f"{perm}Peak9821_windows.bedgraph", f"{perm}Peak9821_randoms.bedgraph",
    f"{perm}Peak9823_windows.bedgraph", f"{perm}Peak9823_randoms.bedgraph",
    f"{perm}Peak9817_windows.bedgraph", f"{perm}Peak9817_randoms.bedgraph",
    f"{perm}AEL4_windows.bedgraph", f"{perm}AEL4_randoms.bedgraph",
    f"{perm}5P3_windows.bedgraph", f"{perm}5P3_randoms.bedgraph"
]

# Significant enhancer regions
overlaps = {
    "sim_A1.0": (13056591, 13057759),
    "sim_st10": (13060273, 13060969),
    "sim_D2.1": (13064360, 13066277),
    "sim_E2.3": (13066277, 13067711),
    "EV_overlap": (13067711, 13068610),
    "sim_VT040842": (13068610, 13069801),
    "sim_1.6MLE": (13070722, 13072353)
}

def scatter_distance():
    """
    Plots motif similarity scores for each window region as a scatterplot.
    """
    for file in bedgraph_files:
        if "_windows.bedgraph" not in file:
            continue  # Skip non-window files

        file_path = os.path.join(data_dir, file)
        df = pd.read_csv(file_path, sep="\t", header=None, names=["chrom", "start", "end", "distance"])
        df["midpoint"] = (df["start"] + df["end"]) / 2
        motif_name = file.split("_")[0]  # Extract the enhancer name from the filename
        
        # Create scatterplot
        plt.figure(figsize=(14, 6))
        plt.scatter(df["midpoint"], df["distance"], color="dodgerblue", alpha=0.7, s=20)
        plt.xlabel("Genomic Midpoint (bp)")
        plt.ylabel(f"{file_place} Distance")
        name = f"{file_place} Distances: {motif_name} vs windows"
        plt.title(name)
        plt.ylim(0, 1)
        plt.xlim(13056267, 13086627)
        plt.grid(True, alpha=0.3)
        save_name = f"{motif_name}_distance_scatter.png"
        plot_path = os.path.join(scatter_dir, save_name)
        plt.savefig(plot_path, dpi=300)
        print(f"Scatterplot saved to {plot_path}")
        
        for start, end in overlaps.values():
            rect = mpatches.Rectangle(
                (start, plt.ylim()[0]),
                end - start,
                plt.ylim()[1] - plt.ylim()[0],
                linewidth=0,
                edgecolor=None,
                facecolor='lightblue',
                alpha=0.2,
                label=None
            )
            plt.gca().add_patch(rect)
        plt.show()
        plt.close()

def visualize():
    """
    Visualizes motif similarity and overlays empirical p-values for each window region.
    """
    for file in bedgraph_files:
        if "_windows.bedgraph" not in file:
            continue  # Skip non-window files
            
        file_path = os.path.join(data_dir, file)
        df = pd.read_csv(file_path, sep="\t", header=None, names=["chrom", "start", "end", file_place])
        
        # Compute the midpoint (average genomic coordinate)
        df["midpoint"] = (df["start"] + df["end"]) / 2
        
        # Sort by chromosome and midpoint
        df.sort_values(["chrom", "midpoint"], inplace=True)
        
        # Trim data of outliers as needed
        q_low = df[file_place].quantile(0)
        q_high = df[file_place].quantile(1)
        df_trimmed = df[(df[file_place] >= q_low) & (df[file_place] <= q_high)]
        avg_y = df_trimmed[file_place].median()
        
        # Extract a name from filename like "rand_IntP2B.fa" -> "IntP2B"
        enhancer_name = file.replace("_randoms.bedgraph", "").replace("_windows.bedgraph","")
        
        # Extract if the file is Random or Window
        types = file.split("_")[-1].replace("windows.bedgraph", "Windows").replace("randoms.bedgraph", "Randoms").replace("RandMerg", "Random (merged)")
        
        pval_file_csv = os.path.join(data_dir, file.replace("windows.bedgraph", "pvals.csv"))
        if os.path.exists(pval_file_csv):
            print(f"Looking for p-value file: {pval_file_csv}")
            pval_df = pd.read_csv(pval_file_csv)
            merged = pd.merge(
                df_trimmed,
                pval_df[["chrom", "start", "end", "raw_pval"]],
                on=["chrom", "start", "end"],
                how="left"
            )
            
            # Create scatterplot
            plt.figure(figsize=(14, 6))
            plt.axhline(y=avg_y, color='black', linestyle='--', linewidth=1, label=f'Average = {avg_y:.3f}')
            scatter = plt.scatter(
                merged["midpoint"],
                merged[file_place],
                c=merged["raw_pval"],
                cmap="plasma",
                alpha=0.9,
                s=20,
                label=None
            )
            plt.colorbar(scatter, label="Empirical p-value")
            
            # Overlay significant points
            sig = merged["raw_pval"] < pval
            plt.scatter(
                merged.loc[sig, "midpoint"],
                merged.loc[sig, file_place],
                color='lime',
                alpha=0.9,
                s=20,
                label=f"Empirical p < {pval}"
            )
            # Overlay points with p-value == 0
            # Shouldn't be any
            sig_zero = merged["raw_pval"] == 0
            plt.scatter(
                merged.loc[sig_zero, "midpoint"],
                merged.loc[sig_zero, file_place],
                color='red',
                alpha=0.9,
                s=20,
                label="p = 0"
            )
        else:
            print(f"Warning: No p-value file found for {file}")
        plt.legend(title="Chromosome")
        plt.xticks(np.linspace(df_trimmed["midpoint"].min(), df["midpoint"].max(), num=10, dtype=int), rotation=45)

        # Draw boxes for each region in overlaps
        for start, end in overlaps.values():
            rect = mpatches.Rectangle(
                (start, plt.ylim()[0]),
                end - start,
                plt.ylim()[1] - plt.ylim()[0],
                linewidth=0,
                edgecolor=None,
                facecolor='lightblue',
                alpha=0.2,
                label=None
            )
            plt.gca().add_patch(rect)
        plt.xlabel("Genomic Midpoint (bp)")
        plt.ylabel(f"{file_place} Distance")
        name = f"Genomic Distribution of {file_place} Distances for {perm} Q={cutoff_qval} P={pval}: {enhancer_name} vs {types}"
        save_name=  f"{file_place}_Q={cutoff_qval}_P={pval}_{enhancer_name}_vs_{types}.png"
        plt.title(name)
        plt.legend(title="Chromosome")
        plt.grid(True)
        plot_path = os.path.join(scatter_dir, save_name)
        plt.savefig(plot_path, dpi=300)
        print(f"Plot saved to {plot_path}")
        plt.show()

# OPTIONAL
def plot_all_and_matching_midpoints(file1, file2, label1="All", label2="Matching", color1="gray", color2="red", output_file=None):
    """
    Plots all points from file1 in gray, and overlays in red those with matching midpoints in file2.
    """
    # Read first file (bedgraph: chrom, start, end, score)
    df1 = pd.read_csv(file1, sep="\t", header=None, names=["chrom", "start", "end", "score"])
    df1["start"] = pd.to_numeric(df1["start"], errors="coerce")
    df1["end"] = pd.to_numeric(df1["end"], errors="coerce")
    df1["score"] = pd.to_numeric(df1["score"], errors="coerce")
    df1["midpoint"] = (df1["start"] + df1["end"]) // 2
    
    # Read second file (can be bedgraph or BED)
    df2 = pd.read_csv(file2, sep="\t", header=None, names=["chrom", "start", "end", "name", "pvalue", "type"])
    df2["start"] = pd.to_numeric(df2["start"], errors="coerce")
    df2["end"] = pd.to_numeric(df2["end"], errors="coerce")
    df2["midpoint"] = (df2["start"] + df2["end"]) // 2
    
    # Find matching midpoints (and chrom)
    match = pd.merge(df1, df2[["chrom", "midpoint"]], on=["chrom", "midpoint"], how="inner")
    
    plt.figure(figsize=(14, 6))
    plt.scatter(df1["midpoint"], df1["score"], color=color1, alpha=0.7, label=label1)
    
    if not match.empty:
        plt.scatter(match["midpoint"], match["score"], color=color2, alpha=0.8, label=label2)
    
    # Draw boxes for each region in overlaps
    for start, end in overlaps.values():
        rect = mpatches.Rectangle(
            (start, plt.ylim()[0]),
            end - start,
            plt.ylim()[1] - plt.ylim()[0],
            linewidth=0,
            edgecolor=None,
            facecolor='lightblue',
            alpha=0.2,
            label=None
        )
        plt.gca().add_patch(rect)
    plt.xlabel("Genomic Midpoint (bp)")
    plt.ylabel("Score")
    plt.title("All points (gray) and matching midpoints (red, boxed overlaps)")
    plt.legend()
    plt.tight_layout()
    if output_file:
        plt.savefig(output_file, dpi=300)
        print(f"Plot saved to {output_file}")
    plt.show()

def pvalue(window_file, random_file, label, output_bedfile=None):
    """
    Computes empirical p-values for each window and applies Benjamini-Hochberg correction.
    Saves results to BED and CSV files.
    """
    if output_bedfile is None:
        output_bedfile = os.path.join(data_dir, f"{label}.bed")
        
    window_path = os.path.join(data_dir, window_file)
    random_path = os.path.join(data_dir, random_file)
    
    win_df = pd.read_csv(window_path, sep="\t", header=None, names=["chrom", "start", "end", "dist"])
    rand_df = pd.read_csv(random_path, sep="\t", header=None, names=["chrom", "start", "end", "dist"])
    win_df.sort_values(["dist"], inplace=True)
    
    # Compute raw empirical p-values
    raw_pvals = []
    coords = []
    for _, row in win_df.iterrows():
        distance = row["dist"]
        
        # Number of random points closer than minimum window value
        k = (rand_df["dist"] <= distance).sum()
        
        # Number of random points
        N = len(rand_df)
        
        # P-value calculation with +1 correction
        p_raw = (k+1) / N if N != 0 else 1.0
        raw_pvals.append(p_raw)
        coords.append((row["chrom"], row["start"], row["end"]))
        
    # Apply Benjamini-Hochberg correction
    pvals_array = np.array(raw_pvals)
    _, pvals_bh, _, _ = multipletests(pvals_array, alpha=0.05, method='fdr_bh')
    
    # Output BED file using raw p-values for all windows
    bed_lines = []
    name_score = 1
    for (chrom, start, end), p_raw in zip(coords, raw_pvals):
        name = f"{label}-{name_score}"
        score = round(p_raw, 6)
        strand = "."
        if p_raw <= pval:
            bed_lines.append([chrom, start, end, name, score, strand])
        name_score += 1
    
    # Save the BED file
    bed_df = pd.DataFrame(bed_lines, columns=["chrom", "start", "end", "name", "score", "strand"])
    bed_df.to_csv(output_bedfile, sep="\t", header=False, index=False)
    print(f"BED file with raw p-values {pval} saved to: {output_bedfile}")
    
    
    # OPTIONAL BH-adjustment
    # Save raw and BH-adjusted p-values to CSV for heatmap usage
    raw_bh_df = pd.DataFrame({
        "chrom": [c[0] for c in coords],
        "start": [c[1] for c in coords],
        "end":   [c[2] for c in coords],
        "raw_pval": raw_pvals,
        "bh_pval": pvals_bh
    })
    raw_bh_path = os.path.join(data_dir, f"{label}_pvals.csv")
    raw_bh_df.to_csv(raw_bh_path, index=False)
    print(f"Raw + BH p-values {pval} saved to: {raw_bh_path}")
    
    
    # Identify BH-significant points (adjusted p < 0.05)
    significant = pvals_bh < 0.05
    
    # Create color array: red for significant, gray for non-significant
    colors = ['crimson' if sig else 'lightgray' for sig in significant]
    
    # Plot BH Corrected ovalues as a scatterplot
    plt.figure(figsize=(6, 5))
    plt.scatter(raw_pvals, pvals_bh, color=colors, alpha=0.6, edgecolor='k', linewidth=0.2)
    plt.plot([0, 1], [0, 1], linestyle="--", color="black", alpha=0.5)
    plt.axhline(0.05, color="blue", linestyle=":", label="BH p = 0.05")
    plt.axvline(0.05, color="green", linestyle=":", label="Raw p = 0.05")
    
    # Labels and formatting
    plt.xlabel("Raw p-value")
    plt.ylabel("BH Adjusted p-value")
    plt.title(f"{label}: Raw vs. BH-Adjusted P-values")
    plt.legend(loc="lower right", fontsize=8)
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def plot_jaccard_peaks(bedgraph_file, value_col="distance", out_dir=None, height=None, distance=10, save_peaks=True):
    """
    Finds and plots peaks in motif similarity scores across the genome.
    Optionally saves peak coordinates to a TSV file.
    """
    df = pd.read_csv(bedgraph_file, sep="\t", header=None, names=["chrom", "start", "end", value_col])
    df["midpoint"] = (df["start"] + df["end"]) // 2
    
    # Find peaks
    inverted = -df[value_col].values
    peaks, properties = find_peaks(inverted, height=.05, distance=distance)
    print(f"Found {len(peaks)} peaks in {os.path.basename(bedgraph_file)}.")
    
    # Create plot
    plt.figure(figsize=(14, 5))
    plt.plot(df["midpoint"], df[value_col], label=value_col)
    plt.scatter(df["midpoint"].iloc[peaks], df[value_col].iloc[peaks], color="red", label="Peaks")
    plt.xlabel("Genomic Midpoint (bp)")
    plt.ylabel(value_col)
    plt.ylim(0)
    plt.title(f"{value_col} across genome with peaks: {os.path.basename(bedgraph_file)}")
    plt.legend()
    plt.tight_layout()
    
    # Draw boxes for each region in overlaps
    for start, end in overlaps.values():
        rect = mpatches.Rectangle(
            (start, plt.ylim()[0]),
            end - start,
            plt.ylim()[1] - plt.ylim()[0],
            linewidth=0,
            edgecolor=None,
            facecolor='lightblue',
            alpha=0.2,
            label=None
        )
        plt.gca().add_patch(rect)
        
    # OPTIONALLY save peaks
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
        plot_path = os.path.join(out_dir, f"{os.path.basename(bedgraph_file).replace('.bedgraph','')}_{value_col}_peaks.png")
        plt.savefig(plot_path, dpi=300)
        print(f"Peak plot saved to {plot_path}")
    plt.show()
    if save_peaks and out_dir:
        peak_df = df.iloc[peaks][["chrom", "start", "end", value_col]]
        peak_file = os.path.join(out_dir, f"{os.path.basename(bedgraph_file).replace('.bedgraph','')}_{value_col}_peaks.tsv")
        peak_df.to_csv(peak_file, sep="\t", index=False)
        print(f"Peak coordinates saved to {peak_file}")
    else:
        print("❌ NOT SAVED")

if __name__ == "__main__":
    # Run analysis and visualizations
    pvalue(f"{perm}st10_windows.bedgraph", f"{perm}st10_randoms.bedgraph", f"{perm}st10")
    pvalue(f"{perm}st10R_windows.bedgraph", f"{perm}st10R_randoms.bedgraph", f"{perm}st10R")
    pvalue(f"{perm}IntP2A_windows.bedgraph", f"{perm}IntP2A_randoms.bedgraph", f"{perm}IntP2A")
    pvalue(f"{perm}IntP2B_windows.bedgraph", f"{perm}IntP2B_randoms.bedgraph", f"{perm}IntP2B")
    pvalue(f"{perm}pCRM_windows.bedgraph", f"{perm}pCRM_randoms.bedgraph", f"{perm}pCRM")
    pvalue(f"{perm}Peak9820_windows.bedgraph", f"{perm}Peak9820_randoms.bedgraph", f"{perm}Peak9820")
    pvalue(f"{perm}Peak9821_windows.bedgraph", f"{perm}Peak9821_randoms.bedgraph", f"{perm}Peak9821")
    pvalue(f"{perm}Peak9823_windows.bedgraph", f"{perm}Peak9823_randoms.bedgraph", f"{perm}Peak9823")
    pvalue(f"{perm}Peak9817_windows.bedgraph", f"{perm}Peak9817_randoms.bedgraph", f"{perm}Peak9817")
    pvalue(f"{perm}AEL4_windows.bedgraph", f"{perm}AEL4_randoms.bedgraph", f"{perm}AEL4")
    pvalue(f"{perm}5P3_windows.bedgraph", f"{perm}5P3_randoms.bedgraph", f"{perm}5P3")
    
    visualize()
    scatter_distance()
    
    for file in bedgraph_files:
        if "_windows.bedgraph" in file:
            plot_jaccard_peaks(os.path.join(data_dir, file), value_col="distance", out_dir=scatter_dir, height=0.5)
    
    print("finished 1")