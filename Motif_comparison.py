# -*- coding: utf-8 -*-
"""
Created on Mon Jun 16 14:01:33 2025

@author: Isabella Buck

This script compares motif p-values and ranks across enhancer regions and reference overlaps.
It generates summary tables, heatmaps, and visualizes enhancer/overlap regions.
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os
import numpy as np
import seaborn as sns

from Motif_dist import data_dir, file_place, pval, cutoff_qval, perm, mut
from Motif_visualizer import overlaps

# Output directories for results and plots
output_dir = rf"C:\Users\izzye\OneDrive - Johns Hopkins\Documents\BRIGHT\motif{mut}\full\{file_place}\csv"
comp_dir = rf"C:\Users\izzye\OneDrive - Johns Hopkins\Documents\BRIGHT\motif{mut}\full\images\overlap"
heat_dir = rf"C:\Users\izzye\OneDrive - Johns Hopkins\Documents\BRIGHT\motif{mut}\full\images\heatmaps"
os.makedirs(comp_dir, exist_ok=True)
os.makedirs(output_dir, exist_ok=True)
os.makedirs(heat_dir, exist_ok=True)

# List of merged BED files to process
files = [f"{perm}st10_merged_sorted.bed",
        f"{perm}st10R_merged_sorted.bed",
        f"{perm}IntP2A_merged_sorted.bed",
        f"{perm}IntP2B_merged_sorted.bed",
        f"{perm}pCRM_merged_sorted.bed",
        f"{perm}Peak9820_merged_sorted.bed",
        f"{perm}Peak9821_merged_sorted.bed",
        f"{perm}Peak9823_merged_sorted.bed",
        f"{perm}Peak9817_merged_sorted.bed",
        f"{perm}AEL4_merged_sorted.bed",
        f"{perm}5P3_merged_sorted.bed"
        ]

def comparison(files):
    """
    Summarizes p-values and ranks for enhancer overlaps, outputs CSV tables and heatmap.
    """
    results = []
    for file in files:
        name = file.split('_')[0]
        file_path = os.path.join(data_dir, file)
        df = pd.read_csv(file_path, sep="\t", header=None, names=["chrom", "start", "end", "pvalue"])
        df["rank"] = df["pvalue"].rank(method="dense", ascending=True)
        label = file.replace("_merged_sorted.bed", "")
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
    output_df = pd.DataFrame(results, columns=["Enhancer", "Sequence", "p-value", "Rank"])
   
    # Group duplicate entries & store lists of values per Enhancer+Sequence
    grouped_df = output_df.groupby(["Enhancer", "Sequence"]).agg(list).reset_index()
    
    # Enforce order for Enhancer & Sequence
    grouped_df["Enhancer"] = pd.Categorical(grouped_df["Enhancer"], categories=[f.replace("_merged_sorted.bed", "") for f in files], ordered=True)
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
    pivot_pval.to_csv(os.path.join(output_dir, f"full-{perm}qval-{cutoff_qval}_pvalues-{pval}_output.csv"))
    pivot_rank.to_csv(os.path.join(output_dir, f"full-{perm}qval-{cutoff_qval}_ranks-{pval}_output.csv"))
        
    # Prepare numeric p-values for heatmap coloring
    pval_numeric = grouped_df.copy()
    pval_numeric["p-value"] = pval_numeric["p-value"].apply(
        lambda x: min(map(float, str(x).split(','))) if x else np.nan
    )
    pivot_pval_numeric = pval_numeric.pivot(index="Enhancer", columns="Sequence", values="p-value")
        
    # Ensure all rows/columns are represented
    pivot_pval_numeric = pivot_pval_numeric.reindex(
        index=[f.replace("_merged_sorted.bed", "") for f in files],
        columns=list(overlaps.keys())
    )
    print("pivot_pval_numeric:")
    print(pivot_pval_numeric)
        
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
    plt.title(f"P-Value-{pval} {perm} Heatmap")
    plt.tight_layout()
    combo_path = os.path.join(heat_dir, f"full-{perm}{file_place}_Qval({cutoff_qval})_Pval({pval})_heatmap.png")
    plt.savefig(combo_path, dpi=300)
    print(f"Combined heatmap saved to: {combo_path}")
    plt.show()

def compare(files):
    """
    Visualizes enhancer regions and reference overlaps as colored boxes.
    """
    fig, ax = plt.subplots(figsize=(12, 4))
    # Plot enhancer regions as blue horizontal boxes
    y_position = 0
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
    ax.set_xlim(13056267, 13086627)
    ax.set_ylim(-1, y_position + len(overlaps) + 1)
    y_positions = [y + 0.25 for y in range(len(files) + len(overlaps))]
    ax.set_yticks(y_positions)
    ax.set_yticklabels([bed_file.replace("_merged_sorted.bed", "") for bed_file in files] + list(overlaps.keys()))
    ax.set_xlabel("Genomic Coordinates")
    ax.set_title(f"Enhancer vs. Overlap Regions-{pval} {perm}")
    plot_path = os.path.join(comp_dir, f"full-{perm}{file_place}_Qval({cutoff_qval})_Pval({pval})_enhancer_overlap_plot_.png")
    plt.savefig(plot_path, dpi=300, bbox_inches="tight")
    print(f"Plot saved to: {plot_path}")
    plt.show()

def annotate_merged_bed_with_motifs(
    merged_bed_file,
    fimo_file,
    output_file,
    cutoff_metric='q-value',
    cutoff_qval=0.05
):
    """
    Annotates merged BED regions with motif hits from a FIMO file.
    Outputs a TSV file listing motifs overlapping each region.
    """
    import pandas as pd
    # Read merged BED file
    try:
        bed_df = pd.read_csv(
            merged_bed_file, sep="\t", header=None,
            names=["chrom", "start", "end", "p-value"]
        )
    except Exception as e:
        print(f"Could not read {merged_bed_file}: {e}")
        return
    if bed_df.empty:
        print(f"Warning: {merged_bed_file} is empty, no regions to annotate.")
        pd.DataFrame(columns=[
            "Chromosome", "Start", "End", "Region_p-value", "Motif_ID", "Motif_Name"
        ]).to_csv(output_file, sep="\t", index=False)
        print(f"Motif-annotated merged BED saved to {output_file}")
        return
            
    # Read FIMO file
    try:
        fimo_df = pd.read_csv(fimo_file, sep='\t', comment='#')
    except Exception as e:
        print(f"Could not read {fimo_file}: {e}")
        return
    if cutoff_metric in fimo_df.columns:
        fimo_df = fimo_df[fimo_df[cutoff_metric] < cutoff_qval]
            
    # Extract chrom, start, stop from sequence_name if not present
    if 'sequence_name' in fimo_df.columns and not {'chrom', 'start', 'stop'}.issubset(fimo_df.columns):
        fimo_df[['chrom', 'start', 'stop']] = fimo_df['sequence_name'].str.split('_', expand=True)
        fimo_df['start'] = fimo_df['start'].astype(int)
        fimo_df['stop'] = fimo_df['stop'].astype(int)
    if not {'chrom', 'start', 'stop', 'motif_id', 'motif_alt_id'}.issubset(fimo_df.columns):
        print(f"Missing columns in {fimo_file}")
        pd.DataFrame(columns=[
            "Chromosome", "Start", "End", "Region_p-value", "Motif_ID", "Motif_Name"
        ]).to_csv(output_file, sep="\t", index=False)
        return
            
    # Ensure types are correct
    bed_df['start'] = bed_df['start'].astype(int)
    bed_df['end'] = bed_df['end'].astype(int)
    fimo_df['start'] = fimo_df['start'].astype(int)
    fimo_df['stop'] = fimo_df['stop'].astype(int)
    rows = []
    for idx, row in bed_df.iterrows():
        chrom = str(row['chrom'])
        start = row['start']
        end = row['end']
            
        # Find motif hits on the same chrom that overlap this region
        hits = fimo_df[
            (fimo_df['chrom'].astype(str) == chrom) &
            (fimo_df['start'] < end) &
            (fimo_df['stop'] > start)
        ]
        motif_ids = hits['motif_id'].astype(str).tolist()
        motif_names = hits['motif_alt_id'].astype(str).tolist()
        if not motif_ids:
            rows.append({
                "Chromosome": chrom,
                "Start": start,
                "End": end,
                "Region_p-value": row["p-value"],
                "Motif_ID": "",
                "Motif_Name": ""
            })
        else:
            for m_id, m_name in zip(motif_ids, motif_names):
                rows.append({
                    "Chromosome": chrom,
                    "Start": start,
                    "End": end,
                    "Region_p-value": row["p-value"],
                    "Motif_ID": m_id,
                    "Motif_Name": m_name
                })
    pd.DataFrame(rows).to_csv(output_file, sep="\t", index=False)
    print(f"Motif-annotated merged BED saved to {output_file}")

def compare_window_and_enhancer_motifs_no_coords(
    merged_with_motifs_file,
    enhancer_fimo_file,
    output_file,
    cutoff_metric='q-value',
    cutoff_qval=0.05
):
    """
    Compares motifs found in merged windows to those in enhancer regions.
    Outputs a CSV file listing shared and unique motifs.
    """
    import pandas as pd
    # Load merged windows with motifs
    df = pd.read_csv(merged_with_motifs_file, sep="\t")
    # Collect motif_id and motif_alt_id pairs in all merged windows
    window_motifs = set()
    for m_id, m_alt in zip(df["Motif_ID"], df["Motif_Name"]):
        if pd.notna(m_id) and m_id:
            for id_val, alt_val in zip(m_id.split(','), m_alt.split(',')):
                window_motifs.add((id_val.strip(), alt_val.strip()))
    # Load enhancer FIMO file and filter
    fimo = pd.read_csv(enhancer_fimo_file, sep="\t", comment="#")
    if cutoff_metric in fimo.columns:
        fimo = fimo[fimo[cutoff_metric] < cutoff_qval]
    enhancer_motifs = set(zip(fimo["motif_id"].astype(str), fimo["motif_alt_id"].astype(str)))
    # Compare
    shared = window_motifs & enhancer_motifs
    only_in_windows = window_motifs - enhancer_motifs
    only_in_enhancer = enhancer_motifs - window_motifs
    # Pad lists for DataFrame
    max_len = max(len(shared), len(only_in_windows), len(only_in_enhancer))
    def pad(lst):
        return list(lst) + [("", "")] * (max_len - len(lst))
    shared_padded = pad(shared)
    only_in_windows_padded = pad(only_in_windows)
    only_in_enhancer_padded = pad(only_in_enhancer)
    # Build DataFrame
    df_out = pd.DataFrame({
        "Shared_Motif_ID": [x[0] for x in shared_padded],
        "Shared_Motif_Alt_ID": [x[1] for x in shared_padded],
        "Only_in_Windows_Motif_ID": [x[0] for x in only_in_windows_padded],
        "Only_in_Windows_Motif_Alt_ID": [x[1] for x in only_in_windows_padded],
        "Only_in_Enhancer_Motif_ID": [x[0] for x in only_in_enhancer_padded],
        "Only_in_Enhancer_Motif_Alt_ID": [x[1] for x in only_in_enhancer_padded],
    })
    df_out.to_csv(output_file, index=False)
    print(f"Motif comparison saved to {output_file}")

# Directories for window and enhancer motif files
wind_dir = rf"C:\Users\izzye\OneDrive - Johns Hopkins\Documents\BRIGHT\motif{mut}\full\windows"
out_dir = rf"C:\Users\izzye\OneDrive - Johns Hopkins\Documents\BRIGHT\motif{mut}\full\{file_place}\comparison"
enhancer_fimo_dir = rf"C:\Users\izzye\OneDrive - Johns Hopkins\Documents\BRIGHT\motif{mut}\full\enhancers"
os.makedirs(wind_dir, exist_ok=True)
os.makedirs(out_dir, exist_ok=True)

# OPTIONAL
# Annotate merged/sorted BED files with motifs from corresponding FIMO files
for bedfile in files:
    base = bedfile.replace("_merged_sorted.bed", "")
    # Remove perm prefix if present
    if base.startswith(perm):
        window_base = base[len(perm):]
    else:
        window_base = base
    fimo_file = os.path.join(wind_dir, f"{window_base}_window_fimo.tsv")
    output_file = os.path.join(out_dir, f"full-{base}-pval={pval}-qval={cutoff_qval}_motifs.tsv")
    merged_bed_path = os.path.join(data_dir, bedfile)
    annotate_merged_bed_with_motifs(
        merged_bed_file=merged_bed_path,
        fimo_file=fimo_file,
        output_file=output_file,
        cutoff_metric='q-value',
        cutoff_qval=cutoff_qval
    )

# Compare motifs in windows and enhancers for selected bases
for base in ["st10", "st10R", "IntP2A", "IntP2B", "pCRM"]:
    merged_with_motifs = os.path.join(out_dir, f"full-{perm}{base}-pval={pval}-qval={cutoff_qval}_motifs.tsv")
    enhancer_fimo = os.path.join(enhancer_fimo_dir, f"{base}_fimo.tsv")
    output_file = os.path.join(out_dir, f"full-{base}-pval={pval}-qval={cutoff_qval}_window_vs_enhancer_motif_comparison.csv")
    compare_window_and_enhancer_motifs_no_coords(
        merged_with_motifs_file=merged_with_motifs,
        enhancer_fimo_file=enhancer_fimo,
        output_file=output_file,
        cutoff_metric='q-value',
        cutoff_qval=cutoff_qval
    )

if __name__ == "__main__":
    compare(files)
    comparison(files)
