# -*- coding: utf-8 -*-
"""
Created on Tue Jun 10 11:36:29 2025

@author: Isabella Buck

This script sorts and merges motif BED files, retaining the minimum p-value for overlapping intervals.
It outputs merged BED files for further comparison.
"""

import sys
import os
import pandas as pd
from Motif_visualizer import data_dir, perm

# List of motif BED files to process
files = [
    f"{perm}st10.bed", 
    f"{perm}st10R.bed", 
    f"{perm}IntP2A.bed",
    f"{perm}IntP2B.bed", 
    f"{perm}pCRM.bed",
    f"{perm}Peak9820.bed",
    f"{perm}Peak9821.bed",
    f"{perm}Peak9823.bed",
    f"{perm}Peak9817.bed",
    f"{perm}AEL4.bed",
    f"{perm}5P3.bed"
]

def sort_and_merge_bed_file(file, data_dir):
    """
    Sorts a BED file by chromosome and start, merges overlapping intervals,
    and retains the minimum p-value for each merged region.
    Saves the result as a new BED file.
    """
    input_path = os.path.join(data_dir, file)
    base = file.split(".bed")[0]
    output_file = os.path.join(data_dir, f"{base}_merged_sorted.bed")
    try:
        # Read BED file (expects columns: chrom, start, end, name, p-value, extra)
        df = pd.read_csv(input_path, sep="\t", header=None, 
                         names=["chrom", "start", "end", "name", "p-value", "extra"])

        # Convert columns to correct types
        df["start"] = df["start"].astype(int)
        df["end"] = df["end"].astype(int)
        df["p-value"] = df["p-value"].astype(float)
        
        # Sort by chromosome & start position
        df.sort_values(by=["chrom", "start"], inplace=True)
        
        # Merge overlapping intervals while retaining min p-value
        merged_intervals = []
        prev_chrom, prev_start, prev_end, prev_pval = None, None, None, None
        for _, row in df.iterrows():
            chrom, start, end, pval = row["chrom"], row["start"], row["end"], row["p-value"]
            
            # If current interval overlaps previous, merge and keep min p-value
            if prev_chrom == chrom and start <= prev_end:
                prev_end = max(prev_end, end)
                prev_pval = min(prev_pval, pval)
            else:
                if prev_chrom is not None:
                    merged_intervals.append([prev_chrom, prev_start, prev_end, prev_pval])
                prev_chrom, prev_start, prev_end, prev_pval = chrom, start, end, pval
       
        # Add the last interval
        if prev_chrom is not None:
            merged_intervals.append([prev_chrom, prev_start, prev_end, prev_pval])
        
        # Save merged BED file (without the name/extra columns)
        merged_df = pd.DataFrame(merged_intervals, columns=["chrom", "start", "end", "p-value"])
        merged_df.to_csv(output_file, sep="\t", index=False, header=False)
        print(f"Processed {file} -> {output_file} (Sorted & Merged, keeping min p-value)")
    except FileNotFoundError:
        print(f"Error: {file} not found", file=sys.stderr)
    except Exception as e:
        print(f"Error processing {file}: {e}", file=sys.stderr)

if __name__ == "__main__":
    # Process each BED file in the list
    for file in files:
        sort_and_merge_bed_file(file, data_dir)
    
    print("finished 2")