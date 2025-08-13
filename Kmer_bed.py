# -*- coding: utf-8 -*-
"""
Created on Tue Jun 10 12:00:00 2025

@author: Isabella Buck

This script concatenates, sorts, and merges BED files, retaining minimum p-values for overlapping intervals.

To use: Configure the parameters directly in the script below, then run with F5 in Spyder.
"""

import sys
import os
import pandas as pd
import subprocess

from Kmer_dist import k, distance, enhancers

# Direct parameter configuration - edit these values as needed
data_dir = rf"C:\Users\izzye\OneDrive - Johns Hopkins\Documents\BRIGHT\seq\{distance}\k={k}"

files = []

# Load bed files
for enhancer in enhancers:
    files.extend([enhancer + ".bed", enhancer + "_rev.bed"])

def concatenate_sort_merge_bed_files(file_list, data_dir):
    """
    Concatenates matching BED file pairs, sorts them, and merges overlapping intervals while retaining min p-values.
    """
    grouped_files = {}

    # Group matching file pairs by prefix
    for file in file_list:
        base_name = file.split("_rev")[0] if "_rev" in file else file.split(".bed")[0]
        grouped_files.setdefault(base_name, []).append(file)
    
    for base, pair in grouped_files.items():
        print(base)
        if len(pair) == 2:  # Only process pairs
            output_file = os.path.join(data_dir, f"{base}_merged_sorted.bed")
            try:
                # Read BED files without identifier column
                df_list = [pd.read_csv(os.path.join(data_dir, f), sep="\t", header=None, 
                                       names=["chrom", "start", "end", "name", "p-value", "extra"])
                           [["chrom", "start", "end", "p-value"]] 
                           for f in pair]  # Dropping the "name" column
                
                merged_df = pd.concat(df_list, ignore_index=True)

                # Convert columns to correct types
                merged_df["start"] = merged_df["start"].astype(int)
                merged_df["end"] = merged_df["end"].astype(int)
                merged_df["p-value"] = merged_df["p-value"].astype(float)

                # Sort by chromosome & start position
                merged_df.sort_values(by=["chrom", "start"], inplace=True)

                # Merge overlapping intervals while retaining min p-value
                merged_intervals = []
                prev_chrom, prev_start, prev_end, prev_pval = None, None, None, None

                for _, row in merged_df.iterrows():
                    chrom, start, end, pval = row["chrom"], row["start"], row["end"], row["p-value"]
                    
                    if prev_chrom == chrom and start <= prev_end:  
                        prev_end = max(prev_end, end)
                        prev_pval = min(prev_pval, pval)
                        
                    else:
                        if prev_chrom is not None:
                            merged_intervals.append([prev_chrom, prev_start, prev_end, prev_pval])
                        prev_chrom, prev_start, prev_end, prev_pval = chrom, start, end, pval

                if prev_chrom is not None:
                    merged_intervals.append([prev_chrom, prev_start, prev_end, prev_pval])

                # Save merged BED file without the name column
                merged_df = pd.DataFrame(merged_intervals, columns=["chrom", "start", "end", "p-value"])
                merged_df.to_csv(output_file, sep="\t", index=False, header=False)
                print(f"Processed {pair} -> {output_file} (Sorted & Merged, keeping min p-value)")
                
            except FileNotFoundError:
                print(f"Error: One or more files in {pair} not found", file=sys.stderr)
            except Exception as e:
                print(f"Error processing {pair}: {e}", file=sys.stderr)

if __name__ == "__main__":
    concatenate_sort_merge_bed_files(files, data_dir)
    results_comp = subprocess.run(["python", "Kmer_comparison.py"], capture_output=True, text=True, check=True)
    
    print(results_comp.stdout)
    print("finished 2")