# -*- coding: utf-8 -*-
"""
Created on Tue Jun 10 11:36:29 2025

@author: Isabella Buck

This script visualizes k-mer distance results from BedGraph files and computes empirical p-values.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
import subprocess

from Kmer_dist import k, distance, pval

# Modify these values as needed
data_dir = rf"C:\Users\izzye\OneDrive - Johns Hopkins\Documents\BRIGHT\seq\{distance}\k={k}"  # Directory containing BedGraph files
show_plots = False    # Set to True to show plots interactively, False to close them automatically

print(data_dir)

# Load BedGraph file
bedgraph_files = [
    "BedResultsWindow_st10.bedgraph",
    "BedResultsRandom_st10.bedgraph",
    "BedResultsWindow_st10_rev.bedgraph",
    "BedResultsRandom_st10_rev.bedgraph",
    "BedResultsWindow_st10R.bedgraph",
    "BedResultsRandom_st10R.bedgraph",
    "BedResultsWindow_st10R_rev.bedgraph",
    "BedResultsRandom_st10R_rev.bedgraph",
    "BedResultsWindow_IntP2A.bedgraph",
    "BedResultsRandom_IntP2A.bedgraph",
    "BedResultsWindow_IntP2A_rev.bedgraph",
    "BedResultsRandom_IntP2A_rev.bedgraph",
    "BedResultsWindow_IntP2B.bedgraph",
    "BedResultsRandom_IntP2B.bedgraph",
    "BedResultsWindow_IntP2B_rev.bedgraph",
    "BedResultsRandom_IntP2B_rev.bedgraph",
    "BedResultsWindow_pCRM.bedgraph",
    "BedResultsRandom_pCRM.bedgraph",
    "BedResultsWindow_pCRM_rev.bedgraph",
    "BedResultsRandom_pCRM_rev.bedgraph",
    "BedResultsWindow_Peak9820.bedgraph",
    "BedResultsRandom_Peak9820.bedgraph",
    "BedResultsWindow_Peak9820_rev.bedgraph",
    "BedResultsRandom_Peak9820_rev.bedgraph",
    "BedResultsWindow_Peak9821.bedgraph",
    "BedResultsRandom_Peak9821.bedgraph",
    "BedResultsWindow_Peak9821_rev.bedgraph",
    "BedResultsRandom_Peak9821_rev.bedgraph",
    "BedResultsWindow_Peak9823.bedgraph",
    "BedResultsRandom_Peak9823.bedgraph",
    "BedResultsWindow_Peak9823_rev.bedgraph",
    "BedResultsRandom_Peak9823_rev.bedgraph",
    "BedResultsWindow_Peak9817.bedgraph",
    "BedResultsRandom_Peak9817.bedgraph",
    "BedResultsWindow_Peak9817_rev.bedgraph",
    "BedResultsRandom_Peak9817_rev.bedgraph",
    "BedResultsWindow_5P3.bedgraph",
    "BedResultsRandom_5P3.bedgraph",
    "BedResultsWindow_5P3_rev.bedgraph",
    "BedResultsRandom_5P3_rev.bedgraph",
    "BedResultsWindow_AEL4.bedgraph",
    "BedResultsRandom_AEL4.bedgraph",
    "BedResultsWindow_AEL4_rev.bedgraph",
    "BedResultsRandom_AEL4_rev.bedgraph"
]

def visualize(bedgraph_file):
    """
    Visualize k-mer distance distributions for each BedGraph file.
    """
    for file in bedgraph_files:
        file_path = os.path.join(data_dir, file)
        df = pd.read_csv(file_path, sep="\t", header=None, names=["chrom", "start", "end", f"{distance}_distance"])

        # Compute the midpoint (average genomic coordinate)
        df["midpoint"] = (df["start"] + df["end"]) / 2

        # Sort by chromosome and midpoint
        df.sort_values(["chrom", "midpoint"], inplace=True)

        # Trim data of outliers if needed
        q_low = df[f"{distance}_distance"].quantile(0)
        q_high = df[f"{distance}_distance"].quantile(1)
        df_trimmed = df[(df[f"{distance}_distance"] >= q_low) & (df[f"{distance}_distance"] <= q_high)]

        avg_y = df_trimmed[f"{distance}_distance"].median()
        print(f"\n{file} with avg {avg_y}")

        # Extract a name from filename like "rand_IntP2B.fa" -> "IntP2B"
        enhancer_name = file.split("_")[1].replace(".sorted.bedgraph", "").replace(".bedgraph", "")

        # Extract if the file is Random or Window
        types = file.split("_")[0].replace("BedResultsWindow", "Window").replace("BedResultsRandom", "Random").replace("RandMerg", "Random (merged)")

        # Create the scatterplot figure
        plt.figure(figsize=(14, 6))

        # Plot the average distance line
        plt.axhline(y=avg_y, color='black', linestyle='--', linewidth=1, label=f'Average = {avg_y:.3f}')

        # Use seaborn to color different chromosomes
        sns.scatterplot(data=df_trimmed, x="midpoint", y=f"{distance}_distance", hue="chrom", palette="mako", alpha=.9)

        # Set x-axis labels explicitly to genomic coordinate values
        plt.xticks(np.linspace(df_trimmed["midpoint"].min(), df["midpoint"].max(), num=10, dtype=int), rotation=45)

        # Label axes and title
        plt.xlabel("Genomic Midpoint (bp)")
        plt.ylabel(f"{distance} Distance")
        plt.ylim(0)
        plt.title(f"Genomic Distribution of {distance} Distances: {enhancer_name} vs {types}")
        plt.legend(title="Chromosome")
        plt.grid(True)

        # Show the plot
        if show_plots:
            plt.show()
        else:
            plt.close()

        # OPTIONALLY Make a histogram of the scatterplot
        sns.histplot(df_trimmed[f"{distance}_distance"], bins=50)
        plt.xlabel(f"{distance} Distance")
        plt.title(f'Distance Distribution: {enhancer_name} vs {types}')
        plt.show()
        plt.close()

def pvalue(window_file, random_file, label, output_bedfile=None):
    """
    Compute p-values for each window and save to BED file.
    """
    if output_bedfile is None:
        output_bedfile = os.path.join(data_dir, f"{label}.bed")

    count = 0
    bound = 500
    partials = 0
    full = 0
    p_cutoff = pval
    o_cutoff = 50
    name_score = 1

    window_path = os.path.join(data_dir, window_file)
    random_path = os.path.join(data_dir, random_file)

    win_df = pd.read_csv(window_path, sep="\t", header=None, names=["chrom", "start", "end", "dist"])
    rand_df = pd.read_csv(random_path, sep="\t", header=None, names=["chrom", "start", "end", "dist"])

    win_df.sort_values(["dist"], inplace=True)

    bed_lines = []
    for _, row in win_df.iterrows():
        distance = row["dist"]
        start = row["start"]
        end = row["end"]
        chrom = row["chrom"]

        k = (rand_df["dist"] <= distance).sum()
        N = len(rand_df)
        p_raw = k/N
        if k < bound:
            count += 1

        name_score_st = str(name_score)
        name = (f"{label}-{name_score_st}")
        score = round(p_raw, 6)
        strand = "."

        # See which seq overlap
        if 13060273 <= start <= 13060969 and 13060273 <= end <= 13060969:
            overlap_e = ((end - start)/(end-start))*100
            overlap_w = ((end - start)/(13060969 - 13060273))*100
            full += 1
            if p_raw <p_cutoff and overlap_e > o_cutoff:
                print(f"\033[44m{name} Full overlap: E in W: {overlap_e:.4}% W in E: {overlap_w:.4}% Pvalue: {p_raw}\033[0m Window dist: {distance}  Start: {start} End: {end} Randoms #: {k} \033[0m")
            else:
                print(f"\033[41m{name} Full overlap: E in W: {overlap_e:.4}% W in E: {overlap_w:.4}% Pvalue: {p_raw} Window dist: {distance}  Start: {start} End: {end} Randoms #: {k} \033[0m")
        elif 13060273 <= start <= 13060969 or 13060273 <= end <= 13060969:
            if 13060273 <= start <= 13060969:
                overlap_e = ((13060969 - start) / (end-start))*100
                overlap_w = ((13060969 - start) / (13060969 - 13060273))*100
            else:
                overlap_e = ((end - 13060273) / (end-start))*100
                overlap_w = ((end - 13060273) / (13060969 - 13060273))*100
            partials +=1
            if p_raw <p_cutoff and overlap_e > o_cutoff:
                print(f"\033[44m{name} Partial overlap: E in W: {overlap_e:.4}% W in E: {overlap_w:.4}% Pvalue: {p_raw} Window dist: {distance}  Start: {start} End: {end} Randoms #: {k} \033[0m")
            elif overlap_e > o_cutoff:
                print(f"\033[41m{name} Partial overlap: E in W: {overlap_e:.4}% W in E: {overlap_w:.4}% Pvalue: {p_raw}\033[0m Window dist: {distance}  Start: {start} End: {end} Randoms #: {k}\033[0m")
            else:
                print(f"{name} Partial overlap: E in W: {overlap_e:.4}% W in E: {overlap_w:.4}% Pvalue: {p_raw} Window dist: {distance}  Start: {start} End: {end} Randoms #: {k}")
        elif p_raw < p_cutoff:
            print(f"{name} No overlap Pvalue: {p_raw} Window dist: {distance}  Start: {start} End: {end} Randoms #: {k}")

        if p_raw < p_cutoff:
            bed_lines.append([chrom, start, end, name, score, strand])
            
        # Increase the score on each iteration
        name_score += 1

    bed_df = pd.DataFrame(bed_lines, columns=["chrom", "start", "end", "name", "score", "strand"])
    bed_df.to_csv(output_bedfile, sep="\t", header=False, index=False)
    print(f"Total - Full Overlaps: {full} Partial Overlaps: {partials}")
    print(f"Number of windows closer than {bound} randoms: {count}")
    print("\n")
    
def compare_distributions(window_file, random_file, label):
    """
    Compare distributions of window and random distances, plot, and compute p-values.
    """
    rand_df = pd.read_csv(os.path.join(data_dir, random_file), sep="\t", header=None, names=["chrom", "start", "end", f"{distance}_distance"])

    stats = rand_df[f"{distance}_distance"].describe()

    plt.figure(figsize=(8, 4))
    sns.boxplot(x=rand_df[f"{distance}_distance"], color="lightblue")

    plt.text(stats["max"]-.8, .25, f"Max: {stats['max']:.2f}", color="black")
    plt.text(stats["75%"]+.2, .4, f"Q3: {stats['75%']:.2f}", color="black")
    plt.text(stats["50%"]-.8, .5, f"Median: {stats['50%']:.2f}", color="blue")
    plt.text(stats["25%"]-1.8, .4, f"Q1: {stats['25%']:.2f}", color="black")
    plt.text(stats["min"]-.8, .25, f"Min: {stats['min']:.2f}", color="black")
    plt.text(stats["mean"]-.8, .45, f"Mean: {stats['mean']:.2f}", color="red")

    plt.xlabel(f"{distance} Distance")
    plt.title(f"Boxplot of Random {distance} Distances {label}")
    plt.close()

    pvalue(window_file, random_file, label)
  
def compare_pvals(pvals):
    """
    Barplot comparing p-values for different enhancer windows.
    """
    labels = ["st10", "st10_rev", "st10R", "st10R_rev", "IntP2A", "IntP2A_rev", "IntP2B", "IntP2B_rev", "pCRM", "pCRM_rev"]

    plt.figure(figsize=(12, 6))
    plt.bar(labels, pvals, color='skyblue')
    plt.axhline(.05, color='blue', linestyle='--', label='p = 0.05')
    plt.axhline(.01, color='black', linestyle='--', label='p = 0.15')
    plt.ylabel('pvalue')
    plt.title('Statistical Significance of Enhancer Distance Distributions')
    plt.xticks(rotation=45)
    plt.ylim(0, .06)
    plt.legend()
    plt.tight_layout()
    plt.close()

if __name__ == "__main__":
    # Run comparison for each pairing
    # Compare window file, random file, name
    compare_distributions("BedResultsWindow_st10.bedgraph", "BedResultsRandom_st10.bedgraph", "st10")
    compare_distributions("BedResultsWindow_st10_rev.bedgraph", 'BedResultsRandom_st10_rev.bedgraph', "st10_rev")
    compare_distributions("BedResultsWindow_st10R.bedgraph", "BedResultsRandom_st10R.bedgraph", "st10R")
    compare_distributions("BedResultsWindow_st10R_rev.bedgraph", "BedResultsRandom_st10R_rev.bedgraph", "st10R_rev")
    compare_distributions("BedResultsWindow_IntP2A.bedgraph", "BedResultsRandom_IntP2A.bedgraph", "IntP2A")
    compare_distributions("BedResultsWindow_IntP2A_rev.bedgraph", "BedResultsRandom_IntP2A_rev.bedgraph", "IntP2A_rev")
    compare_distributions("BedResultsWindow_IntP2B.bedgraph", "BedResultsRandom_IntP2B.bedgraph", "IntP2B")
    compare_distributions("BedResultsWindow_IntP2B_rev.bedgraph", "BedResultsRandom_IntP2B_rev.bedgraph", "IntP2B_rev")
    compare_distributions("BedResultsWindow_pCRM.bedgraph", "BedResultsRandom_pCRM.bedgraph", "pCRM")
    compare_distributions("BedResultsWindow_pCRM_rev.bedgraph", "BedResultsRandom_pCRM_rev.bedgraph", "pCRM_rev")
    compare_distributions("BedResultsWindow_Peak9820.bedgraph", "BedResultsRandom_Peak9820.bedgraph", "Peak9820")
    compare_distributions("BedResultsWindow_Peak9820_rev.bedgraph", "BedResultsRandom_Peak9820_rev.bedgraph", "Peak9820_rev")
    compare_distributions("BedResultsWindow_Peak9821.bedgraph", "BedResultsRandom_Peak9821.bedgraph", "Peak9821")
    compare_distributions("BedResultsWindow_Peak9821_rev.bedgraph", "BedResultsRandom_Peak9821_rev.bedgraph", "Peak9821_rev")
    
    visualize(bedgraph_files)
    
    results_bed = subprocess.run(["python", "Kmer_bed.py"])
    print(results_bed.stdout)
    print("finished 1")

# OPTIONALLY
#compare_pvals(pvals)
