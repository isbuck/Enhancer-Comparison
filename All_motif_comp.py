# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 11:27:36 2025

@author: Isabella Buck

All-vs-All Motif Comparison for Windows
Uses FIMO files and directories from Motif_Comp.py

This script computes all-vs-all motif distance matrices (Jaccard or Canberra) for genomic windows,
performs clustering, and visualizes results as heatmaps, dendrograms, and cluster scatterplots.
"""

import pandas as pd
import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform
import matplotlib.patches as mpatches
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

# Choose 'jaccard' or 'canberra'
file_place = "jaccard"  
cutoff_qval = .05

window_dir = r"C:\Users\izzye\OneDrive - Johns Hopkins\Documents\BRIGHT\motif\full\windows"
output_dir = rf"C:\Users\izzye\OneDrive - Johns Hopkins\Documents\BRIGHT\motif\full\{file_place}\full_comp"

os.makedirs(window_dir, exist_ok=True)
os.makedirs(output_dir, exist_ok=True)

enhancers = ['IntP2B', 'IntP2A', 'st10', 'st10R', 'pCRM', 'Peak9820',
             'Peak9821', 'AEL4', '5P3', 'Peak9817', 'Peak9823']

window_files = []
for enhancer in enhancers:
    window_files.append(enhancer + "_window_fimo.tsv")

def extract_motif_sets(fimo_file, cutoff_metric='q-value', cutoff_qval=1):
    """
    Extracts motif sets for each region from a FIMO file, filtering by cutoff.
    Returns a dict: {region: set of motif_ids}
    """
    try:
        df = pd.read_csv(fimo_file, sep='\t', comment='#')
        if cutoff_metric in df.columns:
            df = df[df[cutoff_metric] < cutoff_qval]
        if df.empty or 'sequence_name' not in df.columns or 'motif_id' not in df.columns:
            return {}
        motif_sets = {
            region: set(group['motif_id'])
            for region, group in df.groupby('sequence_name')
        }
        return motif_sets
    except Exception as e:
        print(f"Error reading {fimo_file}: {e}")
        return {}

def extract_motif_vectors(fimo_file, cutoff_metric='q-value', cutoff_qval=1):
    """
    Extracts motif score vectors for each region from a FIMO file, filtering by cutoff.
    Returns a dict: {region: np.array of motif scores}
    """
    try:
        df = pd.read_csv(fimo_file, sep='\t', comment='#')
        if cutoff_metric in df.columns:
            df = df[df[cutoff_metric] < cutoff_qval]
        if df.empty or 'sequence_name' not in df.columns or 'motif_id' not in df.columns or 'score' not in df.columns:
            return {}
        grouped = df.groupby(['sequence_name', 'motif_id'], as_index=False)['score'].sum()
        matrix = grouped.pivot(index='sequence_name', columns='motif_id', values='score').fillna(0)
        matrix = matrix.sort_index().sort_index(axis=1)
        motif_vectors = {region: matrix.loc[region].values for region in matrix.index}
        return motif_vectors
    except Exception as e:
        print(f"Error reading {fimo_file}: {e}")
        return {}

def jaccard_distance(set1, set2):
    """
    Computes Jaccard distance between two sets.
    """
    if not set1 and not set2:
        return 0.0
    return 1 - len(set1 & set2) / len(set1 | set2)

def canberra_distance(vec1, vec2):
    """
    Computes Canberra distance between two vectors.
    """
    from scipy.spatial.distance import canberra
    return canberra(vec1, vec2)

def all_vs_all_comparison(base, metric):
    """
    For a given window base name, computes all-vs-all motif distance matrix,
    saves CSV, plots heatmap, performs clustering, and visualizes clusters.
    """
    fimo_file = os.path.join(window_dir, f"{base}_window_fimo.tsv")
    print(f"Processing {fimo_file} for all-vs-all {metric} comparison...")

    # Compute distance matrix
    # Jaccard distance
    if metric == 'jaccard':
        motif_sets = extract_motif_sets(fimo_file, cutoff_metric='q-value', cutoff_qval=cutoff_qval)
        regions = list(motif_sets.keys())
        matrix = np.zeros((len(regions), len(regions)))
        for i, reg_i in enumerate(regions):
            for j, reg_j in enumerate(regions):
                matrix[i, j] = jaccard_distance(motif_sets[reg_i], motif_sets[reg_j])
    
    # Canberra distance
    else:  
        motif_vectors = extract_motif_vectors(fimo_file, cutoff_metric='q-value', cutoff_qval=cutoff_qval)
        regions = list(motif_vectors.keys())
        matrix = np.zeros((len(regions), len(regions)))
        for i, reg_i in enumerate(regions):
            for j, reg_j in enumerate(regions):
                matrix[i, j] = canberra_distance(motif_vectors[reg_i], motif_vectors[reg_j])

    df_matrix = pd.DataFrame(matrix, index=regions, columns=regions)
    out_csv = os.path.join(output_dir, f"{base}_all_vs_all_{metric}_distance.csv")
    df_matrix.to_csv(out_csv)
    print(f"Saved all-vs-all {metric} distance matrix to {out_csv}")

    # Plot heatmap
    plt.figure(figsize=(10, 8))
    cbar_label = f"{metric.capitalize()} Distance"
    sns.heatmap(df_matrix, cmap="viridis", cbar=True, cbar_kws={'label': cbar_label})

    # Extract start coordinates from region names for axis ticks
    coords = []
    for r in regions:
        coords.append(int(r.split('_')[1])) 

    # Choose evenly spaced ticks for clarity
    n_ticks = 6
    tick_indices = np.linspace(0, len(coords)-1, n_ticks, dtype=int)
    tick_labels = [f"{coords[i]:,}" for i in tick_indices]

    plt.xticks(tick_indices, tick_labels, rotation=45)
    plt.yticks(tick_indices, tick_labels, rotation=0)
    plt.xlabel("Genomic coordinate (window start)")
    plt.ylabel("Genomic coordinate (window start)")
    plt.title(f"{base} All-vs-All {metric.capitalize()} Motif Distance")
    plt.tight_layout()
    heatmap_png = os.path.join(output_dir, f"{base}_all_vs_all_{metric}_heatmap.png")
    plt.savefig(heatmap_png, dpi=300)
    plt.close()
    print(f"Saved heatmap to {heatmap_png}")

    # Clustering and dendrogram
    try:
        condensed = squareform(df_matrix.values)
        Z = linkage(condensed, method='average')
        plt.figure(figsize=(12, 5))
        dendrogram(Z, labels=[str(c) for c in coords])
        plt.xlabel("Genomic coordinate (window start)")
        plt.ylabel("Distance")
        plt.title(f"{base} Window Motif Clustering Dendrogram")

        # Highlight overlap regions on dendrogram labels
        overlaps = {
            "sim_A1.0": (13056591, 13057759),
            "sim_st10": (13060273, 13060969),
            "sim_D2.1": (13064360, 13066277),
            "sim_E2.3": (13066277, 13067711),
            "EV_overlap": (13067711, 13068610),
            "sim_VT040842": (13068610, 13069801),
            "sim_1.6MLE": (13070722, 13072353)
        }
        ax = plt.gca()
        xlabs = ax.get_xmajorticklabels()
        for lbl in xlabs:
            try:
                coord = int(lbl.get_text())
                for start, end in overlaps.values():
                    if start <= coord <= end:
                        lbl.set_color('deepskyblue')
                        break
            except Exception:
                pass
        plt.tight_layout()
        dendro_png = os.path.join(output_dir, f"{base}_all_vs_all_{metric}_dendrogram.png")
        plt.savefig(dendro_png, dpi=300)
        plt.close()
        print(f"Saved dendrogram to {dendro_png}")
    except Exception as e:
        print(f"Clustering failed for {base}: {e}")

    # K-means clustering for Canberra
    if metric == 'canberra':
        motif_vectors = extract_motif_vectors(fimo_file, cutoff_metric='q-value', cutoff_qval=cutoff_qval)
        regions = list(motif_vectors.keys())
        X = np.array([motif_vectors[r] for r in regions])
        
        # Number of clusters
        k = 4  
        kmeans = KMeans(n_clusters=k, random_state=42).fit(X)
        labels = kmeans.labels_
        region_clusters = dict(zip(regions, labels))
        print(f"K-means clusters for {base}:")
        print(region_clusters)
        cluster_df = pd.DataFrame({'region': regions, 'cluster': labels})
        cluster_csv = os.path.join(output_dir, f"{base}_k={k}kmeans_clusters.csv")
        cluster_df.to_csv(cluster_csv, index=False)
        print(f"Saved k-means cluster assignments to {cluster_csv}")

        # Plot clusters along genome
        coords = [int(r.split('_')[1]) for r in regions]
        plt.figure(figsize=(12, 3))
        plt.scatter(coords, labels, c=labels, cmap='tab10', s=60)
        plt.xlabel("Genomic coordinate (window start)")
        plt.ylabel("K-means Cluster")
        plt.title(f"{base} Canberra k={k}K-means Window Clusters Along Genome")
        plt.yticks(range(k))
        plt.tight_layout()
        cluster_png = os.path.join(output_dir, f"{base}_canberra_k={k}kmeans_clusters_along_genome.png")
        plt.savefig(cluster_png, dpi=300)
        plt.close()
        print(f"Saved Canberra k-means cluster scatterplot to {cluster_png}")

    # Jaccard clustering and visualization
    if metric == 'jaccard':
        try:
            condensed = squareform(df_matrix.values)
            Z = linkage(condensed, method='average')
            
            # Number of clusters
            k = 5
            cluster_labels = fcluster(Z, k, criterion='maxclust')
            cluster_df = pd.DataFrame({'region': regions, 'cluster': cluster_labels})
            cluster_csv = os.path.join(output_dir, f"{base}_k={k}jaccard_clusters.csv")
            cluster_df.to_csv(cluster_csv, index=False)
            print(f"Saved Jaccard cluster assignments to {cluster_csv}")

            # Plot clusters along genome
            plt.figure(figsize=(12, 3))
            plt.scatter(coords, cluster_labels, c=cluster_labels, cmap='tab10', s=60)
            plt.xlabel("Genomic coordinate (window start)")
            plt.ylabel("Cluster")
            plt.title(f"{base} k={k}Jaccard Window Clusters Along Genome")
            plt.yticks(range(1, k+1))
            plt.tight_layout()
            
            # Add rectangles for specific genomic regions
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
            cluster_png = os.path.join(output_dir, f"{base}_k={k}jaccard_clusters_along_genome.png")
            plt.savefig(cluster_png, dpi=300)
            plt.close()

            # Silhouette analysis
            max_k = k
            sil_scores = []
            for test_k in range(2, max_k+1):
                test_labels = fcluster(Z, test_k, criterion='maxclust')
                sil = silhouette_score(df_matrix.values, test_labels, metric='precomputed')
                sil_scores.append(sil)
            plt.figure()
            plt.plot(range(2, max_k+1), sil_scores, marker='o')
            plt.xlabel('Number of clusters (k)')
            plt.ylabel('Silhouette Score')
            plt.title(f'Silhouette Analysis for {base} (Jaccard)')
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f"{base}_jaccardk={k}_silhouette.png"))
            plt.close()

            # Plot average distance by genome position and color by cluster
            avg_distances = []
            for i in range(len(df_matrix)):
                avg = (df_matrix.iloc[i, :].sum() - df_matrix.iloc[i, i]) / (len(df_matrix) - 1)
                avg_distances.append(avg)
            plt.figure(figsize=(12, 4))
            plt.scatter(coords, avg_distances, c=cluster_labels, cmap='tab10', s=60)
            plt.xlabel("Genomic coordinate (window start)")
            plt.ylabel("Average Jaccard Distance to Other Windows")
            plt.title(f"{base} Jaccard Distance by Genome Position (Colored by Cluster)")
            plt.colorbar(label="Cluster")
            plt.tight_layout()
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
            scatter_dist_png = os.path.join(output_dir, f"{base}_jaccardk={k}_avg_distance_by_genome_colored_by_cluster.png")
            plt.savefig(scatter_dist_png, dpi=300)
            plt.close()

        except Exception as e:
            print(f"Jaccard clustering failed for {base}: {e}")

if __name__ == "__main__":
    metric = file_place  # 'jaccard' or 'canberra'
    for file in window_files:
        base = file.replace("_window_fimo.tsv", "")
        all_vs_all_comparison(base, metric=file_place)