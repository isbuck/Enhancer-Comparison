import pandas as pd
import os
from scipy.stats import kendalltau
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.patches as mpatches
import seaborn as sns
import numpy as np

# Configuration parameters
cut = 2
cutoff_qval = .05
pval = .01
kd_cut = .05
ful_mid = "full"
metric = "jaccard"
coords = (13056267, 13086627)

def extract_ordered_motifs(fimo_file, cutoff_metric='q-value', cutoff_qval=1):
    """
    Extract and order motifs from a FIMO output file based on genomic position.
    """
    try:
        df = pd.read_csv(fimo_file, sep='\t', comment='#')
        if df.empty: 
            print(f"XX {fimo_file} is empty. No motif hits.")
            return []
        if cutoff_metric in df.columns:
            df = df[df[cutoff_metric] < cutoff_qval]
        print(f"{fimo_file}: {len(df)} motif hits after q-value < {cutoff_qval}")
        
        def get_abs_start(row):
            try:
                chrom, region_start, region_stop = row['sequence_name'].split('_')
                region_start = int(region_start)
                abs_start = region_start + int(row['start'])
                return abs_start
            except Exception:
                return None
        
        df['abs_start'] = df.apply(get_abs_start, axis=1)
        df = df.sort_values(by='abs_start')
        return list(df['motif_id'])
    except pd.errors.EmptyDataError:
        print(f"XX {fimo_file} is empty or has no columns.")
        return []
    
def extract_ordered_motifs_from_df(df, cutoff_metric='q-value', cutoff_qval=1):
    """
    Extract and order motifs from a pandas DataFrame containing FIMO results.
    """
    if cutoff_metric in df.columns:
        df = df[df[cutoff_metric] < cutoff_qval]
    df = df.copy()
    
    def get_abs_start(row):
        try:
            chrom, region_start, region_stop = row['sequence_name'].split('_')
            region_start = int(region_start)
            abs_start = region_start + int(row['start'])
            return abs_start
        except Exception:
            return None
    
    df['abs_start'] = df.apply(get_abs_start, axis=1)
    df = df.sort_values(by='abs_start')
    return list(df['motif_id'])

def motif_order_kd(list1, list2, min_shared=None):
    """
    Compute Kendall's tau correlation between two ordered motif lists.
    """
    if min_shared is None:
        min_shared = cut
    
    # Find shared motifs
    shared = [m for m in list1 if m in list2]
    shared = list(dict.fromkeys(shared))  # Remove duplicates, preserve order
    print(f"Shared motifs: {len(shared)}")
    
    if len(shared) < min_shared:
        return None  # Not enough motifs for ranking
    
    ranks1 = [list1.index(m) for m in shared]
    ranks2 = [list2.index(m) for m in shared]
    tau, _ = kendalltau(ranks1, ranks2)
    return tau

def compare_enhancer_to_windows(enhancer_file, window_file, enhancer_dir, window_dir, output_file, 
                               cutoff_metric='q-value', cutoff_qval=0.05):
    """
    Compare motif ordering between an enhancer and sliding window regions.
    """
    enhancer_path = os.path.join(enhancer_dir, enhancer_file)
    enhancer_motifs = extract_ordered_motifs(enhancer_path, cutoff_metric, cutoff_qval)
    
    window_path = os.path.join(window_dir, window_file)
    df_win = pd.read_csv(window_path, sep='\t', comment='#')
    results = []
    
    for region_id, group in df_win.groupby("sequence_name"):
        window_motifs = extract_ordered_motifs_from_df(group, cutoff_metric, cutoff_qval)
        shared = [m for m in enhancer_motifs if m in window_motifs]
        shared = list(dict.fromkeys(shared))  # Remove duplicates, preserve order
        kd_score = motif_order_kd(enhancer_motifs, window_motifs)
        motif_hit_count = len(shared)  # Total motif hits in this window region
        weighted_kd = kd_score * motif_hit_count if kd_score is not None else None
        results.append([
            region_id,
            kd_score,
            motif_hit_count,
            weighted_kd
        ])
    
    df_out = pd.DataFrame(results, columns=[
        'Window_Region', 'Kd_score', 'Num_Shared_Motifs', 'Weighted_Kd'
    ])
    df_out.to_csv(output_file, sep='\t', index=False)
    print(f"Kd scores saved to {output_file}")

def annotate_bed_with_kd(bed_file, kd_file, output_file, min_pval=1e-4):
    """
    Annotate BED file regions with Kendall's tau correlation scores.
    """
    bed_df = pd.read_csv(bed_file, sep="\t", header=None, names=["chrom", "start", "end", "name", "pvalue", "extra"])
    kd_df = pd.read_csv(kd_file, sep="\t")

    if "chrom" not in kd_df.columns:
        kd_df[["chrom", "start", "end"]] = kd_df["Window_Region"].str.split("_", expand=True)
        kd_df["start"] = kd_df["start"].astype(int)
        kd_df["end"] = kd_df["end"].astype(int)
    
    print(bed_df["chrom"])
    print(kd_df["chrom"])
    bed_df["chrom"] = bed_df["chrom"].astype(str)
    kd_df["chrom"] = kd_df["chrom"].astype(str)
    bed_df["start"] = bed_df["start"].astype(int)
    kd_df["start"] = kd_df["start"].astype(int)
    bed_df["end"] = bed_df["end"].astype(int)
    kd_df["end"] = kd_df["end"].astype(int)
    print(bed_df["chrom"])
    print(kd_df["chrom"])
    
    merged = pd.merge(
        bed_df,
        kd_df[["chrom", "start", "end", "Kd_score", "Weighted_Kd", "Kd_score_norm", "Weighted_Kd_norm"]],
        on=["chrom", "start", "end"],
        how="left"
    )    
    merged["pvalue"] = pd.to_numeric(merged["pvalue"], errors="coerce")
    merged["pvalue"] = merged["pvalue"].clip(lower=min_pval)
    merged["score"] = merged["Weighted_Kd"] * (-np.log10(merged["pvalue"]))

    # Normalize scores
    score_min = merged["score"].min()
    score_max = merged["score"].max()
    if score_max > score_min:
        merged["score_norm"] = (merged["score"] - score_min) / (score_max - score_min)
    else:
        merged["score_norm"] = merged["score"]

    # Only keep selected columns in the output
    columns_to_keep = ["chrom", "start", "end", "pvalue", "Kd_score", "Weighted_Kd", "Kd_score_norm", "Weighted_Kd_norm","score", "score_norm"]
    merged[columns_to_keep].to_csv(output_file, sep="\t", index=False, header=True)
    print(f"Annotated BED file saved to {output_file}")
    
    # Print and save top 5 by score
    merged["score"] = pd.to_numeric(merged["score"], errors="coerce")
    top5 = merged[merged["score"].notna()].sort_values("score", ascending=False).head(5)
    print(f"\nTop 5 windows by score for {os.path.basename(output_file)}:")
    print(top5[["chrom", "start", "end", "score", "score_norm"]])

    txt_output = output_file.replace(".bed", "_top5.txt").replace(".tsv", "_top5.txt")
    top5[["chrom", "start", "end", "score"]].to_csv(txt_output, sep="\t", index=False)
    print(f"Top 5 windows written to {txt_output}")

def annotate_merged(bed_file, kd_file, output_file, overlaps=None, plot_out_dir=None, min_pval=1e-4):
    """
    Annotate merged BED regions with best Kd scores from overlapping windows.
    """
    bed_df = pd.read_csv(bed_file, sep="\t", header=None, names=["chrom", "start", "end", "pvalue"])
    kd_df = pd.read_csv(kd_file, sep="\t")

    bed_df["chrom"] = bed_df["chrom"].astype(str)
    kd_df["chrom"] = kd_df["chrom"].astype(str)
    bed_df["start"] = bed_df["start"].astype(int)
    kd_df["start"] = kd_df["start"].astype(int)
    bed_df["end"] = bed_df["end"].astype(int)
    kd_df["end"] = kd_df["end"].astype(int)

    if "chrom" not in kd_df.columns:
        kd_df[["chrom", "start", "end"]] = kd_df["Window_Region"].str.split("_", expand=True)
        kd_df["start"] = kd_df["start"].astype(int)
        kd_df["end"] = kd_df["end"].astype(int)

    # For each merged region, find overlapping windows and select the best score
    best_rows = []
    for _, row in bed_df.iterrows():
        chrom, start, end, pvalue = row["chrom"], row["start"], row["end"], row["pvalue"]
        overlaps_df = kd_df[
            (kd_df["chrom"] == chrom) &
            (kd_df["start"] < end) &
            (kd_df["end"] > start)
        ]
        if not overlaps_df.empty:
            valid = overlaps_df[overlaps_df["Weighted_Kd"].notna()]
            if not valid.empty:
                best = valid.loc[valid["Weighted_Kd"].idxmax()]
                weighted_kd = best["Weighted_Kd"]
                try:
                    if pvalue <= 0:
                        pvalue = min_pval
                    logpval = -np.log10(pvalue)
                except Exception:
                    logpval = 0
                try:
                    score = float(weighted_kd) * logpval if weighted_kd != "NA" and logpval != 0 else "NA"
                except Exception:
                    score = "NA"
                best_rows.append([
                    chrom, start, end, pvalue,
                    best["Kd_score"], weighted_kd, best["Kd_score_norm"], best["Weighted_Kd_norm"], score
                ])
            else:
                best_rows.append([chrom, start, end, pvalue, "NA", "NA", "NA", "NA", "NA"])
        else:
            best_rows.append([chrom, start, end, pvalue, "NA", "NA", "NA", "NA", "NA"])

    columns_to_keep = ["chrom", "start", "end", "pvalue", "Kd_score", "Weighted_Kd", "Kd_score_norm", "Weighted_Kd_norm", "logpval_x_weightedkd"]
    df = pd.DataFrame(best_rows, columns=columns_to_keep)

    excel_output = output_file.replace(".bed", ".xlsx").replace(".tsv", ".xlsx")
    df.to_excel(excel_output, index=False)
    df.to_csv(output_file, sep="\t", index=False, header=True)
    print(f"Annotated merged BED file saved to {output_file}")

def plot_kd_boxes_and_heatmap(annotated_bed_files, overlaps, out_dir, value_col="Weighted_Kd", 
                              kd_percentile=None, title_suffix=""):
    """
    Create comprehensive visualization of Kd scores with box plots and heatmaps.
    """
    if kd_percentile is None:
        kd_percentile = kd_cut
    
    os.makedirs(out_dir, exist_ok=True)
    
    # Box Plot
    fig, ax = plt.subplots(figsize=(12, 4))
    y_position = 0
    labels = []

    for bed_file in annotated_bed_files:
        # Get base name to find the corresponding kd file
        base = os.path.basename(bed_file).replace("_merged_sorted.bed", "").replace("_annotated.bed","")
        print(base)
        kd_file = os.path.join(
            os.path.dirname(bed_file), f"{base}_kd_scores_with_norm.tsv"
        )
        print(kd_file)
        
        # Read merged regions
        merged_df = pd.read_csv(bed_file, sep="\t")
        
        # Read kd windows
        if os.path.exists(kd_file):
            kd_df = pd.read_csv(kd_file, sep="\t")
            if "chrom" not in kd_df.columns:
                kd_df[["chrom", "start", "end"]] = kd_df["Window_Region"].str.split("_", expand=True)
                kd_df["start"] = kd_df["start"].astype(int)
                kd_df["end"] = kd_df["end"].astype(int)
        else:
            kd_df = None

        # Plot merged regions (colored)
        for _, row in merged_df.iterrows():
            region_start = row["start"]
            region_end = row["end"]
            color = "blue"
            alpha = 0.6
            try:
                norm_val = float(row["Weighted_Kd_norm"])
                if not np.isnan(norm_val):
                    color = plt.cm.coolwarm(norm_val)
                    alpha = 0.8
            except Exception:
                pass
            ax.add_patch(
                patches.Rectangle(
                    (region_start, y_position),
                    region_end - region_start,
                    0.5,
                    color=color,
                    alpha=alpha,
                )
            )
            ax.text(
                (region_start + region_end) / 2,
                y_position + 0.25,
                f"{float(row[value_col]):.2f}" if row[value_col] not in ["NA", "", None, np.nan] else "",
                ha="center",
                fontsize=6,
                color="black"
            )

        # Plot only the top x% Kd windows NOT overlapping any merged region (gray)
        if kd_df is not None and not kd_df.empty:
            kd_df_valid = kd_df[kd_df["Weighted_Kd"].notna()]
            if not kd_df_valid.empty:
                cutoff = kd_df_valid["Weighted_Kd"].quantile(1-kd_percentile)
                top_kd_df = kd_df_valid[kd_df_valid["Weighted_Kd"] >= cutoff]
                for _, win in top_kd_df.iterrows():
                    win_start = int(win["start"])
                    win_end = int(win["end"])
                    # Check overlap with any merged region
                    overlaps_any = any(
                        (win_start < row["end"]) and (win_end > row["start"])
                        for _, row in merged_df.iterrows()
                    )
                    if not overlaps_any:
                        ax.add_patch(
                            patches.Rectangle(
                                (win_start, y_position),
                                win_end - win_start,
                                0.5,
                                color="gray",
                                alpha=0.5,
                                edgecolor="black",
                                hatch="//"
                            )
                        )
                        ax.text(
                            (win_start + win_end) / 2,
                            y_position + 0.25,
                            f"{float(win['Weighted_Kd']):.2f}",
                            ha="center",
                            fontsize=6,
                            color="black"
                        )
                        
        labels.append(base)
        y_position += 1

    for i, (name, (start, end)) in enumerate(overlaps.items()):
        y_top = y_position + i + 0.5
        ax.add_patch(patches.Rectangle((start, y_top - 0.5), end - start, 0.5, color="purple", alpha=0.6))
        ax.plot([start, start], [y_top, -1], linestyle="dashed", color="navy", alpha=0.9)
        ax.plot([end, end], [y_top, -1], linestyle="dashed", color="navy", alpha=0.9)
        ax.text(start + 100, y_top, name, fontsize=6, color="indigo", weight="bold")
    
    ax.set_xlim(coords)
    ax.set_ylim(-1, y_position + len(overlaps) + 1)
    ax.set_yticks([y + 0.25 for y in range(len(labels) + len(overlaps))])
    ax.set_yticklabels(labels + list(overlaps.keys()))
    ax.set_xlabel("Genomic Coordinates")
    ax.set_title(f"{ful_mid} enhancer vs. Overlap Regions (by {value_col}) [cut ={cut}, qval={cutoff_qval}, pval = {pval}] {title_suffix}")
    plot_path = os.path.join(out_dir, f"kd_enhancer_overlap_plot_{value_col}{ful_mid}_[cut={cut}_qval={cutoff_qval}_pval={pval}].png")
    
    # Legend patches
    kd_patch = mpatches.Patch(color=plt.cm.coolwarm(0.8), label="Merged region (colored by Weighted_Kd)")
    gray_patch = mpatches.Patch(facecolor="gray", edgecolor="black", hatch="//", label=f"Top {(kd_percentile*100)}% Kd window (not merged)")
    purple_patch = mpatches.Patch(color="purple", alpha=0.6, label="Reference overlap region")

    plt.legend(handles=[kd_patch, gray_patch, purple_patch], loc="upper right", fontsize=8, frameon=True)
    plt.savefig(plot_path, dpi=300, bbox_inches="tight")
    print(f"Plot saved to: {plot_path}")
    plt.show()

    # Heatmap
    heatmap_data = []
    heatmap_ranks = []
    total_vals_list = []
    for bed_file in annotated_bed_files:
        base = os.path.basename(bed_file).replace("_merged_sorted.bed", "")
        kd_file = os.path.join(os.path.dirname(bed_file), f"{base}_kd_scores_with_norm.tsv")
        merged_df = pd.read_csv(bed_file, sep="\t")
        kd_df = pd.read_csv(kd_file, sep="\t")
        all_vals = kd_df["Weighted_Kd"].replace("NA", np.nan).astype(float)
        all_vals = all_vals[~np.isnan(all_vals)]
        total_vals_list.append(len(all_vals))
        
        from scipy.stats import rankdata
        all_ranks = {}
        if len(all_vals) > 0:
            ranks = rankdata(-all_vals, method='min')
            for v, r in zip(all_vals, ranks):
                all_ranks[v] = int(r)
        
        row_vals = []
        row_ranks = []
        for name, (start_ref, end_ref) in overlaps.items():
            overlaps_df = merged_df[(merged_df["start"] < end_ref) & (merged_df["end"] > start_ref)]
            if not overlaps_df.empty:
                val = overlaps_df["Weighted_Kd"].astype(float).max()
            else:
                val = np.nan
            row_vals.append(val)
            row_ranks.append(all_ranks.get(val, ""))
        heatmap_data.append(row_vals)
        heatmap_ranks.append(row_ranks)

    heatmap_array = np.array(heatmap_data)
    heatmap_ranks = np.array(heatmap_ranks)

    annot = np.empty_like(heatmap_array, dtype=object)
    for i in range(heatmap_array.shape[0]):
        for j in range(heatmap_array.shape[1]):
            val = heatmap_array[i, j]
            rank = heatmap_ranks[i, j]
            total_in_file = total_vals_list[i]
            if np.isnan(val):
                annot[i, j] = ""
            else:
                annot[i, j] = f"{val:.2f}\n({rank}/{total_in_file})"
    
    sns.heatmap(
        heatmap_array,
        annot=annot,
        fmt="",
        cmap="flare",
        xticklabels=list(overlaps.keys()),
        yticklabels=labels,
        cbar_kws={'label': value_col}
    )
    plt.title(f"{ful_mid} Heatmap of {value_col} in Overlap Regions [cut ={cut}, qval={cutoff_qval}, pval = {pval}]")
    plt.tight_layout()
    heatmap_path = os.path.join(out_dir, f"kd_heatmap_{value_col}_{ful_mid}_[cut={cut}_qval={cutoff_qval}_pval={pval}].png")
    plt.savefig(heatmap_path, dpi=300)
    print(f"Heatmap saved to: {heatmap_path}")
    plt.show()

def plot_toppct_merged_boxes_and_heatmap(annotated_bed_files, overlaps, out_dir, value_col="Weighted_Kd", 
                                          percentile_cutoff=None):
    """
    Create visualization showing only top percentile merged regions by absolute Weighted_Kd values.
    """
    if percentile_cutoff is None:
        percentile_cutoff = kd_cut
    
    os.makedirs(out_dir, exist_ok=True)
    fig, ax = plt.subplots(figsize=(12, 4))
    y_position = 0
    labels = []

    for bed_file in annotated_bed_files:
        base = os.path.basename(bed_file).replace("_merged_sorted.bed", "")
        merged_df = pd.read_csv(bed_file, sep="\t")
        merged_df["Abs_Weighted_Kd"] = np.abs(merged_df["Weighted_Kd"].replace("NA", np.nan).astype(float))

        # Select top percentile merged regions by abs(Weighted_Kd)
        if not merged_df["Abs_Weighted_Kd"].isna().all():
            cutoff = merged_df["Abs_Weighted_Kd"].quantile(1-percentile_cutoff)
            merged_top = merged_df[merged_df["Abs_Weighted_Kd"] >= cutoff]
        else:
            merged_top = merged_df.iloc[[]]

        # Normalize for coloring
        if not merged_top.empty:
            min_abs = merged_top["Abs_Weighted_Kd"].min()
            max_abs = merged_top["Abs_Weighted_Kd"].max()
            merged_top["Abs_Weighted_Kd_norm"] = (merged_top["Abs_Weighted_Kd"] - min_abs) / (max_abs - min_abs) if max_abs > min_abs else merged_top["Abs_Weighted_Kd"]
        else:
            merged_top["Abs_Weighted_Kd_norm"] = np.nan

        # Plot only top percentile merged regions (colored)
        for _, row in merged_top.iterrows():
            region_start = row["start"]
            region_end = row["end"]
            color = "blue"
            alpha = 0.6
            try:
                norm_val = float(row["Abs_Weighted_Kd_norm"])
                if not np.isnan(norm_val):
                    color = plt.cm.coolwarm(norm_val)
                    alpha = 0.8
            except Exception:
                pass
            ax.add_patch(
                mpatches.Rectangle(
                    (region_start, y_position),
                    region_end - region_start,
                    0.5,
                    color=color,
                    alpha=alpha,
                )
            )
            ax.text(
                (region_start + region_end) / 2,
                y_position + 0.25,
                f"{float(row['Abs_Weighted_Kd']):.2f}" if not np.isnan(row['Abs_Weighted_Kd']) else "",
                ha="center",
                fontsize=6,
                color="black"
            )

        labels.append(base)
        y_position += 1

    # Plot overlap regions (purple)
    for i, (name, (start, end)) in enumerate(overlaps.items()):
        y_top = y_position + i + 0.5
        ax.add_patch(mpatches.Rectangle((start, y_top - 0.5), end - start, 0.5, color="purple", alpha=0.6))
        ax.plot([start, start], [y_top, -1], linestyle="dashed", color="navy", alpha=0.9)
        ax.plot([end, end], [y_top, -1], linestyle="dashed", color="navy", alpha=0.9)
        ax.text(start + 100, y_top, name, fontsize=6, color="indigo", weight="bold")
    
    ax.set_xlim(coords)
    ax.set_ylim(-1, y_position + len(overlaps) + 1)
    ax.set_yticks([y + 0.25 for y in range(len(labels) + len(overlaps))])
    ax.set_yticklabels(labels + list(overlaps.keys()))
    ax.set_xlabel("Genomic Coordinates")
    ax.set_title(f"{ful_mid} enhancer vs. Overlap Regions (top {percentile_cutoff}% Abs_Weighted_Kd merged) [cut ={cut}, qval={cutoff_qval}, pval = {pval}]")
    plot_path = os.path.join(out_dir, f"kd_enhancer_overlap_plot_top{int(percentile_cutoff*100)}pct_Abs_Weighted_Kd{ful_mid}_[cut={cut}_qval={cutoff_qval}_pval={pval}].png")

    # Legend patches
    kd_patch = mpatches.Patch(color=plt.cm.coolwarm(0.8), label=f"Top {percentile_cutoff}% merged region (by Abs_Weighted_Kd)")
    purple_patch = mpatches.Patch(color="purple", alpha=0.6, label="Reference overlap region")
    plt.legend(handles=[kd_patch, purple_patch], loc="upper right", fontsize=8, frameon=True)
    plt.savefig(plot_path, dpi=300, bbox_inches="tight")
    print(f"Plot saved to: {plot_path}")
    plt.show()

def plot_logpval_x_weightedkd_all(annotated_bed_files, overlaps, out_dir, title_suffix=""):
    """
    Create visualization of combined log p-value × weighted Kd scores.
    """
    os.makedirs(out_dir, exist_ok=True)
    fig, ax = plt.subplots(figsize=(14, 6))
    y_position = 0
    labels = []

    for bed_file in annotated_bed_files:
        base = os.path.basename(bed_file).replace(f"_p={pval}_merged_sorted.bed", "").split("_")[2]
        print(base)
        df = pd.read_csv(bed_file, sep="\t")
        if "logpval_x_weightedkd" not in df.columns:
            print(f"Skipping {bed_file}: no 'logpval_x_weightedkd' column.")
            continue

        for _, row in df.iterrows():
            region_start = row["start"]
            region_end = row["end"]
            score = row["logpval_x_weightedkd"]
            # Color by sign: blue for positive, red for negative, gray for nan/zero
            if pd.isna(score):
                color = "gray"
            else:
                color = plt.cm.coolwarm(score)
            ax.add_patch(
                mpatches.Rectangle(
                    (region_start, y_position),
                    region_end - region_start,
                    0.5,
                    color=color,
                    alpha=0.7,
                )
            )
            ax.text(
                (region_start + region_end) / 2,
                y_position + 0.25,
                f"{score:.2f}" if not pd.isna(score) else "",
                ha="center",
                fontsize=6,
                color="black"
            )
        labels.append(base)
        y_position += 1

    # Plot overlap regions (purple)
    for i, (name, (start, end)) in enumerate(overlaps.items()):
        y_top = y_position + i + 0.5
        ax.add_patch(mpatches.Rectangle((start, y_top - 0.5), end - start, 0.5, color="purple", alpha=0.6))
        ax.plot([start, start], [y_top, -1], linestyle="dashed", color="navy", alpha=0.9)
        ax.plot([end, end], [y_top, -1], linestyle="dashed", color="navy", alpha=0.9)
        ax.text(start + 100, y_top, name, fontsize=6, color="indigo", weight="bold")
    
    ax.set_xlim(coords)
    ax.set_ylim(-1, y_position + len(overlaps) + 1)
    ax.set_yticks([y + 0.25 for y in range(len(labels) + len(overlaps))])
    ax.set_yticklabels(labels + list(overlaps.keys()))
    ax.set_xlabel("Genomic Coordinates")
    ax.set_title(f"{ful_mid} enhancer vs. Overlap Regions (by logpval_x_weightedkd) [cut ={cut}, qval={cutoff_qval}, pval = {pval}] {title_suffix}")
    plot_path = os.path.join(
        out_dir,
        f"kd_enhancer_overlap_plot_logpval_x_weightedkd{ful_mid}_[cut={cut}_qval={cutoff_qval}_pval={pval}].png"
    )

    # Legend
    blue_patch = mpatches.Patch(color="blue", label="logpval_x_weightedkd < 0")
    red_patch = mpatches.Patch(color="red", label="logpval_x_weightedkd > 0")
    gray_patch = mpatches.Patch(color="gray", label="logpval_x_weightedkd = 0 or NA")
    purple_patch = mpatches.Patch(color="purple", alpha=0.6, label="Reference overlap regions")
    plt.legend(handles=[red_patch, blue_patch, gray_patch, purple_patch], loc="upper right", fontsize=8, frameon=True)
    plt.savefig(plot_path, dpi=300, bbox_inches="tight")
    print(f"Combined score plot saved to: {plot_path}")
    plt.show()

if __name__ == "__main__":
    enhancers = [
        "st10", "st10R", "IntP2A", "IntP2B", "pCRM", "Peak9820",
        "Peak9821", "Peak9823", "Peak9817", "AEL4", "5P3"]
    
    enhancer_files = []
    window_files = []
    for enhancer in enhancers:
        enh_file = enhancer + "_fimo.tsv"
        enhancer_files .append(enh_file)
        
        wind_file = enhancer + "_window_fimo.tsv"
        window_files.append(wind_file)

    # overlaps
    overlaps = {
        "sim_A1.0": (13056591, 13057759),
        "sim_st10": (13060273, 13060969),
        "sim_D2.1": (13064360, 13066277),
        "sim_E2.3": (13066277, 13067711),
        "EV_overlap": (13067711, 13068610),
        "sim_VT040842": (13068610, 13069801),
        "sim_1.6MLE": (13070722, 13072353)
    }

    enhancer_dir = rf"C:\Users\izzye\OneDrive - Johns Hopkins\Documents\BRIGHT\motif\{ful_mid}\enhancers"
    window_dir = rf"C:\Users\izzye\OneDrive - Johns Hopkins\Documents\BRIGHT\motif\{ful_mid}\windows"
    output_dir = rf"C:\Users\izzye\OneDrive - Johns Hopkins\Documents\BRIGHT\motif\{ful_mid}\kd_scores"
    plot_out_dir = rf"C:\Users\izzye\OneDrive - Johns Hopkins\Documents\BRIGHT\motif\{ful_mid}\kd_plots"
    bed_dir = rf"C:\Users\izzye\OneDrive - Johns Hopkins\Documents\BRIGHT\motif\{ful_mid}\{metric}\qval={cutoff_qval}\pval={pval}\bedgraph"
    kd_dir = rf"C:\Users\izzye\OneDrive - Johns Hopkins\Documents\BRIGHT\motif\{ful_mid}\kd_scores"
    save_dir = rf"C:\Users\izzye\OneDrive - Johns Hopkins\Documents\BRIGHT\motif\{ful_mid}\kd_scores"
    os.makedirs(save_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)

    for enhancer_file, window_file in zip(enhancer_files, window_files):
        output_file = os.path.join(output_dir, f"{enhancer_file.replace('_fimo.tsv', '')}_kd_scores.tsv")
        compare_enhancer_to_windows(enhancer_file, window_file, enhancer_dir, window_dir, output_file)

    # After all Kd score files are created, plot them
    kd_score_files = [os.path.join(output_dir, f"{e.replace('_fimo.tsv', '')}_kd_scores.tsv") for e in enhancer_files]
    scatter_dir = r"C:\Users\izzye\OneDrive - Johns Hopkins\Documents\BRIGHT\motif\{ful_mid}\kd_scatterplots"

    for bed_file in os.listdir(bed_dir):
        if not (bed_file.endswith(".bed") or bed_file.endswith(".bedgraph")):
            continue

        if "_merged" in bed_file:
            base = bed_file.split("_")[0]
            
            # Adjust this pattern to match your kd file naming
            kd_file = os.path.join(kd_dir, f"{base}_kd_scores_with_norm.tsv")
            print(f"looking for {kd_file}")
            bed_path = os.path.join(bed_dir, bed_file)
            output_path = os.path.join(save_dir, f"{ful_mid}_{metric}_{base}_q={cutoff_qval}_p={pval}merged_sorted.bed")
            if os.path.exists(kd_file):
                annotate_merged(bed_path, kd_file, output_path, overlaps=overlaps, plot_out_dir=plot_out_dir)
            else:
                print(f"No Kd file found for {bed_file}")

        else:
            print(bed_file)
            base = bed_file.split(".")[0]
            print(base)
            
            # Adjust pattern to match kd file naming
            kd_file = os.path.join(kd_dir, f"{base}_kd_scores_with_norm.tsv")
            print(f"looking for {kd_file}")
            bed_path = os.path.join(bed_dir, bed_file)
            output_path = os.path.join(save_dir, f"{base}_annotated.bed")
            if os.path.exists(kd_file):
                annotate_bed_with_kd(bed_path, kd_file, output_path)
            else:
                print(f"No Kd file found for {bed_file}")

    # List of annotated BED files
    annotated_bed_files = [
        os.path.join(save_dir, f"{base}_merged_sorted.bed")
        for base in enhancers
        if os.path.exists(os.path.join(save_dir, f"{base}_merged_sorted.bed"))
    ]
    
    ann_bed_files = [
        os.path.join(save_dir, f"{ful_mid}_{metric}_{base}_q={cutoff_qval}_p={pval}merged_sorted.bed")
        for base in enhancers
        if os.path.exists(os.path.join(save_dir, f"{ful_mid}_{metric}_{base}_q={cutoff_qval}_p={pval}merged_sorted.bed"))
    ]
    
    # Call plotting functions
    plot_kd_boxes_and_heatmap(annotated_bed_files, overlaps, plot_out_dir, value_col="Weighted_Kd")
    plot_toppct_merged_boxes_and_heatmap(annotated_bed_files, overlaps, plot_out_dir, value_col="Weighted_Kd")

    # List of all window Kd score files (with norm)
    kd_score_files = [
        os.path.join(output_dir, f"{base}_annotated.bed")
        for base in enhancers
        if os.path.exists(os.path.join(output_dir, f"{base}_annotated.bed"))
    ]

plot_logpval_x_weightedkd_all(ann_bed_files, overlaps, plot_out_dir)


