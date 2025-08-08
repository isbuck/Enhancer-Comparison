# Overview
Various other scripts used.

## Shuffle
Randomly shuffles enhancer base pairs by x%.

## Motif_Ordering
- extract_ordered_motifs: Extract and order motifs from a FIMO output file based on genomic position
- extract_ordered_motifs_from_df: Extract and order motifs from a pandas DataFrame containing FIMO results 
- motif_order_kd: Compute Kendall's tau correlation between two ordered motif lists
- compare_enhancer_to_windows: Compare motif ordering between an enhancer and window regions
- annotate_bed_with_kd: Annotate BED file regions with Kendall's tau correlation scores
- annotate_merged: Annotate merged BED regions with best Kd scores from overlapping windows
- plot_kd_boxes_and_heatmap: Create comprehensive visualization of Kd scores with box plots and heatmaps
- plot_top10pct_merged_boxes_and_heatmap: Create visualization showing only top percentile merged regions by absolute Weighted_Kd values
- plot_logpval_x_weightedkd_all: Create visualization of combined log p-value × weighted Kd scores





