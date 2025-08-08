# Motif Based Approach

## Overview
The motif based approach compares windows from a locus of a genome, or randomly selected windows from the entire genome to a given enhancer sequence in order to determine with locus window sequence most likely correlates to that enhancer based on motif content.

The following scripts may be ran in order: Motif_dist -> Motif_visualizer -> Motif_bed -> Motif_comparison

## Approach
### Motif_dist
Takes tsv output FIMo files from windows in a locus of a genome and randomly selected windows from the entire genome. In our case we used the Drosophila genome and the Drosophila sim locus. These windows are then analyzed based on a their motif content for a distance metric.

Locus windows and random windows have their motifs extracted which are then filtered based on a qvalue cutoff. We used a cutoff of qvalue = .05 for _FlyFactorSurvey_ output files and no cutoff (qvalue = 1) for _MidlineSpecific_ output files.

Each window is then scored based on its motif contnent in relation to its respective enhancer by either Jaccard or Summed Canberra distance. Jaccard distance is 1 - intersection/union. Summed Canberra distance calculates the canberra distance between a window files and enhancer file after building a summed score matrix based on the motif scores.

A bedGRAPH file is output for each Locus Window file or Random Window file containing the chromose, window start, window end, score for each window in that file.

There is an optional build in for permuted files.

### Motif_Visualizer
Takes the bedfile output from Motif_dist as input.

Visualizes a two scatterplots of the genomic distribution of the windows for their score. One with significant pvalues overlayed.

Optionally plots jaccard peaks to help better visualize scatterplots.

Caclulates the empirical pvalue for each locus window to enhancer as opposed to the random windows to enhancer. Filters by a sepecified pvalue cutoff. Appends results to a bed file with the enhancer name (ex: st10, st10R). 

Optionally applys the Benjamin-Hochberg correction to the pvalues and saves raw and corrected pvalues to csv file. Plots BH adjusted points as a scatterplot.

Optionally creates a bar plot comparing pvalues for different enhancer windows.

### Motif_bed
Takes the bed output from Motif_visualizer as input.

Sorts and merges and overlapping windows in files into file with chromose, window start, window end, pvalue. Output bed files of results.

### Motif_comparison
Takes the bed output from Motif_bed as input.

Outputs overlap plot based on reference overlap regions (significant comparison enhancers).

Saves pvalue output to csv file and creates a scaled heatmap based on the reference overlap regions and the rank of each window hit from all the hits.

Optionally annotates merged bed files with motif hits.

Optionally compares motif hits in merged windows to those in enhancer regions and outputs a csv file of shared and unique hits.
