# k-mer Based Approach

## Overview
The k-mer based approach compares windows from a locus of a genome, or randomley selected windows from the entire genome to a given enhancer sequence in order to determine with locus window sequence most likely correlates to that enhancer.

The following scripts may be ran in order:
Kmer_dist -> Kmer_visualizer -> Kmer_bed -> Kmer_comparison

## Approach
### Kmer_dist:
Takes FASTA sequence files from windows in a locus of a genome and randomly selected windows from the entire genome. In our case we used the _Drosophila_ genome and the _Drosophila sim locus_. These windows are then analyzed based on a k-mer and distance metric. 

In our case k-mer values 2-6 were used and the following Alfpy word-based distance metrics: 
canberra, braycurtis, chebyshev, csv, euclid_norm, euclid_squared, google, jsd, kld, lcc, manhattan, minkowski

Locus windows and random windows are run in batches and scored based on their k-mer similarty to their respective enhancer. A bedGRAPH file is output for each Locus Window file or Random Window file and reverse complement containing the chromose, window start, window end, score for each window in that file.

### Kmer_visualizer:
Takes the bedGRAPH output from the Kmer_dist as input.

Visualizes a scatterplot of the genomic distribution of the windows for the score and k-mer metrics. Optionally creates a histogram of the scatterplot.

Compares distributions of window and random distances through a boxplot.

Caclulates the empirical pvalue for each locus window to enhancer as opposed to the random windows to enhancer. Filters by a sepecified pvalue cutoff. Appends results to a bed file with the enhancer name (ex: st10, st10_rev). Output how much overlap there is in each window.

Optionally creates a bar plot comparing pvalues for different enhancer windows.

### Kmer_bed:
Takes the bed output from Kmer_visualizer as input.

Groups files by their enhancer (i.e. st10 and st10_rev). Sorts and merges and overlapping windows in grouped files into one file with chromose, window start, window end, pvalue. Output bed files of results.

### Kmer_comparison:
Takes the bed output from Kmer_bed as input.

Outputs overlap plot based on reference overlap regions (significant comparison enhancers). 

Saves pvalue output to csv file and creates a scaled heatmap based on the reference overlap regions and the rank of each window hit from all the hits.
