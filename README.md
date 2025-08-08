# k-mer Based Approach

## Overview
The k-mer based approach compares windows from a locus of a genome, or randomley selected windows from the entire genome to a given enhancer sequence in order to determine with locus window sequence most likely correlates to that enhancer.

The following scripts may be ran in order:
Kmer_dist -> Kmer_visualizer -> Kmer_bed -> Kmer_comparison

## Approach
### K_mer_dist:
Takes FASTA sequence files from windows in a locus of a genome and randomly selected windows from the entire genome. In our case we used the _Drosophila_ genome and the _Drosophila sim locus_. These windows are then analyzed based on a k-mer and distance metric. 

In our case k-mer values 2-6 were used and the following Alfpy word-based distance metrics: 
canberra, braycurtis, chebyshev, csv, euclid_norm, euclid_squared, google, jsd, kld, lcc, manhattan, minkowski
