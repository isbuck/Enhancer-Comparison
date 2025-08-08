# k-mer Based Approach

## Overview
The k-mer based approach compares windows from a locus of a genome, or randomley selected windows from the entire genome to a given enhancer sequence in order to determine with locus window sequence most likely correlates to that enhancer.

## Approach
takes FASTA sequence files from windows in a locus of a genome and randomly selected windows from the entire genome. In our case we used the _Drosophila_ genome and the _Drosophila sim locus_. These windows are then analyzed for their k-mers in the range of 2-6 (any value can be input). 
