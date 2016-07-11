# CodeShare
This repository  contains code files I have been generating to learn BioInformatics.
Motif.py, Replication.py files contains functions to find motifs and patterns in a DNA sequence. It includes greedyMotif, RandomMotifSearch, GibbsSampler algorithms etc. The main aim of code is to find origin of replication "oriC" which in most cases can be found with pattern detection in sequence.

Network.h and Network.cpp files describe creation of weighted graph (gene expression) and find out trans associated gene based on CNV value. After getting all output data k-means clustering is done to get cancer sub types. Equations for finding trans associated genes are taken from a research paper.
Here we create a network of 13038 proteins. Weight of edges is decided by gene expression values of 997 tumor pateints. Vertex weight is combination of its own copy number variation value and its neighour's weight.  
