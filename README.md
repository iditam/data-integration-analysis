# Data integration analysis
Code scripts used for the analysis for the paper
"Multi-omics data integration analysis identifies the spliceosome as a key regulator of DNA double-strand break repair"

## CladePP.R
This script implements the main methods of the cladePP algorithm.
It takes as the following inputs:
clusters - a rds file with a named list of hclust object representing hierarchial clustering of the NPP matrix, one object per clade. The names of the list elements should correspond to the clades.
genelist - A txt file with the gene symbols of the genes of interest.
outfile - Path to file to which the output shall be written.The output is a table indicating for each gene with which gold standard genes it is clustered,
in which clade and its MRS score.

## run_aggregated_classifier.R
The classifier code, input: Gs genes, input data matrix. output: table of all genes with classifier scores

## utils.R
Functions of ratio score computing, used by the run_aggregated_classifier.R script

## module_prediction.ipynb
the code with all the predictions etc

## draw_networkx_for_prediction.ipynb
creates edges for network figures
