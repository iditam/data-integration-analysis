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
The code of the Naïve Bayesian classifier

Input: Gold standard genes, input data matrix

Output: Table of all genes with classifier scores

## utils.R
Functions of ratio score computing, used by the run_aggregated_classifier.R script

## module_prediction.ipynb
Script for characterizing the top classifier hits.
Using XGBoost algorithm for module assignment, LIME explainer to explain the most confident predictions for each functional module and SHAP method to estimate overall feature importance.

Input files: 

classifier_464_genes.csv - list of 464 classifier genes

HRR_gold_standard_genes_and_modules.csv - list of 78 GS genes with modules

string_fields.csv - data from STRING

S4_classifier_extra_fields.csv

MRS: 6 files of "HPA_coexpression_*.csv" and 6 files of "CladePP_*.csv" (one for each module)

