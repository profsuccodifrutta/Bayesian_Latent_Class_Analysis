# Bayesian Latent Class Analysis
This repository contains code and final report for a unsupervised latent clustering problem (labels not observed) involving only cathegorical variables. In this project the model was used to identify
latent patterns in marketing data for a telecommunication company operating trough a physical retail network

## Applications
The model can be applied whenever you are dealing with only cathegorical data, especially meaningful for marketing segmentation tasks, healthcare and epydemiology, political sciences and voting  behaviour.
Thanks to its probabilistic foundation interpretability is very high. 
In general BLCA is a valuable tool both for explanotary and confirmatory analysis.

## Content
* full_conditional_functions: contains R code for the sampling from the full conditional distribution of the 3 parameter of the model
* gibbs_sampler_and_posterior_analysis: contains implementation of the gibbs sampler, diagnostics, and posterior analysis of all parameters
* number_of_cluster_selection: contain frequentist approach based on integrated complete likelihood for selecting the number of cluster. To enrich the analysis one could add a prior on the number of cluster, to have a full bayesian approach
* PPC: contains implementatio of posterior predictive checks, crucial for checking model adequacy with respect to the observed data
