# README
Analysis scripts for **"Differential and co-expression modeling reveal ACE and perceived stress trajectory coordinate human sperm sncRNA"**.

This repository contains the R scripts used to generate the results and figures presented in the associated article.

Data Inputs

The scripts require the following input data:
    **Expression Data:** sncRNA expression normalized via the TMM method.
    **List of dynamic sncRNAs**: sncRNA set selected according to the methodology described in the article.
    **List of dynamic PSS Samples**: Sample subsets identified based on the Perceived Stress Scale (PSS) trajectories dynamics, described in the article.
  

DREAM_models.R - differential expression analysis and volcano plot
WGCNA_coexpression.R - sncRNA co-expression analysis, co-expression modules, etc.
WGCNA_preservation.R - sncRNA co-expression modules preservation test.
DINGO.R - Differential sncRNA co-expression analysis. 
