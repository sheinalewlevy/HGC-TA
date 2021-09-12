# HGC-TA
This is the data and code for: "Cross-cultural variation in child and adolescent time allocation to work and play in twelve hunter-gatherer and mixed-subsistence societies."

Datasets are the following:
- 12 datasets for 11 hunter-gatherer and mixed subsistence societies. The datasets are named as **Society_Author**. Note that data for Tsimane should be requested directly from Jonathan Stieglitz (jonathan.stieglitz@gmail.com)
- **group_variables.csv**: The society-level variables used in the analysis
- **dataset.csv**: The processed dataset used in the analysis

Data analysis R scripts are the following:
- **0_datapreparation.R**: Code for preparing the dataset for analysis, in-text Figure 2, and additional checks. This script can be used alongside the society-specific datasets. 
- **1_ecologicalvariables**.R: Code used to generate the climatic and risk variables. Note that we are sharing this script for transparency purposes only--we have not supplied the Environmental_variables file because this file contains GPS points for each field site. In order to protect the privacy of participants, who for the most part come from small communities with few members, we have chosen not to publish these GPS points.
- **2_models.R**: Code used in the stan models. This code can be used in conjunction with dataset.csv. Note that results will not match those in the paper without Tsimane data. Note as well that models take weeks to fit.
- **3_figures.R**: Code used to investigate posterior distributions and make in-text and supplementary figures

Model outputs are the following:
- **postX.rda**: Posterior distributions from each model. X=model numbers from 2_models.R
- **WAIC.rda**: WAIC values for the main models in the text
- **TAX.rda**: Model estimates. X=model numbers from 2_models.R
