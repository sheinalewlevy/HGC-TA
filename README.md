# HGC-TA
This is the data and code for: "Socioecology shapes child and adolescent activity budgets in twelve hunter-gatherer and mixed-subsistence forager societies"

Datasets are the following:
- 12 datasets for 11 hunter-gatherer and mixed subsistence societies. The datasets are named as **Society_Author**. Note that data for Tsimane should be requested directly from Jonathan Stieglitz (jonathan.stieglitz@gmail.com)
- **group_variables.csv**: The society-level variables used in the analysis
- **dataset.csv**: The processed dataset used in the analysis
- **TAX.csv**: Model estimates. X=model numbers from 2_models.R

Data analysis R scripts are the following:
- **0_datapreparation.R**: Code for preparing the dataset for analysis, in-text Figure 2, and additional checks. This script can be used alongside the society-specific datasets. 
- **1_ecologicalvariables.R**: Code used to generate the climatic and risk variables. Note that we are sharing this script for transparency purposes only--we have not supplied the Environmental_variables file because this file contains GPS points for each field site. In order to protect the privacy of participants, who for the most part come from small communities with few members, we have chosen not to publish these GPS points.
- **2_models.R**: Code used in the stan models. This code can be used in conjunction with dataset.csv. Note that results will not match those in the paper without Tsimane data. Note as well that models take weeks to fit.
- **3_figures.R**: Code used to investigate posterior distributions and make in-text and supplementary figures. This could cal be used in conjunction with dataset.csv and all .rda files (see below). Note that results will not match those in the paper without Tsimane data. 

Due to file size limits, model outputs are [available for download as a zip file here](https://www.dropbox.com/s/chs9eyejzoernv2/post_rv.zip?dl=0) or upon request from Sheina Lew-Levy (sheinalewlevy@gmail.com). Model outputs are:
- **postX.rda**: Posterior distributions from each model. X=model numbers from 2_models.R
- **TAwaic.rda**: WAIC values for the main models in the text

*Update 02/Feb/2022*

During the revision stage, model numbers were changed throughout to match the updated order presented in the text. These changes were made to all models in **2_models.R** and are reflected in **TAX.csv**, **3_figures.R**, **postX.rda** and **TAwaic.rda**.
