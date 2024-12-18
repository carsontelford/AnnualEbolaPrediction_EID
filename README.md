# AnnualEbolaPrediction_EID
 Repo for replicating analysis to estimate annual odds of ebolavirus spillover in equatorial Africa.

 This repository contains 3 R script files and 1 csv with cleaned data to run the analysis of ebolavirus spillover data from 2001-2022. 
 1. P1 Absence point generation: Script to generate 10,000 absence points throughout equatorial Africa with weighting based on log human population count, and merging with data on ebolavirus spillover events from 2001-2022.
 2. P2 Model fit and CV: Script to analyze the cleaned dataset titled tdf_clean. This dataset contains covariate values for the spillover events that occurred from 2001-2022 in the year they occurred, as well as covariate values for the absence locations in the year they were randomly assigned. Cross validation is performed using model ensembles. We run 100 models, each model randomly sampling 50 absence locations per presence location. AUC, sensitivity, and specificity are calculated.
 3. P3 Prediction grid: Script to fit the ensemble of models to prediction grids for 2021 and 2022 to estimate the odds of ebolavirus spillover in each prediction grid cell in equatorial Africa.
 4. tdf_clean: post processed dataset contained covariate values surrounding each presence and absence location, with rows for 24 spillover events and 10000 random absence locations.

Analyses outlined in scripts P1 and P3 require large files that cannot be accomodated on this repository, however R code is provided to outline the processes. tdf_clean contains postprocessed data with covariate values extracted surrounding each presence/absence location, therefore the complete analysis in script P2 can be done with this dataset. 

Additional scripts saved as .txt files contain java script to extract covariates within buffers in google earth engine.
 
