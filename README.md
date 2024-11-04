# spbrom

The R package *spbrom* is the companion package for the paper Castillo-Mateo et al. "???". The package includes scripts and functions to reproduce all the results from the paper. 

## Workflow of the paper

* The data that contains basic information for the weather stations and that will be updated in subsequent scripts is in `data\stations.rda`
* The scripts to reproduce the results from the paper are in `inst\scripts\`:
  + `01_1_download_Tmax_Tmin_data.R` downloads the daily maximum temperature data from ECA
  + `01_2_grid_dist.R` builds a 25 x 25 km grid for peninsular Spain and obtains the distance to the coast for each grid cell and weather station
  + `01_3_EDA.R` includes code to reproduce some figures for the exploratory data analysis
  + `02_1_data.R` builds the data frame that will be used to fit the models
  + `02_2_KFCV.R` obtains the metrics using 10-fold cross-validation for each model
  + `03_1_fitting.R` implements model fitting for the full model
  + `03_2_checks.R` obtains DIC for all models and does the posterior predictive checks (adequacy) for the full model
  + `04_1_results_parameters.R` uses the full model parameters to do posterior inference
  + `04_2_results_inference` obtains model-based tools using posterior predictive samples from the full model
