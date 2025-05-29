# spbrom

The R package *spbrom* is the companion package for the paper Castillo-Mateo et al. "Joint space-time modelling for upper daily maximum and minimum temperature record-breaking". The package includes scripts and functions to reproduce all the results from the paper. 

## Installation

You can install version 0.0.1 from the provided folder
```s
if (!require("remotes")) install.packages("remotes")
remotes::install_local("reproducibility_materials.zip")
```

You can install the **development** version from
[GitHub](https://github.com/JorgeCastilloMateo/spbrom)
```s
if (!require("remotes")) install.packages("remotes")
remotes::install_github("JorgeCastilloMateo/spbrom")
```

## Workflow of the paper

* The data that contains basic information for the weather stations and that will be updated in subsequent scripts is in `data\stations.rda`
* The scripts to reproduce in order the results from the paper are in `inst\scripts\`:
  + `000_downloadECAD.R` downloads the daily maximum and minimum temperature data from ECA
  + `001_grid.R` builds a 25 x 25 km grid for peninsular Spain and obtains the elevation and distance to the coast for each grid cell and weather station
  + `002_eda.R` includes code to reproduce the figures for the exploratory data analysis
  + `003_data.R` builds the data frame that will be used to fit the models
  + `004_logSDx.R` interpolation of the log-standard deviation of daily maximum temperatures for the anisotropic models
  + `005_KFCV.R` obtains the metrics using 10-fold cross-validation for each model
  + `006_model_fitting.R` implements model fitting for each model
  + `007_traceplots.R` convergence traceplots for the full model
  + `008_DIC.R` computes the DIC for each fitted model
  + `009_checking.R` posterior predictive checks (adequacy) for the full model
  + `010_results_parameters.R` uses the full model parameters to do posterior inference
  + `011_results_inference` obtains model-based tools using posterior predictive samples from the full model
  + `MVN` is a C++ implementtion to sample efficiently from the multivariate normal distribution
