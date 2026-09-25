## A Bayesian hierarchical model for preferentially sampled spatial extremes data: an application to air quality.

This repository contains the minimal code and data to fit the proposed models and to (partially) reproduce the analysis performed in A Bayesian hierarchical model for preferentially sampled spatial extremes data: an application to air quality by Bianca Brusco and Veronica J. Berrocal. 


### Contents: Functions and Model Fitting
 
The folder `R` contains the code to fit the Extrmal Preferential model proposed in the paper, as well as the comparison Baseline model. In particular

  - The Extremal Preferential Model  (Spatial Only) : `R/epsm.R`. 
  - The Extremal Preferential Model  (Spatio-Temporal) : `R/epsm-st.R`. 
  - The Extremal Baseline Model  (Spatial Only) : `R/ebm.R`. 
  - The Extremal Baseline Model  (Spatio-Temporal) : `R/ebm-st.R`  

### Contents: Maximum PM2.5 in California Analysis
  
The folder `data-analysis` contains the data and the code to fit the models applied to the PM2.5 Pollution data in california, as well as code to generate predictive distributions and perform the attainment analysis presented in the paper.

  - `data` folder with the clean data file used for the analysis.  
  - `R/fit-models.R` fits the model of interest.  
  - `R/mod1-posteriorpreds.R` : generate predictive distributions for preferential model.  
  - `R/modb1-base-posteriorpreds.R` : generate predictive distributions for baseline model.  
  - `R/predictive-utils.R`: functions to generate posterior predictive/ posterior distributions.  
  - `R/attainment-analysis.R` : reproduce NAAQS attainment analysis   
  - `results` : a thinned version of the MCMC used in the analysis.   


### Examples

- [Fitting models Example Vignette](https://github.com/biabbiassago/epsm/blob/main/vignettes/example-fits.md)
- `inst/scripts/example.R`
  
  
 
 
