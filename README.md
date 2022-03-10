## About this repository
This repository contains an R script for a latent class model to estimate sensitivity and specificity of multiple repeated diagnostic tests, as it was used for the publication by Wichert et al. ... The R script consists mainly of the JAGS model and can be run in R via runjags or rjags, for example.

## About the model
This latent class model was proposed by Wang and Hanson in 2019 (doi: [10.1002/sim.8114](https://doi.org/10.1002/sim.8114)). Our implementation builds on code provided by the same author for a related model without repeated measurements (Wang et al. 2020, doi: [10.1177/0962280219852649](https://doi.org/10.1177/0962280219852649)). We adapted the code to repeated tests in accordance with the definition by Wang &amp; Hanson (2019). To estimate sensitivity and specificity of test combinations, we further expanded the model to parallel testing, as described in Wang &amp; Hanson (2019). 

## Usage and arguments
In order to run this model in R successfully, the following 
