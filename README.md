## About this repository
This repository contains an R script for a latent class model to estimate sensitivity and specificity of multiple repeated diagnostic tests, as it was used for the publication by Wichert et al. ... The R script consists mainly of the JAGS model and can be run in R via runjags or rjags, for example.

## About the model
This latent class model was proposed by Wang and Hanson in 2019 to estimate sensitivity and specificity of multiple diagnostic tests, which are applied repeatedly to the same subjects, if a gold standard is not or only partially available (doi: [10.1002/sim.8114](https://doi.org/10.1002/sim.8114)). Our implementation builds on code for a related model without repeated measurements provided by the same author (Wang, Lin &amp; Nelson 2020, doi: [10.1177/0962280219852649](https://doi.org/10.1177/0962280219852649), model M2). We adapted the code to repeated tests in accordance with the definition by Wang &amp; Hanson (2019). To estimate sensitivity and specificity of test combinations, we further expanded the model to parallel testing, as described in Wang &amp; Hanson (2019). 

## Usage and arguments
The following arguments have to be specified:
<table>
  <tr>
    <td>K</td>
    <td>number of diagnostic tests</td>
  </tr>
  <tr>
    <td>J</td>
    <td>number of time points</td>
  </tr>
  <tr>
    <td>N</td>
    <td>number of subjects</td>
  </tr>
  <tr>
    <td>x</td>
    <td>array of dimension N &#10799; K &#10799; J containing the binary test results (0: negative, 1: positive)</td>
  </tr>
  <tr>
    <td>ref</td>
    <td>vector of length N indicating the true disease status of each subject (-1: unknown, 0: negative, 1: positive)</td>
  </tr>
  <tr>
    <td>z</td>
    <td>vector of length N containing zeros (for Poisson zero trick)</td>
  </tr>
  <tr>
    <td>ncomb</td>
    <td>number of test combinations to be considered for parallel testing</td>
  </tr>
  <tr>
    <td>comb</td>
    <td>matrix of dimension ncomb &#10799; K. Each row codes a selection of the given binary tests using zeros for ‘not included’ and ones for ‘included’, e.g., if Test 1 and 3 out of three tests should be combined into a parallel test, this would be given as '1, 0, 1'.</td>
  </tr>
  <tr>
    <td>z2</td>
    <td>matrix of dimension ncomb &#10799; J containing zeros (Poisson zero trick for sensitivity of parallel tests)</td>
  </tr>
  <tr>
    <td>z3</td>
    <td>matrix of dimension ncomb &#10799; J containing zeros (Poisson zero trick for specificity of parallel tests)</td>
  </tr>
</table> 
