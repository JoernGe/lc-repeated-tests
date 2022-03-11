## About this repository
This repository contains an R script for a latent class model to estimate sensitivity and specificity of multiple repeated diagnostic tests, as it was used for the publication by Wichert et al. ... The R script consists mainly of the JAGS model and can be run in R via <code>runjags</code> or <code>rjags</code>, for example.

## About the model
This latent class model was proposed by Wang and Hanson in 2019 to estimate sensitivity and specificity of multiple diagnostic tests, which are applied repeatedly to the same subjects, if a gold standard is not or only partially available (doi: [10.1002/sim.8114](https://doi.org/10.1002/sim.8114)). 

Our implementation builds on code for a related model without repeated measurements (Wang, Lin &amp; Nelson 2020, doi: [10.1177/0962280219852649](https://doi.org/10.1177/0962280219852649), model M2). We adapted the code to repeated tests in accordance with the definition by Wang &amp; Hanson (2019). To estimate sensitivity and specificity of test combinations, we further expanded the model to parallel testing, as described in Wang &amp; Hanson (2019). 

Parallel testing means that multiple tests are applied to the same subject and if at least one test gives a positive result, the subject is diagnosed as positive. Thus, for a negative test result, all included tests have to be negative. Such a combination of test results can include multiple tests and/or multiple time points. In our implementation, the number of time points included in a parallel test can not be specified, but our code returns the sensitivity and specificity of a parallel test for 1, 2, ..., J time points, if J gives the maximum number of time points of the study. This behaviour was hard-coded for our purposes to inspect the diagnostic accuracy of parallel tests with increasing numbers of applications.

## Usage and arguments
Download the R script into your working directory and run it once in R to create <code>repeated_measurements_parallel_tests.bug</code>. This file is used as input for <code>runjags</code> or <code>rjags</code>.

The following arguments have to be specified in R:
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
    <td>matrix of dimension ncomb &#10799; K. Each row codes a combination of the diagnostic tests using zeros for ‘not included’ and ones for ‘included’, e.g., if Test 1 and 3 out of three tests should be combined into a parallel test, this would be given as '1, 0, 1'.</td>
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

For example, if a study comprised 50 participants with unkown disease status, which were tested at 10 time points using three diagnostic tests, and the diagnostic accuracy of repeated applications of Test 1 should be evaluated over time, the arguments could be specified in R in the following way:

```r
comb_mat <- matrix(c(1, 0, 0), nrow = 1, ncol = 3)

list_data <- list("K" = 3,
                  "J" = 10,
                  "N" = 50,
                  "x" = data_array,
                  "ref" = rep(-1, 50),
                  "z" = rep(0, 50),
                  "ncomb" = nrow(comb_mat),
                  "comb" = comb_mat,
                  "z2" = matrix(0, nrow = nrow(comb_mat), ncol = 10),
                  "z3" = matrix(0, nrow = nrow(comb_mat), ncol = 10))
```
We assume here that the test results are stored in <code>data_array</code>. If parallel tests are not required, parameter <code>ncomb</code> can simply be set to zero to suppress calculations for parallel tests.

The following variables of our model may be included in the list of monitored variables:
<table>
  <tr>
    <td>se</td>
    <td>estimated sensitivity of each test (vector of length K)</td>
  </tr>
  <tr>
    <td>sp</td>
    <td>estimated specificity of each test (vector of length K)</td>
  </tr>
  <tr>
    <td>rp</td>
    <td>estimated correlation between repeated measurements in the diseased population for each test (vector of length K)</td>
  </tr>
  <tr>
    <td>rn</td>
    <td>estimated correlation between repeated measurements in the non-diseased population for each test (vector of length K)</td>
  </tr>
  <tr>
    <td>cop</td>
    <td>estimated pairwise correlation between tests in the diseased population (matrix of dimension K &#10799; K, here a strictly upper triangular matrix, i.e. values on the main diagonal and below are set to zero to omit duplications)</td>
  </tr>
  <tr>
    <td>con</td>
    <td>estimated pairwise correlation between tests in the non-diseased population (matrix of dimension K &#10799; K, a strictly upper triangular matrix like cop)</td>
  </tr>
  <tr>
    <td>pi</td>
    <td>estimated disease prevalence</td>
  </tr>
  <tr>
    <td>prob</td>
    <td>‘likelihood contribution’ of each subject (vector of length N)</td>
  </tr>
  <tr>
    <td>omega1</td>
    <td>mode of beta distribution of test sensitivities (see Wang, Lin &amp; Nelson (2020))</td>
  </tr>
  <tr>
    <td>omega2</td>
    <td>mode of beta distribution of test specificities (see Wang, Lin &amp; Nelson (2020))</td>
  </tr>
  <tr>
    <td>kappa1</td>
    <td>spread of beta distribution of test sensitivities (see Wang, Lin &amp; Nelson (2020))</td>
  </tr>
  <tr>
    <td>kappa2</td>
    <td>spread of beta distribution of test specificities (see Wang, Lin &amp; Nelson (2020))</td>
  </tr>
  <tr>
    <td>par_se</td>
    <td>estimated sensitivity of each parallel test (matrix of dimension ncomb &#10799; J, first index: row index of parallel test in matrix comb, second index:
number of repeated measurements included in the parallel test)</td>
  </tr>
  <tr>
    <td>par_sp</td>
    <td>estimated specificity of each parallel test (matrix of dimension ncomb &#10799; J, see par_se)</td>
  </tr>
</table>

For example, the model could be run using <code>runjags</code> like this:

```r
library(runjags)

out <- run.jags(model = "repeated_measurements_parallel_tests.bug",
                data = list_data, 
                monitor = c("se", "sp", "kappa1", "kappa2", "omega1", "omega2", 
                            "rp", "rn", "cop", "con", "pi",
                            "par_se", "par_sp", "prob"),
                n.chains = 3)
```
\- see the reference manual and vignettes for <code>runjags</code> ([here](https://cran.r-project.org/web/packages/runjags/index.html)) or <code>rjags</code> ([here](https://cran.r-project.org/web/packages/rjags/index.html)) for details on running JAGS models in R.
