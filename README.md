# Assessing-seasonality-and-the-role-of-its-potential-drivers-in-environmental-epidemiology-a-tutoria

This repository stores the updated R code to reproduce the analyisis presented in the article:

Madaniyazi, Lina, et al. "Assessing seasonality and the role of its potential drivers in environmental epidemiology: a tutorial." International Journal of Epidemiology (2022): dyac115. DOI: https://doi.org/10.1093/ije/dyac115. 

# Data
The example dataset are adapted from the original data available from [https://github.com/gasparrini/2014_armstrong_BMCmrm_Codedata](url). The R script main_analysis.R in the folder includes the code for loading the dataset from the original sources. 

# R code
The five R scripts reproducdes all the steps of the analysis and the full results. Specifically:

cyclic.R is the R function for the cyclic spline function to assess seasonality.

findmin.R is the R function to estimate the trough of seasonality.

findmax.R is the R function to estimate the peak of seasonality.

attrs.R is the R function to estimate attributable fraction of seasonality.

main_analysis.R is the codes to reproduce the examples.
