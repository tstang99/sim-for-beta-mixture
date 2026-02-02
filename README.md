# Simulation code for the paper "Improving parameter estimation and order selection for beta mixture models"

This repository contains the R code to reproduce the simulation and application results 
in *“Improving parameter estimation and order selection for beta mixture models.”* 
The main script sources helper routines in `functions/` and 
(i) runs the three-component beta-mixture simulation comparing PMLE, MLE, SAM, and 
SAR using Wasserstein distance, Jensen–Shannon divergence, and RMSLE; 
(ii) benchmarks computation time across sample sizes/model orders 
and compares traditional order selection vs. BMR; and 
(iii) reproduces the DNA methylation application.