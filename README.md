# Penalized generalized estimating equations for relative risk regression 
  
## Table of contents
   * [How to cite?](#how-to-cite)
   * [Contents overview](#contents-overview)
   * [Illustrative analysis](#illustrative-analysis)

## How to cite?

See our manuscript on [bioRxiv](XXX).

# Contents overview

The repository containts code for simulation of correlated binary data as described in [Quaqish (2003)](https://doi.org/10.1093/biomet/90.2.455), and fitting (i) log-link generalized estimating equations with identity variance and unknown dispersion (RR-GEE) and (ii) a penalized version (RR-PGEE), with the gradient of the Jeffreys-prior added to the estimating equations. The latter ensures finiteness of the estimates when boundary estimates occur.  

## Illustrative analysis
Our manuscript compares the two log-link GEEs through extensive simulations and the illustrative example here includes the simulation of 100 data sets from the base simulation scenario (as referred to in the manuscript). The code includes the following steps
  * simulating the data through specifying the simulation parameters such as number of subjects, number of time points, within-cluster correlation, etc, using the code as part of the R package "binarySimCLF_1.0.tar" customized to simulate log-link binary data; 
  * fitting RR-GEE and RR-PGEE using the code provided by us in functions xxx;
  * obtaining the simulation summaries as discussed in the [manusctipt](XXX) in Tables 2, 3, and 4;
  * replicating Figure 1 of the manuscript. 

Note that the minimal example here is for lower number of data sets (100, not 1000) compared to the manuscript to ensure the code can run for a couple of minutes and the code included is meant to replicate only one simulation scenario, i.e. one line of each of the tables mentioned. 
