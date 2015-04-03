# genesigprecision
R Code for reproducing analysis for "Measuring the contribution of genomic predictors through estimator precision gain"

### R packages
Required packages:
BatchJobs, rpart, rattle, xtable 

### File Description
* functions.R contains the code necessary for estimator adjustment, as describe in the "Methods" section of the paper.
* genesigprecision_data.Rda contains the 296-patient dataset referred to in the paper.
* par_fun.R contains the parallelized simulation function that computes the adjusted estimators 100 times. This is run 100 times each via BatchJobs to attain 10,000 total simulations.
* genesigprecision_run.R contains code to run the analysis and reproduce the results in the paper.
* .BatchJobs.R and simple.tmpl are configuration files necessary for running the parallelized job on a cluster.

### Instructions

NOTE: This analysis relies on a Sun Grid Engine cluster. We provide the configuration files and code necessary to reproduce this analysis
on an SGE platform. Please examine the code in genesigprecision_run.R before you run it because it contains the commands to create a registry
and submit SGE batch jobs via the BatchJobs R package.

1. Install the required packages.
2. Download all files in this github repository and keep them together in a local directory on your cluster.
3. Run the code from genesigprecision_run.R.
	* Please keep an eye for "submitJobs" commands in this file - these submit the batch jobs and are at the mercy of the cluster.