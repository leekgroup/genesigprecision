#'
#' Code to recreate results from Supplement for "Measuring the contribution of genomic predictors through estimator precision gain"
#'

# load necessary libraries and dataset
load("genesigprecision_data.Rda")
library(Biobase)
library(BatchJobs)
library(xtable)
library(genefu)
source("par_fun.R")
source("par_fun_internal.R")
source("par_fun_double.R")
source("par_fun_fewer.R")
source("make_table.R")

# We extend the permuted example concept to examine what happens when we double the sample size

reg <- makeRegistry(id="full_permute_double",seed=10284)
ids <- batchMap(reg, par_fun_double, 31389:31488, more.args=list(rand=T))
done <- submitJobs(reg, wait=function(retries) 100, max.retries=10)

waitForJobs(reg)

y <- loadResults(reg)

make_table(y)

# We extend the permuted example concept to examine what happens when we use fewer covariates

reg <- makeRegistry(id="full_permute_fewer",seed=10284)
ids <- batchMap(reg, par_fun_fewer, 31389:31488, more.args=list(rand=T))
done <- submitJobs(reg, wait=function(retries) 100, max.retries=10)

waitForJobs(reg)

y <- loadResults(reg)

make_table(y)

# Generate histogram figure of variability of gain approximations (moved to supplement)

id <- letters[1:10]
seed <- seq(31389, 31389 + 102*9, 102)
input <- data.frame("id"=id,"seed"=seed, stringsAsFactors=F)
input <- lapply(1:10, function(x) input[x,])

# NOTE: This is going to kick off 10 sets of 100 jobs each (1000 jobs).
# This may take a while or you may run into job limits. Please be
# careful about executing this code!

reg <- makeRegistry(id="full_resample",seed=10284)
ids <- batchMap(reg, par_fun_internal, input)
done <- submitJobs(reg, wait=function(retries) 100, max.retries=10)

waitForJobs(reg)

y <- loadResults(reg)
vars <- lapply(y, function(x){apply(x,2,var)})
g_clin <- lapply(vars, function(x){(x[4]-x[6])/x[4]})
g_cpg <- lapply(vars, function(x){(x[10]-x[12])/x[10]})

par(mfrow=c(1,2), mar = c(5,4.5,4,2))
hist(unlist(g_clin)*100, main="Distribution of Percent Gain, Clinical Only", xlab="% Gain Due to Clinical Factors",cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
hist(unlist(g_cpg)*100, main="Distribution of Percent Gain, Clinical + Genomic", xlab="% Gain Due to Clinical + Genomic Factors",cex.lab=1.5, cex.axis=1.5, cex.main=1.5)