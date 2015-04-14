#'
#' Code to reproduce analysis from "Measuring the contribution of genomic predictors through estimator precision gain"
#'
#' This script contains code necessary to reproduce figures and tables presented
#' in the stated paper.
#'
#' Please ensure that this code is run from the same directory that
#' contains the accompanying files available on github.
#'
#' Please also run this code carefully. There are commands to kick off
#' jobs in an SGE cluster environment. We've put in comments where these
#' occur, but please be mindful of any "submitJobs" commands.
#'

# Load necessary libraries

library(rpart)
library(BatchJobs)
library(xtable)
load("genesigprecision_data.Rda")

source("par_fun.R")

# This block will kick off 100 jobs on your cluster
reg <- makeRegistry(id="full",seed=10284)
ids <- batchMap(reg, par_fun, 10289:10388)
done <- submitJobs(reg, wait=function(retries) 100, max.retries=10)

waitForJobs(reg) # This will wait until all jobs return with results and
		     # provides a progress bar

y <- loadResults(reg) # Aggregate results

outmat <- do.call(rbind, y)
means <- matrix(colMeans(outmat), 4, 3, byrow=T)
vars <- apply(outmat, 2, var)

col <- rot <- vector("numeric", 4) # These are the col and rot approximations
for(i in 1:4){
	idx <- i + 2*(i-1)
	rot[i] <- (vars[idx] - vars[idx + 1])/vars[idx]
	col[i] <- (vars[idx] - vars[idx + 2])/vars[idx]
}

vars <- matrix(vars, 4, 3, byrow=T)

final <- cbind(means[,1], vars[,1], means[,2], vars[,2], means[,3], vars[,3], rot, col)

xtable(final)

# Generate histogram figure of variability of gain approximations

source("par_fun_internal.R")

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


# Here, we'll do the permuted example to show no gain when Y and W are uncorrelated.

# This will kick off 100 jobs on your cluster.
reg <- makeRegistry(id="full_permute",seed=10284)
ids <- batchMap(reg, par_fun, 31389:31488, more.args=list(rand=T))
done <- submitJobs(reg, wait=function(retries) 100, max.retries=10)

waitForJobs(reg)

y <- loadResults(reg)

outmat <- do.call(rbind, y)
means <- matrix(colMeans(outmat), 4, 3, byrow=T)
vars <- apply(outmat, 2, var)
col <- rot <- vector("numeric", 4)
for(i in 1:4){
	idx <- i + 2*(i-1)
	rot[i] <- (vars[idx] - vars[idx + 1])/vars[idx]
	col[i] <- (vars[idx] - vars[idx + 2])/vars[idx]
}

vars <- matrix(vars, 4, 3, byrow=T)

final <- cbind(means[,1], vars[,1], means[,2], vars[,2], means[,3], vars[,3], rot, col)

xtable(final, digits=5)