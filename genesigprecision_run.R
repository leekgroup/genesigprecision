#'
#' Code to reproduce analysis from "Measuring the contribution of genomic predictors through estimator precision gain"
#'
#' This script contains code necessary to reproduce figures and tables presented
#' in the stated paper.
#'
#' Please ensure that this code is run from the same directory that
#' contains the accompanying files available on github.

# Load necessary libraries

library(rpart)
library(BatchJobs)
library(xtable)
load("genesigprecision_data.Rda")

source("par_fun.R")

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

clin_out <- cpg_out <- vector("numeric", 100)
for(i in 1:100){
	submat <- outmat[-sample(1:nrow(outmat), nrow(outmat)/4, replace=F),]
	tmp <- apply(submat, 2, var)
	clin_out[i] <- (tmp[4] - tmp[6])/tmp[4]
	cpg_out[i] <- (tmp[10] - tmp[12])/tmp[12]	
}

par(mfrow=c(1,2))
hist(clin_out*100, main="Distribution of Percent Gain, Clinical Only", xlab="% Gain Due to Clinical Factors")
hist(cpg_out*100, main="Distribution of Percent Gain, Clinical + Genomic", xlab="% Gain Due to Clinical + Genomic Factors")

# Here, we'll do the permuted example to show no gain when Y and W are uncorrelated.

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