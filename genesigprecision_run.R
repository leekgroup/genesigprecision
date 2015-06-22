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

# Load necessary libraries and data
library(Biobase)
library(BatchJobs)
library(xtable)
library(genefu)
library(rpart)
load("genesigprecision_data.Rda")

# Load a series of parallelized functions for different aspects of the simulation
source("par_fun.R")
source("par_fun_eset.R")
source("make_table.R")

# This block will kick off 100 jobs on your cluster
reg <- makeRegistry(id="full",seed=10284)
ids <- batchMap(reg, par_fun, 10289:10388)
done <- submitJobs(reg, wait=function(retries) 100, max.retries=10)

waitForJobs(reg) # This will wait until all jobs return with results and
		     # provides a progress bar

y <- loadResults(reg) # Aggregate results

make_table(y)

# This function prepares data for simulation run
# It relabels the relevant data columns and makes MammaPrint predictions using the genefu package
prep_eset_sim <- function(dat){
	# Need annotation to make MammaPrint predictions via genefu package
	colnames(fData(dat)) <- c("EntrezGene.ID", "Symbol")
	pd <- pData(dat)
	pd$mammaprint <- gene70(t(exprs(dat)), fData(dat), do.mapping=T)$risk
	# Indicators for grade
	pd$g2ind <- ifelse(pd$grade == 2, 1, 0)
	pd$g3ind <- ifelse(pd$grade == 3, 1, 0)
	pd$y <- pd$surv5
	pd
}

# Dataset GSE19615
dat_19615 <- prep_eset_sim(GSE19615)

# This block will kick off 100 jobs on your cluster
reg <- makeRegistry(id="GSE19615",seed=10284)
ids <- batchMap(reg, par_fun_eset, 10289:10388, more.args=list(pd=dat_19615))
done <- submitJobs(reg, wait=function(retries) 100, max.retries=10)

waitForJobs(reg) # This will wait until all jobs return with results and
		     # provides a progress bar

y <- loadResults(reg) # Aggregate results

make_table(y)

# Dataset GSE11121
dat_11121 <- prep_eset_sim(GSE11121)

# This block will kick off 100 jobs on your cluster
reg <- makeRegistry(id="GSE11121",seed=10284)
ids <- batchMap(reg, par_fun_eset, 10289:10388, more.args=list(pd=dat_11121))
done <- submitJobs(reg, wait=function(retries) 100, max.retries=10)

waitForJobs(reg) # This will wait until all jobs return with results and
		     # provides a progress bar

y <- loadResults(reg) # Aggregate results

make_table(y)

# Dataset GSE7390
dat_7390 <- prep_eset_sim(GSE7390)
dat_7390 <- dat_7390[,-which(is.na(pData(dat_7390)$grade))]

# This block will kick off 100 jobs on your cluster
reg <- makeRegistry(id="GSE7390",seed=10284)
ids <- batchMap(reg, par_fun_eset, 10289:10388, more.args=list(pd=dat_7390))
done <- submitJobs(reg, wait=function(retries) 100, max.retries=10)

waitForJobs(reg) # This will wait until all jobs return with results and
		     # provides a progress bar

y <- loadResults(reg) # Aggregate results

make_table(y)

# Here, we'll do the permuted example to show no gain when Y and W are uncorrelated.

# This will kick off 100 jobs on your cluster.
reg <- makeRegistry(id="full_permute",seed=10284)
ids <- batchMap(reg, par_fun, 31389:31488, more.args=list(rand=T))
done <- submitJobs(reg, wait=function(retries) 100, max.retries=10)

waitForJobs(reg)

y <- loadResults(reg)

make_table(y)