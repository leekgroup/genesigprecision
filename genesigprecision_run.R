#'
#' Code to reproduce analysis from "Measuring the contribution of genomic predictors through estimator precision gain"
#'
#' This script contains code necessary to reproduce figures and tables presented
#' in the stated paper.
#'
#' Please ensure that this code is run from the same directory that
#' contains the accompanying files available on github.

# Load necessary libraries/data
library(rpart)
library(rattle)
library(BatchJobs)
library(xtable)
load("genesigprecision_data.Rda")
set.seed(42341)

# First, we generate an example decision tree- Figure 1 in the paper
sampfun <- function(pd){ pd[sample(nrow(pd), replace=TRUE),]}
pdtmp <- sampfun(pd) # We sample the data as per the scheme selected
mod_clin_no_er <- rpart(y~Characteristics.Age + Characteristics.TumorSize + g2ind + g3ind, data=pdtmp)
mod_clin_no_er <- prune(mod_clin_no_er, cp=mod_clin_no_er$cptable[which.min(mod_clin_no_er$cptable[,"xerror"]),"CP"])
drawTreeNodes(mod_clin_no_er)

# Next, we provide code to run the two simulations described- the complete label
# permutation (Table 2) and the resampling simulation (Table 3). We also 
# generate sample size calculations for Table 4.

# Note: Cluster behavior can be unpredictable and thus all desired jobs
# may not finish on a particular run. Due to the use of the same list of 
# random seeds, if all jobs finish then the table results will match those
# presented in the paper.

source("par_fun.R")

# This section will produce the permuted sampling results

reg <- makeRegistry(id="mrun_rand")
ids <- batchMap(reg, fun=par_fun, 40193:40292, more.args=list(rand=T))

# NOTE: This is the command that kicks off the batch jobs. You will need to wait until it is finished to load all of the results.
done <- submitJobs(reg, wait=function(retries) 100, max.retries=10)

# (not run) You may use this line to wait until the system has run all your
# jobs before proceeding. Note that this may get stuck if a job fails or times out for
# a cluster-related reason.

# while(showStatus(reg)$done < showStatus(reg)$n){Sys.sleep(5)}

y <- loadResults(reg)

outmat <- do.call(rbind, y)
means <- matrix(colMeans(outmat), 5, 3, byrow=T)
vars <- apply(outmat, 2, var)
col <- rot <- vector("numeric", 5)
for(i in 1:5){
	idx <- i + 2*(i-1)
	rot[i] <- (vars[idx] - vars[idx + 1])/vars[idx]
	col[i] <- (vars[idx] - vars[idx + 2])/vars[idx]
}

vars <- matrix(vars, 5, 3, byrow=T)

final <- cbind(means[,1], vars[,1], means[,2], vars[,2], means[,3], vars[,3], rot, col)
colnames(final) <- c("mu_una", "s_una", "mu_rot", "s_rot", "mu_col", "s_col", "g_rot", "g_col")
rownames(final) <- c("no_age", "age", "gen", "c_g", "cg")

xtable(final, digits=5)

# This section will produce the resampling simulation results

reg <- makeRegistry(id="mrun")
ids <- batchMap(reg, fun=par_fun, 40193:40292)

# NOTE: This is the command that kicks off the batch jobs. You will need to wait until it is finished to load all of the results.
done <- submitJobs(reg, wait=function(retries) 100, max.retries=10)

# (not run) You may use this line to wait until the system has run all your
# jobs before proceeding. Note that this may get stuck if a job fails or times out for
# a cluster-related reason.

# while(showStatus(reg)$done < showStatus(reg)$n){Sys.sleep(5)}

y <- loadResults(reg)
outmat <- do.call(rbind, y)
means <- matrix(colMeans(outmat), 5, 3, byrow=T)
vars <- apply(outmat, 2, var)
col <- rot <- vector("numeric", 5)
for(i in 1:5){
	idx <- i + 2*(i-1)
	rot[i] <- (vars[idx] - vars[idx + 1])/vars[idx]
	col[i] <- (vars[idx] - vars[idx + 2])/vars[idx]
}

vars <- matrix(vars, 5, 3, byrow=T)

final <- cbind(means[,1], vars[,1], means[,2], vars[,2], means[,3], vars[,3], rot, col)
colnames(final) <- c("mu_una", "s_una", "mu_rot", "s_rot", "mu_col", "s_col", "g_rot", "g_col")
rownames(final) <- c("no_age", "age", "gen", "c_g", "cg")

xtable(final, digits=5)

1 - 1/(1+rot)
1 - 1/(1+col)

# This section will produce the resampling simulation results without decision-tree modeling of the covariates (for comparison)

source("par_fun_nomod.R")

reg <- makeRegistry(id="mrun_nomod")
ids <- batchMap(reg, fun=par_fun_nomod, 40193:40292)

# NOTE: This is the command that kicks off the batch jobs. You will need to wait until it is finished to load all of the results.
done <- submitJobs(reg, wait=function(retries) 100, max.retries=10)

# (not run) You may use this line to wait until the system has run all your
# jobs before proceeding. Note that this may get stuck if a job fails or times out for
# a cluster-related reason.

# while(showStatus(reg)$done < showStatus(reg)$n){Sys.sleep(5)}

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
colnames(final) <- c("mu_una", "s_una", "mu_rot", "s_rot", "mu_col", "s_col", "g_rot", "g_col")
rownames(final) <- c("no_age", "age", "gen", "cg")

xtable(final, digits=5)