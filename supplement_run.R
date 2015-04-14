#'
#' Code to recreate results from Supplement for "Measuring the contribution of genomic predictors through estimator precision gain"
#'

# load necessary libraries and dataset
load("genesig_supplement_data.Rda")
library(Biobase)
library(BatchJobs)
library(xtable)
library(genefu)
source("par_fun_eset.R")

# This function runs the same simulation routine for a given dataset
run_eset_sim <- function(dat, id){

	# Need annotation to make MammaPrint predictions via genefu package
	colnames(fData(dat)) <- c("EntrezGene.ID", "Symbol")
	pd <- pData(dat)
	pd$mammaprint <- gene70(t(exprs(dat)), fData(dat), do.mapping=T)$risk
	# Indicators for grade
	pd$g2ind <- ifelse(pd$grade == 2, 1, 0)
	pd$g3ind <- ifelse(pd$grade == 3, 1, 0)
	pd$y <- pd$surv5

	# Parallelized run of 10,000 simulations
	reg <- makeRegistry(id=id,seed=10284)
	ids <- batchMap(reg, par_fun_eset, 10734:10833, more.args=list(pd=pd))
	done <- submitJobs(reg, wait=function(retries) 100, max.retries=10)

	# To not step on the registry's toes while it's initializing
      Sys.sleep(5)

	# Wait until jobs are done
      jobs_done <- FALSE
      while(!jobs_done){
		jobs_done <- waitForJobs(reg,progressbar=FALSE)
      }

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

	cbind(means[,1], vars[,1], means[,2], vars[,2], means[,3], vars[,3], rot, col)
}

dat <- GSE19615

final_eset_19615 <- run_eset_sim(dat, "eset_19615")

dat <- GSE11121

final_eset_11121 <- run_eset_sim(dat, "eset_11121")

dat <- GSE7390
# Need to drop 2 patients missing grade
dat <- dat[,-which(is.na(pData(dat)$grade))]

final_eset_7390 <- run_eset_sim(dat, "eset_7390")

xtable(final_eset_19615, digits=5)
xtable(final_eset_11121, digits=5)
xtable(final_eset_7390, digits=5)




